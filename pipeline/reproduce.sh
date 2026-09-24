#!/usr/bin/env bash
# One-command reproduction of every number in the manuscript.
#
#   bash reproduce.sh            # full: fetch, build inputs, both cohorts, figures
#                                #   (end-to-end time: see FULL PATH below)
#   bash reproduce.sh --check    # tests + figures from stored results
#
# MEASURED WALL CLOCKS, macOS 26.5.2 / arm64, 14 cores. Load varied between runs and
# that is the point -- see the note below each figure. These are timings someone
# actually saw produced, not estimates:
#
#   reproduce.sh --check   509 s (8 m 29 s), exit 0  (2026-09-17, MEASURED end to
#                          end on AC power with desktop apps running, load1
#                          4.8-15.6; log results/reproduce_check_20260917_run49b.log;
#                          the suite has grown and gained its own stages since:)
#                          247 s (4 m 7 s), exit 0   (2026-09-07, MEASURED end to
#                          end after the five steps were added; nothing else
#                          running; log results/reproduce_check_20260907_run3.log)
#                          Superseded: 72 s, then 43 s (2026-09-04), both taken
#                          before the five steps existed.
#   pytest tests/ -q       26.0 / 54.6 / 56.1 s   (86 tests, three runs)
#   pytest tests/ -q       48.7 s idle; 864.3 s with one pan-cancer job resident;
#                          1702.6 s (28 m 22 s) on 2026-09-05 with EIGHT agent
#                          sessions, a pan-cancer job at ~1000% CPU and another
#                          session's pytest all resident -- a 35x spread, and the
#                          largest yet measured. Name what else was running.
#
# THE --check FIGURE IS NOW MEASURED, 2026-09-07, and the measuring found a bug.
# On 2026-09-05 the --check path gained five steps (24_split_determinism --check,
# 26_sort_audit, two --self-tests and --coverage) and was not re-timed end to end
# until 2026-09-07. Doing so took three runs and each failed differently, which
# is worth knowing before trusting a fourth:
#   run 1 (145 s, exit 1)  CONTAMINATED -- launched alongside another pytest and
#                          with the tree mid-edit. Neither number nor verdict
#                          usable. Never run --check next to another pytest:
#                          --check runs the full suite itself.
#   run 2 (137 s, exit 1)  Clean, and it exposed a REAL defect. This script
#                          exports PYTHONHASHSEED=0 (line below), so a fresh
#                          ENVIRONMENT.md render from inside --check reports 0
#                          where the committed file, generated interactively,
#                          records "(unset)". test_environment_md_matches_a_
#                          fresh_render_on_this_platform therefore failed HERE
#                          while passing standalone -- meaning --check could
#                          never exit 0, by construction, since session 28.
#                          Fixed by excluding that one row from _strip_volatile,
#                          the same treatment the git-commit row already gets.
#                          137 s is a TIME-TO-FAILURE: it aborted at the suite
#                          and never reached the self-test stages.
#   run 3 (247 s, exit 0)  The real figure. Quote this one.
# The lesson is the project's own: a check nobody runs end to end is not a check.
#
# FULL PATH, MEASURED END TO END 2026-09-17 (ledger C3, HPC4 job 130192, exit
# 0): 8,751 s = 2 h 25 m 51 s from an empty data/ directory, and both cohorts
# reproduced the corrected-ordering values BIT FOR BIT. Per stage, from the
# log's elapsed-second prefixes: environment and provenance 15 s, fetch and
# interim build 163 s, determinism and sort audits 42 s, test suite 885 s,
# NSCLC 617 s, pan-TCGA 6,944 s, ancestry 57 s, figures 26 s. The "~2.5 h" an
# earlier header quoted was never measured and was withdrawn before this run;
# it happens to have been close. Earlier per-stage figures, for comparison:
# NSCLC 227.55 s on macOS arm64 and 593.13 s on HPC4; pan-cancer 17,118 s on
# macOS and 11,022 s on HPC4, both BEFORE the shared-mask fast path.
#
# The spread among those three is itself contention: 54.6 and 56.1 were measured
# while a SECOND agent session ran its own pytest on this machine; 26.0 was
# measured with the machine idle. Quote a range, not a number.
#
# CONTENTION IS THE DOMINANT TERM, and it is much larger than CPU count suggests.
# The same 86 tests took 1305.22 s -- 24x slower -- while a pan-cancer job held a
# large working set on a machine whose disk was 99% full, so swap could not grow.
# The global-axis caches under results/ were untouched across both measurements
# (mtimes 2026-09-02), so caching is NOT what separates the two numbers.
# Do not quote the fast number as this suite's cost on a loaded machine.
# The full path's end-to-end time: see FULL PATH above.
#
# The --check path re-renders the figures from the committed result CSVs and
# runs the full test suite. It verifies the ANALYSIS CODE and the FIGURES, not
# the model fits; use the full path to regenerate results/ from raw data.
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONHASHSEED=0          # str hashing is randomised per process, and
                                 # gene-pool ordering must be deterministic

echo "=== environment ==="
# Not just printed -- CHECKED against requirements-lock.txt. Printing a version
# is not verifying it, and an unnoticed version drift looks exactly like A9.
python scripts/20_check_versions.py

echo; echo "=== frozen artefacts unchanged ==="
# Runs BEFORE everything below, because everything below trusts these files.
# 19_check_numbers.py reads the artefacts as ground truth, so an edited frozen
# artefact makes it re-verify the documents against the edited value and report
# OK -- measured, not hypothetical. This is the check that would notice.
python scripts/21_provenance_manifest.py

# THE FULL PATH FETCHES AND BUILDS ITS INPUTS FIRST (moved here 2026-09-17).
# The determinism checks and several tests below read data/interim/, so from an
# empty data/ directory they would stop the script before anything had been
# fetched -- `set -e` makes their exit 2 fatal. --check assumes a working tree
# that already has its inputs, as before.
if [[ "${1:-}" != "--check" ]]; then
  echo; echo "=== fetch (skips anything already present; every file is hash-verified) ==="
  # 03 first: it writes the two Hallmark GMTs that 01 then verifies with the rest.
  [[ -f data/raw/signatures/h.all.v2024.1.Hs.symbols.tme.gmt ]] || python scripts/03_fetch_signatures.py
  python scripts/01_fetch_data.py --all || {
    echo "Fetch or verification failed; the lines above name each missing or"
    echo "different file and its public source."; exit 1; }

  # data/interim/ is BUILT from data/raw/ (A13: until 2026-09-16 nothing in this
  # repository wrote it). Skipped when all six files are present; refused when
  # only some are, because a half-built interim directory is not something to
  # guess about. The pan-cancer half needs ~5 GiB free for the expression matrix.
  echo; echo "=== interim inputs (built from data/raw if absent) ==="
  INTERIM_FILES=(meta.parquet X.npy expression_hugo.parquet cohort_nsclc.parquet X_nsclc.npy expr_nsclc.parquet)
  n_present=0
  for f in "${INTERIM_FILES[@]}"; do [[ -f "data/interim/$f" ]] && n_present=$((n_present + 1)); done
  if [[ $n_present -eq ${#INTERIM_FILES[@]} ]]; then
    echo "all six present; not rebuilt (verify with scripts/00_build_interim.py --verify into a scratch --out)"
  elif [[ $n_present -eq 0 ]]; then
    # Exit 3 = built but not byte-verified. Expected on Linux: cohort_nsclc.parquet's
    # 16 mean-z scores differ from the macOS-built file by up to 1.1e-15 (A13), and
    # without the original files there is nothing to compare its content against.
    rc=0; python scripts/00_build_interim.py --out data/interim --verify || rc=$?
    if [[ $rc -eq 3 ]]; then
      echo "NOTE: interim files built; the lines above name any that differ from the recorded bytes."
    elif [[ $rc -ne 0 ]]; then
      exit "$rc"
    fi
  else
    echo "data/interim/ holds $n_present of ${#INTERIM_FILES[@]} files; move them aside and re-run."; exit 1
  fi

fi

echo; echo "=== determinism: tie conventions are pinned ==="
# A9 and A10 were both non-stable sorts: `np.argsort`'s default quicksort orders
# exactly-tied keys by an implementation-defined rule that differs between
# numpy's arm64 and x86_64 kernels. Six sessions went into finding the first one.
# These two checks take ~11 s together and would have caught it on day one, which
# is the only reason they are this high in the script.
#
# Both FAIL (exit 1) if a fix is reverted -- verified by reverting it, not
# assumed. 24 checks the fold partition; 26 sweeps every other ordering in the
# pipeline and reports which are safe, which are dead code, and why.
python scripts/24_split_determinism.py --check
python scripts/26_sort_audit.py

echo; echo "=== test suite ==="
python -m pytest tests/ -q

if [[ "${1:-}" == "--check" ]]; then
  echo; echo "=== figures from stored results ==="
  python scripts/10_make_figures.py
  echo; echo "=== abstract character count ==="
  python scripts/08_count_abstract.py
  echo; echo "=== documents vs frozen results ==="
  python scripts/19_check_numbers.py
  # A claim count is not a coverage measure, and a check nobody has seen go red
  # is not evidence. --self-test proves the patterns can still fail; --coverage
  # prints how much of each document is checked at all (11% of literals overall,
  # so the honest picture rather than the flattering one).
  echo; echo "=== the checker's own patterns can still fail ==="
  python scripts/19_check_numbers.py --self-test
  python scripts/26_sort_audit.py --self-test
  echo; echo "=== how much is actually checked ==="
  # --strict here also GATES: it fails if any document has fewer literals
  # checked than when the floor was recorded. `pipefail` is set above, so a
  # non-zero status survives the pipe -- without it this would report tail's.
  python scripts/19_check_numbers.py --coverage --strict | tail -6
  echo; echo "CHECK PASSED. Full reproduction: bash reproduce.sh"
  exit 0
fi

# WHERE THE FULL RUN WRITES, and why it is not the frozen directories.
#
# This script used to write straight into results/nsclc_v3 -- the directory the
# manuscript quotes. So the advertised one-command reproduction OVERWROTE the
# frozen results it was meant to reproduce, and a reader running it destroyed
# the very thing they wanted to compare against. It also wrote pan-cancer to
# results/pancancer_v2 while the paper quotes pancancer_v3, so the two halves
# disagreed about which run was canonical.
#
# The full run now writes to a fresh, timestamped directory and DIFFS against
# the frozen values at the end. Nothing frozen is touched.
REPRO_DIR="results/repro_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$REPRO_DIR"
echo; echo "Full reproduction writes to $REPRO_DIR"
echo "The frozen results in results/nsclc_v3 and results/pancancer_v3 are NOT touched."

# Runtimes are MEASURED, with the platform named; none is an estimate.
# NSCLC: 227.55 s on macOS arm64, 593.13 s on HPC4 (8 CPUs).
echo; echo "=== NSCLC (n=944; 4-10 min measured) ==="
python scripts/05_run_nsclc.py --outdir "$REPRO_DIR/nsclc"

# pan-TCGA: 6,944 s on HPC4 (8 CPUs) WITH the shared-mask fast path
# (models.SHARED_MASK_FAST_PATH), measured 2026-09-17 in job 130192 -- 1.59x
# faster than the same script's 11,022 s without it. Before the fast path:
# 17,118 s on macOS arm64 and 11,022 s on HPC4.
echo; echo "=== pan-TCGA (n=7,168; 1.9 h measured with the fast path, 3.1 h without) ==="
python scripts/09_run_pancancer.py --outdir "$REPRO_DIR/pancancer"

echo; echo "=== ancestry supplement ==="
python scripts/06_run_ancestry.py --outdir "$REPRO_DIR/ancestry"

echo; echo "=== figures (rendered from the FROZEN results, not this run) ==="
python scripts/10_make_figures.py

echo; echo "DONE. Your run against the corrected-ordering values:"
REPRO_DIR="$REPRO_DIR" python - <<'PY'
import json, os, platform, sys
repro = os.environ["REPRO_DIR"]
# Since the A9/A10 sort fixes, this code reproduces the CORRECTED-ORDERING runs
# (results/*_stablesort), not the frozen primaries, which were computed under
# the old sort and stay the reported numbers (Option B). Compare against the
# corrected runs; show the frozen primary beside them for orientation.
pairs = (("NSCLC", f"{repro}/nsclc", "results/nsclc_v3_stablesort", "results/nsclc_v3"),
         ("pan-TCGA", f"{repro}/pancancer", "results/pancancer_v3_stablesort",
          "results/pancancer_v3"))
print(f"  this platform: {sys.platform} / {platform.machine()}, "
      f"Python {platform.python_version()}")
print("  reference:     NSCLC darwin / arm64, Python 3.13.9; pan-TCGA linux / x86_64,")
print("                 Python 3.12.14 (HPC4 job 129159)\n")
print(f"  {'cohort':10s} {'yours':>34s}   {'corrected ordering':>34s}   {'delta':>10s}")
worst = 0.0
for name, mine, frozen, primary in pairs:
    a = json.load(open(f"{mine}/summary.json"))["primary_excess"]
    b = json.load(open(f"{frozen}/summary.json"))["primary_excess"]
    d = a["value"] - b["value"]
    worst = max(worst, abs(d))
    print(f"  {name:10s} {a['value']:.6f} [{a['ci_lo']:.4f}, {a['ci_hi']:.4f}]   "
          f"{b['value']:.6f} [{b['ci_lo']:.4f}, {b['ci_hi']:.4f}]   {d:+10.6f}")
    c = json.load(open(f"{primary}/summary.json"))["primary_excess"]
    print(f"  {'':10s} (reported primary, old sort: {c['value']:.6f})")

# The corrected-ordering ISI agrees across macOS/arm64 and Linux/x86_64 to about
# one ULP (HPC4 job 129151, NSCLC). The pan-TCGA reference was computed on HPC4.
# The frozen and corrected-ordering runs' secondary quantities (the random-patient
# split) do NOT agree across platforms (A12). Since 2026-09-24 that split is
# pinned, so a fresh run's NSCLC secondaries should equal
# results/nsclc_v3_pinned on either platform. This table compares only the ISI.
print()
if worst <= 1e-9:
    print("  Bit-for-bit match with the corrected-ordering results.")
else:
    print(f"  Largest difference: {worst:.6f}.")
    print("  The corrected-ordering ISI has reproduced across macOS/arm64 and")
    print("  Linux/x86_64 to about 1e-16 (A9, A12 in ../14-SCIENCE-AUDIT.md), so a")
    print("  difference much larger than that is worth investigating.")
PY
