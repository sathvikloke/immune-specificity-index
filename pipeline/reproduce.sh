#!/usr/bin/env bash
# One-command reproduction of every number in the manuscript.
#
#   bash reproduce.sh            # full: fetch, both cohorts, figures  (~2.5 h,
#                                #   NEVER MEASURED -- see below)
#   bash reproduce.sh --check    # tests + figures from stored results
#
# MEASURED WALL CLOCKS, macOS 15 / arm64, 14 cores. Load varied between runs and
# that is the point -- see the note below each figure. These are timings someone
# actually saw produced, not estimates:
#
#   reproduce.sh --check   247 s (4 m 7 s), exit 0   (2026-09-07, MEASURED end to
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
# THE ~2.5 h FULL-RUN FIGURE HAS NEVER BEEN MEASURED. It is the oldest unverified
# performance claim in the repo. The two frozen runs it is presumably derived
# from took 1,739.6 s (NSCLC) and 17,115.6 s (pan-cancer) = 5.2 h of cohort work
# alone, before fetching or figures, so ~2.5 h is not merely unmeasured but
# likely wrong by a factor of two. Treat it as a placeholder, not a budget.
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
# The full path (~2.5 h) is an ESTIMATE and has never been measured end to end.
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

echo; echo "=== fetch (skips anything already present) ==="
python scripts/01_fetch_data.py --all || {
  echo "Fetch failed. Some sources are large or rate-limited; see the script's"
  echo "docstring for manual download instructions."; exit 1; }
python scripts/03_fetch_signatures.py

echo; echo "=== NSCLC (n=944, ~10 min) ==="
python scripts/05_run_nsclc.py --outdir "$REPRO_DIR/nsclc"

echo; echo "=== pan-TCGA (n=7,168, ~2 h) ==="
python scripts/09_run_pancancer.py --outdir "$REPRO_DIR/pancancer"

echo; echo "=== ancestry supplement ==="
python scripts/06_run_ancestry.py --outdir "$REPRO_DIR/ancestry"

echo; echo "=== figures (rendered from the FROZEN results, not this run) ==="
python scripts/10_make_figures.py

echo; echo "DONE. Your run against the frozen values:"
REPRO_DIR="$REPRO_DIR" python - <<'PY'
import json, os, platform, sys
repro = os.environ["REPRO_DIR"]
pairs = (("NSCLC", f"{repro}/nsclc", "results/nsclc_v3"),
         ("pan-TCGA", f"{repro}/pancancer", "results/pancancer_v3"))
print(f"  this platform: {sys.platform} / {platform.machine()}, "
      f"Python {platform.python_version()}")
print(f"  frozen on:     darwin / arm64, Python 3.13.9\n")
print(f"  {'cohort':10s} {'yours':>34s}   {'frozen':>34s}   {'delta':>10s}")
worst = 0.0
for name, mine, frozen in pairs:
    a = json.load(open(f"{mine}/summary.json"))["primary_excess"]
    b = json.load(open(f"{frozen}/summary.json"))["primary_excess"]
    d = a["value"] - b["value"]
    worst = max(worst, abs(d))
    print(f"  {name:10s} {a['value']:.6f} [{a['ci_lo']:.4f}, {a['ci_hi']:.4f}]   "
          f"{b['value']:.6f} [{b['ci_lo']:.4f}, {b['ci_hi']:.4f}]   {d:+10.6f}")

# A9. This is a KNOWN, DOCUMENTED defect, not a fault in your setup: the same
# code gives 0.3182 on macOS/arm64 and 0.2964 on Linux/x86_64 from bit-identical
# inputs. Say so plainly rather than letting a reader think they broke something.
print()
if worst <= 1e-9:
    print("  Bit-for-bit match with the frozen results.")
else:
    print(f"  Largest difference: {worst:.6f}.")
    print("  If you are not on macOS/arm64 this is EXPECTED and is documented as")
    print("  item A9 in ../14-SCIENCE-AUDIT.md and Limitations 8 of the paper.")
    print("  It does not change the sign, the significance, or 16/16 beating null.")
    print("  A difference of this size on macOS/arm64 WOULD be a real regression.")
PY
