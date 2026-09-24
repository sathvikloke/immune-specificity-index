#!/usr/bin/env python3
"""Count an AACR abstract against the 2,600-character limit. Fails loudly if over.

    python scripts/08_count_abstract.py ../07-ABSTRACT-DRAFT.md

THE RULE, verbatim from the AACR 2026 Call for Abstracts (p.26):
  "The combined length of the abstract body, title, and tables may not exceed
   2,600 characters, not including spaces and the author string. Tables count for
   800 characters against the limit. Submission cannot be completed for abstracts
   that exceed this limit."

So: spaces are FREE, the author string is FREE, a table costs a flat 800, and
going over is hard-blocked at submission rather than trimmed. Counting by words
or by len() including spaces will both mislead — len() including spaces
overcounts by roughly 17%, which is enough to make a fitting abstract look
over-length and provoke needless cutting.

This reads the TITLE and BODY sections of the draft (everything between the
`## TITLE` / `## BODY` headings and the next `---`), so the notes-to-self in the
same file are not counted.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

LIMIT = 2600
TABLE_COST = 800
NEXTGEN_LIMIT = 8000
DEFAULT_PATH = (Path(__file__).resolve().parents[2] / '07-ABSTRACT-DRAFT.md')


def extract(md: str) -> tuple[str, str]:
    """Pull the TITLE and BODY sections out of the draft."""
    def section(name: str) -> str:
        m = re.search(rf"^## {name}\s*\n(.*?)(?=^## |\Z)", md, re.S | re.M)
        return m.group(1).strip() if m else ""

    body = section("BODY")
    # Drop a trailing horizontal rule. The section runs to the next `## `
    # heading, and the draft closes the body with a `---` rule before its notes;
    # the dry run says to paste "between ## BODY and the next ---", so the rule
    # is never pasted. Until 2026-09-17 it was counted: 3 characters, which made
    # every reported total 3 too high (2,515 for a pasted 2,512).
    body = re.sub(r"\n\s*(?:-{3,}|\*{3,}|_{3,})\s*\Z", "", body).strip()
    # Drop the markdown bold markers used for the required element labels; they
    # are formatting for reading the draft, not characters in the submission.
    body = body.replace("**", "")
    return section("TITLE").replace("\n", " ").strip(), body


def count(text: str) -> int:
    """Characters excluding ALL whitespace, per the CFA rule."""
    return len(re.sub(r"\s", "", text))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", type=Path, nargs="?", default=DEFAULT_PATH)
    ap.add_argument("--tables", type=int, default=0, help="number of tables included")
    ap.add_argument("--nextgen", action="store_true",
                    help="score against the NextGen Stars rule instead: 8,000 "
                         "characters INCLUDING spaces (the opposite convention "
                         "to the regular abstract, which is why this is a flag "
                         "rather than a second script)")
    args = ap.parse_args()

    if args.nextgen:
        path = (args.path if args.path != DEFAULT_PATH
                else DEFAULT_PATH.parent / "13-NEXTGEN-EXTENDED-ABSTRACT.md")
        title, body = extract(path.read_text())
        if not title or not body:
            print(f"Could not find '## TITLE' and '## BODY' in {path}")
            return 1
        total = len(title) + len(body)          # spaces COUNT here
        print("  NextGen Stars rule: 8,000 characters INCLUDING spaces")
        print(f"  title    {len(title):5d}")
        print(f"  body     {len(body):5d}")
        print(f"  {'-' * 16}")
        print(f"  TOTAL    {total:5d} / {NEXTGEN_LIMIT}   ({total / NEXTGEN_LIMIT:.0%})")
        print(f"  headroom {NEXTGEN_LIMIT - total:5d}   ({len(body.split())} words)")
        if total > NEXTGEN_LIMIT:
            print(f"\n  OVER by {total - NEXTGEN_LIMIT}.")
            return 1
        print("\n  Within the limit.")
        return 0

    # A raw FileNotFoundError traceback here is not a diagnosis. `reproduce.sh
    # --check` calls this script, and the PUBLIC SNAPSHOT deliberately does not
    # stage `07-ABSTRACT-DRAFT.md` -- it is a working document with budget notes
    # and open TODOs, not a deliverable. So every reader running the deposit's
    # advertised entry point reaches this line with the file absent. Measured
    # 2026-09-07 inside a freshly staged snapshot: it crashed with a traceback.
    if not args.path.is_file():
        print(f"MISSING {args.path} -- this script counts an abstract draft, "
              "and that draft is not part of the public snapshot.")
        print("  It is a working document, excluded on purpose. If you are "
              "reading this from the deposit, nothing is wrong with the "
              "deposit; this script simply has no input here.")
        return 2

    title, body = extract(args.path.read_text())
    if not title or not body:
        print(f"Could not find '## TITLE' and '## BODY' sections in {args.path}")
        return 1

    n_title, n_body = count(title), count(body)
    n_tables = args.tables * TABLE_COST
    total = n_title + n_body + n_tables

    print(f"  title    {n_title:5d}")
    print(f"  body     {n_body:5d}")
    if args.tables:
        print(f"  tables   {n_tables:5d}  ({args.tables} x {TABLE_COST})")
    print(f"  {'-' * 16}")
    print(f"  TOTAL    {total:5d} / {LIMIT}   ({total / LIMIT:.0%})")
    print(f"  headroom {LIMIT - total:5d}")

    # For sanity: what the naive counts would have said.
    naive = len(title) + len(body) + n_tables
    print(f"\n  (len() INCLUDING spaces would read {naive} — {naive / max(total,1) - 1:+.0%}; "
          "spaces are free under the CFA rule)")
    print(f"  (body is {len(body.split())} words)")

    if total > LIMIT:
        print(f"\n  OVER by {total - LIMIT}. Submission would be hard-blocked.")
        return 1
    print("\n  Within the limit.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
