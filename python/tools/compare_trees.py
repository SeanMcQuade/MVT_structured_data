#!/usr/bin/env python3
"""Compare two result trees file by file, and say which fields moved.

Answers "what actually changed?" between two builds of the same product - a
rebuild after a code change, or the same day built on two platforms. Digests
alone say *whether* files differ; the ``--fields`` pass says *which numbers*
differ and by how much, which is the part that decides whether a change is
acceptable.

Two passes, cheap first:

``md5``  (default)
    Hash the corresponding files in both trees and report matches, differences,
    and files present in only one. Streams in chunks, so a 700 MB segment costs
    no more memory than a small one.

``fields``
    For the files that differ, decode both and walk the trajectory records to
    find which keys differ, how many records each key differs in, and the
    largest absolute and relative gap. This decodes whole documents, so it is
    memory-hungry (a released ``full`` segment is ~700 MB on disk and several GB
    parsed) - it only runs on files md5 already flagged.

Examples
--------
    # what did the quadrature fix change in `full`?
    python tools/compare_trees.py ../../results_groundtruth ../../results \\
        --product full --day 17

    # same, and name the fields
    python tools/compare_trees.py ../../results_groundtruth ../../results \\
        --product full --day 17 --fields

    # two platforms, whole product
    python tools/compare_trees.py /Volumes/ExtremePro/results ../../results \\
        --product full
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CHUNK = 8 << 20


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def segment_files(root: Path, product: str, days: Iterable[int]) -> List[Path]:
    """Every segment file of `product` under `root`, as tree-relative paths.

    Matches the pipeline's own naming rather than a bare ``*.json`` glob: the
    dataset_info.json sidecar lives at the product root, and globbing it in is
    exactly the bug that broke the micro stage once already.
    """
    if product == "gps":
        # gps has no day folders; select by filename so --day still applies.
        return sorted(
            p.relative_to(root)
            for day in days
            for p in (root / "gps").glob(f"CIRCLES_GPS_*2022-11-{day}.json"))

    found: List[Path] = []
    for day in days:
        day_dir = root / product / f"2022-11-{day}"
        if not day_dir.is_dir():
            continue
        found.extend(sorted(p.relative_to(root) for p in day_dir.glob("I-24*.json")))
    return found


def compare_digests(a_root: Path, b_root: Path, rel_paths: List[Path]) -> Tuple[List[Path], List[Path]]:
    differing: List[Path] = []
    only_one: List[Path] = []
    for i, rel in enumerate(rel_paths, 1):
        a, b = a_root / rel, b_root / rel
        if not a.is_file() or not b.is_file():
            only_one.append(rel)
            print(f"  [{i}/{len(rel_paths)}] MISSING on one side  {rel}")
            continue
        # Size is a free pre-filter: different sizes cannot have equal digests.
        if a.stat().st_size != b.stat().st_size or md5(a) != md5(b):
            differing.append(rel)
            print(f"  [{i}/{len(rel_paths)}] DIFFERS  {rel}")
        else:
            print(f"  [{i}/{len(rel_paths)}] same     {rel}")
        sys.stdout.flush()
    return differing, only_one


def _numeric_gap(x, y) -> Tuple[float, float] | None:
    """Absolute and relative gap, or None when the two are equivalent.

    NaN equals NaN here: MATLAB writes NaN as JSON null and two nulls in the
    same slot are not a difference. The relative term is reported against the
    larger magnitude, and suppressed near zero, because dividing by a
    near-zero denominator manufactures huge ratios out of negligible gaps -
    that mistake produced a bogus "3e-3 relative drift" earlier in this work.
    """
    if isinstance(x, bool) or isinstance(y, bool):
        return None if x == y else (1.0, 1.0)
    if not isinstance(x, (int, float)) or not isinstance(y, (int, float)):
        return None if x == y else (float("nan"), float("nan"))
    if x is None or y is None:
        return None if x == y else (float("nan"), float("nan"))
    if math.isnan(x) and math.isnan(y):
        return None
    if x == y:
        return None
    gap = abs(x - y)
    scale = max(abs(x), abs(y))
    rel = gap / scale if scale > 1e-12 else float("nan")
    return gap, rel


def compare_fields(a_path: Path, b_path: Path) -> Dict[str, Dict[str, float]]:
    """Per-key difference summary between two trajectory documents."""
    with a_path.open() as fh:
        a_doc = json.load(fh)
    with b_path.open() as fh:
        b_doc = json.load(fh)
    if len(a_doc) != len(b_doc):
        return {"<record count>": {"records": abs(len(a_doc) - len(b_doc)),
                                   "max_abs": float("nan"), "max_rel": float("nan")}}

    summary: Dict[str, Dict[str, float]] = {}
    for a_rec, b_rec in zip(a_doc, b_doc):
        for key in set(a_rec) | set(b_rec):
            av, bv = a_rec.get(key), b_rec.get(key)
            pairs = zip(av, bv) if isinstance(av, list) and isinstance(bv, list) else [(av, bv)]
            worst_abs = worst_rel = 0.0
            hit = False
            for x, y in pairs:
                gap = _numeric_gap(x, y)
                if gap is None:
                    continue
                hit = True
                if not math.isnan(gap[0]):
                    worst_abs = max(worst_abs, gap[0])
                if not math.isnan(gap[1]):
                    worst_rel = max(worst_rel, gap[1])
            if hit:
                entry = summary.setdefault(key, {"records": 0, "max_abs": 0.0, "max_rel": 0.0})
                entry["records"] += 1
                entry["max_abs"] = max(entry["max_abs"], worst_abs)
                entry["max_rel"] = max(entry["max_rel"], worst_rel)
    return summary


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare two result trees and report what differs.")
    parser.add_argument("tree_a", type=Path, help="baseline tree (e.g. results_groundtruth)")
    parser.add_argument("tree_b", type=Path, help="tree to compare against it")
    parser.add_argument("--product", default="slim",
                        choices=["slim", "full", "gps"], help="which product (default slim)")
    parser.add_argument("--day", type=int, action="append", choices=[16, 17, 18],
                        help="day to compare; repeatable (default all three)")
    parser.add_argument("--fields", action="store_true",
                        help="for differing files, name the fields that moved "
                             "(decodes whole documents; memory-hungry)")
    args = parser.parse_args(argv)

    days = args.day or [16, 17, 18]
    a_root, b_root = args.tree_a.resolve(), args.tree_b.resolve()
    for root in (a_root, b_root):
        if not root.is_dir():
            parser.error(f"not a directory: {root}")

    rel_paths = segment_files(a_root, args.product, days)
    if not rel_paths:
        print(f"no {args.product} files for day(s) {days} under {a_root}")
        return 2

    print(f"A = {a_root}")
    print(f"B = {b_root}")
    print(f"comparing {len(rel_paths)} {args.product} file(s), day(s) {days}\n")

    differing, only_one = compare_digests(a_root, b_root, rel_paths)
    same = len(rel_paths) - len(differing) - len(only_one)
    print(f"\nidentical={same}  differ={len(differing)}  missing-on-one-side={len(only_one)}")

    if differing and args.fields:
        print("\nfield-level differences:")
        for rel in differing:
            print(f"\n  {rel}")
            summary = compare_fields(a_root / rel, b_root / rel)
            if not summary:
                print("    (digests differ but no field difference found - "
                      "formatting or key order)")
            for key, stat in sorted(summary.items(), key=lambda kv: -kv[1]["records"]):
                print(f"    {key:<48} records={stat['records']:<6} "
                      f"max_abs={stat['max_abs']:.6g}  max_rel={stat['max_rel']:.3g}")

    return 1 if (differing or only_one) else 0


if __name__ == "__main__":
    raise SystemExit(main())
