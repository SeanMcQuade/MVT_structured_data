#!/usr/bin/env python3
"""Check that mvtpy reproduces MATLAB's JSON output byte for byte.

Two modes, both run against files the MATLAB pipeline already produced:

``tokens``  (default, memory-light)
    Scan the file for numeric literals, parse each one, re-encode it with
    :func:`mvtpy.matjson.encode_number`, and compare the text. This exercises
    the number formatter over millions of real values without decoding the
    whole document, which matters because a released segment is ~900 MB.

``file``
    Decode a (small) JSON document and re-encode the whole thing, comparing
    bytes. Use it on fixtures and on the prefix extracted by ``--records``.

Examples
--------
    python tools/json_parity.py ../../results/slim/2022-11-17/I-24MOTION_2022-11-17_07-59-59.json
    python tools/json_parity.py --mode file --records 200 <same file>
    python tools/json_parity.py --mode tokens --limit 5000000 <file>
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import matjson  # noqa: E402

# Numeric literals, ignoring anything inside a string.
_STRING = re.compile(r'"(?:[^"\\]|\\.)*"')
_NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def check_tokens(path: Path, limit: int, chunk_bytes: int) -> int:
    """Compare every numeric literal against our re-encoding of it."""
    text = path.read_text(encoding="utf-8", errors="strict") if chunk_bytes <= 0 else _read_prefix(path, chunk_bytes)
    text = _STRING.sub('""', text)

    checked = 0
    mismatches = 0
    for match in _NUMBER.finditer(text):
        token = match.group(0)
        ours = matjson.encode_number(float(token))
        if ours != token:
            mismatches += 1
            if mismatches <= 20:
                print(f"  mismatch at offset {match.start()}: MATLAB {token!r} != mvtpy {ours!r}")
        checked += 1
        if limit and checked >= limit:
            break

    print(f"{path.name}: {checked} numeric tokens checked, {mismatches} mismatched")
    return mismatches


def check_file(path: Path, records: int) -> int:
    """Decode and re-encode a document (or its first `records` elements)."""
    text = path.read_text(encoding="utf-8")
    if records:
        text = extract_prefix_records(text, records)

    decoded = matjson.loads(text)
    ours = matjson.dumps(decoded)

    if ours == text:
        print(f"{path.name}: byte-identical after decode/encode ({len(text)} bytes)")
        return 0

    index = _first_difference(text, ours)
    print(f"{path.name}: MISMATCH at byte {index}")
    print(f"  MATLAB: ...{text[max(0, index - 60):index + 60]}...")
    print(f"  mvtpy : ...{ours[max(0, index - 60):index + 60]}...")
    return 1


def extract_prefix_records(text: str, records: int) -> str:
    """Return a valid JSON array holding the first `records` top-level elements.

    Used to carve a small, self-consistent fixture out of a multi-hundred-MB
    released segment without decoding the whole file.
    """
    if not text.startswith("["):
        raise ValueError("expected a top-level JSON array")

    depth = 0
    count = 0
    in_string = False
    escaped = False

    for index, char in enumerate(text):
        if in_string:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                in_string = False
            continue
        if char == '"':
            in_string = True
        elif char in "[{":
            depth += 1
        elif char in "]}":
            depth -= 1
            if depth == 1:
                count += 1
                if count >= records:
                    return text[: index + 1] + "]"
            elif depth == 0:
                return text[: index + 1]
    raise ValueError("unterminated JSON array")


def _read_prefix(path: Path, chunk_bytes: int) -> str:
    with path.open("rb") as handle:
        raw = handle.read(chunk_bytes)
    # Drop a possibly truncated trailing token.
    return raw.decode("utf-8", errors="ignore").rsplit(",", 1)[0]


def _first_difference(left: str, right: str) -> int:
    for index, (a, b) in enumerate(zip(left, right)):
        if a != b:
            return index
    return min(len(left), len(right))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+", type=Path, help="MATLAB-produced JSON files")
    parser.add_argument("--mode", choices=("tokens", "file"), default="tokens")
    parser.add_argument("--limit", type=int, default=0,
                        help="stop after this many numeric tokens (0 = all)")
    parser.add_argument("--chunk-bytes", type=int, default=64 * 1024 * 1024,
                        help="bytes to read in tokens mode (0 = whole file)")
    parser.add_argument("--records", type=int, default=0,
                        help="in file mode, only compare the first N array elements")
    args = parser.parse_args()

    failures = 0
    for path in args.paths:
        if args.mode == "tokens":
            failures += check_tokens(path, args.limit, args.chunk_bytes)
        else:
            failures += check_file(path, args.records)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
