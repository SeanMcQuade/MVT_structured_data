"""Streaming reader for the raw I-24 MOTION segment files.

A raw 10-minute segment is a single JSON array of trajectory objects, roughly
2 GB on disk. Decoding one with ``json.loads`` needs many gigabytes of RAM, and
the lane-identification stage has to traverse the file twice: once to estimate
the driving line (a whole-file statistic) and once to assign lanes per
trajectory. This module yields one decoded trajectory at a time instead, so
memory stays bounded by the largest single record.

The scanner tracks nesting depth using only the structural characters, skipping
string contents (and escapes within them), then hands each complete top-level
element to ``json.loads``.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Callable, Iterator, Optional

from .matjson import loads as _loads

__all__ = ["iter_trajectories", "count_trajectories", "DEFAULT_CHUNK_BYTES"]

#: Read size per I/O call. Large enough to amortize syscalls, small enough that
#: the working set stays modest.
DEFAULT_CHUNK_BYTES = 8 * 1024 * 1024

_STRUCTURAL = re.compile(r'["\\{}\[\]]')


def iter_trajectories(
    path: "str | Path",
    chunk_bytes: int = DEFAULT_CHUNK_BYTES,
    decoder: Optional[Callable[[str], Any]] = None,
) -> Iterator[Any]:
    """Yield each top-level element of a JSON array file, one at a time.

    Parameters
    ----------
    path:
        File containing a single JSON array (the raw MOTION segment format).
    chunk_bytes:
        Bytes per read.
    decoder:
        Callable applied to each element's text. Defaults to
        :func:`mvtpy.matjson.loads`, which decodes every number as a float so
        that ``-0`` survives.

    Notes
    -----
    Elements are decoded lazily, so a caller that only needs a few fields can
    stop early without paying for the rest of the file.
    """
    decode = decoder or _loads
    path = Path(path)

    buffer = ""
    buffer_start = 0      # absolute-ish offset of buffer[0] within the scan
    depth = 0
    in_string = False
    escaped = False
    array_opened = False  # the file's own enclosing '[' is not an element
    element_start: Optional[int] = None
    position = 0          # scan position within `buffer`

    with path.open("r", encoding="utf-8") as handle:
        while True:
            chunk = handle.read(chunk_bytes)
            if not chunk:
                break
            buffer += chunk

            for match in _STRUCTURAL.finditer(buffer, position):
                index = match.start()
                char = match.group(0)

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
                elif char == "[" and not array_opened and depth == 0:
                    # Enclosing array of the file itself; elements start next.
                    array_opened = True
                elif char in "{[":
                    if depth == 0:
                        element_start = index
                    depth += 1
                elif char in "}]":
                    if depth == 0:
                        # The array's own closing bracket.
                        continue
                    depth -= 1
                    if depth == 0 and element_start is not None:
                        yield decode(buffer[element_start:index + 1])
                        element_start = None

            position = len(buffer)

            # Drop everything before the element currently being accumulated,
            # so the buffer never grows past one record plus a chunk.
            if depth == 0 and not in_string:
                buffer_start += len(buffer)
                buffer = ""
                position = 0
            elif element_start is not None and element_start > 0:
                buffer = buffer[element_start:]
                position -= element_start
                buffer_start += element_start
                element_start = 0

    if depth != 0:
        raise ValueError(f"unterminated JSON structure in {path}")


def count_trajectories(path: "str | Path", chunk_bytes: int = DEFAULT_CHUNK_BYTES) -> int:
    """Count top-level elements without fully decoding them."""
    return sum(1 for _ in iter_trajectories(path, chunk_bytes, decoder=lambda text: None))
