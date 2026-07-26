"""Map each raw MOTION segment to the processed file it produces.

Port of ``Scripts/+mvt/manifest.m``, ``segmentName.m`` and
``rawFirstTimestamp.m``.

Processed filenames are content-derived: a stage only learns that
``..._fri_0_00.json`` becomes ``I-24MOTION_2022-11-18_05-59-59.json`` after
reading the first trajectory's ``first_timestamp``. The build system needs that
mapping *before* running anything, so that freshness is decided without the
multi-gigabyte decode, and so that per-segment work units can be scheduled.

The timestamp is found with a prefix scan rather than a full parse, and the
whole mapping is cached under ``results/.mvt/manifests``.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Optional
from zoneinfo import ZoneInfo

__all__ = ["Segment", "segment_name", "raw_first_timestamp", "manifest",
           "raw_dir", "DAY_ABBREV"]

_CENTRAL = ZoneInfo("America/Chicago")

#: The raw files embed the weekday of the run, and the stage globs match on it.
DAY_ABBREV = {16: "wed", 17: "thu", 18: "fri"}

_FIRST_TIMESTAMP = re.compile(r'"first_timestamp"\s*:\s*([-+0-9.eE]+)')

#: Prefix sizes tried in order when scanning for the first timestamp.
_PREFIX_BYTES = (4_000_000, 64_000_000)


@dataclass(frozen=True)
class Segment:
    """One raw 10-minute segment and the processed file it becomes."""

    seq: int                 # 0..23, parsed from the raw name
    raw_name: str
    raw_path: Path
    first_timestamp: float
    output_name: str

    def output_path(self, slim_dir: Path) -> Path:
        return Path(slim_dir) / self.output_name


def segment_name(first_timestamp: float) -> str:
    """Processed filename for a segment, from its first trajectory's timestamp.

    Nashville local time at second resolution, e.g.
    ``I-24MOTION_2022-11-18_05-59-59.json``. Released filenames must not
    change, so this is the single place the conversion happens - as in MATLAB,
    where both the stage and the manifest call ``mvt.segmentName``.
    """
    stamp = datetime.fromtimestamp(first_timestamp, _CENTRAL)
    return f"I-24MOTION_{stamp.strftime('%Y-%m-%d_%H-%M-%S')}.json"


def raw_first_timestamp(raw_file: "str | Path") -> float:
    """First trajectory's ``first_timestamp``, without decoding the whole file.

    Scans a prefix for the field; falls back to widening the prefix, and only
    then to a streaming decode of the first element.
    """
    raw_file = Path(raw_file)
    if not raw_file.is_file():
        raise FileNotFoundError(f"No such file: {raw_file}")

    for prefix_bytes in _PREFIX_BYTES:
        with raw_file.open("r", encoding="utf-8", errors="ignore") as handle:
            chunk = handle.read(prefix_bytes)
        match = _FIRST_TIMESTAMP.search(chunk)
        if match:
            return float(match.group(1))

    from .rawio import iter_trajectories

    for record in iter_trajectories(raw_file):
        return float(record["first_timestamp"])
    raise ValueError(f"No first_timestamp field found in {raw_file}")


def raw_dir(data_dir: "str | Path", day: int) -> Path:
    return Path(data_dir) / "i24motion" / f"2022-11-{day}"


def manifest(data_dir: "str | Path", results_dir: "str | Path", day: int,
             force: bool = False, log=None) -> List[Segment]:
    """Raw segment -> processed name for one day, cached under results/.mvt.

    The cache is reused while the raw listing (names, sizes, modification
    times) is unchanged, so the prefix scans happen once.
    """
    log = log or (lambda *_: None)
    listing = _raw_listing(data_dir, day)
    if not listing:
        raise FileNotFoundError(
            f"No raw MOTION files matching *_{DAY_ABBREV[day]}_0_*.json in "
            f"{raw_dir(data_dir, day)}")

    cache = Path(results_dir) / ".mvt" / "manifests" / f"segments_2022-11-{day}.json"
    fingerprint = [[path.name, path.stat().st_size, int(path.stat().st_mtime)]
                   for path in listing]
    if not force and cache.is_file():
        try:
            cached = json.loads(cache.read_text(encoding="utf-8"))
            if cached.get("fingerprint") == fingerprint:
                return [Segment(seq=entry["seq"], raw_name=entry["raw_name"],
                                raw_path=Path(entry["raw_path"]),
                                first_timestamp=entry["first_timestamp"],
                                output_name=entry["output_name"])
                        for entry in cached["segments"]]
        except (ValueError, KeyError):
            pass          # unreadable cache is not fatal; rebuild it

    segments = []
    for path in listing:
        timestamp = raw_first_timestamp(path)
        segments.append(Segment(seq=_sequence(path.name), raw_name=path.name,
                                raw_path=path, first_timestamp=timestamp,
                                output_name=segment_name(timestamp)))
        log(f"    manifest: {path.name} -> {segments[-1].output_name}")

    try:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(json.dumps(
            {"fingerprint": fingerprint,
             "segments": [{"seq": s.seq, "raw_name": s.raw_name,
                           "raw_path": str(s.raw_path),
                           "first_timestamp": s.first_timestamp,
                           "output_name": s.output_name} for s in segments]},
            indent=1), encoding="utf-8")
    except OSError:
        pass              # caching is an optimization, never a requirement
    return segments


def _raw_listing(data_dir: "str | Path", day: int) -> List[Path]:
    directory = raw_dir(data_dir, day)
    if not directory.is_dir():
        raise FileNotFoundError(f"Raw MOTION folder does not exist: {directory}")
    # Skip macOS AppleDouble siblings (._<uuid>__fri_0_00.json): they match the
    # glob but are not data. The MATLAB stages filter these too.
    return sorted(path for path in directory.glob(f"*_{DAY_ABBREV[day]}_0_*.json")
                  if path.is_file() and not path.name.startswith("."))


def _sequence(raw_name: str) -> int:
    """The 0..23 index encoded at the end of a raw segment name."""
    match = re.search(r"_(\d+)\.json$", raw_name)
    return int(match.group(1)) if match else -1
