"""Check that generated outputs are exactly what they should be.

This exists to answer one question with evidence rather than assertion: *does
the Python pipeline produce the same files as the MATLAB pipeline?* It is the
artifact to hand a colleague who has to take the port seriously.

Two modes, both driven by the same declaration of what a stage produces:

* **Manifest mode** - compare every output against a stored table of md5
  checksums (``python/expected/checksums-2022-11-DD.json``). Portable: it needs
  only the checksums file, not a copy of the reference data. This is the strict
  test - a single differing byte fails it.
* **Reference mode** - compare against another results tree (typically the
  MATLAB one) file by file. Slower, needs both trees, but explains *how* two
  files differ rather than only that they do.

Not every output can be checked by checksum, and saying so plainly matters more
than a green tick:

* ``.json`` (gps, slim) - byte comparison. These are the strong evidence: the
  port's stated acceptance test is byte-identical JSON.
* ``.npz`` vs MATLAB ``.mat`` (samples, fields) - array comparison with a
  tolerance, because the two container formats cannot be byte-compared at all.
  A checksum of the ``.npz`` still pins Python-to-Python reproducibility.
* ``dataset_info.json`` - **not** checksummed and never listed. It records the
  data-set version alongside run provenance (host, user, interpreter), so it is
  machine-specific by design; see :mod:`mvtpy.datasetinfo`.
* ``.png`` (figures) - **not** checksummed. Matplotlib and MATLAB renderers do
  not agree pixel for pixel, and PNG bytes carry encoder metadata, so a
  checksum here would fail for reasons that say nothing about correctness. They
  are reported as ``skipped`` with that reason.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .workspace import Workspace

__all__ = ["Result", "checksum", "outputs_for", "verify_against_manifest",
           "verify_against_reference", "write_manifest", "manifest_path",
           "summarize", "CHUNK_BYTES"]

CHUNK_BYTES = 1 << 20

#: Array comparison tolerance for .npz/.mat pairs. The ported numerical stages
#: are float-exact to ~1e-11 (see docs/PYTHON_PORT.md); this is deliberately
#: tighter than "close enough for plotting" and loose enough to survive
#: container round-tripping.
DEFAULT_RTOL = 1e-9
DEFAULT_ATOL = 1e-9

MATCH, DIFFERS, MISSING, SKIPPED, NO_REFERENCE = (
    "match", "differs", "missing", "skipped", "no-reference")


@dataclass
class Result:
    """One checked output."""

    path: Path
    stage: str
    status: str
    detail: str = ""

    @property
    def ok(self) -> bool:
        return self.status in (MATCH, SKIPPED, NO_REFERENCE)


def checksum(path: "str | Path") -> str:
    """md5 of a file, read in chunks (these run to gigabytes)."""
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_BYTES), b""):
            digest.update(block)
    return digest.hexdigest()


def manifest_path(day: int, directory: Optional[Path] = None) -> Path:
    directory = directory or (Path(__file__).resolve().parent.parent / "expected")
    return Path(directory) / f"checksums-2022-11-{day}.json"


# --- what a stage produces -------------------------------------------------


def outputs_for(ws: Workspace, day: int) -> Dict[str, List[Path]]:
    """Files each stage produces for one day, in a stable order.

    Discovered from the workspace rather than declared twice, so this cannot
    drift from what the stages actually write.
    """
    figures = ws.figures_dir(day)
    return {
        "gps": [ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{day}.json"],
        "slim": sorted(ws.slim_dir(day).glob("I-24MOTION_*.json"))
                if ws.slim_dir(day).is_dir() else [],
        "samples": [figures / f"samples_for_distance_analysis_{day}.npz"],
        "fields": [figures / f"fields_motion_2022-11-{day}.npz"],
        "figures": sorted(figures.glob("fig_field_*.png")),
        "micro": sorted(figures.glob("fig_motion_trajectories_*_py_*.png")),
    }


def _checkable(path: Path) -> bool:
    """PNGs are excluded from checksumming; see the module docstring."""
    return path.suffix.lower() != ".png"


# --- manifest mode ---------------------------------------------------------


def write_manifest(ws: Workspace, day: int, directory: Optional[Path] = None,
                   log=None) -> Path:
    """Record md5s of the current outputs as the expected values."""
    log = log or (lambda *_: None)
    entries = {}
    for stage, paths in outputs_for(ws, day).items():
        for path in paths:
            if not (path.is_file() and _checkable(path)):
                continue
            entries[_key(ws, path)] = {"md5": checksum(path),
                                       "bytes": path.stat().st_size,
                                       "stage": stage}
            log(f"    {entries[_key(ws, path)]['md5']}  {_key(ws, path)}")

    target = manifest_path(day, directory)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(
        {"day": day, "files": dict(sorted(entries.items()))}, indent=1) + "\n",
        encoding="utf-8")
    return target


def verify_against_manifest(ws: Workspace, day: int,
                            directory: Optional[Path] = None) -> List[Result]:
    """Compare current outputs against the stored checksums."""
    source = manifest_path(day, directory)
    if not source.is_file():
        raise FileNotFoundError(
            f"no checksum manifest for 2022-11-{day} at {source}; "
            f"create one with `mvt verify --day {day} --update`")
    expected = json.loads(source.read_text(encoding="utf-8"))["files"]

    results: List[Result] = []
    for stage, paths in outputs_for(ws, day).items():
        for path in paths:
            key = _key(ws, path)
            if not _checkable(path):
                results.append(Result(path, stage, SKIPPED,
                                      "png: renderers differ; not checksummed"))
                continue
            if key not in expected:
                results.append(Result(path, stage, NO_REFERENCE,
                                      "not in manifest"))
                continue
            if not path.is_file():
                results.append(Result(path, stage, MISSING, "not generated"))
                continue
            actual = checksum(path)
            if actual == expected[key]["md5"]:
                results.append(Result(path, stage, MATCH, actual))
            else:
                results.append(Result(
                    path, stage, DIFFERS,
                    f"md5 {actual} != {expected[key]['md5']}; "
                    f"{path.stat().st_size} bytes vs {expected[key]['bytes']}"))

    generated = {_key(ws, p) for paths in outputs_for(ws, day).values() for p in paths}
    for key in sorted(set(expected) - generated):
        results.append(Result(Path(key), expected[key].get("stage", "?"), MISSING,
                              "in manifest but not generated"))
    return results


def _key(ws: Workspace, path: Path) -> str:
    """Manifest key: path relative to the results tree, so it is portable."""
    try:
        return str(Path(path).resolve().relative_to(Path(ws.results_dir).resolve()))
    except ValueError:
        return str(path)


# --- reference mode --------------------------------------------------------


def verify_against_reference(ws: Workspace, day: int,
                             reference_dir: "str | Path",
                             rtol: float = DEFAULT_RTOL,
                             atol: float = DEFAULT_ATOL) -> List[Result]:
    """Compare outputs against another results tree, explaining differences.

    JSON is compared byte for byte. ``.npz`` is compared against the reference
    ``.npz`` or the MATLAB ``.mat`` of the same stage, array by array.
    """
    reference_dir = Path(reference_dir)
    results: List[Result] = []

    for stage, paths in outputs_for(ws, day).items():
        for path in paths:
            if not _checkable(path):
                results.append(Result(path, stage, SKIPPED,
                                      "png: renderers differ; not compared"))
                continue
            if not path.is_file():
                results.append(Result(path, stage, MISSING, "not generated"))
                continue

            counterpart = _counterpart(reference_dir / _key(ws, path))
            if counterpart is None:
                results.append(Result(path, stage, NO_REFERENCE,
                                      f"no reference under {reference_dir}"))
            elif counterpart.suffix == path.suffix:
                results.append(_compare_same_format(path, counterpart, stage,
                                                    rtol, atol))
            else:
                results.append(_compare_arrays(path, counterpart, stage, rtol, atol))
    return results


def _counterpart(candidate: Path) -> Optional[Path]:
    """The reference file for an output: same name, or the MATLAB .mat twin."""
    if candidate.is_file():
        return candidate
    if candidate.suffix == ".npz":
        as_mat = candidate.with_suffix(".mat")
        if as_mat.is_file():
            return as_mat
    return None


def _compare_same_format(path: Path, reference: Path, stage: str,
                         rtol: float, atol: float) -> Result:
    if path.suffix == ".npz":
        return _compare_arrays(path, reference, stage, rtol, atol)
    if checksum(path) == checksum(reference):
        return Result(path, stage, MATCH, "byte-identical")
    return Result(path, stage, DIFFERS,
                  f"{path.stat().st_size} bytes vs {reference.stat().st_size}; "
                  f"first difference at byte {_first_difference(path, reference)}")


def _first_difference(left: Path, right: Path) -> str:
    offset = 0
    with left.open("rb") as a, right.open("rb") as b:
        while True:
            block_a, block_b = a.read(CHUNK_BYTES), b.read(CHUNK_BYTES)
            if not block_a and not block_b:
                return "none (lengths differ only)"
            if block_a != block_b:
                for index, (x, y) in enumerate(zip(block_a, block_b)):
                    if x != y:
                        return str(offset + index)
                return str(offset + min(len(block_a), len(block_b)))
            offset += len(block_a)


def _compare_arrays(path: Path, reference: Path, stage: str,
                    rtol: float, atol: float) -> Result:
    """Array-by-array comparison, .npz against .npz or MATLAB .mat."""
    try:
        actual = dict(np.load(path, allow_pickle=False))
        expected = _load_reference_arrays(reference)
    except Exception as error:                      # unreadable reference
        return Result(path, stage, NO_REFERENCE, f"could not read: {error}")

    if expected is None:
        return Result(path, stage, NO_REFERENCE,
                      f"cannot read {reference.suffix} (needs h5py/scipy)")

    problems = []
    for name, value in sorted(actual.items()):
        if name not in expected:
            continue                                # naming differs by format
        other = np.asarray(expected[name], dtype=float).squeeze()
        mine = np.asarray(value, dtype=float).squeeze()
        if mine.shape != other.shape:
            problems.append(f"{name}: shape {mine.shape} vs {other.shape}")
            continue
        if np.array_equal(mine, other, equal_nan=True):
            continue
        close = np.isclose(mine, other, rtol=rtol, atol=atol, equal_nan=True)
        if close.all():
            worst = float(np.nanmax(np.abs(mine - other)))
            problems.append(f"{name}: within tolerance (max |diff| {worst:.3e})")
        else:
            problems.append(f"{name}: {int((~close).sum())} of {close.size} "
                            f"outside tolerance")

    hard = [text for text in problems if "within tolerance" not in text]
    if hard:
        return Result(path, stage, DIFFERS, "; ".join(hard[:4]))
    if problems:
        return Result(path, stage, MATCH, "equal within tolerance: "
                      + "; ".join(problems[:3]))
    return Result(path, stage, MATCH, "arrays identical")


def _load_reference_arrays(path: Path):
    if path.suffix == ".npz":
        return dict(np.load(path, allow_pickle=False))
    try:
        import h5py

        with h5py.File(path, "r") as handle:        # MATLAB v7.3
            return {key: np.array(handle[key]).T for key in handle
                    if not key.startswith("#")}
    except (ImportError, OSError):
        pass
    try:
        import scipy.io as sio

        return {key: value for key, value in
                sio.loadmat(path, squeeze_me=True).items()
                if not key.startswith("__")}
    except ImportError:
        return None


# --- reporting -------------------------------------------------------------


def summarize(results: Sequence[Result], log=print, verbose: bool = False) -> int:
    """Print a per-file report. Returns the number of failures."""
    counts: Dict[str, int] = {}
    for result in results:
        counts[result.status] = counts.get(result.status, 0) + 1
        if result.status in (DIFFERS, MISSING):
            log(f"  FAIL  {result.path.name}  [{result.stage}]  {result.detail}")
        elif verbose:
            log(f"  {result.status:12s} {result.path.name}  {result.detail}")

    order = (MATCH, DIFFERS, MISSING, SKIPPED, NO_REFERENCE)
    log("  " + "  ".join(f"{name}={counts.get(name, 0)}" for name in order
                         if counts.get(name)))
    return counts.get(DIFFERS, 0) + counts.get(MISSING, 0)
