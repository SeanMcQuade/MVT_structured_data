"""Write the ``dataset_info.json`` sidecar that says what a product folder holds.

Python counterpart of ``Scripts/+mvt/writeDatasetInfo.m``, emitting the same
fields so a folder is self-describing regardless of which implementation built
it. Two copies of ``slim`` at 2.1 and 2.1.1 differ only in the last decimal of a
few fuel totals and are otherwise indistinguishable, so without this a recipient
cannot tell them apart.

The file deliberately mixes reproducible facts (version, product, copyright)
with run provenance (host, user, interpreter), so it is **not** byte-identical
across machines and is excluded from the checksum manifests - see
:mod:`mvtpy.verify`, whose ``outputs_for`` never lists it.
"""

from __future__ import annotations

import getpass
import json
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

__all__ = ["write", "build_info", "FILENAME"]

FILENAME = "dataset_info.json"

_SCHEME = ("MAJOR.MINOR.PATCH; PATCH keeps every field, meaning and format "
           "identical and changes only the last decimal of a small number of "
           "values. See docs/DATA_CHANGELOG.md.")


def build_info(data_dir: Path, product: str, day: Optional[int] = None) -> dict:
    """Assemble the sidecar contents, describing the data in ``data_dir``."""
    # Imported here, not at module scope: __init__ defines DATA_VERSION after
    # it imports this module, so a top-level import would be circular.
    from . import DATA_VERSION

    data_dir = Path(data_dir)
    files = [p for p in data_dir.iterdir()
             if p.is_file() and not p.name.startswith(".") and p.name != FILENAME]
    newest = max((p.stat().st_mtime for p in files), default=None)

    info = {
        "dataset": "CIRCLES MegaVanderTest (MVT) derived data",
        "product": product,
        "data_version": DATA_VERSION,
        "version_scheme": _SCHEME,
        "upstream_motion_data": ("I-24 MOTION, versioned independently by the "
                                 "observatory and unmodified by this pipeline"),
        "copyright": "(C) 2026 CIRCLES Consortium",
        "license": "BSD-3-Clause",
        "documentation": "docs/DATA_CHANGELOG.md, docs/REPRODUCIBLE_QUADRATURE.md",
    }
    if day is not None:
        info["day"] = f"2022-11-{day}"
    info["files"] = len(files)
    # From the newest output, not the clock: re-running on an unchanged product
    # must not invent a new timestamp.
    if newest is not None:
        info["generated_utc"] = (datetime.fromtimestamp(newest, timezone.utc)
                                 .strftime("%Y-%m-%dT%H:%M:%SZ"))
    info["generated_by"] = _provenance()
    return info


def write(data_dir, product: str, day: Optional[int] = None,
          info_dir=None) -> Optional[Path]:
    """Describe the data in ``data_dir``, writing the sidecar into ``info_dir``.

    ``info_dir`` defaults to the parent of ``data_dir``: the sidecar must not
    sit beside the data, because any consumer globbing ``*.json`` in a folder of
    trajectory JSON will swallow it - which is exactly how the MATLAB micro
    stage broke. Pass ``info_dir=data_dir`` for a product whose data is not in a
    per-day subfolder.

    Writes nothing for a folder with no outputs, rather than claiming one was
    produced. Failure to write is never fatal to a pipeline run.
    """
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return None
    info_dir = Path(info_dir) if info_dir is not None else data_dir.parent
    info = build_info(data_dir, product, day)
    if not info.get("files"):
        return None
    info_dir.mkdir(parents=True, exist_ok=True)
    target = info_dir / FILENAME
    try:
        target.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    except OSError:
        return None
    return target


def _provenance() -> dict:
    """Who and what produced this. Best effort; a missing field is never fatal."""
    info = {
        "implementation": "Python (mvtpy)",
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "host": "",
        "user": "",
        "code_commit": "",
    }
    try:
        info["host"] = socket.gethostname()
    except OSError:
        pass
    try:
        info["user"] = getpass.getuser()
    except Exception:                                    # noqa: BLE001
        pass
    try:
        repo = Path(__file__).resolve().parents[2]
        sha = subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        if sha.returncode == 0:
            info["code_commit"] = sha.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return info
