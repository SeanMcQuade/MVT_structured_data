"""Locate the data and results trees, mirroring the MATLAB ``mvt.paths``.

The Python port follows the same layout convention as the pipeline: the code
lives in ``MVT_structured_data/python`` and, by default, ``data`` and
``results`` are siblings of ``MVT_structured_data``. Any of those can be
overridden by argument or environment variable, which is what lets the Docker
image write into a mounted host directory.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

__all__ = ["Workspace"]


@dataclass(frozen=True)
class Workspace:
    """Resolved data/results locations for one run.

    Resolution order for each of data/results (first hit wins):
      1. an explicit path passed in,
      2. the ``MVT_DATA_DIR`` / ``MVT_RESULTS_DIR`` environment variable,
      3. the default sibling of the repository (``<workspace>/data`` etc.).
    """

    data_dir: Path
    results_dir: Path

    @classmethod
    def resolve(cls, data_dir: Optional[str] = None,
                results_dir: Optional[str] = None) -> "Workspace":
        # python/mvtpy/workspace.py -> mvtpy -> python -> repo -> workspace root
        repo_root = Path(__file__).resolve().parents[2]
        workspace_root = repo_root.parent

        data = _pick(data_dir, "MVT_DATA_DIR", workspace_root / "data")
        results = _pick(results_dir, "MVT_RESULTS_DIR", workspace_root / "results")
        return cls(data_dir=data, results_dir=results)

    # Convenience locations, matching the MATLAB layout ------------------------
    def cars_dir(self) -> Path:
        return self.data_dir / "cars"

    def motion_dir(self, day: int) -> Path:
        return self.data_dir / "i24motion" / f"2022-11-{day}"

    def gps_dir(self) -> Path:
        return self.results_dir / "gps"

    def slim_dir(self, day: int) -> Path:
        return self.results_dir / "slim" / f"2022-11-{day}"

    def figures_dir(self, day: int) -> Path:
        return self.results_dir / "figures" / f"2022-11-{day}"


def _pick(explicit: Optional[str], env_var: str, default: Path) -> Path:
    if explicit:
        return Path(explicit).expanduser().resolve()
    from_env = os.environ.get(env_var)
    if from_env:
        return Path(from_env).expanduser().resolve()
    return default.resolve()
