"""A make for the Python pipeline: dependency graph, staleness, and -j.

Port of the Makefile plus ``Scripts/+mvt/isStale.m``, in pure Python and with
no external build tool. It answers the same two questions make does - *what is
out of date* and *what can run at once* - over the ported stages.

The unit of work
----------------
Every buildable thing is a :class:`Unit`: a stage, a day, and for ``slim`` a
single 10-minute segment. Segment-level units are what make ``-j`` behave the
way you would want without any shard arithmetic: with 3 days there are 3 gps
units and 72 slim units, so ``-j10`` fills all ten workers across days *and*
segments by itself, and ``-j3`` does not sit idle behind one long day. The
scheduler simply runs any ready unit when a slot frees.

``slim`` still accepts an explicit ``--shard k/N`` for parity with the MATLAB
``Shard [k N]`` option and for driving it by hand.

Staleness
---------
The classic makefile rule, matching ``mvt.isStale``: a target is stale when an
output is missing, when a data input is newer, or when the *code* that produces
it is newer. Comparisons carry a one-second tolerance because filesystems
disagree about sub-second mtimes. A unit is also stale when anything it depends
on is going to rebuild, so a single changed module propagates forward.

Each unit runs in its own subprocess. That mirrors the MATLAB pipeline's
process-level parallelism (there is no Parallel Computing Toolbox anywhere in
it), and it means a stage's peak memory is returned to the OS when it finishes
rather than accumulating across a run.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from .workspace import Workspace

__all__ = ["Unit", "STAGES", "plan", "execute", "is_stale", "DEFAULT_TARGETS"]

#: Seconds of slack when comparing mtimes (APFS/network volumes disagree on
#: sub-second resolution).
TOLERANCE = 1.0

#: Stage order, and which stages a default `build` runs.
DEFAULT_TARGETS = ("gps", "slim", "samples", "fields", "figures", "micro")

#: Code each stage depends on, relative to the mvtpy package. Editing one of
#: these makes that stage's outputs stale, exactly as the Makefile's
#: `$(SRC_stage)` prerequisites do.
STAGE_SOURCES = {
    "gps": ("gpsassemble.py", "gpsruns.py", "gpsmatch.py", "matjson.py",
            "matround.py", "rawio.py", "workspace.py"),
    "slim": ("slim.py", "lanes.py", "kinematics.py", "fuel.py", "avdist.py",
             "matjson.py", "matround.py", "rawio.py", "segments.py", "workspace.py"),
    "samples": ("samples.py", "rawio.py", "workspace.py"),
    "fields": ("fields.py", "rawio.py", "workspace.py"),
    "figures": ("plotting.py", "avanalysis.py", "fields.py", "rawio.py",
                "workspace.py"),
    "micro": ("microplot.py", "rawio.py", "workspace.py"),
}

STAGES = tuple(STAGE_SOURCES)


@dataclass
class Unit:
    """One schedulable piece of work."""

    stage: str
    day: int
    key: str
    outputs: List[Path]
    inputs: List[Path]
    argv: List[str]
    deps: List[str] = field(default_factory=list)
    label: str = ""

    def sources(self) -> List[Path]:
        package = Path(__file__).resolve().parent
        return [package / name for name in STAGE_SOURCES[self.stage]]


# --- staleness ------------------------------------------------------------


def is_stale(outputs: Sequence[Path], inputs: Sequence[Path],
             sources: Sequence[Path], force: bool = False):
    """Make's freshness test. Returns ``(stale, reason)``.

    Compares the *oldest* output against the *newest* prerequisite, so a
    partially rebuilt target is correctly reported as stale. Missing
    prerequisites are ignored here - the stage itself raises the more
    informative error when a required input is absent.
    """
    if force:
        return True, "forced rebuild"
    if not outputs:
        return True, "no outputs declared"

    missing = [path for path in outputs if not Path(path).exists()]
    if missing:
        return True, f"missing output {_short(missing[0])}"

    oldest_output = min((Path(p).stat().st_mtime, p) for p in outputs)
    prerequisites = [Path(p) for p in list(inputs) + list(sources)]
    existing = [(p.stat().st_mtime, p) for p in prerequisites if p.exists()]
    if not existing:
        return False, "up to date"

    newest_input = max(existing)
    if newest_input[0] > oldest_output[0] + TOLERANCE:
        kind = "source" if Path(newest_input[1]) in set(map(Path, sources)) else "input"
        return True, (f"{kind} {_short(newest_input[1])} is newer than "
                      f"output {_short(oldest_output[1])}")
    return False, "up to date"


def _short(path) -> str:
    path = Path(path)
    return "/".join(path.parts[-3:]) if len(path.parts) > 3 else str(path)


# --- the graph ------------------------------------------------------------


def plan(ws: Workspace, days: Sequence[int], targets: Sequence[str] = DEFAULT_TARGETS,
         force: bool = False, log=None) -> List[Unit]:
    """Build the unit graph for the requested days and stages.

    Only ``slim`` fans out per segment; every other stage is one unit per day.
    Units are returned in dependency order.
    """
    log = log or (lambda *_: None)
    targets = [stage for stage in DEFAULT_TARGETS if stage in set(targets)]
    units: List[Unit] = []
    base = _base_argv(ws)

    for day in days:
        slim_dir = ws.slim_dir(day)
        figures = ws.figures_dir(day)
        gps_file = ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{day}.json"
        gps_key = f"gps/{day}"

        if "gps" in targets:
            units.append(Unit(
                stage="gps", day=day, key=gps_key,
                outputs=[gps_file],
                inputs=_listing(ws.cars_dir()) + _listing(_raw_dir(ws, day)),
                argv=base + ["gps", "--day", str(day)],
                label=f"gps 2022-11-{day}"))

        slim_keys: List[str] = []
        if "slim" in targets:
            from . import segments as segment_map

            for segment in segment_map.manifest(ws.data_dir, ws.results_dir, day, log=log):
                key = f"slim/{day}/{segment.seq:02d}"
                slim_keys.append(key)
                units.append(Unit(
                    stage="slim", day=day, key=key,
                    outputs=[segment.output_path(slim_dir)],
                    inputs=[segment.raw_path, gps_file],
                    argv=base + ["slim", "--day", str(day),
                                 "--segment", str(segment.seq)],
                    deps=[gps_key] if "gps" in targets else [],
                    label=f"slim 2022-11-{day} #{segment.seq:02d}"))

        # Stages that consume the whole day's slim tree. Their data inputs are
        # resolved at plan time; on a cold tree that list is empty, which is
        # fine - the slim units they depend on will have produced it by then.
        day_inputs = _listing(slim_dir, "I-24MOTION_*.json")
        for stage, outputs in (
            ("samples", [figures / f"samples_for_distance_analysis_{day}.npz"]),
            ("fields", [figures / f"fields_motion_2022-11-{day}.npz"]),
        ):
            if stage in targets:
                units.append(Unit(
                    stage=stage, day=day, key=f"{stage}/{day}",
                    outputs=outputs, inputs=day_inputs,
                    argv=base + [stage, "--day", str(day)],
                    deps=list(slim_keys), label=f"{stage} 2022-11-{day}"))

        if "figures" in targets:
            units.append(Unit(
                stage="figures", day=day, key=f"figures/{day}",
                outputs=[figures / f"fig_field_2022-11-{day}_{name}.png"
                         for name in ("Rho", "Q", "F", "U", "Phi", "Psi")],
                inputs=day_inputs + [gps_file],
                argv=base + ["figures", "--day", str(day)],
                deps=list(slim_keys) + ([f"samples/{day}"] if "samples" in targets else []),
                label=f"figures 2022-11-{day}"))

        if "micro" in targets:
            stem = (f"fig_motion_trajectories_"
                    f"{datetime(2022, 11, day).strftime('%Y%m%d')}_west_laneall_py")
            units.append(Unit(
                stage="micro", day=day, key=f"micro/{day}",
                outputs=[figures / f"{stem}_{suffix}.png"
                         for suffix in ("lowres", "zoomwin", "zoom")],
                inputs=day_inputs,
                argv=base + ["micro", "--day", str(day)],
                deps=list(slim_keys), label=f"micro 2022-11-{day}"))

    _mark_stale(units, force)
    return units


def _base_argv(ws: Workspace) -> List[str]:
    return [sys.executable, "-m", "mvtpy",
            "--data-dir", str(ws.data_dir), "--results-dir", str(ws.results_dir)]


def _raw_dir(ws: Workspace, day: int) -> Path:
    return ws.data_dir / "i24motion" / f"2022-11-{day}"


def _listing(directory: Path, pattern: str = "*") -> List[Path]:
    if not directory.is_dir():
        return []
    return sorted(path for path in directory.glob(pattern)
                  if path.is_file() and not path.name.startswith("."))


def _mark_stale(units: Sequence[Unit], force: bool) -> None:
    """Annotate each unit with `stale`/`reason`, propagating along dependencies."""
    rebuilding = set()
    for unit in units:
        stale, reason = is_stale(unit.outputs, unit.inputs, unit.sources(), force)
        if not stale:
            upstream = [key for key in unit.deps if key in rebuilding]
            if upstream:
                stale, reason = True, f"dependency {upstream[0]} is rebuilding"
        unit.stale = stale          # type: ignore[attr-defined]
        unit.reason = reason        # type: ignore[attr-defined]
        if stale:
            rebuilding.add(unit.key)


# --- the scheduler --------------------------------------------------------


def execute(units: Sequence[Unit], jobs: int = 1, dry_run: bool = False,
            log: Optional[Callable[[str], None]] = None,
            keep_going: bool = False, progress=None) -> int:
    """Run every stale unit, up to `jobs` at once, respecting dependencies.

    Returns the number of failed units. Like make, a failure stops new work
    from being scheduled unless `keep_going`; units already running finish.

    `progress` is a :class:`mvtpy.progress.Progress`. By default one is created
    that draws a live bar on a terminal and degrades to one durable line per
    event when output is piped to a file - which is how these runs are usually
    captured.
    """
    log = log or print
    pending = {unit.key: unit for unit in units if getattr(unit, "stale", True)}
    fresh = [unit for unit in units if not getattr(unit, "stale", True)]

    if dry_run:
        for unit in fresh:
            log(f"[mvt] up to date: {unit.label}")
        for unit in units:
            if unit.key in pending:
                log(f"[mvt] would build {unit.label}: {getattr(unit, 'reason', '')}")
        log(f"[mvt] {len(pending)} unit(s) to build, {len(fresh)} up to date")
        return 0

    if progress is None:
        from .progress import Progress
        progress = Progress(total=len(pending), jobs=jobs, fresh=len(fresh))
    for unit in fresh:
        progress.info(f"[mvt] up to date: {unit.label}")

    done = set(unit.key for unit in fresh)
    failed: set = set()
    running: Dict[str, tuple] = {}
    started = time.time()
    total = len(pending)
    completed = 0

    def ready_unit() -> Optional[Unit]:
        for unit in pending.values():
            if any(dep in failed for dep in unit.deps):
                continue
            # A dependency still queued or in flight blocks this unit; one that
            # is not in the plan at all (not requested this run) does not.
            if any(dep in pending or dep in running for dep in unit.deps):
                continue
            return unit
        return None

    while pending or running:
        while len(running) < max(1, jobs) and not (failed and not keep_going):
            unit = ready_unit()
            if unit is None:
                break
            del pending[unit.key]
            progress.started(unit.key, unit.label, getattr(unit, "reason", ""))
            # Each child writes to its own file rather than a pipe: with many
            # concurrent jobs a full pipe buffer would block a child that
            # nobody is reading yet.
            sink = tempfile.TemporaryFile(mode="w+")
            running[unit.key] = (unit, subprocess.Popen(
                unit.argv, stdout=sink, stderr=subprocess.STDOUT, text=True),
                time.time(), sink)

        if not running:
            break

        # Reap whichever children have finished, then wait a beat.
        finished = [key for key, (_, process, _, _) in running.items()
                    if process.poll() is not None]
        if not finished:
            progress.tick()             # refresh elapsed times while we wait
            time.sleep(0.2)
            continue

        for key in finished:
            unit, process, start, sink = running.pop(key)
            completed += 1
            detail = ""
            if process.returncode == 0:
                done.add(key)
            else:
                failed.add(key)
                sink.seek(0)
                detail = "\n".join(sink.read().rstrip().splitlines()[-15:])
            sink.close()
            progress.finished(key, unit.label, process.returncode == 0,
                              time.time() - start, detail)

    # Anything still queued was blocked by a failure (or by a failed dependency).
    skipped = len(pending)
    if skipped:
        progress.info(f"[mvt] {skipped} unit(s) not attempted after failure")
    from .progress import format_duration
    progress.close(f"[mvt] {completed - len(failed)} built, {len(failed)} failed, "
                   f"{len(fresh)} up to date in "
                   f"{format_duration(time.time() - started)}")
    return len(failed) + skipped
