"""Tests for the pure-Python make: staleness, the graph, and the -j scheduler.

These use synthetic units whose `argv` is a trivial Python one-liner, so the
scheduler is exercised without running the real pipeline.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import build  # noqa: E402
from mvtpy.workspace import Workspace  # noqa: E402


def touch(path: Path, age: float = 0.0) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x", encoding="utf-8")
    if age:
        stamp = time.time() - age
        os.utime(path, (stamp, stamp))
    return path


def unit(key, outputs, argv, deps=(), stage="micro", day=18) -> build.Unit:
    return build.Unit(stage=stage, day=day, key=key, outputs=list(outputs),
                      inputs=[], argv=list(argv), deps=list(deps), label=key)


def writer(path: Path, seconds: float = 0.0, code: int = 0):
    """argv for a child that creates `path`, optionally slowly, then exits."""
    return [sys.executable, "-c",
            f"import time,pathlib,sys; time.sleep({seconds}); "
            f"pathlib.Path({str(path)!r}).write_text('done'); sys.exit({code})"]


# --- staleness ------------------------------------------------------------


def test_missing_output_is_stale(tmp_path):
    stale, reason = build.is_stale([tmp_path / "out"], [], [])
    assert stale and "missing output" in reason


def test_no_declared_outputs_is_stale():
    assert build.is_stale([], [], [])[0]


def test_up_to_date_when_output_is_newest(tmp_path):
    source = touch(tmp_path / "code.py", age=100)
    data = touch(tmp_path / "in.json", age=100)
    out = touch(tmp_path / "out.png")
    stale, reason = build.is_stale([out], [data], [source])
    assert not stale and reason == "up to date"


def test_newer_input_is_stale(tmp_path):
    out = touch(tmp_path / "out.png", age=100)
    data = touch(tmp_path / "in.json")
    stale, reason = build.is_stale([out], [data], [])
    assert stale and "input" in reason


def test_newer_source_is_stale_and_says_so(tmp_path):
    """Editing a stage must rebuild its outputs - the bug the MATLAB
    skip-if-exists guards had."""
    out = touch(tmp_path / "out.png", age=100)
    source = touch(tmp_path / "code.py")
    stale, reason = build.is_stale([out], [], [source])
    assert stale and reason.startswith("source ")


def test_oldest_output_decides_so_partial_rebuilds_are_stale(tmp_path):
    data = touch(tmp_path / "in.json", age=50)
    fresh = touch(tmp_path / "a.png")
    old = touch(tmp_path / "b.png", age=100)
    assert build.is_stale([fresh, old], [data], [])[0]


def test_force_overrides_everything(tmp_path):
    out = touch(tmp_path / "out.png")
    assert build.is_stale([out], [], [], force=True)[0]


def test_missing_prerequisites_are_ignored(tmp_path):
    """The stage itself raises the better error for an absent input."""
    out = touch(tmp_path / "out.png")
    assert not build.is_stale([out], [tmp_path / "gone.json"], [])[0]


# --- the graph ------------------------------------------------------------


def test_plan_fans_slim_out_per_segment_and_wires_dependencies(tmp_path, monkeypatch):
    from mvtpy import segments as segment_map

    fake = [segment_map.Segment(seq=i, raw_name=f"r{i}.json",
                                raw_path=tmp_path / f"r{i}.json",
                                first_timestamp=0.0,
                                output_name=f"I-24MOTION_2022-11-18_0{i}.json")
            for i in range(3)]
    monkeypatch.setattr(segment_map, "manifest", lambda *a, **k: fake)

    ws = Workspace(data_dir=tmp_path / "data", results_dir=tmp_path / "results")
    units = build.plan(ws, [18])
    by_key = {u.key: u for u in units}

    assert sum(1 for u in units if u.stage == "slim") == 3
    assert by_key["slim/18/00"].deps == ["gps/18"]
    # Day-level stages wait for every segment of their day.
    for stage in ("samples", "fields", "micro"):
        assert set(by_key[f"{stage}/18"].deps) == {"slim/18/00", "slim/18/01", "slim/18/02"}
    assert "samples/18" in by_key["figures/18"].deps


def test_plan_respects_target_selection(tmp_path, monkeypatch):
    from mvtpy import segments as segment_map
    monkeypatch.setattr(segment_map, "manifest", lambda *a, **k: [])

    ws = Workspace(data_dir=tmp_path / "data", results_dir=tmp_path / "results")
    units = build.plan(ws, [18], targets=["fields"])
    assert [u.stage for u in units] == ["fields"]


def test_staleness_propagates_to_dependents(tmp_path):
    """A rebuilt dependency makes its dependents stale even if their own
    outputs look current."""
    upstream = unit("up", [tmp_path / "missing.json"], ["true"])
    downstream = unit("down", [touch(tmp_path / "done.png")], ["true"], deps=["up"])
    build._mark_stale([upstream, downstream], force=False)
    assert upstream.stale and downstream.stale
    assert "dependency up" in downstream.reason


# --- the scheduler --------------------------------------------------------


def test_dry_run_builds_nothing(tmp_path, capsys):
    out = tmp_path / "out.txt"
    units = [unit("a", [out], writer(out))]
    build._mark_stale(units, force=False)
    assert build.execute(units, dry_run=True) == 0
    assert not out.exists()
    assert "would build" in capsys.readouterr().out


def test_executes_and_reports_up_to_date_units(tmp_path, capsys):
    made = tmp_path / "made.txt"
    already = touch(tmp_path / "already.txt")
    units = [unit("a", [made], writer(made)), unit("b", [already], writer(already))]
    build._mark_stale(units, force=False)

    assert build.execute(units) == 0
    assert made.read_text() == "done"
    assert "up to date: b" in capsys.readouterr().out


def test_dependencies_run_in_order(tmp_path):
    first, second = tmp_path / "1.txt", tmp_path / "2.txt"
    # The second unit exits non-zero unless `first` finished before it started,
    # so a scheduler that ignored deps would fail this rather than pass it.
    guard = (f"import pathlib, sys\n"
             f"if not pathlib.Path({str(first)!r}).exists():\n"
             f"    sys.exit(3)\n"
             f"pathlib.Path({str(second)!r}).write_text('done')\n")
    units = [
        unit("first", [first], writer(first, seconds=0.5)),
        unit("second", [second], [sys.executable, "-c", guard], deps=["first"]),
    ]
    build._mark_stale(units, force=False)
    assert build.execute(units, jobs=4) == 0
    assert second.read_text() == "done"


def test_jobs_actually_run_concurrently(tmp_path):
    """-j N must overlap independent units, not serialize them."""
    outs = [tmp_path / f"{i}.txt" for i in range(4)]
    units = [unit(f"u{i}", [out], writer(out, seconds=0.6)) for i, out in enumerate(outs)]
    build._mark_stale(units, force=False)

    started = time.time()
    assert build.execute(units, jobs=4) == 0
    elapsed = time.time() - started
    assert all(out.exists() for out in outs)
    assert elapsed < 4 * 0.6, f"4 units of 0.6s took {elapsed:.1f}s; not parallel"


def test_failure_is_reported_and_blocks_dependents(tmp_path, capsys):
    bad, downstream = tmp_path / "bad.txt", tmp_path / "down.txt"
    units = [
        unit("bad", [bad], [sys.executable, "-c", "import sys; sys.exit(2)"]),
        unit("down", [downstream], writer(downstream), deps=["bad"]),
    ]
    build._mark_stale(units, force=False)

    assert build.execute(units) > 0
    assert not downstream.exists()
    output = capsys.readouterr().out
    assert "FAILED bad" in output
    assert "not attempted" in output


def test_keep_going_still_builds_unrelated_units(tmp_path):
    bad, other = tmp_path / "bad.txt", tmp_path / "other.txt"
    units = [
        unit("bad", [bad], [sys.executable, "-c", "import sys; sys.exit(2)"]),
        unit("other", [other], writer(other)),
    ]
    build._mark_stale(units, force=False)

    assert build.execute(units, jobs=2, keep_going=True) == 1
    assert other.exists()
