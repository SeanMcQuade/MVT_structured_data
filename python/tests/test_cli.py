"""Tests for workspace resolution and the CLI argument handling."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import cli  # noqa: E402
from mvtpy.workspace import Workspace  # noqa: E402


def test_workspace_defaults_to_sibling_layout():
    ws = Workspace.resolve()
    repo = Path(__file__).resolve().parents[2]
    assert ws.data_dir == (repo.parent / "data").resolve()
    assert ws.results_dir == (repo.parent / "results").resolve()


def test_workspace_explicit_overrides(tmp_path):
    ws = Workspace.resolve(str(tmp_path / "d"), str(tmp_path / "r"))
    assert ws.data_dir == (tmp_path / "d").resolve()
    assert ws.results_dir == (tmp_path / "r").resolve()


def test_workspace_env_overrides(monkeypatch, tmp_path):
    monkeypatch.setenv("MVT_DATA_DIR", str(tmp_path / "envdata"))
    monkeypatch.setenv("MVT_RESULTS_DIR", str(tmp_path / "envresults"))
    ws = Workspace.resolve()
    assert ws.data_dir == (tmp_path / "envdata").resolve()
    assert ws.results_dir == (tmp_path / "envresults").resolve()


def test_workspace_derived_locations():
    ws = Workspace.resolve("/d", "/r")
    assert ws.slim_dir(16) == Path("/r/slim/2022-11-16")
    assert ws.motion_dir(17) == Path("/d/i24motion/2022-11-17")
    assert ws.gps_dir() == Path("/r/gps")
    assert ws.figures_dir(18) == Path("/r/figures/2022-11-18")


def test_cli_requires_a_stage():
    with pytest.raises(SystemExit):
        cli.main([])


def test_cli_path_flags_work_before_and_after_stage(tmp_path):
    # Missing data should raise a clear SystemExit, not an argparse error,
    # proving the flags are parsed in both positions.
    for argv in (["--results-dir", str(tmp_path), "samples", "--day", "16"],
                 ["samples", "--day", "16", "--results-dir", str(tmp_path)]):
        with pytest.raises(SystemExit) as excinfo:
            cli.main(argv)
        assert "not found" in str(excinfo.value)


def test_cli_rejects_unknown_day():
    with pytest.raises(SystemExit):
        cli.main(["samples", "--day", "15"])
