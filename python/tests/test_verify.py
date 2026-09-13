"""Tests for the output verifier.

This is the piece that has to be trustworthy when someone asks "is the Python
output really the same?", so the tests care most about it *failing* when it
should: a changed byte, a missing file, a shape mismatch. A verifier that
reports success too easily is worse than none.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import verify  # noqa: E402
from mvtpy.workspace import Workspace  # noqa: E402

DAY = 18


def build_tree(root: Path, gps_text="gps-payload", slim_text="slim-payload",
               arrays=None) -> Workspace:
    results = root / "results"
    (results / "gps").mkdir(parents=True)
    (results / "gps" / f"CIRCLES_GPS_10Hz_2022-11-{DAY}.json").write_text(
        gps_text, encoding="utf-8")

    slim = results / "slim" / f"2022-11-{DAY}"
    slim.mkdir(parents=True)
    (slim / f"I-24MOTION_2022-11-{DAY}_05-59-59.json").write_text(
        slim_text, encoding="utf-8")

    figures = results / "figures" / f"2022-11-{DAY}"
    figures.mkdir(parents=True)
    arrays = arrays if arrays is not None else {"t": np.arange(5.0)}
    np.savez(figures / f"fields_motion_2022-11-{DAY}.npz", **arrays)
    np.savez(figures / f"samples_for_distance_analysis_{DAY}.npz", **arrays)
    (figures / f"fig_field_2022-11-{DAY}_Rho.png").write_bytes(b"\x89PNG-not-really")
    return Workspace(data_dir=root / "data", results_dir=results)


def test_checksum_matches_hashlib(tmp_path):
    import hashlib

    path = tmp_path / "f.bin"
    path.write_bytes(b"x" * 3_000_000)          # spans several read chunks
    assert verify.checksum(path) == hashlib.md5(b"x" * 3_000_000).hexdigest()


# --- manifest mode --------------------------------------------------------


def test_manifest_round_trip_verifies(tmp_path):
    ws = build_tree(tmp_path)
    verify.write_manifest(ws, DAY, tmp_path / "expected")
    results = verify.verify_against_manifest(ws, DAY, tmp_path / "expected")
    assert verify.summarize(results, log=lambda *_: None) == 0
    assert any(r.status == verify.MATCH for r in results)


def test_a_single_changed_byte_fails(tmp_path):
    ws = build_tree(tmp_path)
    verify.write_manifest(ws, DAY, tmp_path / "expected")
    target = ws.slim_dir(DAY) / f"I-24MOTION_2022-11-{DAY}_05-59-59.json"
    target.write_text("slim-payloae", encoding="utf-8")     # one byte differs

    results = verify.verify_against_manifest(ws, DAY, tmp_path / "expected")
    assert verify.summarize(results, log=lambda *_: None) == 1
    failed = [r for r in results if r.status == verify.DIFFERS]
    assert len(failed) == 1 and "md5" in failed[0].detail


def test_a_deleted_output_fails(tmp_path):
    ws = build_tree(tmp_path)
    verify.write_manifest(ws, DAY, tmp_path / "expected")
    (ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{DAY}.json").unlink()

    results = verify.verify_against_manifest(ws, DAY, tmp_path / "expected")
    assert verify.summarize(results, log=lambda *_: None) == 1
    assert any(r.status == verify.MISSING for r in results)


def test_pngs_are_skipped_not_silently_passed(tmp_path):
    """Renderers differ, so PNGs are excluded - but visibly, with a reason."""
    ws = build_tree(tmp_path)
    verify.write_manifest(ws, DAY, tmp_path / "expected")
    results = verify.verify_against_manifest(ws, DAY, tmp_path / "expected")

    skipped = [r for r in results if r.status == verify.SKIPPED]
    assert skipped and all(r.path.suffix == ".png" for r in skipped)
    assert "png" in skipped[0].detail
    # and no png is in the manifest at all
    manifest = json.loads((tmp_path / "expected" /
                           f"checksums-2022-11-{DAY}.json").read_text())
    assert not any(key.endswith(".png") for key in manifest["files"])


def test_missing_manifest_says_how_to_make_one(tmp_path):
    ws = build_tree(tmp_path)
    with pytest.raises(FileNotFoundError, match="--update"):
        verify.verify_against_manifest(ws, DAY, tmp_path / "expected")


def test_manifest_keys_are_relative_so_the_tree_can_move(tmp_path):
    ws = build_tree(tmp_path)
    verify.write_manifest(ws, DAY, tmp_path / "expected")
    manifest = json.loads((tmp_path / "expected" /
                           f"checksums-2022-11-{DAY}.json").read_text())
    assert all(not Path(key).is_absolute() for key in manifest["files"])
    assert f"gps/CIRCLES_GPS_10Hz_2022-11-{DAY}.json" in manifest["files"]


# --- reference mode -------------------------------------------------------


def test_reference_mode_reports_identical_trees(tmp_path):
    ours = build_tree(tmp_path / "a")
    theirs = build_tree(tmp_path / "b")
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)
    assert verify.summarize(results, log=lambda *_: None) == 0


def test_reference_mode_locates_the_first_differing_byte(tmp_path):
    ours = build_tree(tmp_path / "a", gps_text="gps-payloadXY")
    theirs = build_tree(tmp_path / "b", gps_text="gps-payloadAB")
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)

    failed = [r for r in results if r.status == verify.DIFFERS]
    assert failed and "first difference at byte 11" in failed[0].detail


def test_reference_mode_flags_array_shape_mismatch(tmp_path):
    ours = build_tree(tmp_path / "a", arrays={"t": np.arange(5.0)})
    theirs = build_tree(tmp_path / "b", arrays={"t": np.arange(6.0)})
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)

    failed = [r for r in results if r.status == verify.DIFFERS]
    assert any("shape" in r.detail for r in failed)


def test_arrays_within_tolerance_pass_but_say_so(tmp_path):
    ours = build_tree(tmp_path / "a", arrays={"t": np.arange(5.0) + 1e-12})
    theirs = build_tree(tmp_path / "b", arrays={"t": np.arange(5.0)})
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)

    npz = [r for r in results if r.path.suffix == ".npz"]
    assert npz and all(r.status == verify.MATCH for r in npz)
    assert any("within tolerance" in r.detail for r in npz)


def test_arrays_outside_tolerance_fail(tmp_path):
    ours = build_tree(tmp_path / "a", arrays={"t": np.arange(5.0) + 1e-3})
    theirs = build_tree(tmp_path / "b", arrays={"t": np.arange(5.0)})
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)
    assert any(r.status == verify.DIFFERS and "outside tolerance" in r.detail
               for r in results)


def test_absent_reference_is_reported_not_counted_as_success(tmp_path):
    ours = build_tree(tmp_path / "a")
    (tmp_path / "b").mkdir()
    results = verify.verify_against_reference(ours, DAY, tmp_path / "b")
    assert all(r.status in (verify.NO_REFERENCE, verify.SKIPPED) for r in results)
    # No reference is not a pass, but it is not a failure either - it must be
    # visible in the summary rather than folded into "match".
    assert not any(r.status == verify.MATCH for r in results)


def test_nan_arrays_compare_equal(tmp_path):
    """Field arrays carry NaN by design; NaN != NaN must not fail the check."""
    values = np.array([1.0, np.nan, 3.0])
    ours = build_tree(tmp_path / "a", arrays={"t": values})
    theirs = build_tree(tmp_path / "b", arrays={"t": values.copy()})
    results = verify.verify_against_reference(ours, DAY, theirs.results_dir)
    assert verify.summarize(results, log=lambda *_: None) == 0
