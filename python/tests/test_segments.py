"""Tests for the raw-segment -> processed-name mapping.

The mapping is load-bearing for the build system: it decides output filenames
*before* anything is decoded, so a wrong name means make rebuilds forever.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import segments  # noqa: E402


#: First trajectory of the first 2022-11-18 raw segment (..._fri_0_00.json).
#: The released file it produces is I-24MOTION_2022-11-18_05-59-59.json.
FIRST_RAW_TIMESTAMP = 1668772799.9833353


def test_segment_name_matches_the_released_convention():
    assert segments.segment_name(FIRST_RAW_TIMESTAMP) == \
        "I-24MOTION_2022-11-18_05-59-59.json"


def test_segment_name_truncates_rather_than_rounds_seconds():
    """This very timestamp is .983 past the second: rounding would name the
    released file 06-00-00 and every rebuild would disagree with the data."""
    assert segments.segment_name(1668772799.9833353).endswith("05-59-59.json")
    assert segments.segment_name(1668772800.0).endswith("06-00-00.json")


def test_segment_name_is_nashville_local_time():
    """Central, not UTC: this instant is 11:59:59 UTC and 05:59:59 local."""
    assert "05-59-59" in segments.segment_name(FIRST_RAW_TIMESTAMP)


def make_raw(directory: Path, name: str, first_timestamp: float) -> Path:
    path = directory / name
    path.write_text(json.dumps([
        {"first_timestamp": first_timestamp, "timestamp": [first_timestamp]}]),
        encoding="utf-8")
    return path


def test_raw_first_timestamp_reads_from_the_prefix(tmp_path):
    path = make_raw(tmp_path, "x__fri_0_00.json", FIRST_RAW_TIMESTAMP)
    assert segments.raw_first_timestamp(path) == pytest.approx(FIRST_RAW_TIMESTAMP)


def test_raw_first_timestamp_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        segments.raw_first_timestamp(tmp_path / "nope.json")


def build_day(tmp_path, count=3):
    raw = tmp_path / "data" / "i24motion" / "2022-11-18"
    raw.mkdir(parents=True)
    for index in range(count):
        make_raw(raw, f"uuid__fri_0_{index:02d}.json", FIRST_RAW_TIMESTAMP + 600 * index)
    return tmp_path / "data", tmp_path / "results"


def test_manifest_maps_every_segment_and_caches(tmp_path):
    data, results = build_day(tmp_path)
    first = segments.manifest(data, results, 18)
    assert [s.seq for s in first] == [0, 1, 2]
    assert first[0].output_name == "I-24MOTION_2022-11-18_05-59-59.json"
    assert first[1].output_name == "I-24MOTION_2022-11-18_06-09-59.json"

    cache = results / ".mvt" / "manifests" / "segments_2022-11-18.json"
    assert cache.is_file()

    # A second call reuses the cache and returns the same mapping.
    assert [s.output_name for s in segments.manifest(data, results, 18)] == \
           [s.output_name for s in first]


def test_manifest_rebuilds_when_the_raw_listing_changes(tmp_path):
    data, results = build_day(tmp_path, count=2)
    segments.manifest(data, results, 18)
    make_raw(data / "i24motion" / "2022-11-18", "uuid__fri_0_02.json",
             FIRST_RAW_TIMESTAMP + 1800)
    assert len(segments.manifest(data, results, 18)) == 3


def test_manifest_survives_a_corrupt_cache(tmp_path):
    data, results = build_day(tmp_path)
    segments.manifest(data, results, 18)
    cache = results / ".mvt" / "manifests" / "segments_2022-11-18.json"
    cache.write_text("{not json", encoding="utf-8")
    assert len(segments.manifest(data, results, 18)) == 3


def test_manifest_skips_appledouble_siblings(tmp_path):
    data, results = build_day(tmp_path)
    (data / "i24motion" / "2022-11-18" / "._uuid__fri_0_00.json").write_text(
        "not data", encoding="utf-8")
    assert len(segments.manifest(data, results, 18)) == 3


def test_manifest_errors_clearly_on_a_missing_or_empty_day(tmp_path):
    data, results = build_day(tmp_path)
    with pytest.raises(FileNotFoundError, match="does not exist"):
        segments.manifest(data, results, 16)

    empty = tmp_path / "data" / "i24motion" / "2022-11-17"
    empty.mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="No raw MOTION files"):
        segments.manifest(data, results, 17)
