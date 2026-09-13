"""Tests for the streaming reader used on multi-gigabyte raw MOTION segments."""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import matjson  # noqa: E402
from mvtpy.rawio import count_trajectories, iter_trajectories  # noqa: E402


@pytest.fixture()
def array_file(tmp_path):
    document = [
        {"_id": {"$oid": "a1"}, "timestamp": [1.5, 2.5], "note": "braces {} and \"quotes\""},
        {"_id": {"$oid": "b2"}, "timestamp": [], "nested": {"deep": [1, [2, 3]]}},
        {"_id": {"$oid": "c3"}, "timestamp": [-0.0], "escape": "back\\slash"},
    ]
    path = tmp_path / "segment.json"
    path.write_text(matjson.dumps(document), encoding="utf-8")
    return path, document


def test_streams_every_element(array_file):
    path, document = array_file
    assert list(iter_trajectories(path)) == document


def test_structural_characters_inside_strings_do_not_confuse_the_scanner(array_file):
    path, document = array_file
    streamed = list(iter_trajectories(path))
    assert streamed[0]["note"] == document[0]["note"]
    assert streamed[2]["escape"] == document[2]["escape"]


def test_small_chunk_sizes_give_the_same_result(array_file):
    path, document = array_file
    assert list(iter_trajectories(path, chunk_bytes=7)) == document


def test_counting_does_not_decode(array_file):
    path, document = array_file
    assert count_trajectories(path) == len(document)


def test_unterminated_array_is_rejected(tmp_path):
    path = tmp_path / "truncated.json"
    path.write_text('[{"a":1},{"b":', encoding="utf-8")
    with pytest.raises(ValueError):
        list(iter_trajectories(path))


def test_iteration_is_lazy_on_a_real_segment(raw_segment_path):
    """Taking a few records must not read the whole multi-GB file."""
    first = list(itertools.islice(iter_trajectories(raw_segment_path), 3))
    assert len(first) == 3
    assert all("timestamp" in record and "x_position" in record for record in first)
