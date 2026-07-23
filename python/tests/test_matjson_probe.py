"""Cross-check the encoder against a real MATLAB ``jsonencode`` probe.

``python/matlab_probes/probe_jsonencode.m`` records how a given MATLAB release
formats a battery of probe values. Those probes cover the cases the released
data cannot pin down: magnitudes between 1e5 and 1e9, values below 1e-4,
infinities, and structural details.

The whole module skips when no probe file is present, so it costs nothing until
someone runs the probe in MATLAB and commits the result.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import matjson  # noqa: E402

PROBE_DIR = Path(__file__).resolve().parents[1] / "matlab_probes"

#: Probe name -> the Python value corresponding to the MATLAB expression.
#: Mirrors the probe table in probe_jsonencode.m.
SCALAR_PROBES = {
    "e5_integral": 1e5,
    "e5_fraction": 123456.789,
    "e6_integral": 1e6,
    "e6_fraction": 1234567.8912,
    "e7_fraction": 12345678.9123,
    "e8_fraction": 123456789.1234,
    "e9_boundary_below": 999999999.9999,
    "e9_boundary_at": 1e9,
    "e9_fraction": 1668599999.9836,
    "e15": 1e15,
    "e16": 1e16,
    "em4": 0.0001,
    "em5": 0.00001,
    "em6": 0.000001,
    "em7_fraction": 1.2345e-07,
    "digits15": 1668599999.9836,
    "digits17": 1668691841.800001,
    "digits_short": 1.5,
    "zero": 0.0,
    "negative_zero": -0.0,
    "negative_fraction": -3.25,
    "nan": math.nan,
    "inf": math.inf,
    "neg_inf": -math.inf,
    "integral_small": 840.0,
    "integral_negative": -1.0,
    "int32_value": 7,
    "logical_true": True,
    "logical_false": False,
}


def _probe_files():
    return sorted(PROBE_DIR.glob("jsonencode_probe_*.json"))


pytestmark = pytest.mark.skipif(
    not _probe_files(),
    reason="no MATLAB jsonencode probe recorded yet (run matlab_probes/probe_jsonencode.m)",
)


@pytest.mark.parametrize("probe_path", _probe_files())
def test_scalar_probes_match(probe_path):
    probe = json.loads(probe_path.read_text())
    mismatches = []
    for name, value in SCALAR_PROBES.items():
        if name not in probe:
            continue
        ours = matjson.dumps(value)
        if ours != probe[name]:
            mismatches.append((name, probe[name], ours))
    assert not mismatches, (
        f"{probe_path.name}: MATLAB vs mvtpy differ for "
        + ", ".join(f"{n}: {m!r} != {o!r}" for n, m, o in mismatches)
    )


@pytest.mark.parametrize("probe_path", _probe_files())
def test_structural_probes_match(probe_path):
    probe = json.loads(probe_path.read_text())

    if "struct_field_order" in probe:
        assert matjson.dumps({"b": 1.0, "a": 2.0}) == probe["struct_field_order"]
    if "empty_numeric" in probe:
        assert matjson.dumps({"x": []}) == probe["empty_numeric"]
    if "two_element_array" in probe:
        assert matjson.dumps({"x": [5.0, 6.0]}) == probe["two_element_array"]
    if "nested_struct" in probe:
        assert matjson.dumps({"a": {"b": 1.0}}) == probe["nested_struct"]
    if "struct_array" in probe:
        assert matjson.dumps([{"a": 1.0}, {"a": 2.0}]) == probe["struct_array"]
