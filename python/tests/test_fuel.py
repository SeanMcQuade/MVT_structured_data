"""Tests for the ported CIRCLES fuel models.

Three layers:

* the constants in ``mvtpy.fuel`` are re-parsed from ``Models/*.m`` and must
  match, so the port cannot silently drift from the MATLAB source;
* the model structure (families, class mapping) matches the MATLAB switch;
* the fuel rates and totals reproduce the released data exactly when evaluated
  on the raw inputs MATLAB used.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mvtpy import fuel, kinematics as kin  # noqa: E402
from mvtpy.matround import round_decimals  # noqa: E402

MODELS_DIR = Path(__file__).resolve().parents[2] / "Models"
CONSTANT_PATTERN = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*([-\d.eE+]+)\s*;", re.M)


def _matlab_constants(model_name: str) -> dict[str, float]:
    source = (MODELS_DIR / f"fuel_model_{model_name}_simplified.m").read_text()
    # Stop before the shared unit-conversion block so gs2kW is excluded.
    header = source.split("% Other constants")[0]
    return {name: float(value) for name, value in CONSTANT_PATTERN.findall(header)}


@pytest.mark.parametrize("model_name", sorted(fuel.FUEL_MODELS))
def test_constants_match_matlab_source(model_name):
    if not MODELS_DIR.is_dir():
        pytest.skip("Models/ not available")

    expected = _matlab_constants(model_name)
    model = fuel.FUEL_MODELS[model_name]

    mismatches = {
        name: (value, getattr(model, name))
        for name, value in expected.items()
        if getattr(model, name, None) != value
    }
    assert not mismatches, f"{model_name}: {mismatches}"


def test_conversion_factor_matches_matlab():
    source = (MODELS_DIR / "fuel_model_midBase_simplified.m").read_text()
    value = float(re.search(r"gs2kW\s*=\s*([\d.]+)", source).group(1))
    assert fuel.GRAMS_PER_SECOND_TO_KW == value


def test_model_families():
    """Light-duty models cut fuel above vc; heavy-duty ones have a linear floor."""
    assert isinstance(fuel.FUEL_MODELS["midBase"], fuel.CarFuelModel)
    assert isinstance(fuel.FUEL_MODELS["Class8Tractor"], fuel.TruckFuelModel)


def test_coarse_class_mapping_matches_matlab_switch():
    assert fuel.model_for_coarse_class(0).name == "midBase"
    assert fuel.model_for_coarse_class(1).name == "midSUV"
    assert fuel.model_for_coarse_class(2).name == "Pickup"
    assert fuel.model_for_coarse_class(3).name == "Pickup"
    assert fuel.model_for_coarse_class(4).name == "Class8Tractor"
    assert fuel.model_for_coarse_class(5).name == "Pickup"
    with pytest.raises(KeyError):
        fuel.model_for_coarse_class(6)   # motorcycle: absent from the data


def test_idle_and_fuel_cut_branches():
    model = fuel.FUEL_MODELS["midBase"]

    idle_rate, _, _ = model(np.array([0.0]), np.array([0.0]), np.array([0.0]))
    assert idle_rate[0] == model.fc_idle

    # Hard braking well above the cut-off speed burns no fuel.
    cut_rate, _, _ = model(np.array([25.0]), np.array([-4.0]), np.array([0.0]))
    assert cut_rate[0] == 0.0

    # The infeasibility flag marks negative speed (2) and impossible accel (1).
    _, _, flag = model(np.array([-1.0, 10.0]), np.array([0.0, 50.0]), np.array([0.0, 0.0]))
    assert list(flag) == [2.0, 1.0]


def test_power_is_rate_times_conversion():
    model = fuel.FUEL_MODELS["Pickup"]
    rate, power, _ = model(np.array([12.0]), np.array([0.3]), np.array([0.01]))
    assert power[0] == rate[0] * fuel.GRAMS_PER_SECOND_TO_KW


def test_fuel_rate_matches_released_data(segment_pairs, grade_map):
    """Rates computed from raw samples must equal MATLAB's, digit for digit."""
    for pair in segment_pairs:
        segment = pair["segment"]
        time = pair["time"] - pair["time"][0]
        distance = kin.FT_TO_METER * np.abs(pair["x_feet"] - pair["x_feet"][0])
        x_meters = kin.FT_TO_METER * (pair["x_feet"] - kin.ORIGIN_X_FEET)

        speed = kin.speed(distance, time)
        acceleration = kin.acceleration(distance, time)
        grade = grade_map(x_meters, pair["direction"])

        model = fuel.FUEL_MODELS[segment["energy_model"]]
        rate, _, _ = model(speed, acceleration, grade, project=True)

        assert np.array_equal(round_decimals(rate, 4),
                              np.asarray(segment["fuel_rate_grams_per_second"]))


def test_energy_model_choice_matches_released_data(segment_pairs):
    for pair in segment_pairs:
        segment = pair["segment"]
        expected = fuel.model_for_coarse_class(segment["coarse_vehicle_class"]).name
        assert segment["energy_model"] == expected


def test_totals_match_released_data(segment_pairs, grade_map):
    """Total fuel, gallons, and fuel economy reproduce the released values."""
    for pair in segment_pairs:
        segment = pair["segment"]
        time = pair["time"] - pair["time"][0]
        distance = kin.FT_TO_METER * np.abs(pair["x_feet"] - pair["x_feet"][0])
        x_meters = kin.FT_TO_METER * (pair["x_feet"] - kin.ORIGIN_X_FEET)

        rate, _, _ = fuel.FUEL_MODELS[segment["energy_model"]](
            kin.speed(distance, time),
            kin.acceleration(distance, time),
            grade_map(x_meters, pair["direction"]),
            project=True,
        )

        total_grams = kin.trapezoid_integral(time, rate)
        assert round_decimals(total_grams, 4) == segment["total_fuel_consumed_grams"]

        gallons = kin.GRAM_TO_GALLON * total_grams
        assert round_decimals(gallons, 4) == segment["total_fuel_consumed_gallons"]

        travelled = abs(distance[-1] - distance[0])
        assert round_decimals(travelled, 4) == segment["total_distance_traversed_meters"]

        if segment["total_fuel_economy_mpg"] is not None and gallons != 0:
            mpg = (travelled * kin.METER_TO_MILE) / gallons
            assert round_decimals(mpg, 4) == segment["total_fuel_economy_mpg"]
