"""CIRCLES simplified fuel-consumption models (port of ``Models/fuel_model_*.m``).

Each model gives instantaneous fuel rate in grams/second from instantaneous
velocity, acceleration, and road grade, as fitted by the CIRCLES energy team
(model version 3.1, 2023-03-10). ``Scripts/generate_data_mvt_slim.m`` selects a
model per trajectory from the MOTION coarse vehicle class and calls it with
``project = true``.

Two functional families appear in ``Models/``:

* **car family** (Compact, Pickup, midBase, midSUV): the fitted rate is capped
  below by ``beta0`` while ``v <= vc``, and fuel is cut entirely above ``vc``
  when the deceleration is stronger than the fitted brake threshold.
* **truck family** (Class4PND, Class8Tractor): the fitted rate is capped below
  by the linear floor ``h0 + h1 v``, with no fuel-cut branch.

The parameter values below are transcribed from the MATLAB sources;
``tests/test_fuel.py`` re-parses those ``.m`` files and asserts the numbers here
still match, so the two cannot drift apart.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np

__all__ = [
    "CarFuelModel",
    "TruckFuelModel",
    "FUEL_MODELS",
    "GRAMS_PER_SECOND_TO_KW",
    "model_for_coarse_class",
    "fuel_rate",
]

#: grams/sec to kW conversion factor used by every model
GRAMS_PER_SECOND_TO_KW = 42.36


@dataclass(frozen=True)
class _BaseFuelModel:
    """Shared fitted terms: feasibility boundary, cruise/accel/grade polynomial."""

    name: str
    fc_idle: float
    b1: float
    b2: float
    b3: float
    b4: float
    b5: float
    b6: float
    C0: float
    C1: float
    C2: float
    C3: float
    p0: float
    p1: float
    p2: float
    q0: float
    q1: float
    z0: float
    z1: float
    z2: float

    def max_feasible_acceleration(self, v: np.ndarray, g: np.ndarray) -> np.ndarray:
        """Upper boundary of the fitted feasibility region (MATLAB ``ma``)."""
        return (np.minimum(self.b1, self.b2 / np.maximum(v, 1e-12) - self.b3 * v**2)
                - np.minimum(self.b4, self.b5 + self.b6 * v) * g)

    def _fitted_rate(self, v: np.ndarray, a: np.ndarray, g: np.ndarray) -> np.ndarray:
        # Acceleration at which the quadratic term q(v) is minimized.
        aplus = np.maximum(
            a,
            -(self.p0 + self.p1 * v + self.p2 * v**2)
            / (2 * (self.q0 + self.q1 * np.maximum(v, 1e-12))),
        )
        return (
            self.C0 + self.C1 * v + self.C2 * v**2 + self.C3 * v**3      # cruising
            + self.p0 * a + self.p1 * a * v + self.p2 * a * v**2        # linear accel
            + self.q0 * aplus**2 + self.q1 * aplus**2 * v               # quadratic accel
            + self.z0 * g + self.z1 * g * v + self.z2 * g * v**2        # road grade
        )

    def __call__(self, v, a=None, g=None, project: bool = True):
        v = np.asarray(v, dtype=float)
        a = np.zeros_like(v) if a is None else np.asarray(a, dtype=float)
        g = np.zeros_like(v) if g is None else np.asarray(g, dtype=float)
        if not (v.shape == a.shape == g.shape):
            raise ValueError("Inputs must be of identical size.")

        ma = self.max_feasible_acceleration(v, g)
        # flag: 2 where v < 0, 1 where v >= 0 but the request is infeasible
        flag = 2.0 * (v < 0) + (v >= 0) * (a > ma)

        if project:
            a = np.minimum(a, ma)
        v = np.maximum(v, 0.0)

        fc = self._apply_floor(self._fitted_rate(v, a, g), v, a, g)

        # Idle rate for a stationary vehicle.
        fc = np.where((v < 0.1) & (np.abs(a) < 0.01), self.fc_idle, fc)

        power = fc * GRAMS_PER_SECOND_TO_KW
        return fc, power, flag

    def _apply_floor(self, fc, v, a, g):  # pragma: no cover - overridden
        raise NotImplementedError


@dataclass(frozen=True)
class CarFuelModel(_BaseFuelModel):
    """Light-duty models: minimum rate below ``vc``, fuel cut on strong braking."""

    vc: float = 0.0
    beta0: float = 0.0
    a0: float = 0.0
    a1: float = 0.0
    a2: float = 0.0
    a3: float = 0.0
    a4: float = 0.0

    def _apply_floor(self, fc, v, a, g):
        fc = np.maximum(fc, (v <= self.vc) * self.beta0)
        a_brake = (self.a0 + self.a1 * v + self.a2 * g
                   + self.a3 * v**2 + self.a4 * v * g)
        return np.where((v > self.vc) & (a <= a_brake), 0.0, fc)


@dataclass(frozen=True)
class TruckFuelModel(_BaseFuelModel):
    """Heavy-duty models: linear lower bound ``h0 + h1 v``, no fuel cut."""

    h0: float = 0.0
    h1: float = 0.0

    def _apply_floor(self, fc, v, a, g):
        return np.maximum(fc, self.h0 + self.h1 * v)


FUEL_MODELS: Dict[str, _BaseFuelModel] = {
    "Compact": CarFuelModel(
        name="Compact", fc_idle=0.0972, b1=3.3605, b2=41.6037, b3=0.00021189,
        b4=8.9362, b5=3.9757, b6=0.24476, C0=0.15918, C1=0.013463, C2=0,
        C3=3.1889e-05, p0=0.047828, p1=0.086975, p2=6.825e-08, q0=0.0025557,
        q1=0.019099, z0=0.13285, z1=0.77984, z2=0.0019733,
        vc=5.04, beta0=0.0972, a0=-0.26981, a1=-0.0023996, a2=-9.0623,
        a3=-0.00029215, a4=-0.011899),
    "midBase": CarFuelModel(
        name="midBase", fc_idle=0.1271, b1=3.9218, b2=48.9891, b3=0.00013964,
        b4=8.9036, b5=6.1892, b6=0.10048, C0=0.19829, C1=0.021122, C2=0,
        C3=2.7801e-05, p0=0.23956, p1=0.0080592, p2=0.0027737, q0=0,
        q1=0.050556, z0=2.5227, z1=0.76464, z2=0.0060213,
        vc=5.07, beta0=0.1271, a0=-0.15742, a1=-0.00037876, a2=-9.1124,
        a3=-0.00022957, a4=-0.011559),
    "midSUV": CarFuelModel(
        name="midSUV", fc_idle=0.1637, b1=3.3377, b2=53.4583, b3=0.00023901,
        b4=9.1847, b5=8.1403, b6=0.034303, C0=0.22498, C1=0.021292, C2=0,
        C3=3.7654e-05, p0=0.17419, p1=0.094617, p2=0.00071347, q0=0,
        q1=0.02884, z0=2.3211, z1=0.74453, z2=0.013073,
        vc=9.16, beta0=0.1637, a0=-0.26854, a1=-0.0015267, a2=-9.4305,
        a3=-0.00032843, a4=-0.0053817),
    "Pickup": CarFuelModel(
        name="Pickup", fc_idle=0.1999, b1=3.0163, b2=60.3816, b3=0.00038328,
        b4=9.1334, b5=8.9129, b6=0.0063029, C0=0.26318, C1=0.023432, C2=0,
        C3=5.5207e-05, p0=0.23805, p1=0.10287, p2=0.0012594, q0=0,
        q1=0.030277, z0=3.766, z1=0.69241, z2=0.023785,
        vc=11, beta0=0.1999, a0=-0.26463, a1=-0.0013822, a2=-9.4942,
        a3=-0.00039814, a4=-0.004408),
    "Class4PND": TruckFuelModel(
        name="Class4PND", fc_idle=0.1923, b1=1.2639, b2=14.9616, b3=0.00049533,
        b4=9.5693, b5=7.5025, b6=0.099425, C0=0.24292, C1=0.03827, C2=0,
        C3=0.00018703, p0=0.65006, p1=0.33381, p2=0.0025519, q0=0.36738,
        q1=0.042944, z0=2.0694, z1=3.7719, z2=0,
        h0=0, h1=0.031719),
    "Class8Tractor": TruckFuelModel(
        name="Class8Tractor", fc_idle=0.2384, b1=2.4229, b2=8.4463,
        b3=0.00026269, b4=9.7395, b5=8.6171, b6=0.15762, C0=0.59448,
        C1=0.082609, C2=0, C3=0.00027278, p0=0.20476, p1=1.1962, p2=0.019119,
        q0=0, q1=0.14424, z0=0.88147, z1=11.1899, z2=0.18836,
        h0=0.49109, h1=0.025152),
}

#: MOTION coarse vehicle class -> energy model, exactly as the switch statement
#: in generate_data_mvt_slim.m assigns it. Class 6 (motorcycle) is absent from
#: the data and unmapped there as well.
COARSE_CLASS_TO_MODEL = {
    0: "midBase",        # sedan
    1: "midSUV",         # midsize
    2: "Pickup",         # van
    3: "Pickup",         # pickup
    4: "Class8Tractor",  # semi
    5: "Pickup",         # truck
}


def model_for_coarse_class(coarse_class: int) -> _BaseFuelModel:
    """Return the fuel model the pipeline uses for a MOTION coarse class."""
    try:
        return FUEL_MODELS[COARSE_CLASS_TO_MODEL[int(coarse_class)]]
    except KeyError as error:
        raise KeyError(
            f"no energy model for coarse_vehicle_class {coarse_class}"
        ) from error


def fuel_rate(model_name: str, v, a, g, project: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate a named model; mirrors ``fuel_model_<name>_simplified(v,a,g,project)``."""
    return FUEL_MODELS[model_name](v, a, g, project=project)
