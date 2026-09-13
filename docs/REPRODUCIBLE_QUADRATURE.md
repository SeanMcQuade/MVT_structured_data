# Making the fuel quadrature reproducible across machines

*Data set version 2.1.1 — see [DATA_CHANGELOG.md](DATA_CHANGELOG.md) for the release note, and this document for the measurements.*

## Why this change

The trajectory fuel totals were computed with `integrate = @(t,v) dot(t(2:end)-t(1:end-1), (v(1:end-1)+v(2:end))/2)`.
MATLAB's `dot` dispatches to whichever BLAS library ships with the platform —
Accelerate on Apple silicon, MKL on Windows — and those libraries accumulate the
sum in different orders. Floating-point addition is not associative, so a
different order gives a different last bit. This is not a bug in either library;
it is ordinary library freedom. But it means **the same script, on the same
input data, produces different `total_fuel_consumed_grams` on a Mac and on a
PC**, which is what we observed independently in two places: comparing MATLAB
runs across machines, and comparing MATLAB against the Python port. Measured
against the exactly-rounded reference on 60 real trajectories, MATLAB's `dot`
deviates on **19 of 60**, by 1–2 units in the last place (~1.6 × 10⁻¹⁶
relative). The released JSON rounds fuel totals to 4 decimals, so a value that
happens to sit within one bit of a rounding boundary lands on different sides on
different machines, and the published file differs byte-for-byte.

The fix replaces `dot` with compensated (Kahan–Babuška–Neumaier) summation in
`mvt.neumaierDot`, selected by `flag_deterministic_quadrature` (default `true`)
in `generate_data_mvt_slim.m` and `generate_data_mvt_full.m`. Compensated
summation is a fixed sequence of IEEE-754 double operations with no library-,
thread-, or platform-dependent freedom, so every machine produces identical
bits — and so does any other language implementing the same loop, which is how
the Python port (`mvtpy.kinematics.neumaier_dot`) now agrees exactly. It is also
**more accurate**, not merely more consistent: it matched the exactly-rounded
sum on 60 of 60 test trajectories, where `dot` matched 41 of 60. The cost is
about 1 s per 10-minute segment.

## Scale of the difference

| Measure | Value |
| --- | --- |
| Relative change in the fuel integral | ≤ 3.7 × 10⁻¹⁶ (a few parts per quadrillion) |
| Raw integrals that change at all | 19 of 60 sampled trajectories |
| **Published 4-decimal values that change** | **≈ 2 in 20,466 trajectories (0.01%)** |
| Size of such a change | 0.0001 g — one unit in the last recorded decimal (0.1 mg) |
| MATLAB vs Python with the flag on | 60 of 60 bit-identical |

The numerical content of the data is unchanged. A 10⁻¹⁶ relative shift is nine
or ten orders of magnitude below any physical uncertainty in the fuel models,
the GPS positions, or the MOTION trajectories, and the affected values move by
0.1 mg on one trajectory in ten thousand, with no systematic sign. No aggregate
in the paper — fleet fuel consumption, energy savings, per-class averages — is
affected at any reported precision. What changes is that the pipeline now
returns the same answer on every machine, which is the property we need for the
released code and data to be independently reproducible.

## What this means operationally

Because `dot` was platform-dependent, any previously generated `slim`/`full`
outputs are specific to the machine that produced them. Regenerate them with the
flag on, and re-derive the checksum manifests
(`python -m mvtpy verify --update`) from the regenerated tree, so the recorded
md5s describe a platform-independent artifact. Setting
`flag_deterministic_quadrature = false` (and `DETERMINISTIC_QUADRATURE = False`
on the Python side) restores the old expression for comparison — though note
that "the old behaviour" is not a single well-defined target, since it differs
between the two BLAS libraries.

Affected fields: `total_fuel_consumed_grams` (slim and full), and in `full` also
`total_fuel_consumed_flat_road_grams`, `total_reference_fuel_consumed_grams`, and
`total_reference_fuel_consumed_flat_road_grams`. Fields derived from them —
`total_fuel_consumed_gallons`, `total_fuel_economy_mpg` — inherit the change.
