# MVT data dictionary

Every field of the data sets the pipeline produces: units, dtype, rounding, and
which stage writes it. This doubles as the specification for the Python port,
which must reproduce these files byte-for-byte.

Conventions used below:

- **Rounding.** Unless noted, every numeric field in the processed MOTION data
  (`slim`/`full`) is rounded to **4 decimals** immediately before encoding.
  The GPS data is not uniformly rounded; see that section.
- **null / empty.** `NaN` encodes as `null`; an empty array encodes as `[]`.
  These are meaningful: `null` marks a sample with no value; `[]` marks a field
  with no values at all (e.g. no control vehicle was ever nearby).
- **Coordinate frame.** Positions are in the I-24 MOTION v2 frame, origin at the
  Mill Creek Bridge (`x = 309804.0625 ft` in the raw data), x increasing
  eastbound. The processed data reports `x_position_meters` in meters from that
  origin; several internal quantities are measured in feet.
- Encoding details (number formatting, field order, whitespace) are in
  `MATLAB_JSON_FORMAT.md`. Field **order** is significant for byte parity and is
  the order given in each table.

---

## Processed MOTION data — `results/slim/` and `results/full/`

One JSON file per 10-minute segment, an array of trajectory records. `slim`
keeps westbound trajectories and the 31 fields the paper uses; `full` keeps both
directions and adds reference-controller and flat-road fields (45 fields).
Produced by `generate_data_mvt_slim.m` / `generate_data_mvt_full.m` (Python:
`mvtpy.slim`).

`length` in a table's dtype column means the value is an array with one entry
per timestamp sample.

### Fields common to slim and full

| # | Field | Unit | Type | Notes |
|---|---|---|---|---|
| 1 | `trajectory_id` | — | object `{x_oid}` | Raw MOTION `_id` with a `-N` suffix per clipped segment (`64888e87…-0`). |
| 2 | `timestamp` | s (POSIX) | array | Sample times, copied from the raw data. |
| 3 | `x_position_meters` | m | array | `0.3048·(x_raw − 309804.0625)`; longitudinal position from the origin. |
| 4 | `y_position_corrected_meters` | m | array | Lateral position after driving-line correction and clipping, in meters. |
| 5 | `coarse_vehicle_class` | code | scalar | 0 sedan, 1 midsize, 2 van, 3 pickup, 4 semi, 5 truck (6 motorcycle, absent). |
| — | `direction` (full only) | ±1 | scalar | −1 westbound, +1 eastbound. Dropped from slim (all −1). |
| 6 | `first_timestamp` | s | scalar | First sample time of the clipped segment. Names the output file. |
| 7 | `last_timestamp` | s | scalar | Last sample time. |
| 8 | `starting_x` | m | scalar | `x_position_meters(1)`. |
| 9 | `ending_x` | m | scalar | `x_position_meters(end)`. |
| 10 | `length` | ft | scalar | Vehicle length, from raw MOTION. Not converted. |
| 11 | `width` | ft | scalar | Vehicle width, from raw MOTION. |
| 12 | `height` | ft | scalar | Vehicle height, from raw MOTION. |
| 13 | `total_distance_traversed_meters` | m | scalar | `abs(x(end) − x(1))` in meters. |
| 14 | `speed_meters_per_second` | m/s | array | Central difference of distance vs time, ends duplicated. |
| 15 | `acceleration_meters_per_second_per_second` | m/s² | array | Second central difference, first/last repeated. |
| 16 | `road_grade_radians` | rad | array | Grade from `Models/Eastbound_grade_fit.csv`; negated for westbound. |
| 17 | `lane_number` | lane | scalar | Assigned lane: 1 (HOV/leftmost) … 4 (rightmost), 0 off highway, 5 ramp. |
| 18 | `energy_model` | name | string | Fuel model used: `midBase`, `midSUV`, `Pickup`, `Class8Tractor` (by class). |
| 19 | `fuel_rate_grams_per_second` | g/s | array | Instantaneous fuel rate from the class fuel model (road grade applied). |
| 20 | `percent_infeasibility` | % | scalar | Fraction of samples the fuel model flagged as dynamically infeasible. |
| 21 | `total_fuel_consumed_grams` | g | scalar | Trapezoidal integral of the fuel rate over time, by compensated summation so the value does not depend on the platform's BLAS (see [DATA_CHANGELOG.md](DATA_CHANGELOG.md), data version 2.1.1). |
| 22 | `total_fuel_consumed_gallons` | gal | scalar | `3.522294e-4 · grams`. |
| 23 | `total_fuel_economy_mpg` | mi/gal | scalar | `(distance·6.213712e-4) / gallons`; `null` when gallons is 0. |

### Distance-to-AV fields (both, in this order)

Downstream distances are **positive**, upstream **negative**. Each is an array
over samples; `null` marks a sample with no qualifying vehicle, and the whole
field is `[]` when none qualifies anywhere in the segment. "Engaged" restricts
to samples where the control vehicle's controller was active.

| # | Field | Unit | Notes |
|---|---|---|---|
| 24 | `downstream_av_id` | id | Nearest control vehicle ahead, per sample. |
| 25 | `distance_to_downstream_av_meters` | m | Distance to it (≥ 0). |
| 26 | `downstream_engaged_av_id` | id | Nearest *engaged* control vehicle ahead. |
| 27 | `distance_to_downstream_engaged_av_meters` | m | Distance to it (≥ 0). |
| 28 | `upstream_engaged_av_id` | id | Nearest *engaged* control vehicle behind. |
| 29 | `distance_to_upstream_engaged_av_meters` | m | Distance to it (≤ 0). |
| 30 | `upstream_av_id` | id | Nearest control vehicle behind. |
| 31 | `distance_to_upstream_av_meters` | m | Distance to it (≤ 0). |

### Additional fields in `full` only

`full` inserts `direction` after `coarse_vehicle_class` and
`reference_a1/a2_…` after `lane_number`, and adds a block of reference and
flat-road energy fields before the distance-to-AV block. The "reference"
quantities come from a reference controller model; the "flat road" quantities
recompute energy with `road_grade = 0`.

| Field | Unit | Notes |
|---|---|---|
| `reference_a1_meters_per_second_per_second` | m/s² | Reference-controller acceleration, model 1. |
| `reference_a2_meters_per_second_per_second` | m/s² | Reference-controller acceleration, model 2. |
| `fuel_rate_flat_road_grams_per_second` | g/s | Fuel rate with grade set to zero. |
| `percent_infeasibility_flat_road` | % | Infeasibility of the flat-road evaluation. |
| `total_fuel_consumed_flat_road_grams` | g | Integral of the flat-road rate. |
| `total_fuel_consumed_flat_road_gallons` | gal | |
| `total_fuel_economy_flat_road_mpg` | mi/gal | |
| `reference_fuel_rate_grams_per_second` | g/s | Fuel rate under the reference controller. |
| `percent_reference_infeasibility` | % | |
| `total_reference_fuel_consumed_grams` | g | |
| `reference_fuel_rate_flat_road_grams_per_second` | g/s | Reference controller, flat road. |
| `percent_reference_infeasibility_flat_road` | % | |
| `total_reference_fuel_consumed_flat_road_grams` | g | |

Only `slim` is used by the paper; `full` is provided for other researchers.

---

## Assembled GPS data — `results/gps/`

One JSON file per day, an array of control-vehicle run records. Produced by
`assemble_data_GPS.m` (Python: `mvtpy.gpsruns`, partial). Field order below.

| Field | Unit | Type | Notes |
|---|---|---|---|
| `av_id` | id | scalar | Control-vehicle number (car N). |
| `assigned_lane` | lane | scalar | Lane assigned to the run. |
| `direction` | ±1 | scalar | −1 westbound, +1 eastbound. |
| `timestamp` | s (POSIX) | array | Resampled to a 10 Hz grid. |
| `latitude` | ° | array | WGS-84 latitude. |
| `longitude` | ° | array | WGS-84 longitude. |
| `x_position` | m | array | Longitudinal position in the MOTION frame. |
| `y_position` | m | array | Lateral position, sign corrected. |
| `controller_engaged` | 0/1 | array | Whether the vehicle's controller was active. |
| `speed` | m/s | array | CAN-reported speed. |
| `is_server_connected` | 0/1 | array | Whether the server ping was live at that time. |
| `first_timestamp` | s | scalar | Start of the run. |
| `last_timestamp` | s | scalar | End of the run. |
| `control_car` | 0/1 | array | Control state merged from the server ping records. |
| `control_last30` | 0/1 | array | Control active at any point in the preceding 30 s. |

Note the two control signals. `controller_engaged` is on-vehicle; `control_car`
is reconstructed from the server pings and is what the distance-to-AV
computation treats as "engaged" when present (falling back to
`controller_engaged` when it is empty).

---

## Macroscopic fields — `results/figures/*/fields_motion_2022-11-DD.mat`

MATLAB `.mat` file (`-v7`, not JSON) written by `generate_macroscopic_fields.m`.
Fields are kernel-weighted averages over a regular time-space grid (5 s × 50 m,
06:00–10:00 local, x ∈ [0, 6500] m, box kernel, ~100 m spatial window).

| Variable | Contents |
|---|---|
| `field` | struct with `value.<F>`, `name.<F>`, `unit.<F>`, `factor.<F>` per field F |
| `t` | time grid (POSIX seconds) |
| `x` | space grid (m) |
| `direction` | −1 westbound |
| `lane` | 0 (all lanes) |

The fields `F` and their definitions:

| Field | Name | Unit | Definition |
|---|---|---|---|
| `Rho` | density | veh/km | kernel-weighted vehicle presence (`factor` 1000) |
| `Q` | flow | veh/h | `Rho`-weighted speed (`factor` 3600) |
| `F` | fuel rate density | g/s | `Rho`-weighted fuel rate |
| `U` | bulk velocity | m/s | `Q / Rho`, NaN where `Rho < 1e-3` |
| `Phi` | bulk fuel rate | g/s | `F / Rho`, NaN where `Rho < 1e-3` |
| `Psi` | bulk fuel consumption | g/m | `F / Q`, NaN where `Q < 1e-2` |

(`full`-style runs also produce the flat-road variants `F0`, `Phi0`, `Psi0`.)
`factor` is the display scaling applied by the plotting code, not baked into
`value`. These fields underlie Figure 3 and SM5.

---

## Sample collections — `results/figures/*/samples_for_distance_analysis_DD.mat`

MATLAB `.mat` file (`-v7.3`) written by `generate_data_samples.m`: flat column
vectors, one entry per qualifying sample, pooled over the day. Every sample is
within 1000 m of an engaged control vehicle. Consumed by `plot_AV_analysis.m`
to produce the fuel-vs-distance results (Figures 2, SM2, SM3).

| Variable | Unit | Notes |
|---|---|---|
| `samples_dist` | m | Signed distance to the engaged AV (− ahead, + behind). |
| `samples_speed` | m/s | Sample speed. |
| `samples_fr` | g/s | Fuel rate. |
| `samples_fcons` | g/m | Fuel consumption. |
| `samples_class` | code | Coarse vehicle class. |
| `samples_xpos` | m | Longitudinal position. |
| `samples_lane` | lane | Lane number. |
| `samples_t` | s | Seconds after 06:00. |
