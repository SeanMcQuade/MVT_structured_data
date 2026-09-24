# MVT pipeline algorithms

For a stage-by-stage walkthrough in structured English, see
[`pseudocode/`](pseudocode/) — one file per script, at a level between this
document and the MATLAB.

What each stage computes, at a level a reviewer can check against the paper and
a reimplementer can follow. Exact field lists and units are in
`DATA_DICTIONARY.md`; number formatting is in `MATLAB_JSON_FORMAT.md`. The
Python port (`python/mvtpy/`) mirrors these steps and is verified byte-for-byte
against the released output, so its modules are cited as executable references.

## Data-flow graph

```
data/cars/*                      ── assemble_data_GPS ──▶ results/gps/CIRCLES_GPS_10Hz_2022-11-DD.json
data/i24motion/2022-11-DD/*    ──┘

results/gps + raw i24motion    ── generate_data_mvt_slim ─▶ results/slim/2022-11-DD/*.json   (24 files)
                               ── generate_data_mvt_full ─▶ results/full/2022-11-DD/*.json   (24 files)

results/slim/2022-11-DD        ── generate_data_samples ──▶ analysis/2022-11-DD/samples_for_distance_analysis_DD.mat
results/slim/2022-11-DD        ── generate_macroscopic_fields ▶ analysis/2022-11-DD/fields_motion_2022-11-DD.mat
raw i24motion                  ── generate_orig_dist_lanes ▶ analysis/2022-11-DD/I-24MOTION_*_orig_dist_lane.mat  (24/day)
slim + the lane sidecars       ── extract_lane_changes_v_dist_to_av ▶ analysis/2022-11-DD/LC_data_DD.mat
results/slim/2022-11-DD        ── relative_speed_histogram ▶ analysis/2022-11-DD/relspeed_data_DD.mat

fields + gps                   ── plot_macroscopic_fields ─▶ fig_field_*.png            (Fig. 3, SM5)
results/slim (or full)         ── plot_microscopic_trajectories ▶ fig_motion_trajectories_*.png
all three days' samples        ── plot_AV_analysis ───────▶ fig_2_*, fig_3_*, fig_SM2_*  (Fig. 2, SM2, SM3)
LC_data + gps                  ── plotting_LC_analysis ───▶ fig_lc_*.png                 (4/day)
relspeed_data                  ── plot_relative_speed ────▶ fig_relspeed_*.pdf           (2/day)
```

The derived `.mat` intermediates live under `results/analysis/`, apart from the
figures built from them, so that the inputs to a figure can be distributed
without the figure. Stages that read only those `.mat` files need no trajectory
data at all.

Coordinate frame throughout: I-24 MOTION v2, origin at the Mill Creek Bridge
(`x = 309804.0625 ft`), x increasing eastbound. `ft2m = 0.3048`,
`m2mi = 6.213712e-4`, `g2gal = 3.522294e-4`.

---

## Stage 1 — `assemble_data_GPS`

Assembles the control-vehicle (AV) GPS/CAN recordings and the server ping
records into one file of testbed *runs* per day, each matched to its I-24 MOTION
trajectory. (Python: `mvtpy.gpsruns`, run-splitting only.)

1. **Parse runs** (`parse_gps_data`). For each vehicle `circles_v2_1_carN.csv`,
   restrict to the day (03:00–18:00 local), then split into runs: a run starts
   where the vehicle moves inward across the roadway edge while inside the
   testbed bounds (`|rcs_y| < 150 ft`, `−3000 < rcs_x < 29100 ft`) and ends
   where it leaves those bounds. Keep runs longer than 60 s and 1.2 km with a
   healthy GPS fix for ≥ 99% of samples. Run direction is the sign of net x
   displacement.
2. **Preprocess** (`preproc_gps`, `sample_10hz`). Resample each run onto a
   uniform 10 Hz grid.
3. **Merge server state** (`get_connection_status`, `get_control_car_status`).
   Join the 1 Hz `veh_ping_202211DD.csv` records by VIN and time to recover the
   server-side control state (`control_car`), server connectivity
   (`is_server_connected`), and a 30-second look-back (`control_last30`).
   `control_car` is the authoritative "engaged" signal; the on-vehicle
   `controller_engaged` is the fallback.
4. **Assign lanes** (GPS-specific `assign_lanes`) from the lateral position.
5. **Match to MOTION** trajectories in the overlapping 10-minute segments and
   write one JSON record per run.

The matching pass accumulates across MOTION segments into a single output, which
is why this stage is not shardable (unlike stage 2).

---

## Stage 2 — `generate_data_mvt_slim` / `generate_data_mvt_full`

Turns each raw 10-minute MOTION segment into a processed segment. `slim` keeps
westbound trajectories; `full` keeps both and adds reference/flat-road energy.
(Python: `mvtpy.lanes`, `mvtpy.kinematics`, `mvtpy.fuel`, `mvtpy.avdist`,
`mvtpy.slim`.)

Per segment:

1. **Driving line** (`assign_lanes`, whole-segment statistic). Sample westbound
   trajectories longer than 5 s at 1 Hz, drop lateral outliers, divide x into
   200 cells, and in each cell with > 20 samples take the circular mean of the
   lateral offsets (mapped onto a unit circle of period = lane width). This is
   the roadway's lateral "wiggle" as a function of x.
2. **Lane assignment**, per trajectory. Correct lateral position by the
   interpolated driving line: `y_corr = Sw·(y − drivingLine(x)) + Cw`
   (westbound `Sw = 0.98, Cw = 1`). Convert to a fractional lane index, clamp to
   `[0, 5]`, and median-filter over a 10-sample window.
3. **Lane-change clipping** (`clip_lane_changes`). Walk each trajectory, emit a
   segment (`<id>-N`) for each stretch that stays within one lane for at least
   0.5 s, and discard the parts spent changing lanes. This is where one raw
   trajectory becomes several released records.
4. **Kinematics** (`mvtpy.kinematics`). Longitudinal position in meters from the
   origin; speed by central difference of distance vs time (ends duplicated);
   acceleration by second central difference (first/last repeated).
5. **Road grade.** Look up `x` in `Models/Eastbound_grade_fit.csv` (piecewise
   linear percent grade, clamped to the mapped range) and take `arcsin`;
   negate for westbound.
6. **Fuel** (`mvtpy.fuel`). Pick the model from the coarse vehicle class
   (0→midBase, 1→midSUV, 2/3/5→Pickup, 4→Class8Tractor) and evaluate the
   instantaneous rate from `(v, a, grade)` with infeasible accelerations
   projected to the feasibility boundary. Two model families:
   - **light-duty** (Compact, Pickup, midBase, midSUV): a fitted
     cruise+accel+grade polynomial, floored at `beta0` below the cut-off speed
     `vc`, with fuel cut entirely above `vc` on hard braking.
   - **heavy-duty** (Class4PND, Class8Tractor): the same polynomial floored by a
     linear term `h0 + h1·v`, no fuel cut.
   Integrate the rate over time (trapezoid) for total grams, then gallons and
   mpg.
7. **Distance to AVs** (`mvtpy.avdist`). For each sample, find the nearest
   control vehicle ahead (downstream, +) and behind (upstream, −), both
   unconditionally and restricted to engaged samples, by interpolating each
   candidate vehicle's position onto the trajectory's timestamps (NaN outside
   the vehicle's own time range). Same lane and direction required.
8. **Round and write.** Round every numeric field to 4 decimals and write the
   record in the fixed field order (`mvtpy.slim.FIELD_ORDER`) via `jsonencode`.

**Validation.** The Python port reproduces a full released segment byte-for-byte
(md5 `3a7bfc02…`, 409 MB); each intermediate quantity is checked against the
released values in `python/tests/`.

---

## Stage 3 — `generate_data_samples`

Pools, over one day, every trajectory sample that lies within 1000 m of an
*engaged* control vehicle, recording distance-to-AV alongside speed, fuel rate,
fuel consumption, class, position, lane, and time (seconds after 06:00). Reads
the slim segments (falls back to full). Output feeds the AV-effect analysis.
Distances are signed: negative ahead of the AV, positive behind.

---

## Stage 4 — `generate_macroscopic_fields`

Converts the microscopic trajectories of one day into macroscopic traffic-state
fields on a regular grid (5 s × 50 m, 06:00–10:00 local, x ∈ [0, 6500] m).

1. Build the (t, x) grid and a box kernel with ~100 m spatial window and time
   window `hx / v_char` (`v_char = 20 m/s`).
2. Stream each segment's trajectories; for every sub-sampled sample, add its
   kernel-weighted contribution within the window to three accumulators:
   density `Rho`, flow `Q` (speed-weighted), and fuel-rate density `F`.
3. Normalize by the kernel factor, then derive:
   - `U = Q / Rho` (bulk velocity), NaN where `Rho < 1e-3`
   - `Phi = F / Rho` (bulk fuel rate), NaN where `Rho < 1e-3`
   - `Psi = F / Q` (bulk fuel consumption), NaN where `Q < 1e-2`
4. Save `field`, `t`, `x`, `direction`, `lane`.

`Rho` is scaled ×1000 (veh/km) and `Q` ×3600 (veh/h) for display via
`field.factor`; the stored `value` is in SI.

---

## Stages 3b–3d — the lane-change and relative-speed intermediates

Three analyses that were run by hand until 2026-09 and are now stages. Each
writes a `.mat` under `results/analysis/` and is consumed by a plotting stage,
so the figures can be rebuilt without the trajectories.

### `generate_orig_dist_lanes` — lane origin and destination

One `.mat` per raw segment, recording for every released trajectory the lane it
came **from** and the lane it went **to**. The slim data reports only the lane a
clipped trajectory was driven in, so the lane change that produced the clip is
not recoverable from it.

The method is the lane identification and clipping of stage 2, re-run with the
origin and destination carried through the clip. Walking a trajectory:

1. Advance to where the lane estimate is stable
   (`|dlane/dt| < MaxLaneChangeRate`, 0.1 lane widths/s).
2. `origin_lane` for the first clipped piece is the rounded lane at that point;
   for every later piece it is the lane of the piece before it.
3. Find the first index where `|lane - tempLane| > LaneChangeThresh`
   (0.5 lane widths). `destination_lane` is `tempLane ± 1` in the direction of
   the change, or `tempLane` itself if the trajectory ends without changing.
4. Emit the piece, skip forward past the change, and repeat.

> **This code is a fork.** `assign_lanes` and `clip_lane_changes` are duplicated
> from `generate_data_mvt_slim.m` and extended in place. The two copies must
> produce the same clipping, because the sidecar is paired with the slim JSON
> **by index** — entry *i* describes released trajectory *i*. Nothing in the
> staleness machinery couples them, so editing one does not invalidate the
> other's outputs. The consumer checks the counts match and errors if they do
> not, which catches the coarse failure but not a subtler divergence.

### `extract_lane_changes_v_dist_to_av` — lane-change events

Turns a day of slim trajectories plus their sidecars into a list of lane-change
events annotated with how far away the nearest control vehicle was.

A trajectory contributes an event if `origin_lane ≠ lane_number` (it merged
**in**, recorded at the first timestamp) or `destination_lane ≠ lane_number` (it
merged **out**, recorded at the last). The signed magnitude is
`lane_number - origin_lane` for merges in and `destination_lane - lane_number`
for merges out, so positive and negative distinguish the direction of the move.

Each event is emitted **once per AV it can be measured against** — once for the
nearest downstream AV and once for the nearest upstream one, when each exists —
and each row also carries the distance to the nearest *engaged* AV if there was
one. So a single lane change can produce two rows, and the counts in
`LC_data_DD.mat` are row counts, not vehicle counts.

### `relative_speed_histogram` — relative speed to the engaged AV

Pools, over a whole day, the rate at which traffic closes on or falls away from
the engaged AV ahead of it. Per trajectory:

1. Keep the timesteps where a downstream engaged AV exists
   (`~isnan(downstream_engaged_av_id)`).
2. Mark the indices where the tracked AV changes identity, or where the distance
   jumps by more than 70 m/s — a discontinuity in *which* AV is being measured,
   not real motion.
3. Differentiate the distance with a centred stencil, made one-sided at exactly
   those marks, so a switch never contributes a spurious speed.
4. Smooth with 40 passes of a 3-point average using the same stencils.
5. Keep samples between `lower_Bnd` = 30 m and `upper_Bnd` = 350 m. Below 30 m
   the distance is dominated by the AV's own footprint; above 350 m the nearest
   AV is barely influencing the vehicle.

The surviving distance/speed pairs for all 24 segments are saved. Splitting the
save from the plotting is what lets the figures be redrawn without the
trajectories.

---

## Stage 5 — plotting

- **`plot_macroscopic_fields`** → **Figure 3, SM5**. One image per field
  (`Rho, Q, F, U, Phi, Psi`), optionally overlaid with the control-vehicle GPS
  traces coloured by whether control was engaged. Per-field colorbar limits are
  fixed so the three days are visually comparable.
- **`plot_microscopic_trajectories`** → supplementary time-space plots. Each
  trajectory is drawn as a patch whose height encodes vehicle length. Builds
  `*_reduced.mat` caches (in `results/.mvt/cache/`) on first run; the
  high-resolution mode is memory-hungry.
- **`plot_AV_analysis`** → **Figure 2, SM2, SM3**. Pools all three days'
  samples, bins them by signed distance to the nearest engaged AV (350 m
  half-window, a masked band ±30 m around the AV, a 1 m exclusion at 0),
  computes mean / median / "effective" fuel consumption per bin and the
  active-AV counts over time (06:45–09:15), and renders the comparison grids.
  These cross-day figures are written to the shared `results/figures/` folder.
- **`plotting_LC_analysis`** → the lane-change figures, four per day, for lanes
  2, 3 and 4 and for engaged AVs only. See below; this is the most involved of
  the plotting stages.
- **`plot_relative_speed`** → the relative-speed histogram for the day, with
  mean, median, standard deviation and interquartile markers, plus (via
  `binned_relative_speed`) the mean and median relative speed in 10 m bins of
  distance behind the AV. Bin edges run from −350 m to 350 m with an extra
  ±1 m (`AVLOCBUFFER`) pair inserted around zero, so samples at the AV's own
  position do not land in a neighbouring bin.

### `plotting_LC_analysis`: exposure, rate, and cumulative excess

A raw count of lane changes at a given distance from an AV says little, because
traffic spends unequal amounts of time at each distance. The stage therefore
normalises by exposure.

**Exposure.** For each lane, take the control vehicles whose modal
`assigned_lane` is that lane, restricted to before 11:00 local. Lay a 1 m grid
over the span of observed merge positions. At each distinct GPS timestamp,
compute the signed distance from every grid point to the nearest AV ahead and
the nearest behind, histogram those distances into the 10 m bins, divide by the
10 grid points per bin, and weight by the GPS timestep in minutes. Summing over
timestamps gives, per bin, the vehicle-minutes spent at that distance from an
AV. Both an all-AV and an engaged-AV series are accumulated; the figures use the
engaged one.

**Rate.** Merge counts per bin divided by that bin's exposure minutes, so the
unit is merges per minute. Left and right moves are summed for the excess
calculation. A dashed **basal rate** is the mean rate over 350–500 m, taken as
the far-field value where the AV is not influencing behaviour.

**The masked band.** Four bins either side of zero (roughly within 35 m) are
drawn faded and excluded from the fit: at those distances the lane-change
detection and the AV distance are least trustworthy.

**Excess and its extrapolation.** `e(x)` is the rate minus the basal rate, over
30–500 m. Because the masked band leaves no data below ~35 m, a straight line is
fitted to the first ten reliable points (about 35–125 m) and evaluated at 5, 15
and 25 m to fill the gap. `g(x)` is then the cumulative sum of `e(x)` over bins,
plotted solid to 350 m and faded beyond it.

> The extrapolation is a modelling choice, not a measurement: the three points
> nearest the AV are inferred from the trend further out, and they enter every
> later value of `g(x)` through the cumulative sum.

---

## Reproducibility notes

- **Determinism.** Every stage is deterministic: no RNG, no wall-clock input.
  Re-running produces byte-identical JSON and content-identical `.mat` (see
  `make verify-full`).
- **Parallelism changes nothing.** Stage 2 and `generate_orig_dist_lanes` shard
  over the 24 segments; `SHARDS=4` was verified to produce byte-identical output
  to a serial run, and the sidecars rebuilt across 6 shards reproduced a serial
  run's `mvt.matHash` exactly.
- **Figure canvases are pinned.** The relative-speed stages set an explicit
  figure size, because MATLAB otherwise derives the default from the screen and
  `exportgraphics` crops the page to the drawn content — the same code then
  produced a different page size under `matlab -batch` than interactively.
- **Units are converted explicitly** with named factors; MOTION positions and
  vehicle dimensions arrive in feet, processed outputs are meters/mph/gallons.
- **The one subtlety in rounding** is that MATLAB's `round(x, 4)` snaps
  near-ties (within one ULP below the midpoint) away from zero; reproducing it
  is what makes byte parity possible. See `MATLAB_JSON_FORMAT.md` and
  `python/mvtpy/matround.py`.
