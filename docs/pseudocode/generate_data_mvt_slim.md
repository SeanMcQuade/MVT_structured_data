# `generate_data_mvt_slim.m` — stage `slim`

Turns each raw 10-minute I-24 MOTION segment into the released westbound data
set: lanes assigned, lane changes clipped into separate trajectories, fuel
estimated per vehicle, and distance to the nearest control vehicle recorded.

**Reads** `data/i24motion/2022-11-DD/*_0_*.json`, `results/gps/*.json`, `Models/`
**Writes** `results/slim/2022-11-DD/I-24MOTION_<date>_<time>.json` (24 per day)
**Runs** once per day, shardable by segment.

This is the stage the released data comes from, and the one whose output is
verified byte-for-byte against the Python port.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `LaneWidth` | 12 ft | one lane, used as the unit for every lateral threshold |
| `Nr_XCells` | 200 | cells along x used to estimate the driving line |
| `MinCellSamples` | 20 | a cell with fewer samples is assumed unshifted |
| `Y_Up_Lim`, `Y_Low_Lim` | 5, 0.5 lane widths | lateral bounds; samples outside are outliers |
| `Sw`, `Cw` | 0.98, 1 | scale and offset applied to corrected westbound y |
| `LaneChangeThresh` | 0.5 lane | lateral movement that counts as a lane change |
| `MaxLaneChangeRate` | 0.1 lane/s | below this the vehicle is holding its lane |
| `MinClipTime` | 0.5 s | shorter clipped pieces are discarded |
| `ChangeBufferThresh` | 0.2 lane | how close to the new lane before the vehicle has settled |
| `originXPosition` | 309804.0625 ft | the MOTION origin (Mill Creek Bridge) |
| `mcDist` | 0.225 mile | offset from that origin to the grade map's origin |
| `flag_deterministic_quadrature` | `true` | compensated summation for the fuel integral |

Conversions are named: `ft2meterFactor`, `meter2mileFactor`, `g2gallonsFactor`.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

SETUP
    segments  ← the manifest: raw file → output name, without decoding anything
    ERROR unless there are 24 segments
    read the road grade fit from Models/Eastbound_grade_fit.csv
    decode the day's GPS file once, for use by every segment
    put Models/ on the path  (the fuel models are dispatched by name)

FOR each segment owned by this shard:        # worker k of N takes k, k+N, ...
    # Freshness is decided from the manifest, BEFORE the multi-GB decode
    IF the output exists and is newer than the raw file and the code,
       AND Force is not set:
        report the skip and CONTINUE
    IF Clean: delete the existing output
    IF DryRun: CONTINUE

    decode the raw segment
    CHECK that the name derived from the decoded data matches what the
          manifest predicted, else ERROR (the manifest is stale)

    KEEP only westbound trajectories (direction < 0)
    DROP the fields the released data does not carry
         (flags, compute_node_id, fragment_ids, merged_ids, configuration_id,
          fine_vehicle_class, x_score, y_score, road_segment_ids)

    dataLanes ← ASSIGN LANES to every trajectory        → assign_lanes
    dataTemp  ← CLIP LANE CHANGES, splitting each
                trajectory into single-lane pieces      → clip_lane_changes
    distToAvs ← DISTANCE TO CONTROL VEHICLES            → calculate_distance_to_avs

    FOR each clipped trajectory:
        COPY THROUGH
            trajectory_id, timestamp, coarse_vehicle_class,
            first_timestamp, last_timestamp, length, width, height,
            lane_number
        POSITION
            x_position_meters ← (x − origin) in metres
            y_position_corrected_meters ← corrected y in metres
            starting_x, ending_x ← its first and last x
        KINEMATICS
            work in metres from the trajectory's own start, and seconds from
            its own first timestamp
            speed        ← central difference of x, one-sided at the ends
            acceleration ← second central difference, ends repeated
            total_distance_traversed_meters ← |x(end) − x(1)|
        ROAD GRADE
            position along the grade map ← x in miles − mcDist
            find the map cell containing it, clamped to the map's ends
            theta ← asin(slope·x + intercept), negated for westbound
        FUEL
            vehicle class → model:
                0 → midBase      1 → midSUV       2 → Pickup
                3 → Pickup       4 → Class8Tractor  5 → Pickup
            fuel_rate ← fuel_model_<class>_simplified(speed, accel, theta)
            percent_infeasibility ← share of timesteps the model could not meet
            total_fuel_consumed_grams ← trapezoidal integral of fuel_rate over t
                # Compensated summation (mvt.neumaierDot) rather than dot(),
                # so the result does not depend on the platform's BLAS.
            total_fuel_consumed_gallons ← grams × g2gallonsFactor
            total_fuel_economy_mpg ← miles travelled ÷ gallons
        DISTANCE TO AVs
            downstream_av_id, distance_to_downstream_av_meters
            downstream_engaged_av_id, distance_to_downstream_engaged_av_meters
            upstream_av_id, distance_to_upstream_av_meters
            upstream_engaged_av_id, distance_to_upstream_engaged_av_meters
            (empty where no such AV was on the road)
        ROUND every numeric field to 4 decimals
            # Once, here, immediately before encoding. Intermediate values stay
            # at full precision. This is what makes byte-identical output
            # reproducible across implementations.

    encode as JSON and write through a temporary name
    release the decoded segment before the next one

take Models/ off the path
```

## `assign_lanes(data, opts)` — lateral position to lane number

The road is not straight in the MOTION frame, so a fixed lateral threshold
would put the same physical lane in different lane numbers at different points
along the road. The fix is to estimate where the middle of a lane actually is
as a function of x — the *driving line* — and measure against that.

```
1. SAMPLE
   FOR each westbound trajectory longer than 5 seconds:
       take its x and y once per second
   DISCARD samples outside 0.5 to 5 lane widths laterally

2. ESTIMATE THE DRIVING LINE
   divide the sampled x range into Nr_XCells (200) cells
   FOR each cell:
       y_cell ← the cell's y samples, shifted down by half a lane width
       IF the cell holds more than MinCellSamples (20):
           # Average on a circle whose circumference is one lane width, so
           # that vehicles in different lanes reinforce rather than cancel:
           # lane position is periodic in y with period LaneWidth.
           p ← mean(cos(2π·y_cell / LaneWidth))
           q ← mean(sin(2π·y_cell / LaneWidth))
           driving_line(cell) ← atan2(q, p) · LaneWidth / 2π
       ELSE:
           driving_line(cell) ← 0        # assume no shift

3. CORRECT AND ASSIGN
   FOR each trajectory:
       y_corr ← Sw · (y − driving_line interpolated at this x) + Cw
       lane_raw ← (|y_corr| − LaneWidth/2) / LaneWidth, clamped to 0…5
       lane ← median filter of lane_raw over a 10-sample window,
              with the ends held
RETURN y_corr and lane per trajectory
```

Lane 1 is the leftmost (HOV) through 4 rightmost; 0 is off the highway and 5 is
an on/off ramp.

## `clip_lane_changes(data, dataLanes, opts)` — one trajectory per lane

A vehicle that changes lane would otherwise appear as one trajectory spanning
two lanes, which makes per-lane statistics meaningless. This cuts it into
pieces, each entirely within one lane, and discards the transitions.

```
FOR each trajectory:
    pointer ← start of the trajectory
    WHILE not at the end:
        # Start where the vehicle is holding its lane
        skip forward to the first point where |d(lane)/dt| < MaxLaneChangeRate
        IF there is no such point: DISCARD the whole trajectory
        current_lane ← round(lane at that point)

        # Find where it leaves that lane
        change_at ← first point where |lane − current_lane| > LaneChangeThresh
        IF there is none: change_at ← the end of the trajectory

        # Walk back from the change to where the vehicle was last settled
        settled_at ← searching backwards from change_at, the last point
                     within ChangeBufferThresh of current_lane and still
                     below MaxLaneChangeRate

        IF the piece from the start to settled_at lasts at least MinClipTime:
            EMIT it as a trajectory of its own:
                same fields, clipped arrays, lane = current_lane
                trajectory id gets a suffix: -0, -1, -2, …
        # else the piece is too short to be useful, and is dropped

        # Skip past the change itself, to where the vehicle has settled in
        # its new lane, and continue from there
        advance the pointer past the transition
        # A move that starts and does not complete — the vehicle drifts and
        # returns — is treated as a false change: clipped, but not the start
        # of a new lane.
```

## `calculate_distance_to_avs(dataTemp, dataGPS)` — proximity to control vehicles

```
FOR each trajectory:
    avs_on_road ← control vehicle runs that
        are in the SAME lane, AND
        travel in the same direction, AND
        overlap this trajectory in time

    FOR each of upstream and downstream:
        FOR each such AV run:
            project the AV's x onto this trajectory's timestamps
            signed_distance ← (AV x − trajectory x) × direction
                # so positive is always downstream, whichever way traffic runs
            keep it only if its sign matches the direction being measured
            engaged ← the AV's control_car state, interpolated to these times
        AT EACH TIMESTEP:
            distance_to_av ← the smallest such distance, and which AV it was
            distance_to_engaged_av ← the same, over engaged AVs only
    RETURN empty where no AV qualified
```

"Engaged" uses `control_car`, which holds through a stop — see
[`assemble_data_GPS.md`](assemble_data_GPS.md).

## Notes for review

* **Order matters.** Lanes are assigned before clipping, because clipping needs
  the lane estimate; distance to AVs is computed after clipping, so it is
  measured per released trajectory rather than per original one.
* **Rounding happens once**, immediately before encoding. This is deliberate
  and is what makes byte-identical reproduction possible; see
  `docs/MATLAB_JSON_FORMAT.md`.
* **The fuel model is dispatched by name** — `eval` on a string built from the
  vehicle class — so a class name with no matching file in `Models/` fails at
  run time, not load time.
* **`full` is the same routine** with eastbound and reference trajectories kept
  and extra fields computed. See
  [`generate_data_mvt_full.md`](generate_data_mvt_full.md).
