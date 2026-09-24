# `assemble_data_GPS.m` — stage `gps`

Combines the control vehicles' own 10 Hz recordings with the server-side ping
log, puts them in the I-24 MOTION coordinate frame, and removes each vehicle's
GPS bias by matching it against the MOTION trajectory of the same vehicle.

**Reads** `data/cars/cars_gps/circles_v2_1_car*.csv`, `data/cars/cars_vins.csv`,
`data/cars/veh_ping_202211DD.csv`, `data/i24motion/2022-11-DD/*_0_*.json`
**Writes** `results/gps/CIRCLES_GPS_10Hz_2022-11-DD.json`
**Runs** once per day. Not shardable — the bias estimate needs every segment.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `maxMatchDist` | 6 m | how close an AV and a MOTION trajectory must be to count as the same vehicle |
| `maxMatchSpdDiff` | 2 m/s | and how close in speed |
| `maxMatchLaneDiff` | 0.5 lane | and how close in lane |
| `minMatchTime` | 3 s | how long that must hold before the stretch is used |
| `originXPosition` | 309804.0625 ft | the MOTION origin (Mill Creek Bridge) |
| `inVehShift` | 7 ft | GPS antenna to rear bumper, subtracted along the direction of travel |
| `maxGpsInactive` | 1% | a run with more dropout than this is discarded |
| `minRunLength` | 1.2 km | shorter runs are discarded |
| `minRunTime` | 60 s | shorter runs are discarded |
| `XLIMS` | −400 to 25 400 ft | testbed bounds; samples outside are clipped |

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE
    validate the day; read options (Force, Clean, DryRun, Verbose)
    resolve paths from the location of the code, not the working directory
    declare the output:  results/gps/CIRCLES_GPS_10Hz_2022-11-DD.json
    declare the inputs:  the car CSVs, the VIN map, the ping log, the raw segments
    IF the output exists, is newer than every input, and is newer than the
       code that produced it, AND Force is not set:
        report why it was skipped and RETURN
    IF Clean: delete the existing output
    IF DryRun: RETURN before doing any work

1. PARSE the per-vehicle CSVs into runs            → parse_gps_data
2. PRE-PROCESS those runs onto a common 10 Hz grid → preproc_gps

3. FIND THE AV ACTIVITY WINDOW
   day_window ← 06:00 to 10:00 local, as POSIX seconds
   first_run ← earliest start among runs ending after 06:00
   last_run  ← latest end among runs starting before 10:00
   segment_range ← the 10-minute segments spanning first_run to last_run
   # Segments outside this window contain no control vehicle, so they are never
   # decoded. Progress is reported as "segments m-n of 24" for that reason.

4. MATCH EACH AV RUN TO ITS MOTION TRAJECTORY
   FOR each raw MOTION segment in segment_range:
       decode the segment
       runs_here ← AV runs whose time span overlaps this segment
       IF none: skip the segment
       drop the MOTION fields matching does not use (ids, dimensions, scores)
       assign a lane to every MOTION trajectory        → assign_lanes

       FOR each AV run in runs_here:
           candidates ← MOTION trajectories that overlap it in time
                        AND travel in the same direction
           FOR each candidate:
               express the candidate's x in metres relative to the MOTION origin
               dist_to_av ← AV x interpolated onto the candidate's timestamps,
                            minus the candidate's x
               lane_diff  ← AV lane minus candidate lane

               # Cheap rejection before the expensive per-timestep walk
               IF mean(|dist_to_av|) > 200 m OR mean(|lane_diff|) >= 2 lanes:
                   skip this candidate

               speed_diff ← d(dist_to_av)/dt, smoothed (Gaussian, 3 samples)

               # Walk the candidate, accumulating time while all three
               # thresholds hold, and emit a stretch when they stop holding
               matching_time ← 0
               FOR each timestep after the first:
                   IF |dist_to_av| <= maxMatchDist
                      AND |lane_diff| <= maxMatchLaneDiff
                      AND speed_diff <= maxMatchSpdDiff:
                       remember where this stretch began
                       matching_time ← matching_time + the timestep
                   ELSE:
                       IF matching_time >= minMatchTime:
                           RECORD (av run index, dist_to_av, timestamps)
                       reset matching_time
               IF matching_time >= minMatchTime at the end:
                   RECORD the final stretch

5. ESTIMATE EACH RUN'S GPS BIAS
   FOR each AV run:
       pool dist_to_av over all its recorded stretches
       bias ← median of that pool, ignoring NaN
       IF the run matched nothing: bias ← 0
   # The median, not the mean: a few bad matches should not move it.

6. ADD SERVER CONNECTION STATE
   read the ping log and the VIN map              → get_connection_status
   FOR each AV run, for each timestep:
       connected ← was there a ping within the last 2 seconds?

7. ASSEMBLE THE OUTPUT
   FOR each AV run:
       av_id, assigned_lane, direction        ← from the VIN map and the run
       timestamp, latitude, longitude         ← rounded to 6 decimals
       x_position ← run x in metres, MINUS the run's bias
       y_position ← run y in metres
       controller_engaged, speed              ← from the CAN record
       is_server_connected                    ← from step 6
       first_timestamp, last_timestamp        ← ends of the run
       control_car, control_last30            → get_control_car_status

8. CLIP TO THE TESTBED
   FOR each run: keep only samples with XLIMS(1) < x < XLIMS(2),
                 then reset first_timestamp and last_timestamp

9. WRITE
   encode as JSON and write through a temporary name, so an interrupted run
   cannot leave a truncated file that later looks complete
```

## `parse_gps_data(folder, day)` — CSV to runs

```
day_window ← 03:00 to 18:00 local, as POSIX seconds
FOR each circles_v2_1_car*.csv:
    car number ← the digits at the end of the filename
    read the CSV, keeping rows inside day_window
    split the day into runs wherever the recording stops and restarts
    FOR each run:
        DISCARD it unless:
            it covers at least minRunLength (1.2 km) along x, AND
            it lasts at least minRunTime (60 s), AND
            no more than maxGpsInactive (1%) of it has a dead GPS signal, AND
            it stays inside the testbed in x and below maxRcsy in y
        direction ← sign of travel along x (−1 westbound, +1 eastbound)
RETURN the surviving runs
```

## `preproc_gps(runs, folder)` — onto a 10 Hz grid

```
FOR each run:
    shift x by inVehShift (7 ft) against the direction of travel,
        moving the reported position from the GPS antenna to the rear bumper
    grid ← the run's span, rounded outward to whole tenths of a second
    resample x, y, latitude, longitude, state_x, state_y and CAN speed
        onto that grid                                  → sample_10hz
    resample control_active the same way, then round back to a logical
    assigned_lane ← the lane this car was assigned, from the VIN map
```

`sample_10hz` is linear interpolation onto the new timestamps, with the ends
held rather than extrapolated.

## `assign_lanes(data)` — lane from lateral position

The same two-step method the `slim` stage uses, and described in full in
[`generate_data_mvt_slim.md`](generate_data_mvt_slim.md): estimate the driving
line as a function of x by averaging lateral position on a unit circle of one
lane width, subtract it, then divide by the lane width and median-filter.

Here it exists only to give the matcher a lane to compare against.

## `get_connection_status(pings, vins)` — server view of each vehicle

```
FOR each vehicle in the VIN map:
    collect that VIN's rows from the ping log
    RETURN their timestamps, so the main routine can ask
           "was this vehicle pinging at time t?"
```

## `get_control_car_status(run)` — what "under control" means

```
control_car ← controller_engaged, with a correction:
    the controller reports disengaged whenever the vehicle stops, even though
    the driver has not taken over, so a stop inside an engaged stretch stays
    engaged
control_last30 ← was control_car true at any point in the preceding 30 seconds
```

`control_last30` is what the downstream stages mean by "an engaged AV": a
vehicle that has just disengaged is still influencing the traffic behind it.

## Notes for review

* **The bias correction is the point of stage 1.** Steps 4 and 5 exist only to
  produce one number per run, subtracted in step 7. The GPS units have a
  standing offset from the MOTION frame; matching each run to the trajectory of
  the same physical vehicle measures it.
* **Two rejection thresholds, for different reasons.** The 200 m / 2-lane test
  in step 4 is a cheap filter to avoid walking obviously unrelated
  trajectories. The `maxMatch*` thresholds decide what actually counts as a
  match. Only the second set affects the result.
* **Steps 4 and 5 cannot be sharded**, because the median in step 5 needs every
  stretch from every segment. Days are independent, so parallelise over days.
