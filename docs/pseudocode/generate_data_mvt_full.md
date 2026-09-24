# `generate_data_mvt_full.m` — stage `full`

The same processing as `slim`, keeping both directions and computing three
extra families of fields: travel direction, a smooth *reference* trajectory for
each vehicle, and a second fuel estimate as if the road were flat.

**Reads** the same inputs as `slim`
**Writes** `results/full/2022-11-DD/I-24MOTION_<date>_<time>.json` (24 per day)
**Runs** once per day, shardable by segment. **Opt-in**: `make all` does not
build it, and nothing downstream reads it.

Read [`generate_data_mvt_slim.md`](generate_data_mvt_slim.md) first — everything
there applies. This file covers only the differences.

## What `full` does differently

```
# 1. It does not filter by direction.
#    slim:  KEEP only westbound (direction < 0)
#    full:  keep everything, and record which way each vehicle was going
data(i).direction ← veh.direction        # −1 westbound, +1 eastbound

# 2. It keeps one more raw field.
#    slim drops 9 fields; full drops the same list minus 'direction'.

# 3. It adds a reference trajectory (see below).

# 4. It adds a flat-road fuel estimate: the same fuel model, called again with
#    the road grade set to zero.
fuel_rate_flat_road_grams_per_second     ← model(speed, accel, theta = 0)
percent_infeasibility_flat_road
total_fuel_consumed_flat_road_grams      ← same trapezoidal integral
total_fuel_consumed_flat_road_gallons
total_fuel_economy_flat_road_mpg

# 5. It runs the fuel model a third time, on the reference trajectory.
reference_fuel_rate_grams_per_second
percent_reference_infeasibility
total_reference_fuel_consumed_grams
reference_fuel_rate_flat_road_grams_per_second
percent_reference_infeasibility_flat_road
total_reference_fuel_consumed_flat_road_grams
```

## The reference trajectory

A hypothetical version of the same trip: same start and end position, same
start and end speed, same duration, but travelled in **two constant-acceleration
phases** instead of whatever the vehicle actually did. It is the baseline the
fuel comparison is made against — what this vehicle would have burned covering
the same ground without the accelerations that traffic imposed on it.

```
GIVEN x and t measured from the trajectory's own start, and v its speed

# Split the trip in half by time and solve for the two accelerations that
# make the reference start and end exactly where the real trajectory did,
# at the same speeds.
t1 ← t(end) / 2                                  # the switch point

a2 ← (x(end) − x(1) − ½·t1·v(1) − v(end)·t(end) + ½·t1·v(end))
     ÷ (−½·t(end)² + ½·t1·t(end))
a1 ← (v(end) + a2·(t1 − t(end)) − v(1)) ÷ t1

vRef ← v(1) + a1·t                      for t <= t1
       v(end) + a2·(t − t(end))         for t >  t1
aRef ← a1 for t <= t1, else a2
xRef ← the corresponding position, integrated piecewise

# A two-phase fit can imply reversing, which is not physical on a motorway.
IF vRef goes negative anywhere:
    # Fall back to three phases: decelerate to a stop, wait, accelerate away.
    tt1 ← (x(end) − x(1)) / mean(v(1), v(end))       # when it stops
    tt2 ← t(end) − tt1                               # when it starts again
    aa1 ← −mean(v(1), v(end)) / (x(end) − x(1)) · v(1)
    aa2 ← +mean(v(1), v(end)) / (x(end) − x(1)) · v(end)
    vRef ← v(1) + aa1·t   for t <= tt1
           0              for tt1 < t < tt2
           v(end) + aa2·(t − t(end))   for t >= tt2
    aRef, xRef ← correspondingly
    reference_a1, reference_a2 ← aa1, aa2

RECORD reference_a1_meters_per_second_per_second
       reference_a2_meters_per_second_per_second
RUN the fuel model on (vRef, aRef, theta) and again on (vRef, aRef, 0)
```

## Notes for review

* **`full` is a dead end by design.** Nothing downstream reads it; it exists so
  the eastbound data and the reference comparison are available, and for the
  eastbound variant of the trajectory plots. `make all` skips it, `mvt.status`
  reports it as opt-in rather than pending, and it is not in the checksum
  manifests — though it has been verified byte-identical to the Python port by
  direct comparison.
* **It is roughly 60% larger than `slim`** (87 GB against 55 GB), from carrying
  both directions and the extra fields.
* **The fuel model runs four times per trajectory** in `full` — real and
  reference, each with and without grade — against once in `slim`. That is most
  of the extra runtime.
