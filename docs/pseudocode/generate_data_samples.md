# `generate_data_samples.m` — stage `samples`

Flattens a day of trajectories into one long list of samples, each one a single
vehicle at a single instant, tagged with how far it was from an engaged control
vehicle. This is the input to the paper's fuel-versus-distance result.

**Reads** `results/slim/2022-11-DD/*.json` (falls back to `full/` if slim has
fewer than 24 segments)
**Writes** `results/analysis/2022-11-DD/samples_for_distance_analysis_DD.mat`
(~1.4 GB per day)
**Runs** once per day. Not sharded.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `Max_Dist` | 1000 m | a sample is kept only if an engaged AV was within this |
| `sixAM18` | 1668772800 | POSIX seconds at 06:00 on 18 Nov 2022, the time origin |

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

LOCATE the day's processed segments:
    prefer results/slim/2022-11-DD
    IF fewer than 24 are there, try results/full/2022-11-DD
    IF still fewer than 24: ERROR

FOR each of the 24 segments:
    decode it
    SKIP the segment if it has no distance-to-engaged-AV field at all
        # no control vehicle was on the road during it
    EXTRACT its samples                          → stats_to_av_dist
    APPEND them to the day's running vectors

SAVE all eight vectors through a temporary name
```

## `stats_to_av_dist(data, day)` — trajectory to samples

```
GIVEN one decoded segment

PRE-ALLOCATE eight vectors, each as long as every timestep in the segment,
    in narrowed types chosen to keep the file manageable:
        speed, fuel rate, distance   double
        x position                   int16   (metres, rounded)
        vehicle class, lane          uint8
        time                         uint16  (seconds after 06:00)

FOR each trajectory:
    SKIP it if it is eastbound
        # 'full' carries a direction field; 'slim' is westbound already

    SKIP it unless an engaged AV was recorded either upstream or downstream

    FOR each timestep in the trajectory:

        IF an engaged AV was DOWNSTREAM and within Max_Dist:
            EMIT a sample:
                dist  ← distance to that AV        (positive)
                speed ← speed at this instant
                fr    ← fuel rate at this instant
                class ← the vehicle's coarse class
                xpos  ← x position, rounded to whole metres
                lane  ← lane number
                t     ← seconds since 06:00 on this test day

        IF an engaged AV was UPSTREAM and within Max_Dist:
            EMIT a second sample, identically, with
                dist ← distance to that AV         (negative)

    # A timestep with an engaged AV both ahead and behind produces TWO samples.

TRIM the vectors back to the number of samples actually written

fcons ← fr ÷ (speed + 1e-6)
    # Instantaneous fuel consumption, grams per metre rather than per second.
    # The epsilon avoids dividing by zero for a stopped vehicle; a sample at
    # rest therefore reports a very large consumption, which the downstream
    # speed filter removes.
```

## The distance sign convention

`samples_dist` is signed, and the sign says where the AV was relative to the
sampled vehicle:

* **positive** — the engaged AV was **downstream** (ahead); the sample is
  *behind* the AV
* **negative** — the engaged AV was **upstream** (behind); the sample is
  *ahead of* the AV

The fuel-versus-distance figures put this on the x-axis, so the AV sits at zero
with traffic behind it to the right.

## The time origin

`samples_t` is seconds after 06:00 **on that day**, computed by taking the known
06:00 on 18 November and stepping back a whole day per day of difference. It
fits in a `uint16`, which caps it at about 18 hours after 06:00 — ample for a
test window that ends mid-morning.

## Notes for review

* **Roughly 89 million samples per day**, which is why the types are narrowed
  and why the file is 1.4 GB. `x` in particular is rounded to whole metres by
  storing it as `int16`.

* **One timestep can yield two samples**, and both are independent rows: the
  same vehicle at the same instant appears once measured against the AV ahead
  and once against the AV behind. Anything that counts rows is counting
  (vehicle, instant, AV) triples, not vehicles.

* **The `full` fallback is a convenience** for a tree where `slim` was never
  built. The samples differ in one way: `full` contains eastbound trajectories,
  which are skipped explicitly, so the result is the same either way.

* **`fcons` divides by speed**, so it is only meaningful where the vehicle was
  moving. `plot_AV_analysis` applies a minimum-speed filter before using it.
