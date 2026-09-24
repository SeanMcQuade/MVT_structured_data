# `relative_speed_histogram.m` — stage `relspeed`

Measures how fast traffic is closing on, or falling away from, the engaged
control vehicle ahead of it, and pools a whole day of those samples.

**Reads** `results/slim/2022-11-DD/*.json`
**Writes** `results/analysis/2022-11-DD/relspeed_data_DD.mat` (~0.7 GB per day)
**Runs** once per day, about 5 minutes. Not sharded.

This is the expensive half of the relative-speed analysis. It saves its result
so that [`plot_relative_speed`](plot_relative_speed.md) can redraw the figures
in seconds, without the trajectories.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `j_start`, `j_end` | 1, 24 | which of the day's segments to pool |
| `edges` | −30 : 0.5 : 30 m/s | histogram bins (used by the plotting stage) |
| `lower_Bnd` | 30 m | nearest distance to an AV that is kept |
| `upper_Bnd` | 350 m | furthest distance that is kept |
| jump threshold | 70 m/s | above this, the distance changed discontinuously |
| smoothing passes | 40 | 3-point moving averages applied to the speed |

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

slim_files ← the day's 24 released segments, from mvt.expectedOutputs
    # Not from the raw manifest: this stage reads only slim, so it must work
    # on a download with no raw data.

FOR j in j_start..j_end:
    decode slim segment j

    FOR each trajectory in the segment:
        t        ← its timestamps
        av_id    ← downstream_engaged_av_id at each timestep
        av_dist  ← distance_to_downstream_engaged_av_meters at each timestep

        KEEP only the timesteps where an engaged AV was present (av_id not NaN)
        SKIP the trajectory if fewer than 2 such timesteps remain

        # ---- Where the measurement is discontinuous ----------------------
        # Two things break the assumption that av_dist is a smooth function of
        # time: the tracked AV being replaced by a different one, and a jump in
        # the reported distance. Differentiating across either would invent a
        # relative speed that never happened.
        av_switch ← timesteps where av_id changes
        big_jump  ← timesteps where |d(av_dist)/dt| > 70 m/s
        breaks    ← av_switch OR big_jump

        # ---- Differentiate with one-sided stencils at the breaks ----------
        build a forward and a backward index for every timestep:
            normally i+1 and i−1, giving a centred difference
            at a break, pulled in to the same side, giving a one-sided one
        rel_speed ← (av_dist[forward] − av_dist[backward])
                    ÷ (t[forward] − t[backward])

        # ---- Smooth ------------------------------------------------------
        REPEAT 40 times:
            rel_speed ← mean of (rel_speed[backward], rel_speed, rel_speed[forward])
            # The same stencils, so smoothing also never crosses a break

        accumulate av_dist and rel_speed for this segment

    # ---- Keep the usable distance band -----------------------------------
    KEEP samples where lower_Bnd < av_dist < upper_Bnd
        # Below 30 m the distance is dominated by the AV's own length and the
        # measurement is unreliable; beyond 350 m the nearest AV is barely
        # influencing the vehicle.
    accumulate them into the day's pool

    record per-segment statistics (mean, median, standard deviation,
    quartiles) — currently unused downstream

SAVE, through a temporary name:
    filtered_dist_all_files   the pooled distances (m)
    filtered_speed_all_files  the pooled relative speeds (m/s)
    lower_Bnd, upper_Bnd      the bounds used
    j_start, j_end            the segment range pooled
    day                       16, 17 or 18
    # The bounds and range travel with the data, because two days' figures are
    # only comparable if they were pooled the same way.

REPORT an Anderson–Darling test on the pooled speeds
    # h = 1 means the distribution is not Gaussian. Informational only;
    # nothing downstream reads it.
```

## The sign convention

`rel_speed` is the rate of change of the **distance to the AV ahead**:

* **positive** — the gap is opening; the vehicle is falling behind the AV
* **negative** — the gap is closing; the vehicle is catching the AV

Day 17's pooled mean is about −0.35 m/s, so traffic in the band was on average
slowly closing on the AV ahead.

## Notes for review

* **The break handling is the substance of this stage.** Without it, every AV
  handover would produce a large spurious relative speed, and those artefacts
  would dominate the tails of the histogram.

* **40 smoothing passes is a lot.** Each is a 3-point average, so the effective
  window is wide; the intent is a slowly varying closing rate rather than
  per-sample noise. The count is not derived from anything and is a reasonable
  thing to question.

* **Only the downstream engaged AV is considered.** Upstream AVs, and
  unengaged ones, are ignored entirely — the question being asked is how
  traffic behaves *behind* an active controller.

* **The saved file is ~0.7 GB per day** because it keeps every surviving
  sample (about 14 million for day 17) rather than a histogram. That is
  deliberate: it lets the bin edges be changed without re-decoding the day.
