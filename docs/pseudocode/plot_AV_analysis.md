# `plot_AV_analysis.m` — stage `av`

Bins every sample by its distance from an engaged control vehicle and plots
fuel consumption against that distance, pooling all three days. This is the
paper's headline result.

**Reads** all three days' `results/analysis/2022-11-DD/samples_for_distance_analysis_DD.mat`
and `fields_motion_2022-11-DD.mat`
**Writes** into `results/figures/` (the shared folder, not a day folder):
`fig_2_fuel_results_*`, `fig_3_fuel_results_*`, `fig_SM2_vehicle_samples_counts_*`,
each as `.png` and `.fig`
**Runs** once, over all three days. Produces Figures 2, SM2 and SM3.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `testDays` | `[16 17 18]` | all three; the argument is accepted but the stage reads all three regardless |
| `MAXDIST` | 350 m | furthest distance binned |
| `AVLOCBUFFER` | 1 m | a narrow bin either side of zero distance |
| `XWINDOW` | 228–7000 m | the stretch of road considered |
| `TWINDOW` | 06:45–09:15 | the time window considered |
| `minSpdRaw` | 0 m/s | minimum speed for a sample to be kept at all |
| `minSpdMean` | 1 m/s | minimum speed for a sample to enter the *mean* |
| `statsToPlotChoices` | effective, mean, median | the three statistics |
| `figRes`, `figScale` | 1600 × 1000 px, 384 DPI | |

## Main routine

```
GIVEN nothing (or a day list, which is accepted but effectively fixed)

STANDARD PREAMBLE (see assemble_data_GPS.md)

BUILD THE BINS
    evenly spaced from −MAXDIST to +MAXDIST, about 10 m apart
    INSERT an extra pair of edges at ±AVLOCBUFFER
        # A narrow bin either side of zero, so samples level with the AV do not
        # bias a neighbouring bin.

FOR each of the three days:
    LOAD samples_for_distance_analysis_DD.mat
    KEEP samples where:
        speed >= minSpdRaw
        AND x is inside XWINDOW
        AND t is inside TWINDOW (06:45–09:15)
    # Fields are cleared as they are read, one at a time: three days of samples
    # will not all fit in memory otherwise.

FOR each day, FOR each distance bin:
    collect the samples whose distance falls in that bin
    COMPUTE THREE STATISTICS:

        effective ← Σ fuel rate ÷ Σ speed
            # Total fuel burned by everything in the bin, divided by total
            # distance covered by everything in the bin. A fleet average:
            # slow vehicles contribute more fuel and less distance, and it
            # weights each sample by how far it actually went.

        mean ← mean of the per-sample fuel consumption,
               over samples with speed >= minSpdMean
            # Per-sample consumption is fuel rate ÷ speed, which diverges as
            # speed approaches zero; the 1 m/s floor is what keeps the mean
            # finite. This is the reason effective and mean differ.

        median ← median of the per-sample fuel consumption
            # No speed floor needed: the median is insensitive to the tail.

    RECORD the number of samples in the bin

PLOT, for each requested statistic:
    fuel consumption against distance to the engaged AV,
    one line per day, in the fixed day colours (blue, red, green)
    with the AV at zero and traffic behind it on the positive side

ALSO PLOT the number of active control vehicles over time, per day,
    reading the GPS runs and the fields file for each day

SAVE each figure as both .png (384 DPI) and .fig
    fig_2_fuel_results_<stats>_<t1>_<t2>
    fig_3_fuel_results_<stats>_<t1>_<t2>
    fig_SM2_vehicle_samples_counts_<stats>_<t1>_<t2>
```

## The three statistics, and why there are three

They answer different questions about the same bin, and they disagree in ways
that matter:

* **effective** — "how much fuel did the traffic in this bin burn per metre it
  travelled?" A ratio of sums, so it is a genuine fleet average and is what a
  fuel-saving claim should be based on.
* **mean** — "what did the average *sample* consume?" A mean of ratios, which
  weights a slow-moving vehicle the same as a fast one despite its covering
  less ground, and which needs a speed floor to stay finite.
* **median** — the same per-sample quantity, robust to the tail, and a check
  that the mean is not being driven by a few extreme samples.

`fig_2` uses *effective*; `fig_3` shows *mean* and *median* alongside it.

## Notes for review

* **This stage reads all three days regardless of what it is asked for.** The
  `testDays` argument is validated and used to label, but the loops that load
  samples and fields are written over days 16–18. A single-day download will
  therefore fail here while every other stage succeeds — which is why the
  README says the cross-day figures need the full set.

* **The time and space windows discard a lot.** `TWINDOW` of 06:45–09:15 is
  narrower than the 06:00–10:00 the data covers, and `XWINDOW` starts 228 m in.
  Both are choices about where the experiment was in steady state, and both
  move the result.

* **`AVLOCBUFFER` creates two very narrow bins** at ±1 m. They hold few samples
  and are noisy; they exist to keep near-zero samples out of the 10 m bins on
  either side rather than to be read directly.

* **Memory is managed by hand.** Each field is cleared as soon as it has been
  filtered, because three days of samples at roughly 89 million per day will
  not otherwise fit.
