# `plotting_LC_analysis.m` — stage `lcplot`

Turns a day of lane-change events into the paper's lane-change figures: how
long traffic was exposed to an engaged AV at each distance, how often it merged
out at that distance, and the cumulative excess of merges over the far-field
baseline.

**Reads** `results/analysis/2022-11-DD/LC_data_DD.mat`, `results/gps/*.json`
**Writes** four PNGs in `results/figures/2022-11-DD/`:
`fig_lc_exposure_*`, `fig_lc_rate_merge_out_*`, `fig_lc_cumulative_excess_*`,
`fig_lc_cumulative_excess_combined_*`
**Runs** once per day. ~20 minutes per day; the exposure loop dominates.

This is the most involved of the plotting stages, and the one whose method most
needs stating.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `lanesToAnalyze` | `[2 3 4]` | one panel per lane; lane 1 (HOV) excluded |
| `binEdges` | −500 : 10 : 500 m | distance bins relative to the AV |
| `max_T` | 11:00 local | AV runs starting later are dropped |
| `grid_res` | 1 m | resolution of the spatial grid used for exposure |
| basal window | 350–500 m | where the far-field merge rate is measured |
| masked band | 4 bins either side of 0 | roughly ±35 m, excluded as unreliable |
| `twoColWidth/Height` | 22.5 × 11.2 in | figure canvas, set explicitly |

Only the **engaged-AV** case is plotted. The all-AV series is still computed and
kept available, but no figure uses it.

## Why normalise at all

A raw count of lane changes at a given distance from an AV says almost nothing,
because traffic does not spend equal time at every distance. Far from an AV
there is much more road and many more vehicle-minutes, so more merges will
happen there whatever the AV is doing. The stage therefore divides counts by
*exposure* — the vehicle-minutes actually spent at that distance — to get a
rate that can be compared across distances.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)
    STOP with an explanation if LC_data_DD.mat or the GPS file is absent

LOAD all_lane_changes_start and all_lane_changes_end
LOAD the day's GPS runs; keep westbound runs starting before 11:00
dt ← the GPS timestep (the mode of its differences), in minutes

CREATE four off-screen figures at an explicit canvas size:
    exposure, merge-out rate, cumulative excess, combined cumulative excess

FOR each lane in {2, 3, 4}:

    merges_in  ← merge-in rows in this lane
    merges_out ← merge-out rows in this lane
    SKIP the lane if it has no merges at all

    # ---- EXPOSURE ----------------------------------------------------------
    # How many vehicle-minutes were spent at each distance from an AV, in this
    # lane, over the day.
    span ← from the earliest to the latest x at which a merge happened here
    grid ← span laid out every grid_res (1 m)
    points_per_bin ← bin width ÷ grid_res          # 10 grid points per 10 m bin

    lane_avs ← GPS runs whose MODAL assigned_lane is this lane
    pool every such run's (timestamp, x, control_car) into one series

    FOR each distinct timestamp in that pool:
        x_all ← positions of every control vehicle at that instant
        x_eng ← those of them whose controller was engaged

        FOR each of (x_all, x_eng):
            # Signed distance from every grid point to the nearest AV ahead,
            # and to the nearest behind
            distances ← grid − each AV position
            nearest_ahead  ← smallest positive distance at each grid point
            nearest_behind ← largest negative distance at each grid point
            counts ← histogram of those distances into binEdges
            exposure(bin) ← exposure(bin) + counts ÷ points_per_bin × dt
            # ÷ points_per_bin converts "grid points" into "bin widths"
            # × dt weights by how long this instant lasted

    treat zero exposure as missing, not as zero

    # ---- RATE --------------------------------------------------------------
    dist_out   ← distance to the nearest engaged AV for each merge-out here
    change_out ← its signed direction
    rate_right ← histogram(dist_out where change_out > 0) ÷ exposure
    rate_left  ← histogram(dist_out where change_out < 0) ÷ exposure
    rate_total ← rate_right + rate_left            # merges per minute

    basal_rate ← mean of rate_total over bins from 350 to 500 m
        # The far field: far enough that the AV is not influencing behaviour.
        # Drawn as a dashed reference line on the figure.

    # ---- THE MASKED BAND ---------------------------------------------------
    # Within about 35 m of the AV, both the lane-change detection and the
    # distance measurement are least trustworthy, so those bins are drawn
    # faded and excluded from the fit below.
    faded ← the 4 bins nearest zero on each side
    solid ← everything else

    # ---- EXCESS AND ITS EXTRAPOLATION --------------------------------------
    restrict to 30–500 m on the downstream side
    e(x) ← rate_total(x) − basal_rate            # excess merges per minute
    treat missing values as zero

    # The masked band leaves no data below ~35 m, so the first three bins are
    # INFERRED rather than measured:
    fit a straight line to the first 10 reliable points (about 35–125 m)
    evaluate it at x = 5, 15, 25 m and prepend those three values to e(x)

    g(x) ← cumulative sum of e(x) over bins
        # "How many excess merges per minute have accumulated by distance x."

    PLOT g(x) solid out to 350 m and faded beyond it
    PLOT the same curve onto the combined figure, one colour per lane

after the lane loop:
    finish the combined figure: 350 m marker, labels, per-lane legend
    resize the combined axes to match one panel of the tiled figures
        # so the two can be placed side by side at the same scale
SAVE all four figures, then close them
```

## Notes for review

* **The extrapolation is a modelling choice, not a measurement.** The three
  values at 5, 15 and 25 m are inferred from the trend further out, and because
  `g(x)` is a cumulative sum they enter **every** later value of the curve. If
  the merge rate very close to an AV departs from the linear trend, the whole
  curve shifts. This is the assumption a reviewer is most likely to probe.

* **Exposure is the expensive part.** The inner loop runs once per distinct GPS
  timestamp — over 128 000 for one lane on one day — and each iteration
  histograms distances from every grid point to every AV. That is essentially
  all of the ~20 minutes per day.

* **"Modal lane" assigns each AV run to one lane for the whole run.** A control
  vehicle that changed lanes contributes its entire run to whichever lane it
  spent most time in, rather than being split. This matters if any run was
  genuinely split between lanes.

* **Lane 1 is excluded.** `lanesToAnalyze` is `[2 3 4]`; the HOV lane is not
  analysed.

* **Only the engaged-AV case is plotted**, though both are computed. Switching
  the figures to all AVs is a one-line change to `titleSuffix`/`distFields`.
