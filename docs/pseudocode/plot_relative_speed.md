# `plot_relative_speed.m` and `binned_relative_speed.m` — stage `relspeedplot`

Draws the two relative-speed figures from the samples
[`relative_speed_histogram`](relative_speed_histogram.md) pooled.

**Reads** `results/analysis/2022-11-DD/relspeed_data_DD.mat`
**Writes** two PDFs in `results/figures/2022-11-DD/`:
`fig_relspeed_histogram_<yyyyMMdd>.pdf`, `fig_relspeed_behind_av_<yyyyMMdd>.pdf`
**Runs** once per day, about 10 seconds.

Splitting this from the pooling is what lets the figures be redrawn without the
51 GB of trajectories.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `edges` | −30 : 0.5 : 30 m/s | histogram bins |
| `figRes` | 700 × 420 px | figure canvas, pinned (see the note below) |
| `xlim` | ±15 m/s | the plotted range of the histogram |
| `MAXDIST` | 350 m | furthest distance bin |
| `AVLOCBUFFER` | 1 m | a narrow bin placed either side of zero distance |
| `fontSize` | 12 pt | |

## `plot_relative_speed` — the histogram

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)
    STOP with an explanation if relspeed_data_DD.mat is absent

LOAD filtered_dist_all_files, filtered_speed_all_files, j_start, j_end

COMPUTE over the pooled speeds:
    mean, median, standard deviation, first and third quartiles

DRAW, on an off-screen figure of pinned size:
    a histogram of the pooled speeds over `edges`
    a vertical line at the mean, full height
    a vertical line at the median, half height
    a horizontal bar spanning mean ± one standard deviation
    a horizontal bar spanning the interquartile range
    a legend quoting the mean and median to three decimals
    title naming the day; x from −15 to +15 m/s
SAVE as vector PDF

HAND the pooled samples to binned_relative_speed for the second figure
```

## `binned_relative_speed` — mean and median by distance

```
GIVEN the pooled distances and speeds, the day, and an output path

BUILD the distance bins:
    evenly spaced from −MAXDIST to +MAXDIST, about 10 m apart
    then INSERT an extra pair of edges at −AVLOCBUFFER and +AVLOCBUFFER
    # A narrow bin either side of zero, so samples taken essentially level with
    # the AV do not fall into a neighbouring bin and bias it.
bin_centres ← the midpoints

FOR each bin:
    collect the pooled speeds whose distance falls in it
    mean_binned(bin)   ← their mean, ignoring NaN
    median_binned(bin) ← their median

DRAW, on an off-screen figure of pinned size:
    mean against bin centre, solid
    median against bin centre, dashed
    a dotted horizontal line at zero relative speed
    a dotted vertical line at zero distance
    the text "behind AV" near the top left
    title naming the day
    x from 0 to the largest bin centre        # only the behind-AV half is shown
SAVE as vector PDF, then close
```

## Notes for review

* **The figure canvas is pinned, and must stay pinned.** MATLAB otherwise takes
  the default figure size from the screen, and `exportgraphics` crops the PDF
  to the drawn content — so the same code produced a different page size under
  `matlab -batch` than in an interactive session. With `figRes` set explicitly,
  batch runs reproduce the published page geometry exactly (444 × 327 pt for
  the histogram, 450 × 331–332 pt for the binned figure, the 1 pt difference
  between days being genuine, from data-driven y-limits).

* **Only the behind-AV half is plotted**, though both halves are computed: the
  x-limit starts at zero. The ahead-of-AV data is in the file and in the bins.

* **The "ahead of AV" label is commented out** in the source, left from when
  both halves were shown.

* **`j_start`/`j_end` are loaded but unused.** They are provenance — which
  segments were pooled — and are quoted only by title variants that are
  currently commented out.
