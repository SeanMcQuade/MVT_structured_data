# `plot_macroscopic_fields.m` — stage `macro`

Draws the traffic-state fields as time × space heatmaps, optionally with the
control vehicles' own traces laid over them.

**Reads** `results/analysis/2022-11-DD/fields_motion_2022-11-DD.mat`,
`results/gps/*.json`
**Writes** six PNGs per day in `results/figures/2022-11-DD/`,
`fig_field_<yyyyMMdd>_<direction>_<lane>_motion_<field>[_av]_nature_large.png`
**Runs** once per day. Produces the paper's Figure 3 and SM5.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `plotFields` | Rho, Q, F, U, Phi, Psi | one figure each |
| `direction`, `lane` | −1, 0 | westbound, all lanes — must match the fields file |
| `flagPlotAvTrajectories` | 1 | overlay the control-vehicle traces |
| `skipTPlot` | 1 | sub-sampling of those traces (1 = every point) |
| `timeZoomWindow` | 06:10–09:50 | the window actually shown |
| `figRes` | 2500 × 800 px | canvas; DPI 384 at save time |
| `colorbarLimit` | per field | fixed upper limits (see below) |

Fixed colour limits: `Rho` 220 veh/km, `Q` 8000 veh/h, `F` 0.25 g/s, `U` 35 m/s,
`Phi` 4.5 g/s, `Psi` 0.22 g/m. Setting one to 0 falls back to the field's 99th
percentile.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

LOAD fields_motion_2022-11-DD.mat        # field, t, x, direction, lane
IF overlaying AV traces:
    LOAD the day's GPS runs
    KEEP those travelling in `direction`

FOR each field in {Rho, Q, F, U, Phi, Psi}:

    values ← field.value transposed, times its display factor
        # e.g. Rho is stored per metre and displayed per km, factor 1000

    # ---- Upsample for a smooth image --------------------------------------
    build a grid 4× finer in both t and x
    interpolate the field onto it                   # interp2, bilinear
    # Cosmetic only: it removes the blockiness of a 5 s × 50 m grid without
    # adding information.

    DRAW the result as an image, x in km against time

    # ---- Colour scale ------------------------------------------------------
    IF this field has a preset limit:
        clamp the colour axis to [0, limit]
            # Fixed across the three days, so the days are visually comparable.
            # This is the reason the limits are hard-coded rather than derived.
    ELSE:
        use [0, the field's 99th percentile]
            # The percentile, not the maximum: a single extreme cell would
            # otherwise compress the whole scale.

    # ---- AV traces ---------------------------------------------------------
    IF overlaying:
        FOR each control-vehicle run:
            plot its x against its time, sub-sampled by skipTPlot
            colour each point by whether the controller was engaged
        # This is what makes the figure show cause and effect: the waves in the
        # field against the vehicles that were damping them.

    # ---- Framing -----------------------------------------------------------
    label the time axis in local time, with n_xticks (16) ticks
    title with the field's name, unit, direction, lane and the date
    restrict the x-axis to timeZoomWindow (06:10–09:50)

    IF saving:
        name ← fig_field_<yyyyMMdd>_<direction>_<lane>_motion_<field>
        IF the AV traces were overlaid: append _av
        append _nature_large
        set the canvas and axes position explicitly, font size 20
        print at 384 DPI
```

## Notes for review

* **The colour limits are the reason these figures are comparable.** Every one
  is drawn to the same scale across all three days, so a reader can see day 17
  being denser than day 16 rather than each day being auto-scaled to fill its
  own range. Changing a limit changes the appearance of all three days and
  should be done for all of them together.

* **The 4× interpolation is cosmetic.** The underlying resolution is the
  averaging window from the `fields` stage — 100 m and 5 s — not the grid and
  not this upsampling.

* **`direction` and `lane` are declared twice**, here and in
  `generate_macroscopic_fields`, and must agree. The fields file records what it
  was built with, and the title is drawn from this script's values, so a
  mismatch would mislabel rather than fail.

* **The figures are 10000 × 3200 px** at 384 DPI. An earlier version produced
  5000 × 1600; `_nature_large` marks the current size, so a figure from the
  older run is identifiable by its dimensions.
