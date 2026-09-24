# `generate_macroscopic_fields.m` — stage `fields`

Turns a day of individual trajectories into continuous traffic-state fields on
a time × space grid: density, flow, fuel rate, and the bulk quantities derived
from them.

**Reads** `results/slim/2022-11-DD/*.json`
**Writes** `results/analysis/2022-11-DD/fields_motion_2022-11-DD.mat` (~15 MB)
**Runs** once per day. Not sharded.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `hx` | 100 m | averaging window half-width in space |
| `v_char` | 20 m/s | characteristic speed; sets `ht = hx / v_char` = 5 s |
| `t_res`, `x_res` | 5 s, 50 m | grid resolution |
| `ax_t` | 06:00–10:00 | time extent of the grid |
| `ax_x` | 0–6500 m | space extent of the grid |
| `direction` | −1 | westbound |
| `lane` | 0 | all lanes |
| `skip_t_compute` | 5 | use every 5th trajectory sample |
| `kernel` | `'box'` | or `'Gaussian'` |
| `created_fields` | Rho, Q, F, U, Phi, Psi | which fields to build |

The grid is 2881 × 131 at these settings.

## The idea

Each trajectory sample is a point measurement — one vehicle, one instant. A
field value at a grid point is the average over a window around it, built by
letting every sample *deposit* its contribution into the grid cells it
overlaps, weighted by how long that sample lasted. Accumulate over every
trajectory, normalise once at the end, and the result is a smooth field.

Density and flow come out of the same accumulation: density weights each
deposit by the sample's duration, flow weights it by duration × speed.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

SET UP THE KERNEL
    ht ← hx / v_char                              # 5 s from 100 m at 20 m/s
    IF box kernel:
        fac ← 1 / (4·ht·hx)                       # uniform weight, normalised by area
    IF Gaussian kernel:
        widths ← ht/√3, hx/√3                     # same standard deviation as the box
        extend the box to 3 standard deviations
        fac ← 1 / (2π·wt·wx)

BUILD THE GRID
    t ← 06:00 to 10:00 every t_res (5 s)
    x ← 0 to 6500 m every x_res (50 m)
    allocate a zero array per field over that grid

ACCUMULATE
FOR each of the day's 24 segments:
    decode it
    trajectories ← those travelling in `direction` and in `lane` (0 = any)

    FOR each such trajectory:
        speed ← central difference of x over t, signed by direction
        sub-sample t, x, speed and fuel rate every skip_t_compute (5th) point
        dt ← the time each retained sample now represents

        FOR each retained sample:
            # Which grid cells this sample can reach
            box ← grid points within ht in time AND hx in space of the sample
            IF box kernel:  G ← 1 over that box
            IF Gaussian:    G ← the Gaussian weight at each grid point in it

            Rho(box) ← Rho(box) + dt · G                  # presence
            Q(box)   ← Q(box)   + dt · speed · G          # flux
            F(box)   ← F(box)   + dt · fuel_rate · G      # fuel per time
            F0(box)  ← F0(box)  + dt · flat_road_rate · G # only if 'full' input

NORMALISE
    multiply Rho, Q, F (and F0) by fac
        # Once at the end rather than per deposit: the same result, far fewer
        # multiplications.

DERIVE
    U   ← Q / Rho          bulk velocity          (NaN where Rho < 1e-3)
    Phi ← F / Rho          bulk fuel rate         (NaN where Rho < 1e-3)
    Psi ← F / Q            bulk fuel consumption  (NaN where Q   < 1e-2)
    # The guards matter: dividing by a near-empty cell produces a huge value
    # that would dominate the colour scale. Such cells are marked missing.

SAVE field, t, x, direction and lane
```

## The fields

| Field | Meaning | Unit | Built from |
|---|---|---|---|
| `Rho` | vehicle density | veh/km | accumulated presence |
| `Q` | flow rate | veh/h | accumulated flux |
| `F` | fuel rate density | g/s | accumulated fuel rate |
| `U` | bulk velocity | m/s | `Q / Rho` |
| `Phi` | bulk fuel rate | g/s | `F / Rho` |
| `Psi` | bulk fuel consumption | g/m | `F / Q` |

`F0`, `Phi0` and `Psi0` are the flat-road counterparts; they need the
`fuel_rate_flat_road_grams_per_second` field, which only `full` carries, and are
not in the default `created_fields`.

The `factor` recorded alongside each field converts it to its display unit —
1000 for `Rho` (per metre to per km), 3600 for `Q` (per second to per hour).

## Notes for review

* **The window is much wider than the grid spacing** — 100 m against 50 m,
  5 s against 5 s — so neighbouring grid points share samples. That is what
  makes the field smooth, and it means the effective resolution is the window,
  not the grid.

* **`ht = hx / v_char` ties the two window dimensions together** through a
  characteristic speed, so that the averaging region roughly follows how far
  traffic moves in the time it spans. Changing `hx` changes both.

* **`skip_t_compute = 5` uses every fifth sample**, cutting the work fivefold.
  Because each retained sample's `dt` is scaled accordingly, the totals are
  unaffected; only the sampling of the deposit pattern is coarser.

* **The NaN guards are what make the figures readable.** Without them, cells at
  the edge of the data — a single vehicle passing an otherwise empty stretch —
  produce enormous bulk values and flatten the colour scale for everything else.

* **The plotting code at the end of this file is disabled**
  (`flag_plot_field = 0`). `plot_macroscopic_fields` is the stage that draws
  these; the block here is for interactive checking.
