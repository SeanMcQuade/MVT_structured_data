# `plot_microscopic_trajectories.m` — stage `micro`

Draws every vehicle as a band in the time–space plane, its width the vehicle's
length and its colour the vehicle's speed, so that traffic waves appear as
diagonal stripes.

**Reads** `results/slim/2022-11-DD/*.json` (or `full/` for eastbound), via
reduced `.mat` caches it builds itself
**Writes** three PNGs per day in `results/figures/2022-11-DD/`:
`fig_motion_trajectories_<yyyyMMdd>_<direction>_<lane>_lowres.png`, plus
`_zoomwin_lowres` and `_zoom_lowres`
**Runs** once per day, about 20 minutes. Memory-hungry.

## Tunables

| Name | Value | Meaning |
|---|---|---|
| `direction`, `lane` | −1, 0 | westbound, all lanes |
| `ax_t`, `ax_x` | 06:00–10:00, 0–6500 m | extent drawn |
| `skip_t_plot` | 50 | sub-sampling along each trajectory |
| `timeZoomWin` | 06:12–06:19 | the zoom window |
| `xZoomWin` | 845–1495 m | its spatial extent |
| `lowResDPI` | 192 | resolution of the saved figures |
| `highResDPI` | 1536 | only if `flag_save_highres` is on (it is off) |
| `flag_reduce_data_files` | 1 | build and use the `.mat` caches |
| `flag_use_int32_vars` | 1 | store patch coordinates as scaled integers |
| `toInt32Factor` | 1000 | that scaling |

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

CHOOSE THE INPUT
    westbound → results/slim/2022-11-DD
    eastbound → results/full/2022-11-DD     # slim has no eastbound data

BUILD OR REUSE THE REDUCED CACHES
    IF fewer than 24 *_reduced.mat exist in results/.mvt/cache/2022-11-DD:
        FOR each released segment:
            decode it
            DROP the ~26 fields this plot does not use — fuel totals, reference
                 trajectories, AV distances, corrected y, dimensions beyond
                 length
            SAVE as <segment>_reduced.mat
        # No information is added: this is a format conversion that skips the
        # JSON decode on later runs. ~4.6 GB per day against 17 GB of JSON.

OPEN one figure; set up the colour scale and colour bar for speed

FOR each of the 24 cached segments:
    load it
    trajectories ← those in `direction` and `lane` (0 = any)
    SORT them by vehicle length, longest first
        # So that a lorry cannot hide a car drawn underneath it: shorter bands
        # are painted last, on top.

    FOR each trajectory:
        speed ← central difference of x over t, signed by direction
        sub-sample to about length/skip_t_plot points, at least 2

        # ---- One closed band per vehicle ----------------------------------
        # Trace the front of the vehicle forwards in time, then the back
        # backwards, giving a closed polygon whose width is the vehicle length.
        patch_t ← [t forwards ; t backwards]
        patch_x ← [x forwards ; x backwards offset by the vehicle length]
        patch_v ← [speed forwards ; speed backwards]

        IF using int32:
            scale by toInt32Factor and cast
            # Halves the memory of the vertex arrays, which is what makes a
            # day's worth of patches fit at all.

        collect the vertices, the colours, and the vertex count

    # ---- Draw the whole segment as ONE patch object ---------------------
    concatenate every trajectory's vertices into one array
    build a face index matrix, one row per trajectory, NaN-padded to the
        longest trajectory
    DRAW a single patch with those faces, coloured by speed
    # One graphics object per segment rather than per trajectory. With 185 000
    # trajectories in a day, one object each exhausted memory and took the
    # machine down; this is the fix.

LABEL the axes in local time, restrict to 06:10–09:50
SAVE the full figure at lowResDPI

IF zooming:
    DRAW the zoom window as a box on the full figure
    SAVE that as *_zoomwin_lowres.png
    THEN set the axis limits to the zoom window, hide the axes, and make the
        plot box fill the figure on a black background
    SAVE that as *_zoom_lowres.png
    # This second file deliberately has no axes, labels or colour bar: it is
    # an INSET, meant to be placed inside the _zoomwin figure. On its own it
    # looks unfinished, and is supposed to.
```

## Notes for review

* **Two figures, one of which is not standalone.** `_zoomwin` shows the whole
  day with a box marking the zoom; `_zoom` is the contents of that box, with no
  framing, to be dropped into it. Judging `_zoom` on its own is a
  misunderstanding that has already happened once.

* **The single-patch-per-segment construction is a memory fix, not a
  style choice.** Drawing one patch object per trajectory — 185 000 of them —
  exhausted RAM and rebooted the machine. See
  `docs/MEMORY_FIX_microscopic_trajectories.md`.

* **The caches are derived and disposable**, under `results/.mvt/cache/`.
  Deleting them costs one rebuild. They are also the reason this stage can run
  without re-reading the JSON, which is most of its speed on later runs.

* **`flag_save_highres` is off.** At 1536 DPI the render is enormous and has
  historically been able to hang the machine; the low-resolution figures are
  what the paper uses.

* **`skip_t_plot = 50` is aggressive.** Each trajectory is drawn with about a
  fiftieth of its points. At the full-day scale this is invisible; in the zoom
  it is the reason the bands have visible corners.
