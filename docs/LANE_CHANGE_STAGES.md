# Lane-change and relative-speed stages

Four analyses that used to be run by hand are now pipeline stages, with the same
staleness checking, sharding and `make` targets as the rest of the pipeline.

| stage | script | reads | writes |
|---|---|---|---|
| `lanes` | `generate_orig_dist_lanes.m` | raw I-24 MOTION segments | 24 × `I-24MOTION_<ts>_orig_dist_lane.mat` per day |
| `lc` | `extract_lane_changes_v_dist_to_av.m` | slim JSON + the sidecars | `LC_data_DD.mat` |
| `lcplot` | `plotting_LC_analysis.m` | `LC_data_DD.mat` + GPS | 4 × `fig_lc_*.png` per day |
| `relspeed` | `relative_speed_histogram.m`, `binned_relative_speed.m` | slim JSON | 2 × `.pdf` per day |

All outputs land under `results/figures/`: the per-day products in
`results/figures/2022-11-DD/`, and the two relative-speed PDFs at the
`results/figures/` root, under the names they have always had.

`lanes` is the only one of the four that reads the raw data. Its outputs are
small (3.6 MB for all three days), so shipping them alongside a results bundle
spares anyone reproducing the lane-change figures a 55 GB raw download.

## Running them

```bash
make lanes-17                  # 24 sidecars for one day
make SHARDS=6 lanes-17         # the same, six MATLAB processes over 24 segments
make lc-17                     # extract the day's lane-change events
make lcplot-17                 # the four lane-change figures
make relspeed-17               # the two relative-speed figures
make lanes lc lcplot relspeed  # all four, all days in DAYS
make status                    # what is stale, and why
```

`lanes` and `lc` are part of `make data`; `lcplot` and `relspeed` are part of
`make figures`, and therefore of `make all`. From inside MATLAB, use
`mvt.build('lcplot', 17)` or let `run_all_scripts` walk the whole graph.

Only `lanes` shards, because it is the only one of the four that works segment
by segment. `lc` deliberately processes a whole day in one process: it decodes
all 24 slim segments and concatenates them.

Rough timings on a 128 GB machine: `lanes` ~25 s per segment (so ~2 min for a
day across 6 shards), `lc` ~5 min per day, `lcplot` ~6.5 min per day,
`relspeed` ~5 min per day.

## Testing this branch

Build into a separate results tree, so nothing published is touched:

```bash
make RESULTS=/path/to/results_test SHARDS=6 lanes-17
make RESULTS=/path/to/results_test lc-17 lcplot-17 relspeed-17
```

To compare against a known-good tree, `.mat` files cannot be byte-compared —
v7 files are gzip streams carrying a creation timestamp, so two saves of
identical data differ on disk. Use the content hash instead:

```matlab
mvt.matHash('/path/to/results_test/figures/2022-11-17/LC_data_17.mat')
mvt.matHash('/path/to/results/figures/2022-11-17/LC_data_17.mat')
```

Equal hashes mean equal contents. `mvt.matCompare` explains a difference when
they are not equal. PNG and PDF output is not checksummed: renderers do not
agree byte for byte.

## Things worth knowing

**The sidecars pair with the slim JSON by index.** `lc` sets
`origin_lane`/`destination_lane` on slim trajectory *i* from sidecar entry *i*.
That only holds while both were built from the same raw segment with the same
clipping, so `lc` errors out if the counts disagree rather than silently
mispairing lanes with trajectories. If it does, rebuild `lanes` for that day.

**The lane-assignment code is forked.** `generate_orig_dist_lanes` carries its
own copy of `assign_lanes` and `clip_lane_changes`, extended to track the origin
and destination lane through the clipping. The originals live in
`generate_data_mvt_slim`. The two must stay in step, and `mvt.sources` cannot
see the coupling, so editing one will not invalidate the other's outputs.
Factoring them into one shared file is the obvious follow-up.

**Figure canvases are pinned.** `relative_speed_histogram` and
`binned_relative_speed` set an explicit `figRes`, because MATLAB otherwise
derives the default figure size from the screen. `exportgraphics` crops to the
drawn content, so the same code produced a differently sized PDF under
`matlab -batch` than in an interactive session. With the canvas pinned, batch
runs reproduce the published page geometry exactly.

**`plotting_LC_analysis` keeps CRLF line endings**, matching the copy it was
imported from, so that edits made elsewhere merge cleanly.
