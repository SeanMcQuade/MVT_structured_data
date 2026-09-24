# Pseudocode descriptions of the pipeline stages

One file per script in the `run_all_scripts` loop, in the order that loop runs
them. Each describes what the script does in structured English, at a level
between the prose in `ALGORITHMS.md` and the MATLAB itself: enough detail to
follow the logic and check it, without reading the code.

| Order | Stage | Script | Pseudocode |
|---|---|---|---|
| 1 | `gps` | `assemble_data_GPS.m` | [assemble_data_GPS.md](assemble_data_GPS.md) |
| 2 | `full` | `generate_data_mvt_full.m` | [generate_data_mvt_full.md](generate_data_mvt_full.md) |
| 3 | `slim` | `generate_data_mvt_slim.m` | [generate_data_mvt_slim.md](generate_data_mvt_slim.md) |
| 4 | `lanes` | `generate_orig_dist_lanes.m` | [generate_orig_dist_lanes.md](generate_orig_dist_lanes.md) |
| 5 | `lc` | `extract_lane_changes_v_dist_to_av.m` | [extract_lane_changes_v_dist_to_av.md](extract_lane_changes_v_dist_to_av.md) |
| 6 | `lcplot` | `plotting_LC_analysis.m` | [plotting_LC_analysis.md](plotting_LC_analysis.md) |
| 7 | `relspeed` | `relative_speed_histogram.m` | [relative_speed_histogram.md](relative_speed_histogram.md) |
| 8 | `relspeedplot` | `plot_relative_speed.m`, `binned_relative_speed.m` | [plot_relative_speed.md](plot_relative_speed.md) |
| 9 | `samples` | `generate_data_samples.m` | [generate_data_samples.md](generate_data_samples.md) |
| 10 | `fields` | `generate_macroscopic_fields.m` | [generate_macroscopic_fields.md](generate_macroscopic_fields.md) |
| 11 | `macro` | `plot_macroscopic_fields.m` | [plot_macroscopic_fields.md](plot_macroscopic_fields.md) |
| 12 | `micro` | `plot_microscopic_trajectories.m` | [plot_microscopic_trajectories.md](plot_microscopic_trajectories.md) |
| 13 | `av` | `plot_AV_analysis.m` | [plot_AV_analysis.md](plot_AV_analysis.md) |

Stages 1–12 run once per day; stage 13 runs once, over all three days.

## Conventions used here

* **Tunable** marks a constant declared at the top of the script and intended
  to be changed. Its value is the committed one.
* Every stage begins with the same harness preamble — validate the day, read
  options, resolve paths, check staleness, honour `Force`/`Clean`/`DryRun`.
  It is written out in full for stage 1 and referred to as **the standard
  preamble** thereafter.
* Loop variables and field names match the MATLAB, so a line here can be found
  in the source.
* Units are named wherever a value carries one. The pipeline mixes feet (as
  I-24 MOTION records them) with metres (as the released data reports them),
  and the conversion points are called out.

## Related documents

* [`ALGORITHMS.md`](../ALGORITHMS.md) — the same material as prose, with the
  data-flow graph and the mapping to the paper's figures.
* [`DATA_FILES.md`](../DATA_FILES.md) — what is in which file.
* [`DATA_DICTIONARY.md`](../DATA_DICTIONARY.md) — every field, with units.
