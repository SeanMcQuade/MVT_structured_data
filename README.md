
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022.  This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

These scripts provide analysis regarding the energy usage of vehicles during the CIRCLES consortium MegaVanderTest field experiments. The field experiments were conducted in Nashville on highway I-24 during the week of November 14th 2022. 

The repository refers to data recorded from the I-24 MOTION observatory, as well as the CIRCLES partially automated vehicles. Those data must be obtained elsewhere, and will be privately stored or available only to CIRCLES Team members until the paper is released.

# Contents

- [Which route do you want?](#which-route-do-you-want)
- [Before you start](#before-you-start)
- [Route 1 — plot the figures](#route-1--plot-the-figures)
- [Route 2 — rebuild the figure data from the trajectories](#route-2--rebuild-the-figure-data-from-the-trajectories)
- [Route 3 — rebuild the trajectories from the raw data](#route-3--rebuild-the-trajectories-from-the-raw-data)
- [If something looks wrong](#if-something-looks-wrong)
- [Documentation](#documentation)
- [Requirements and demo](#requirements-and-demo)
- [About the data](#about-the-data)
- [Websites](#websites)

# Which route do you want?

Each route is a superset of the one before it. Pick by how far back you want to
start, not by which figures you want.

| | You want to | You download | You run |
|---|---|---|---|
| **Route 1** | plot the figures | the analysis files + GPS, **~5 GB** | `make figures` |
| **Route 2** | rebuild the analysis files from the processed trajectories, then plot | \+ `slim/`, **51 GB** | `make` |
| **Route 3** | rebuild the processed trajectories from the raw recordings, then everything above | \+ `data/`, **59 GB** | `make rebuild` |

Route 1 produces every figure except the microscopic trajectory plots, which
need the trajectories themselves.

# Before you start

The repository must sit **next to** the folders you download:

```
<some folder>/
  MVT_structured_data/   <- this repository
  results/               <- routes 1 and 2 download into here
  data/                  <- route 3 adds this
```

Check the layout before running anything:

```bash
cd MVT_structured_data
./check_data.sh
```

You need MATLAB (R2025a or b; R2022b also works) with the Statistics and
Machine Learning Toolbox. No Parallel Computing Toolbox is used anywhere. A
Python implementation of the released-data stages lives in `python/`; see
[`python/README.md`](python/README.md).

# Route 1 — plot the figures

Download into `results/`:

| What | Into | Size |
|---|---|---|
| the three GPS files | `results/gps/` | 821 MB |
| `fields_motion_2022-11-DD.mat` | `results/analysis/2022-11-DD/` | 46 MB |
| `LC_data_DD.mat` | `results/analysis/2022-11-DD/` | 40 MB |
| `relspeed_data_DD.mat` | `results/analysis/2022-11-DD/` | 690 MB |
| `samples_for_distance_analysis_DD.mat` | `results/analysis/2022-11-DD/` | 3.9 GB |

Then:

```bash
make figures
```

About 10 minutes per day. The figures land in `results/figures/`.

Two notes. The cross-day fuel figures need **all three days** present, so a
single-day download will build everything else and then stop on that one; add
`-k` to let the rest finish. And if you skip the 3.9 GB `samples_*.mat`, the
download drops to ~1.6 GB and you lose only those cross-day figures.

# Route 2 — rebuild the figure data from the trajectories

Add `results/slim/` (51 GB), the processed westbound trajectories, and the
analysis files are regenerated rather than downloaded:

```bash
make
```

Allow a few hours for all three days. This also produces the microscopic
trajectory plots, which route 1 cannot.

**One thing must still be downloaded.** The lane sidecars
(`*_orig_dist_lane.mat`, 3.6 MB for all three days) are derived from the *raw*
recordings, not from `slim/`, because the processed trajectories no longer carry
the lane a vehicle came from. Keep those, or keep the `LC_data_DD.mat` built
from them, or the lane-change figures cannot be made. Everything else on this
route rebuilds from `slim/`.

# Route 3 — rebuild the trajectories from the raw data

Add `data/` (59 GB), the raw I-24 MOTION segments and control-vehicle
recordings, and everything is rebuilt from them:

```bash
make rebuild
```

Hours, not minutes. To spread the heavy stages across processes, sized by
memory rather than cores:

```bash
make rebuild SHARDS=6
```

To build into a fresh tree instead of in place:

```bash
make rebuild RESULTS=/path/to/results_new
```

# If something looks wrong

`make status` says what every stage would do and why, and builds nothing. Start
there.

A stage whose inputs are missing names the file it wanted; check it against the
download table for your route.

If a stage says a script **"is newer than"** a file you downloaded, that is the
one confusing case. Freshness is decided by modification time, and an archive
that preserves timestamps can arrive older than a repository you just cloned.
Nothing is wrong with the data — you downloaded these files rather than building
them, so say so:

```bash
make accept-verified     # checks the published checksums, then marks them current
```

This only saves time; it is not required.

For the full list of stages, targets and options, see
[`docs/MAKE_TARGETS.md`](docs/MAKE_TARGETS.md).

# Documentation

| Document | What it covers |
| --- | --- |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | The data-flow graph and what each stage computes, with the mapping to the paper's figures. Start here. |
| [`docs/DATA_FILES.md`](docs/DATA_FILES.md) | What is in which file: an inventory of every input and output, the columns of each, and the shared vocabularies (vehicle class, direction, lane, coordinate frames). |
| [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) | Every field of every output: units, dtype, rounding, null/empty semantics, and which stage writes it. |
| [`docs/MATLAB_JSON_FORMAT.md`](docs/MATLAB_JSON_FORMAT.md) | How MATLAB's `jsonencode` formats numbers, needed for byte-identical output. |
| [`docs/DATA_CHANGELOG.md`](docs/DATA_CHANGELOG.md) | Versions of the released data set, what changed in each, and how to tell which version a folder holds. |
| [`docs/REPRODUCIBLE_QUADRATURE.md`](docs/REPRODUCIBLE_QUADRATURE.md) | Why the fuel quadrature was made platform-independent, and the size of the difference. |
| [`docs/MAKE_TARGETS.md`](docs/MAKE_TARGETS.md) | Every stage and every `make` target, the options they take, sharding, and the other ways to drive the pipeline. |
| [`docs/LANE_CHANGE_STAGES.md`](docs/LANE_CHANGE_STAGES.md) | The lane-change and relative-speed stages: how to run and test them. |
| [`docs/DATA_RELEASE.md`](docs/DATA_RELEASE.md) | What to publish and how to lay it out, across the raw, derived and analysis layers. |
| [`docs/PYTHON_PORT.md`](docs/PYTHON_PORT.md) | Status of the Python implementation, what is and is not ported, and how its parity is verified. |
| [`python/README.md`](python/README.md) | The Python port: install (venv), run (CLI), Docker, and test. |

This code produces data set version **2.2** (`mvt.dataVersion`). Each product
folder carries a `dataset_info.json` sidecar recording that version alongside
the copyright, the licence, and how the folder was produced — MATLAB release or
Python version, platform, host, user and code commit — so a folder can be
identified on its own. See
[`docs/DATA_CHANGELOG.md`](docs/DATA_CHANGELOG.md).

The pipeline runs two ways: the **MATLAB** scripts under `Scripts/` (the
reference implementation, driven by `make` or `run_all_scripts`), and the
**Python** port under `python/` (`mvtpy`, runnable from a venv or a container).
This README covers the MATLAB side; see `python/README.md` for the Python side.


# Requirements and demo

To carry out the analysis, the following software requirements are needed:

#### Software Dependencies:

- Install [MATLAB](https://mathworks.com/)
- Only `Matlab` must be installed as a required toolbox
- MATLAB v2025a and v2025b have both been tested
- Tested on Mac Sequoia 15.6, Windows 11

#### Installation Guide
- Install MATLAB (approximately 20 minutes, depending on download speeds)
- Fetch the data (will be made publicly available upon publication). The extraction should be to a hard drive with at least 700GB of available space, to ensure enough space to generate additional files.

#### Demo

##### Instructions to run the demo:

- Navigate to the folder `MVT_structured_data/Scripts`.
- Run the file `run_all_scripts.m`

The execution of this file will take significant time. It will reproduce the style of plot for the data provided. 

- Generate structured collections of data samples for this day from the `results/` folder
- Generate the macroscopic fields (from the provided data)
- Plot the macroscopic fields (from the provided data) with the GPS data from the cars (across the entire day)
- Plot the microscopic fields of all cars (from above)
- Carry out the analysis of fuel usage and comparisons (from the provided data)
- Plot the fuel usage and comparisons (from above)

A subset of these plots are synthesized in `results/figures/2022-11-17` in the released data set (available upon conclusion of peer review).

##### Expected output: 

Resulting files should be generated in the following folders for day 2022-11-17

###### Found in `results/figures/2022-11-17`
```
fields_motion_2022-11-17.mat
fig_2_fuel_results_effective_645_915.fig
fig_2_fuel_results_effective_645_915.png
fig_3_fuel_results_effective_mean_median_645_915.fig
fig_3_fuel_results_effective_mean_median_645_915.png
fig_field_20221117_west_laneall_motion_F_av_nature_large.png
fig_field_20221117_west_laneall_motion_Phi_av_nature_large.png
fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png
fig_field_20221117_west_laneall_motion_Q_av_nature_large.png
fig_field_20221117_west_laneall_motion_Rho_av_nature_large.png
fig_field_20221117_west_laneall_motion_U_av_nature_large.png
fig_motion_trajectories_20221117_west_laneall_lowres.png
fig_motion_trajectories_20221117_west_laneall_zoom_lowres.png
fig_motion_trajectories_20221117_west_laneall_zoomwin_lowres.png
fig_SM2_vehicle_samples_counts_effective_mean_median_645_915.fig
fig_SM2_vehicle_samples_counts_effective_mean_median_645_915.png
samples_for_distance_analysis_17.mat
```

An example image that should be produced represents the macroscopic flow rates, with overlay  is `results/figures/2022-11-17/fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png` which describes the bulk fuel consumption with overlay of GPS data from our control cars, indicating when their control was active (or not) during their drives.

![Bulk Fuel results (partial) with only a subset of data from 2022-11-17](../results/figures/2022-11-17/fig_field_20221117_west_laneall_motion_Psi_av_nature_large.png)

**Note** there may be minor errors or warnings thrown, since the data pipeline is intended to reproduce exact figures with comparable max/min values and colors across multiple plots. If only one day with a subset of data into that anaysis, is included, plots may have

- what seems to be missing large portions on the left/right of presented data
- what seems to indicate that axes are 'zoomed out'

These are artifacts of axes bounds that are normalized across multiple days for comparison, and thus the appropriate approach is to open those .fig files as plots with MATLAB, and then zoom in to explore.

###### Found in `results/slim/2022-11-17`

This will include a single .mat file that represents binary data in MATLAB format for quick review and analysis in subsequent plots.

```
I-24MOTION_2022-11-17_07-59-59_reduced.mat
```

###### Notes on warning messages

*Note* several output messages will show in the MATLAB window that show warnings for additional legend entries that are not used. This is due to only a subset of data being shown for the demonstration data.

```
Warning: Ignoring extra legend entries. 
> In legend>process_inputs (line 575)
In legend>make_legend (line 294)
In legend (line 245)
In plot_AV_analysis>plot_one_sided (line 473)
In plot_AV_analysis (line 310)
In reproduce_plots (line 27) 
```

##### Expected Runtime

It should take approximately 4-5 minutes or faster to regenerate plots if all intermediate data already are downloaded, and the `slim/` folder exists.



# About the data

### About data from I-24 MOTION

The `results/slim` data include trajectories recorded and processed by the I-24 MOTION observatory. Primarily the data from I-24 MOTION are comprised of position and speed trajectories from all vehicles detected on the roadway by that observatory. The `results/slim` data align the location and state of each of the CIRCLES control cars includes GPS location data as well as data from on-board vehicle informatics recorded by custom hardware installed with a Raspberri Pi. Those data include information regarding to indicate the control and system state of the vehicle. 

### About data from GPS and the CIRCLES Cars
Two kinds of data were collected from the CIRCLES Cars, and found in the `results/gps` folder in the data archive. Each day provides a separate standalone file for the cars that ran on that day. Information includes GPS data, aligned with timeseries information from data collected directly from the car.

#### Data collection from CIRCLES CARS with team-installed GPS Sensors
Information from team-installed GPS sensors was collected at 10-Hz. These data include the position and speed of the vehicle at each sample point.

#### Data Collection from CIRCLES Cars with team-designed on-board data collection
Information from team-installed computers that interface with the Controller Area Network (CAN) were critical to sensing and control of the experiment cars. These data are aligned with the raw GPS information to provide the state of the vehicle at that time (speed, assigned lane of travel, desired cruise control set point, etc.). 


# Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)

