
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022.  This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

These scripts provide analysis regarding the energy usage of vehicles during the CIRCLES consortium MegaVanderTest field experiments. The field experiments were conducted in Nashville on highway I-24 during the week of November 14th 2022. 

The repository refers to data recorded from the I-24 MOTION observatory, as well as the CIRCLES partially automated vehicles. Those data must be obtained elsewhere, and will be privately stored or available only to CIRCLES Team members until the paper is released.

# Contents
- [Quick start](#quick-start)
- [Running the pipeline](#running-the-pipeline)
- [Documentation](#documentation)
- [Scripts to generate the integrated data set](#scripts-to-generate-the-integrated-data-set)
- [Scripts to plot figures from the article](#scripts-to-plot-figures-from-the-article)
- [Websites](#websites)

# Quick start

## 1. Put the data in place

The pipeline expects this repository to sit **next to** the `data/` and
`results/` folders:

```
<some folder>/
  MVT_structured_data/   <- this repository
  data/                  <- raw inputs (cars/, i24motion/)   — for bootstrapping
  results/               <- processed inputs and outputs (slim/, gps/, figures/)
```

The data are distributed separately, in three routes of increasing size:

| route | download | gets you |
|---|---|---|
| 1 | **~5 GB** — GPS + the derived `.mat` analysis files | every figure but the trajectory plots, with no trajectory data at all |
| 2 | **+51 GB** — the slim trajectories | the trajectory plots, and regenerates the `.mat` rather than downloading it |
| 3 | **+59 GB** — the raw data | rebuilds `gps`, `slim` and everything downstream |

Pick one in [What to download](#what-to-download) before going further. Route 1
is ~1.6 GB if you skip the cross-day fuel figures.

Once you have them, check the layout is correct:

```bash
./check_data.sh          # reports what's present and which workflows can run
```

It tells you whether you can **plot/analyze** (needs `results/slim` +
`results/gps`) and/or **bootstrap from raw data** (needs `data/cars` +
`data/i24motion`). Note that it does not yet know about the smallest download
route: `make figures-from-mat` needs only `results/gps` and the `.mat`
intermediates, so `check_data.sh` will report that plotting is unavailable even
when that route will work. If the folders are elsewhere, point it at them:

```bash
MVT_DATA_DIR=/abs/data MVT_RESULTS_DIR=/abs/results ./check_data.sh
```

## 2. Run it
These examples use 3 concurrent processes. If your machine has a lot of ram, you can bump this up to 8 or 10.

### MATLAB:
```matlab
cd MVT_structured_data/Scripts
make all Workers 3
```

### MATLAB (headless) from command prompt:
```bash
cd MVT_structured_data
make -j3 all
```

### Python from command prompt:
```bash
cd MVT_structured_data/python
./setup_venv.sh --full
source .venv/bin/activate
mvt build -j 3
```

#### MATLAB: more details

```matlab
cd MVT_structured_data/Scripts
make all Workers 3
```

That is the whole thing. It builds every stage, for all three days, in the right
order, with the heavy stage spread across 3 processes. Work that is
already up to date is skipped, so an interrupted run picks up where it left off.

Everything else — choosing the worker count, building one day or one stage,
writing results elsewhere, the Unix `make`, and the Python port — is under
[Running the pipeline](#running-the-pipeline).

## Documentation

| Document | What it covers |
| --- | --- |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | The data-flow graph and what each stage computes, with the mapping to the paper's figures. Start here. |
| [`docs/DATA_FILES.md`](docs/DATA_FILES.md) | What is in which file: an inventory of every input and output, the columns of each, and the shared vocabularies (vehicle class, direction, lane, coordinate frames). |
| [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) | Every field of every output: units, dtype, rounding, null/empty semantics, and which stage writes it. |
| [`docs/MATLAB_JSON_FORMAT.md`](docs/MATLAB_JSON_FORMAT.md) | How MATLAB's `jsonencode` formats numbers, needed for byte-identical output. |
| [`docs/DATA_CHANGELOG.md`](docs/DATA_CHANGELOG.md) | Versions of the released data set, what changed in each, and how to tell which version a folder holds. |
| [`docs/REPRODUCIBLE_QUADRATURE.md`](docs/REPRODUCIBLE_QUADRATURE.md) | Why the fuel quadrature was made platform-independent, and the size of the difference. |
| [`docs/PYTHON_PORT.md`](docs/PYTHON_PORT.md) | Status of the Python implementation and how its parity is verified. |
| [`python/README.md`](python/README.md) | The Python port: install (venv), run (CLI), Docker, and test. |

This code produces data set version **2.1.1** (`mvt.dataVersion`). Each product
folder carries a `dataset_info.json` sidecar recording that version alongside
the copyright, the licence, and how the folder was produced — MATLAB release or
Python version, platform, host, user and code commit — so a folder can be
identified on its own. See
[`docs/DATA_CHANGELOG.md`](docs/DATA_CHANGELOG.md).

The pipeline runs two ways: the **MATLAB** scripts under `Scripts/` (the
reference implementation, driven by `make` or `run_all_scripts`), and the
**Python** port under `python/` (`mvtpy`, runnable from a venv or a container).
This README covers the MATLAB side; see `python/README.md` for the Python side.

## Running the pipeline

The pipeline is incremental: every stage checks whether its outputs are
missing, older than their inputs, or older than the code that produced them,
and does nothing otherwise. Editing a script therefore causes exactly the
affected outputs to rebuild — the older behavior ("output file exists, skip")
meant a code change silently produced nothing until you deleted files by hand.

### Stages

| stage | what it makes | reads |
|---|---|---|
| `gps` | `results/gps/CIRCLES_GPS_10Hz_*.json` | raw car CSVs + MOTION segments |
| `slim` | released westbound trajectories | raw MOTION segments + gps |
| `full` | adds eastbound and reference trajectories | the same (opt-in; not built by `all`) |
| `lanes` | 24 origin/destination-lane sidecars per day | raw MOTION segments |
| `lc` | `LC_data_DD.mat`, the day's lane-change events | slim + the sidecars |
| `samples` | `samples_for_distance_analysis_DD.mat` | slim |
| `fields` | `fields_motion_2022-11-DD.mat` | slim |
| `macro` | macroscopic field figures | `fields_*.mat` + gps |
| `lcplot` | four lane-change figures per day | `LC_data_DD.mat` + gps |
| `av` | cross-day fuel and sample-count figures | `samples_*.mat`, all days |
| `micro` | trajectory figures | reduced plotting caches, derived from slim |
| `relspeed` | two relative-speed figures per day | slim |

`lanes`, `slim` and `full` shard by segment; the rest ignore `Workers`.

### Targets

```matlab
make                          % everything out of date, all three days
make all Workers 6            % the same, with slim and lanes across 6 processes
make data                     % gps, slim, lanes, lc, samples, fields
make figures                  % macro, lcplot, av, micro, relspeed
make figures-from-mat         % only the figures that need no slim tree
make figures-from-slim        % the rest: micro, relspeed
make slim                     % one stage
make slim Days 18             % one stage, one day
make status                   % what is stale, and why; builds nothing
make config                   % resolved paths, days, workers, data version
```

`make` builds `figures-from-mat` before `figures-from-slim`, so a results-only
download produces everything it can before anything reaches for the 51 GB slim
tree. See [What to download](#what-to-download).

Options may follow any target: `Workers`, `Days`, `Force`, `Clean`, `DryRun`,
`Verbose`, `SettleSeconds` (see `Scripts/+mvt/options.m`). Command syntax is
fine — `make slim Days 18 Force true` — as is function syntax,
`make('slim', 'Days', 18)`.

```matlab
make all DryRun true          % plan the whole run, write nothing
make slim Force true          % rebuild regardless of timestamps
```

### If a run stops early

`make` prints the days it is building and any `MVT_*` variables in effect
before it starts:

```
[make] environment: MVT_DAYS = 16
[make] target 'all': gps slim samples fields macro micro av
[make] days 16, 3 worker(s) for sharded stages
```

`setenv` persists for the whole MATLAB session, so an `MVT_DAYS` left over from
an earlier experiment quietly narrows every later build. `make config` shows
what was resolved:

```matlab
make config          % days, paths, workers, data version
setenv('MVT_DAYS','')   % clear an unwanted override
```

Every build writes a timestamped log:

```
[make] logging to <results>/.mvt/logs/make-all-20260728_091706.log
```

It captures the environment, the stages as they run, and — importantly — the
full error report when something fails, which in `-batch` goes to stderr and
would otherwise be missing from the log that exists to explain it. Sharded
stages additionally leave one log per worker in the same folder. `make all Log
false` turns it off.

By default a failing stage stops the run, so one bad day abandons the rest. To
attempt everything and see all the failures together:

```matlab
make all Workers 3 KeepGoing true
```

Each failure is reported as it happens and listed again at the end, and the run
still errors afterwards so it cannot be mistaken for success.

### How many workers

`Workers` is the only setting that turns anything concurrent. Nothing in this
pipeline uses `parfor`, `opts.UseParfor` is declared but never read by any
stage, and no toolbox is required: parallelism comes from running several
MATLAB processes over disjoint shards of a day, which is what
`mvt.runShards` does and what the Makefile does.

**Size it by memory, not cores.** Each worker holds a decoded 10-minute segment
and peaks at several GB, so 6 workers want roughly 30 GB. On a 32 GB machine use
4–6; 12 workers on a 128 GB machine cut a day of `slim` from about 40 minutes to
4. Only `slim` and `full` shard by segment — the other stages ignore `Workers`
rather than doing the same job N times.

Progress is reported as workers finish, and each logs to
`<results>/.mvt/logs/<stage>-<day>-shard<k>of<N>.log`:

```
[mvt] slim 2022-11-18: launching 6 MATLAB workers
[mvt] slim 2022-11-18 workers | [##########..........] 3/6 | 3m15s | ~3m15s left
[mvt] slim 2022-11-18: 6/6 shards ok in 245 s
```

### Checking the outputs

After a build, confirm the bytes match the reference:

```matlab
make all Workers 6
mvt.verify                    % all three days, against python/expected/
mvt.verify Days 18            % one day
mvt.verify Verbose true       % list every file, not just failures
```

It reads the same manifests the Python tool uses, so both implementations check
against one set of expected values, and it needs no Python. Only the JSON
products (`gps`, `slim`) are byte-comparable and therefore checked; figures and
`.mat` files are skipped, and so is `dataset_info.json`, which records host and
user by design.

Each manifest records the data set version it was built from. If that differs
from `mvt.dataVersion()` the check warns first, because the differences are then
expected rather than a defect — comparing a 2.1.1 build against a 2.1 manifest
reports the fuel totals that the deterministic quadrature moved.

### Results somewhere else

Only needed for a non-standard layout — a different disk, or keeping runs side
by side. Paths otherwise resolve from the location of the code:

```matlab
setenv('MVT_DATA_DIR',    'D:\mvt-nature\data')
setenv('MVT_RESULTS_DIR', 'D:\mvt-nature\results-pc')
make all Workers 6
```

### Other ways to run

**The original MATLAB drivers**, unchanged:

```matlab
run_all_scripts                        % all three days
run_all_scripts('Days', 17)            % one day
generate_data_mvt_slim(17)             % a single stage, unchanged call style
mvt.build('slim', 17)                  % the uniform stage entry point
mvt.runShards('slim', 17, 6)           % that stage, across 6 processes
```

`make all` is equivalent to:

```matlab
for day = [16 17 18]
    mvt.build('gps', day)
    mvt.runShards('slim', day, 6)
    for stage = ["samples" "fields" "macro" "micro"]
        mvt.build(char(stage), day)
    end
end
mvt.build('av', [])
```
**Using the make.m file in MATLAB**

`make` here is `Scripts/make.m`, not the Unix tool — no `make`, no toolbox and
no environment variables are needed, on any platform including Windows. The `cd`
is only so MATLAB can find it.

Two commands worth running first:

```matlab
make config     % the paths it resolved, and the data set version
make status     % what is stale, and why; builds nothing
```


**The Unix Makefile** (macOS/Linux; drives MATLAB headlessly):

```bash
cd MVT_structured_data
make status                     # what would run
make -j3 all                    # everything, three days in parallel
```
If you want to run only one stage (slim) or on one day (slim-17):
```bash
make SHARDS=4 slim-17           # one stage, one day, 4 processes
```
To watch what is happening in another tab:
```bash
make watch                      # live progress, from a second terminal
```
To verify the output files match expected checksums:
```bash
make verify                     # check outputs against expected checksums
```

Set `MATLAB=/path/to/matlab` if MATLAB is not at the default macOS location.

**The Python port**, which needs no MATLAB — see
[`python/README.md`](python/README.md):

```bash
cd MVT_structured_data/python
./setup_venv.sh --full
source .venv/bin/activate
mvt build -j 8                  # everything stale, all three days
mvt verify                      # check against the expected checksums
```

### 1. System Requirements

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


## About CIRCLES data analyzed with this software

### About data from I-24 MOTION

The `results/slim` data include trajectories recorded and processed by the I-24 MOTION observatory. Primarily the data from I-24 MOTION are comprised of position and speed trajectories from all vehicles detected on the roadway by that observatory. The `results/slim` data align the location and state of each of the CIRCLES control cars includes GPS location data as well as data from on-board vehicle informatics recorded by custom hardware installed with a Raspberri Pi. Those data include information regarding to indicate the control and system state of the vehicle. 

### About data from GPS and the CIRCLES Cars
Two kinds of data were collected from the CIRCLES Cars, and found in the `results/gps` folder in the data archive. Each day provides a separate standalone file for the cars that ran on that day. Information includes GPS data, aligned with timeseries information from data collected directly from the car.

#### Data collection from CIRCLES CARS with team-installed GPS Sensors
Information from team-installed GPS sensors was collected at 10-Hz. These data include the position and speed of the vehicle at each sample point.

#### Data Collection from CIRCLES Cars with team-designed on-board data collection
Information from team-installed computers that interface with the Controller Area Network (CAN) were critical to sensing and control of the experiment cars. These data are aligned with the raw GPS information to provide the state of the vehicle at that time (speed, assigned lane of travel, desired cruise control set point, etc.). 

## What to download

Three routes, easiest first. Each is a superset of the one before it, so pick by
how much you want to rebuild rather than by which figures you want. Every route
assumes this repository sits **next to** the folder you download into:

```
<some folder>/
  MVT_structured_data/   <- this repository
  results/               <- routes 1 and 2 download into here
  data/                  <- route 3 adds this
```

Run everything from `MVT_structured_data/` (for `make`) or
`MVT_structured_data/Scripts/` (for MATLAB).

### Route 1: figures from the analysis files (~5 GB)

The cheapest way to reproduce figures. No trajectory files at all: the stages
read the derived `.mat` intermediates directly.

Download into `results/`:

| what | into | size | needed by |
|---|---|---|---|
| the three GPS files | `results/gps/` | 821 MB | `macro`, `lcplot` |
| `fields_motion_2022-11-DD.mat` | `results/analysis/2022-11-DD/` | 46 MB | `macro` |
| `LC_data_DD.mat` | `results/analysis/2022-11-DD/` | 40 MB | `lcplot` |
| `relspeed_data_DD.mat` | `results/analysis/2022-11-DD/` | 690 MB | `relspeedplot` |
| `samples_for_distance_analysis_DD.mat` | `results/analysis/2022-11-DD/` | 3.9 GB | `av` |

Then:

```bash
cd MVT_structured_data
make figures-from-mat
```

That produces, into `results/figures/`: six macroscopic field figures and four
lane-change figures per day, two relative-speed figures per day, and three
cross-day fuel and sample-count figures. Roughly 10 minutes per day, most of it
`lcplot`.

This route has been tested end to end from a clean tree: 5.0 GB of downloads
produced 15 figures without ever creating a `results/slim` folder.

Two things to know:

* **The cross-day fuel figures (`av`) need all three days.** That stage reads
  every day's `samples_*.mat` and `fields_motion_*.mat` regardless of what you
  pass as `DAYS`, so a single-day download will build `macro`, `lcplot` and
  `relspeedplot` and then fail on `av`. Use `make -k` to let the rest finish, or
  download all three days.
* Leaving out `samples_for_distance_analysis_DD.mat` drops the download to
  **~1.6 GB** and costs you only those cross-day figures.

The two figures this route cannot make are the microscopic trajectory plots
(`micro`), which need the trajectories themselves. `make figures-from-mat`
deliberately excludes them, so it will not go looking for data you have not
downloaded. Plain `make` would attempt them after everything else.

### Route 2: add the slim trajectories (+51 GB)

Adds `results/slim/`, the released westbound trajectories. The trajectory
figures become available, and most of the `.mat` files above can be regenerated
instead of downloaded.

```bash
make            # regenerates what it can, then every figure
make micro      # or just the trajectory figures
```

**One exception, worth knowing before you delete anything.** `fields_*.mat`,
`samples_*.mat` and `relspeed_data_*.mat` all rebuild from `slim/`. The lane
sidecars do **not**: they are derived from the **raw** MOTION segments, because
the released `slim` data no longer carries the lane a trajectory came from. So
this route needs one of:

* the 24 lane sidecars per day (**3.6 MB** for all three days), from which
  `LC_data_DD.mat` rebuilds in about five minutes; or
* the `LC_data_DD.mat` files themselves (**40 MB**), if you would rather not
  spend the five minutes.

Without either, `lcplot` has nothing to build from and will say so. Both paths
are tested: a tree holding only `slim/`, `gps/` and the analysis files produced
all 15 of one day's figures.

`micro` derives ~4.6 GB per day of reduced plotting caches from slim on its
first run, under `results/.mvt/cache/`; later runs reuse them. Expect that
first run to be slow.

### Route 3: add the raw data (+59 GB)

Adds `data/`, the raw I-24 MOTION segments and car CSVs, and rebuilds
everything from them: `gps`, `slim`, the lane sidecars, then the intermediates
and the figures.

```bash
make rebuild                                   # in place
make RESULTS=/path/to/results_new rebuild      # into a fresh tree instead
```

Hours, not minutes. Spread the sharded stages across processes with `SHARDS`,
and size that by memory rather than cores — each process holds a decoded
segment, several GB at peak:

```bash
make SHARDS=6 slim-17     # one day, six MATLAB processes over 24 segments
```

`full` (the eastbound and reference trajectories) is opt-in even here; add
`make full` if you want it.

### Which figures come from where

| figure | stage | reads | route |
|---|---|---|---|
| macroscopic fields | `macro` | `fields_*.mat` + gps | 1 |
| lane-change exposure, rate, cumulative excess | `lcplot` | `LC_data_DD.mat` + gps | 1 |
| relative speed to the nearest AV | `relspeedplot` | `relspeed_data_DD.mat` | 1 |
| fuel results, sample counts | `av` | `samples_*.mat` | 1 |
| microscopic trajectories | `micro` | slim, via the reduced caches | 2 |

`LC_data_DD.mat` is the one intermediate that route 2 cannot regenerate: the
lane sidecars behind it come from the raw segments, so rebuilding it needs
route 3. Download it (or the sidecars) on any route below that.

### If something looks wrong

`make status` says what each stage would do and why, and builds nothing. Start
there:

```bash
make status
```

**A stage whose inputs are absent** names the file it wanted. Check it against
the download table for your route.

**A stage that says a script "is newer than" a file you downloaded** is the one
confusing case. Freshness is decided by modification time, and an archive that
preserves timestamps can arrive older than the code you just cloned. The stage
is not telling you the data is wrong — only that it cannot prove it is current.
You downloaded these outputs rather than building them, so say so:

```bash
make accept-verified    # checks the published checksums, then marks them current
make accept             # same, without the checksum pass
```

This is worth doing on routes 2 and 3, though it is no longer load-bearing:
stages depend on the files they read rather than on the stages that made them,
so a downloaded tree builds without it. What `accept` saves is the time spent
regenerating intermediates you already downloaded — about five minutes for
`LC_data`, longer for `samples`.

On a machine without GNU make, every target above works from MATLAB too:

```matlab
cd MVT_structured_data/Scripts
make figures-from-mat
make status
```

## Data Install

Download the data to your computer, it will be called either `data` or `results`. In the same folder that contains the data or results, clone this repository.

## Reproduce Plots: Scripts to plot figures from the article.

Ensure you have downloaded all the data, and stored it according to the structure described in the next step. 

### Step 0: Correct folder structure for reproducing plots. 

This git repository should be a *sibling* folder to the `results` folder. 

```
cd MVT_structured_data
ls ..
```
This command should show you 

```
results/
MVT_structured_data
```

You may also have folders like `data` if you have downloaded the bootstrapping/base data.

If you want to confirm what it should look like from the 'parent' directory:

```
[Parent]
  - MVT_structured_data
  - results/
  | - slim
    | - 2022-11-16
      | - I-24MOTION_2022-11-16_05-59-59.json
      | - I-24MOTION_2022-11-16_06-09-59.json
      | - ...
    | - 2022-11-17
      | - ...
    | - 2022-11-18
      | - ...
  | - gps
    | - CIRCLES_GPS_10Hz_2022-11-16.json
    | - CIRCLES_GPS_10Hz_2022-11-17.json
    | - CIRCLES_GPS_10Hz_2022-11-18.json
  | - README.md
```

Once you run the scripts, additional folders and files will be produced in the `results` folder.

### Step 1: Reproduce the plots
To reproduce all the plots for all days, simply run 

```
cd MVT_structured_data/Scripts
```
and in that folder in matlab, run

```
reproduce_plots.m
```

This will run the below scripts for all days as well.

Running `Scripts\plot_AV_analysis.m` generates the results from the article (figure 2, figure SM2, and Figure SM3).

Running``"Scripts\plot_macroscopic_fields.m` to generate the macroscopic fields figures from the article (figure 3, and SM 5, as well as additional fields).

### Step 2: Examine figures

The outputs in `results/figures` provide reproductions of the figures used in the main graphics in the paper.

## Advanced Bootstrap
How to generate the integrated data set. (advanced only)

***Note*** This step requires a different dataset to begin, and generates the `slim` and `full` results that are part of the release. These steps perform alignment of vehicle and I-24 MOTION data, from base files from I-24 MOTION and original GPS and vehicle CAN data files that are assigned by each car.

This step is ADVANCED and OPTIONAL and is included mainly to provide the algorithmic insights to anyone interested in how those data are produced and aligned.

### Step 0 (Advanced): Correct folder structure for bootstrap data synthesis. 

This git repository should be a *sibling* folder to the data folder. 

```
cd MVT_structured_data
ls ..
```
This command should show you 

```
data/
MVT_structured_data
```

A folder named `results` will be created as part of these scripts.

If you want to confirm what it should look like from the 'parent' directory:

```
[Parent]
  - data/
  | - cars
    | - cars_gps
      | - circles_v2_1_car1.csv
      | - circles_v2_1_car2.csv
      | - circles_v2_1_car3.csv
      ...
    | - cars_vins.csv
    | - veh_ping_ ...
    | - ...
  | - i24motion
    | - 2022-11-16
      | - 64888dc ....wed_0_00.json
      | - 64888dc ....wed_0_01.json
      | - 64888dc ....wed_0_02.json
      | - ...
    | - 2022-11-17
      | - ...
    | - 2022-11-18
      | - ...
  | - README.md
  - MVT_structured_data
  - results
```

### Step 1 (Advanced): run scripts to generate the data

### Step 1: run the stages

`make all` runs everything below in order, for all three days. The individual
scripts are listed so they can be run one at a time; each takes a day
(`16`, `17`, or `18`) and writes under `results/`.

| | Script | Produces |
|---|---|---|
| 1 | `assemble_data_GPS(day)` | `results/gps/CIRCLES_GPS_10Hz_2022-11-DD.json` |
| 2 | `generate_data_mvt_slim(day)` | `results/slim/2022-11-DD/` — 24 segment files |
| 3 | `generate_data_samples(day)` | `samples_for_distance_analysis_DD.mat` |
| 4 | `generate_macroscopic_fields(day)` | `fields_motion_2022-11-DD.mat` |
| 5 | `plot_macroscopic_fields(day)` | figure 3 and SM5 |
| 5b | `plot_microscopic_trajectories(day)` | trajectory time-space figures |
| 6 | `plot_AV_analysis()` | figures 2, SM2, SM3 (reads all three days) |

#### `slim` and `full`

`generate_data_mvt_full` produces a second, larger variant that additionally
carries eastbound and reference trajectories plus flat-fuel and direction
fields. **It is not built by `all`, in either the MATLAB or the Python
implementation, and nothing downstream reads it.** Stages 3–6 all read `slim`,
which is the released data set.

Build it only if you need eastbound plots — `direction = 1` in
`plot_microscopic_trajectories.m` switches to it:

```
make full                      # MATLAB; or: make full Days 17
mvt full --day 17              # Python; or: mvt build --target full -j 6
```

Both implementations produce it, and their output is byte-identical.

`make status` reports `full` as `opt-in` rather than `BUILD` while it has never
been built, so a missing `full` tree does not read as pending work. Once you
have built it, it is reported like any other stage.

## Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)
