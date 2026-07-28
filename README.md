
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022. 
This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the 
installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

# Contents
- [Quick start](#quick-start)
- [Running the pipeline](#running-the-pipeline)
- [Documentation](#documentation)
- [Generate integrated data.](#tag1)
- [Plot data and results.](#tag2)
- [Websites](#tag3)

## Quick start

### 1. Put the data in place

The pipeline expects this repository to sit **next to** the `data/` and
`results/` folders:

```
<some folder>/
  MVT_structured_data/   <- this repository
  data/                  <- raw inputs (cars/, i24motion/)   — for bootstrapping
  results/               <- processed inputs and outputs (slim/, gps/, figures/)
```

The data are distributed separately (see the download instructions further
down). Once you have them, check the layout is correct:

```bash
./check_data.sh          # reports what's present and which workflows can run
```

It tells you whether you can **plot/analyze** (needs `results/slim` + `results/gps`)
and/or **bootstrap from raw data** (needs `data/cars` + `data/i24motion`). If the
folders are elsewhere, point it at them:

```bash
MVT_DATA_DIR=/abs/data MVT_RESULTS_DIR=/abs/results ./check_data.sh
```

### 2. Run — pick one

**A. From within MATLAB** (interactive; the reference implementation)

```matlab
cd MVT_structured_data/Scripts
run_all_scripts                 % all three days, skipping work already done
run_all_scripts('Days', 17)     % one day
mvt.status                      % what is stale and why (nothing runs)
```

**B. MATLAB from the command line** (headless, incremental, parallel)

```bash
cd MVT_structured_data
make check-data                 # same layout check as above
make status                     # what would run
make -j3 all                    # everything, three days in parallel
make figures                    # just the plots (needs processed data)
make SHARDS=4 slim-17           # one stage, one day, 4 processes
```

Set `MATLAB=/path/to/matlab` if MATLAB is not at the default macOS location.

**C. Python from the command line** (no MATLAB needed; see `python/README.md`)

```bash
cd MVT_structured_data/python
./setup_venv.sh --full          # create .venv and install (numpy, pandas, matplotlib, ...)
source .venv/bin/activate
mvt fields  --day 16            # macroscopic fields  -> .npz
mvt figures --day 16            # field heatmaps + AV fuel curves -> .png
mvt all     --day 16            # gps, samples, fields, figures
```

The Python port also runs in Docker against a mounted data directory — see
[`python/README.md`](python/README.md).

## Documentation

| Document | What it covers |
| --- | --- |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | The data-flow graph and what each stage computes, with the mapping to the paper's figures. Start here. |
| [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) | Every field of every output: units, dtype, rounding, null/empty semantics, and which stage writes it. |
| [`docs/MATLAB_JSON_FORMAT.md`](docs/MATLAB_JSON_FORMAT.md) | How MATLAB's `jsonencode` formats numbers, needed for byte-identical output. |
| [`docs/PYTHON_PORT.md`](docs/PYTHON_PORT.md) | Status of the Python implementation and how its parity is verified. |
| [`python/README.md`](python/README.md) | The Python port: install (venv), run (CLI), Docker, and test. |

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

### From MATLAB

One command builds everything, from the `Scripts` folder:

```matlab
make                    % every stage, all three days, skipping fresh work
make all Workers 6      % the same, with slim split across 6 MATLAB processes
make status             % what is stale, and why; builds nothing
```

`make` is `Scripts/make.m` — it needs no Unix `make`, no toolbox, and no
environment variables on the standard layout. It works on every platform; see
[Building on Windows](#on-windows-entirely-from-matlab) for the
concurrency details, which matter most where the Makefile cannot run.

The original entry points are unchanged:

```matlab
run_all_scripts                        % all three days, skipping fresh work
run_all_scripts('Days', 17)            % one day
run_all_scripts('Force', true)         % rebuild everything
run_all_scripts('DryRun', true)        % show what would run, change nothing
generate_data_mvt_slim(17)             % single stage, unchanged call style
generate_data_mvt_slim(17, 'Force', true)
mvt.status('Days', 17)                 % what is stale, and why
```

Options (see `Scripts/+mvt/options.m`): `Force`, `Clean`, `DryRun`, `Verbose`,
`Shard`, `Days`, `SettleSeconds`. They may be given as name/value pairs or as an
options struct.

### On Windows, entirely from MATLAB

Two lines. No `make`, no toolbox, no environment variables:

```matlab
cd C:\path\to\MVT_structured_data\Scripts
make all Workers 6
```

That builds every stage, for all three days, in the right order, with `slim`
split across 6 MATLAB processes. Stages that are already up to date are skipped,
so re-running after an interruption picks up where it left off.

`make` here is `Scripts/make.m`, not the Unix tool. The `cd` is only so MATLAB
can find it; paths are resolved from the location of the code itself, so with
the standard layout — this repository sitting beside `data/` and `results/` —
the defaults are already right.

Check before committing to a long run:

```matlab
make config          % the paths it resolved, and the data set version
make status          % what is stale, and why; builds nothing
make all DryRun true % plan the whole run, write nothing
```

Other useful forms:

```matlab
make all Days 18 Workers 6    % one day
make slim Days 18 Workers 6   % one stage
make figures                  % just macro and micro
make slim Force true          % rebuild regardless of timestamps
```

#### Results somewhere else

Only needed for a non-standard layout — a different disk, or keeping runs side
by side:

```matlab
setenv('MVT_DATA_DIR',    'D:\mvt-nature\data')
setenv('MVT_RESULTS_DIR', 'D:\mvt-nature\results-pc')
make all Workers 6
```

#### Driving the stages yourself

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

#### Why `runShards`, and how many workers

`mvt.build` is **sequential**. Nothing in this pipeline uses `parfor`,
`opts.UseParfor` is declared but never read by any stage, and no toolbox is
required. All parallelism — in the Makefile and here — comes from running
several MATLAB processes over disjoint shards of the same day.
`mvt.runShards(stage, day, N)` launches N background MATLAB processes, each
taking segments `k, k+N, k+2N, …`, then waits for them all and raises if any
failed. Progress is reported as workers finish:

```
[mvt] slim 2022-11-18: launching 6 MATLAB workers
[mvt] slim 2022-11-18 workers | [##########..........] 3/6 | 3m15s elapsed | ~3m15s left
[mvt] slim 2022-11-18: 6/6 shards ok in 245 s
```

Each worker logs to `<results>/.mvt/logs/<stage>-<day>-shard<k>of<N>.log`.

**Choose N by memory, not by cores.** Each worker holds a decoded 10-minute
segment and peaks at several GB, so 6 workers wants roughly 30 GB. On a 32 GB
machine use 4–6; 12 workers on a 128 GB machine cut a day of `slim` from about
40 minutes to 4.

Only `slim` and `full` split by segment. `runShards` refuses more than one
worker for the other stages rather than silently doing the same job N times:

```matlab
mvt.runShards('fields', 16, 4)   % errors: mvt:runShards:notSharded
```

Verified on 2022-11-18: 12 workers, 24 segments, 12/12 succeeded, and all 24
outputs byte-identical to the same day built sequentially.

## Scripts to generate the integrated data set.

### Step 0: Correct folder structure before you begin. 

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


### Step 1: run some scripts to generate the data

You have a couple options here.

---
Run "Scripts\generate_data_mvt_full" or "Scripts\generate_data_mvt_slim.":This load files from `../data/cars/gps` and saves assembled GPS.json into Data_2022-11-??__I24_Base for each day. 

Step 2: Once the slimmed or full data is generated, run "Scripts\generate_data_samples.m." It produces mat-files in the Folder called 'Data_Analysis,' one for each day.

Step 3: Run 'Scripts\generate_macroscopic_fields.m.' This will save macoscopic fields data to the "Data\Data_for_Figures" folder as a matlab file ".mat."


## Scripts to plot figures from the article.
Step 4: Run "Scripts\plot_AV_stats.m' to generate the results from the article (figure 2, figure SM2, and Figure SM3).

Step 5: Run "Scripts\plot_macroscopic_stats.m' to generate the macroscopic fields figures from the article (figure 3, and SM 5, as well as additional fields).

## Websites
[Visit the CIRCLES consortium website](https://circles-consortium.github.io/)

[This repository is licensed under the BSD3 license](https://opensource.org/license/bsd-3-clause)
