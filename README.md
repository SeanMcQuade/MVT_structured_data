
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022. 
This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the 
installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

# Contents
- [Quick start](#quick-start)
- [Running the pipeline](#running-the-pipeline)
- [Documentation](#documentation)
- [Scripts to generate the integrated data set](#scripts-to-generate-the-integrated-data-set)
- [Scripts to plot figures from the article](#scripts-to-plot-figures-from-the-article)
- [Websites](#websites)

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

### 2. Run it

From MATLAB:

```matlab
cd MVT_structured_data/Scripts
make all Workers 6
```

That is the whole thing. It builds every stage, for all three days, in the right
order, with the heavy stage spread across 6 MATLAB processes. Work that is
already up to date is skipped, so an interrupted run picks up where it left off.

`make` here is `Scripts/make.m`, not the Unix tool — no `make`, no toolbox and
no environment variables are needed, on any platform including Windows. The `cd`
is only so MATLAB can find it.

Two commands worth running first:

```matlab
make config     % the paths it resolved, and the data set version
make status     % what is stale, and why; builds nothing
```

Everything else — choosing the worker count, building one day or one stage,
writing results elsewhere, the Unix `make`, and the Python port — is under
[Running the pipeline](#running-the-pipeline).

## Documentation

| Document | What it covers |
| --- | --- |
| [`docs/ALGORITHMS.md`](docs/ALGORITHMS.md) | The data-flow graph and what each stage computes, with the mapping to the paper's figures. Start here. |
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

### Targets

```matlab
make                          % everything out of date, all three days
make all Workers 6            % the same, with slim across 6 processes
make data                     % gps, slim, samples, fields
make figures                  % macro, micro
make slim                     % one stage
make slim Days 18             % one stage, one day
make status                   % what is stale, and why; builds nothing
make config                   % resolved paths, days, workers, data version
```

Options may follow any target: `Workers`, `Days`, `Force`, `Clean`, `DryRun`,
`Verbose`, `SettleSeconds` (see `Scripts/+mvt/options.m`). Command syntax is
fine — `make slim Days 18 Force true` — as is function syntax,
`make('slim', 'Days', 18)`.

```matlab
make all DryRun true          % plan the whole run, write nothing
make slim Force true          % rebuild regardless of timestamps
```

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

**The Unix Makefile** (macOS/Linux; drives MATLAB headlessly):

```bash
cd MVT_structured_data
make status                     # what would run
make -j3 all                    # everything, three days in parallel
make SHARDS=4 slim-17           # one stage, one day, 4 processes
make watch                      # live progress, from a second terminal
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
