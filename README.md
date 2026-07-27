
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

Run from the `Scripts` folder as before:

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

No `make` is required. Set the two locations (or let them default to `data/` and
`results/` beside this repository), then run the stages in order:

```matlab
setenv('MVT_DATA_DIR',    'D:\mvt-nature\data')
setenv('MVT_RESULTS_DIR', 'D:\mvt-nature\results-pc')
cd('C:\path\to\MVT_structured_data\Scripts')

mvt.build('gps',     16)      % ~12 min/day
mvt.build('slim',    16)      % ~40 min/day sequentially - see below
mvt.build('samples', 16)
mvt.build('fields',  16)
mvt.build('macro',   16)
mvt.build('micro',   16)
mvt.status('Days', 16)        % what is stale, and why
```

**Concurrency: `mvt.build` on its own is sequential.** Nothing in this pipeline
uses `parfor`; `opts.UseParfor` is declared but never read by any stage, and the
Parallel Computing Toolbox is not required anywhere. All the parallelism in the
Makefile comes from launching *separate MATLAB processes*. `mvt.runShards` does
the same thing from inside MATLAB, so the `make -j8` approach works on Windows:

```matlab
mvt.runShards('slim', 16, 8)   % 8 concurrent MATLAB processes over the 24 segments
mvt.runShards('full', 16, 8)
```

Each worker takes segments `k, k+N, k+2N, …`, logs to
`<results>/.mvt/logs/<stage>-<day>-shard<k>of<N>.log`, and reports its exit
status; the caller blocks until all workers finish and raises if any failed.

Only `slim` and `full` shard by segment, so those are the only stages that gain
from this — `runShards` refuses more than one worker for the others rather than
silently doing the same work N times. **Size the worker count by memory, not
cores**: each process holds a decoded segment and peaks at several GB.

For a whole day, the sequential stages plus a sharded `slim` is the useful
shape:

```matlab
for day = [16 17 18]
    mvt.build('gps', day)
    mvt.runShards('slim', day, 8)
    for stage = ["samples" "fields" "macro" "micro"]
        mvt.build(char(stage), day)
    end
end
```

Note that `gps` must finish before `slim` for the same day (slim reads the
assembled GPS), and `slim` before the rest.

### With make (parallel)

`make` encodes the same dependency graph and runs independent work
concurrently, using ordinary MATLAB batch processes — no Parallel Computing
Toolbox required.

```bash
make config              # show resolved paths and settings
make status              # what is stale, and why
make -n all              # dry run: print the plan
make -j3 all             # everything, three days in parallel
make SHARDS=4 slim-17    # one day, four processes over its 24 segments
make FORCE=1 fields-17   # rebuild even if outputs look fresh
make clean-figures-17    # narrow, explicit cleanup
```

Set `MATLAB=/path/to/matlab` if MATLAB is not at the default macOS location.

### Writing results somewhere else (for comparison)

`RESULTS=` (or the `MVT_RESULTS_DIR` environment variable) sends output to a
different tree, so a new run can be compared against an existing one without
overwriting it:

```bash
make RESULTS=$PWD/../results_verify SHARDS=1 slim-16
md5 ../results/slim/2022-11-16/I-24MOTION_2022-11-16_05-59-59.json \
    ../results_verify/slim/2022-11-16/I-24MOTION_2022-11-16_05-59-59.json
```

`MVT_DATA_DIR` does the same for the raw inputs.

### Accepting existing outputs after a no-op code change

Because staleness is decided from timestamps, editing a stage marks its outputs
for rebuild even when the edit cannot change results. Once you have confirmed
that (by rebuilding into a separate tree and comparing checksums, as above):

```bash
make accept              # all days; mvt.accept('DryRun', true) to preview
make DAYS=17 accept      # one day
```

This only updates timestamps of files that already exist; it never writes data.

Choose the number of concurrent processes by memory, not cores: each worker
peaks in the multi-GB range while decoding a raw segment and encoding its
output. Two to four workers is a reasonable start on a 32 GB machine.

Measured on an Apple Silicon Mac (137 GB RAM), rebuilding one day of `slim`:

| Run | Wall time | Per segment |
| --- | --- | --- |
| serial | ~24 min (est. from 59 s/segment) | 59 s |
| `SHARDS=4` | 10 min 43 s at 402% CPU | ~100 s |

That is a 2.2x speedup rather than 4x: the stage is dominated by reading a
~2 GB raw segment and writing a ~800 MB result, so workers contend for I/O and
memory bandwidth rather than CPU. Shards stayed well balanced (583-634 s each),
which is what the interleaved assignment is for. All 24 outputs were
byte-identical to a serial reference run.

### Tests

Two levels, both driven from `make`:

```bash
make test           # fast suite: ~40 tests, seconds, no data tree required
make verify-full    # check real outputs against the recorded manifests
```

`make test` runs `tests/` through MATLAB's `runtests`: options parsing and
sharding, the staleness rule that decides every rebuild, path and output-name
conventions, and the `.mat` comparators. It works on a machine with no data,
using temporary files.

`make verify-full` compares the actual outputs against `tests/manifests/`,
which record what a known-good run produced: md5 for the JSON files (which are
deterministic), `mvt.matHash` content hashes for the `.mat` files (whose bytes
are not comparable — v7 stores a gzip timestamp, v7.3 is HDF5), and sizes for
figures, reported but never failed since pixels vary with renderer and MATLAB
release.

To record manifests from a tree you trust:

```bash
make RESULTS=$PWD/../results_groundtruth manifests
```

The manifests are small (about 10 KB per day) and belong in version control.

### Derived caches and generated state

Build state lives under `results/.mvt/`:

| Path | Contents |
| --- | --- |
| `results/.mvt/manifests/` | raw segment → output filename maps |
| `results/.mvt/deps/` | cached per-stage code dependency lists |
| `results/.mvt/stamps/` | make stamp files |
| `results/.mvt/cache/` | `*_reduced.mat` plotting caches |

The `*_reduced.mat` plotting caches used to be written next to the released
JSON in `results/slim/<day>/`. They now live in `results/.mvt/cache/<day>/` so
that `results/slim` contains only published artifacts. If you already have the
old caches, move them once with:

```bash
make migrate-cache
```

Everything under `results/.mvt/` is derived and safe to delete; the next run
rebuilds it.

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
