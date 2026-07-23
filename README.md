
These scripts and data are from the CIRCLES consortium MegaVanderTest experiment that was conducted in Nashville on highway I-24 during the week of November 14th 2022. 
This repository contains data recorded from the partially automated vehicles, called the GPS data (this includes GPS location data recoreded by the 
installed raspberri pi and several flags to indicate the state of the vehicle) and data recorded by the I-24 MOTION observatory.

# Contents
- [Running the pipeline](#running-the-pipeline)
- [Generate integrated data.](#tag1)
- [Plot data and results.](#tag2)
- [Websites](#tag3)

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
