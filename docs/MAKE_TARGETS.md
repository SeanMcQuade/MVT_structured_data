# Make targets and stages

The README covers the three routes most people want. This is the full
reference: every stage, every target, and the options that apply to them.

Two front ends drive the same stages. `make` from the repository root uses GNU
make; `make` from `Scripts/` inside MATLAB is `Scripts/make.m`, for machines
without GNU make. Targets and options are the same in both, with the usual
difference in syntax:

```bash
make slim-17 SHARDS=6          # GNU make, from the repository root
```
```matlab
make slim Days 17 Workers 6    % MATLAB, from Scripts/
```

## Stages

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


## Targets

```matlab
make                          % everything out of date, all three days
make all Workers 6            % the same, with slim and lanes across 6 processes
make data                     % gps, slim, lanes, lc, relspeed, samples, fields
make figures                  % the figures built from the .mat analysis files:
                              %   macro, lcplot, av, relspeedplot
make micro                    % the trajectory plots (these read slim)
make slim                     % one stage
make slim Days 18             % one stage, one day
make rebuild                  % force everything, from the raw data up
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


## How many workers

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


## Results somewhere else

Only needed for a non-standard layout — a different disk, or keeping runs side
by side. Paths otherwise resolve from the location of the code:

```matlab
setenv('MVT_DATA_DIR',    'D:\mvt-nature\data')
setenv('MVT_RESULTS_DIR', 'D:\mvt-nature\results-pc')
make all Workers 6
```


## Checking the outputs

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


## If a run stops early

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


## Other ways to run

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

## Target names that changed

`figures` used to mean every figure. It now means the figures that build from
the `.mat` analysis files, which is what a figures-only download can make; the
trajectory plots are `micro`, and `make all` is `figures` plus `micro`. The
interim names `figures-from-mat` and `figures-from-slim` still work and are
aliases for `figures` and `micro`.
