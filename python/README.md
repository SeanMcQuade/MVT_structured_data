# mvtpy — Python implementation of the MVT pipeline

A Python port of the CIRCLES MegaVanderTest data pipeline. It reproduces the
MATLAB stages that turn raw I-24 MOTION and control-vehicle data into the
released data sets and the paper's figures, and is verified against the MATLAB
output (byte-identical for the JSON, exact or float-exact for the `.mat`
producers). See [`../docs/PYTHON_PORT.md`](../docs/PYTHON_PORT.md) for the
stage-by-stage verification status and
[`../docs/ALGORITHMS.md`](../docs/ALGORITHMS.md) for what each stage computes.

The package expects the standard workspace layout — `data/` and `results/` as
siblings of `MVT_structured_data/` — but every location is overridable, so it
runs equally well against mounted directories in a container.

## Install (virtual environment)

From this `python/` directory:

```bash
./setup_venv.sh          # core: numpy, pandas
./setup_venv.sh --full   # + matplotlib, h5py, scipy (figures, reading .mat)
./setup_venv.sh --dev    # + pytest and the full extras (to run the tests)

source .venv/bin/activate
python -m mvtpy --help
```

Or by hand:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[full]"     # or just "." for the core stages
```

Python 3.9+ is required. The core install (numpy + pandas) is enough to run the
`gps`, `samples`, and `fields` stages; `full` adds matplotlib for figures and
h5py/scipy for reading MATLAB `.mat` reference files.

Before running, confirm the data are laid out where the pipeline expects
(this repository as a sibling of `data/` and `results/`):

```bash
../check_data.sh          # or: MVT_DATA_DIR=/abs/data MVT_RESULTS_DIR=/abs/results ../check_data.sh
```

## Run

The CLI (`mvt`, or `python -m mvtpy`) runs one stage for one or more days:

```bash
mvt gps      --day 16          # assemble the control-vehicle GPS JSON
mvt slim     --day 16          # released slim segments from raw MOTION
mvt samples  --day 16          # distance-to-AV samples  -> .npz
mvt fields   --day 16          # macroscopic fields       -> .npz
mvt figures  --day 16          # field heatmaps + AV fuel curves -> .png
mvt micro    --day 16          # microscopic trajectory time-space figures -> .png
mvt all      --day 16          # all six above, in order

mvt samples --day 16 --day 17 --day 18   # several days
mvt fields                                # all three days (default)
```

### `mvt verify` — proving the outputs are the expected ones

```bash
mvt verify --day 18                       # against the stored checksums
mvt verify --day 18 --results-dir /tmp/x  # verify a tree built elsewhere
mvt verify --day 18 --reference ../results   # explain how it differs
mvt verify --day 18 --update              # re-record the expected checksums
```

or `make verify` / `make verify-against RESULTS=…` / `make verify-update`.

The expected md5s live in `python/expected/checksums-2022-11-DD.json` and were
**derived from the MATLAB outputs**, so a pass means a pipeline reproduced the
released files byte for byte. `mvt verify` exits non-zero on any mismatch, so it
works as a test.

What can and cannot be checked this way, because the distinction matters more
than a green tick:

| Output | Check | Why |
| --- | --- | --- |
| `gps`, `slim` (`.json`) | **md5, byte-exact** | the port's stated acceptance test |
| `samples`, `fields` (`.npz`) | array comparison vs the MATLAB `.mat` | the container formats cannot be byte-compared at all; the `.npz` md5 still pins Python-to-Python reproducibility |
| figures (`.png`) | **skipped**, and reported as skipped | matplotlib and MATLAB renderers do not agree pixel for pixel, and PNG bytes carry encoder metadata — a checksum would fail for reasons unrelated to correctness |

Current status: **75 of 75 byte-comparable files match** — the 3 assembled GPS
files and 72 slim segments, for all three days, verified against manifests
derived from a MATLAB build. The same manifests are matched by MATLAB on macOS
and on Windows, so two independent implementations on two operating systems
produce identical bytes. See
[../docs/PYTHON_PORT.md](../docs/PYTHON_PORT.md) for how that was reached and
[../docs/DATA_CHANGELOG.md](../docs/DATA_CHANGELOG.md) for the data set version.

### `mvt build` — the make

`build` is a pure-Python equivalent of the MATLAB `Makefile`: it plans the whole
dependency graph, rebuilds only what is stale, and runs units in parallel.

```bash
mvt build -j 8                  # everything stale, all three days
mvt build --day 18 -n           # dry run: what is stale, and why
mvt build --target slim -j 6    # one stage
mvt build --day 18 -B           # force, ignore timestamps
mvt build -j 8 -k               # keep going past a failure
```

A target is stale when an output is missing, when a data input is newer, or
when the **code** that produces it is newer — the classic makefile rule, matching
`Scripts/+mvt/isStale.m`. Staleness also propagates: if `slim` is rebuilding,
everything downstream of it is stale too.

`-j` needs no shard arithmetic. Every stage is one unit per day except `slim`,
which fans out to **one unit per 10-minute segment** — so a three-day build is 3
gps units and 72 slim units, and `-j10` fills all ten workers across days *and*
segments by itself. Each unit is its own process (no threads, matching the
MATLAB pipeline's process-level parallelism), so peak memory is returned to the
OS as units finish. **Cap `-j` by memory, not cores**: `slim` peaks near 5 GB per
process, so `-j8` wants ~40 GB.

`slim` also accepts explicit `--segment N` and `--shard K/N` for driving it by
hand, mirroring the MATLAB `Shard [k N]` option.

#### Progress display

On a terminal, `build` draws a live block — a bar, one line per running unit
with its elapsed time, and an ETA from observed unit durations:

```
[████████████░░░░░░░░░░░░░░░░░░] 12/29  41%  8m12s elapsed  ~8m58s left
  ⠦ slim 2022-11-18 #06  3m21s
  ⠦ slim 2022-11-18 #07  2m48s
  12 built, 15 queued
```

When output is **redirected**, it degrades to one durable line per event with no
escape codes, so a `tee`d log stays greppable and survives an interrupted run.
`--progress bar|plain|auto` overrides the choice, and `MVT_NO_PROGRESS=1` or
`TERM=dumb` forces plain.

Locations default to the sibling `data/` and `results/` folders and can be
pointed anywhere:

```bash
mvt samples --day 16 --data-dir /abs/data --results-dir /abs/results
# or via the environment (what the container uses)
MVT_DATA_DIR=/abs/data MVT_RESULTS_DIR=/abs/results mvt fields --day 16
```

Inputs and outputs, by stage (day `DD`):

| Stage | Reads | Writes |
| --- | --- | --- |
| `gps` | `data/cars/*`, `data/i24motion/2022-11-DD/` | `results/gps/CIRCLES_GPS_10Hz_2022-11-DD.json` |
| `slim` | `data/i24motion/2022-11-DD/`, the GPS JSON, `Models/Eastbound_grade_fit.csv` | `results/slim/2022-11-DD/I-24MOTION_*.json` (24 files) |
| `samples` | `results/slim/2022-11-DD/` | `results/figures/2022-11-DD/samples_for_distance_analysis_DD.npz` |
| `fields` | `results/slim/2022-11-DD/` | `results/figures/2022-11-DD/fields_motion_2022-11-DD.npz` |
| `figures` | `results/slim/`, `results/gps/`, the samples `.npz` | `results/figures/2022-11-DD/fig_*.png` |
| `micro` | `results/slim/2022-11-DD/` | `results/figures/2022-11-DD/fig_motion_trajectories_*_py_{lowres,zoomwin,zoom}.png` |

The `.npz` outputs carry the same arrays as the MATLAB `.mat` (load with
`numpy.load`).

`gps` is the heavy stage (~12 min/day; it streams the day's raw MOTION segments
for the matching pass); `samples` and `fields` are ~4 min/day (dominated by
decoding the slim JSON).

`micro` is stage 5b, the one the MATLAB pipeline could not finish without
exhausting memory. It is kept out of `all` because a cold run reads the whole
15 GB/day of slim JSON (~2 min/segment); each segment's triangulation is then
cached in `results/.mvt/cache/2022-11-DD/*_micro.npz`, and later runs reuse it
in about a second per segment. `--no-cache` forces a re-read, `--skip-t N`
changes the plotting sub-sample (default 50; smaller is denser and slower).
Output names carry a `_py` marker so they never overwrite the MATLAB figures of
the same stage, which land in the same folder.

The memory behavior is the point of this stage, so it is worth stating: one
matplotlib artist is created **per segment file**, not per trajectory (~24 for a
day, against ~370,000 patches in the original MATLAB), and the artist is a
`TriMesh` holding flat arrays rather than a `PolyCollection`, which would build
one `Path` object per polygon. Peak RSS is well under 1 GB for a day.

## Run in Docker (mounted directories)

The image bundles only the code; the data and results stay on the host and are
mounted, so outputs land directly in your `results/` folder.

```bash
# build (from this python/ directory)
docker build -t mvtpy .

# run a stage, mounting data read-only and results writable
docker run --rm \
  -v /abs/path/to/data:/work/data:ro \
  -v /abs/path/to/results:/work/results \
  mvtpy samples --day 16
```

The image sets `MVT_DATA_DIR=/work/data` and `MVT_RESULTS_DIR=/work/results`, so
no path flags are needed once the volumes are mounted. A ready-made service is
in `docker-compose.yml`:

```bash
DATA_DIR=/abs/data RESULTS_DIR=/abs/results \
  docker compose run --rm mvtpy all --day 16
```

## Test

```bash
pip install -e ".[dev]"
python -m pytest tests -q
```

The suite has two layers: fast unit tests (no data needed) and parity tests that
compare against the released `data/`+`results/` when present, and that skip
cleanly otherwise. A few full-scale parity tests are opt-in via `MVT_RUN_SLOW=1`.

## Package layout

| Module | Role |
| --- | --- |
| `matjson`, `matround` | MATLAB-compatible JSON encoding and rounding |
| `rawio` | streaming reader for the multi-GB raw MOTION JSON |
| `kinematics`, `fuel`, `lanes`, `avdist`, `slim` | stage 2: the `slim` JSON producer |
| `gpsruns`, `gpsassemble`, `gpsmatch` | stage 1: control-vehicle GPS assembly |
| `samples` | stage 3: distance-to-AV samples |
| `fields` | stage 4: macroscopic fields |
| `avanalysis` | stage 6 core: fuel-vs-distance binning |
| `plotting` | matplotlib figures (matching colors/layout, not pixel-exact) |
| `microplot` | stage 5b: microscopic trajectory time-space figures, batched per file |
| `segments` | raw segment → processed filename, resolved before any decode |
| `build` | the make: dependency graph, staleness, `-j` scheduler |
| `workspace`, `cli` | path resolution and the command-line entry point |
