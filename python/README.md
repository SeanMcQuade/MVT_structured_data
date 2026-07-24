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

## Run

The CLI (`mvt`, or `python -m mvtpy`) runs one stage for one or more days:

```bash
mvt gps      --day 16          # assemble the control-vehicle GPS JSON
mvt samples  --day 16          # distance-to-AV samples  -> .npz
mvt fields   --day 16          # macroscopic fields       -> .npz
mvt figures  --day 16          # field heatmaps + AV fuel curves -> .png
mvt all      --day 16          # the four above, in order

mvt samples --day 16 --day 17 --day 18   # several days
mvt fields                                # all three days (default)
```

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
| `samples` | `results/slim/2022-11-DD/` | `results/figures/2022-11-DD/samples_for_distance_analysis_DD.npz` |
| `fields` | `results/slim/2022-11-DD/` | `results/figures/2022-11-DD/fields_motion_2022-11-DD.npz` |
| `figures` | `results/slim/`, `results/gps/`, the samples `.npz` | `results/figures/2022-11-DD/fig_*.png` |

The `.npz` outputs carry the same arrays as the MATLAB `.mat` (load with
`numpy.load`). Producing the `slim` JSON itself from raw MOTION is implemented in
`mvtpy.slim` but not wired into the CLI, since it needs the assembled GPS first
and writes ~0.9 GB per segment — add it if you need the full bootstrap from raw
data.

`gps` is the heavy stage (~12 min/day; it streams the day's raw MOTION segments
for the matching pass); `samples` and `fields` are ~4 min/day (dominated by
decoding the slim JSON).

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
| `workspace`, `cli` | path resolution and the command-line entry point |
