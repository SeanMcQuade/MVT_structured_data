"""Command-line entry point for the Python MVT pipeline.

Runs the ported stages against a data tree and writes outputs into a results
tree, so the whole analysis can run from a shell or inside a container:

    python -m mvtpy gps      --day 16
    python -m mvtpy samples  --day 16
    python -m mvtpy fields   --day 16
    python -m mvtpy figures  --day 16
    python -m mvtpy all      --day 16

Locations default to the sibling ``data``/``results`` folders of the repository
and can be pointed elsewhere with ``--data-dir`` / ``--results-dir`` or the
``MVT_DATA_DIR`` / ``MVT_RESULTS_DIR`` environment variables - which is how the
Docker image writes into a mounted host directory.

Notes
-----
* ``gps`` reads the raw vehicle CSVs and MOTION segments and writes the
  assembled GPS JSON. It is the heavy stage (the MOTION-matching pass).
* ``samples`` and ``fields`` read the released ``slim`` JSON and write ``.npz``
  (portable) next to where the MATLAB ``.mat`` would go.
* ``figures`` renders the field heatmaps and the AV fuel curves (needs
  matplotlib).
* Producing the ``slim`` JSON itself from raw MOTION is available in
  ``mvtpy.slim`` but not wired here, as it needs the assembled GPS first and
  writes ~0.9 GB per segment; add it if you need the full bootstrap.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from .workspace import Workspace

DAYS = (16, 17, 18)


def main(argv=None) -> int:
    # Shared options accepted both before and after the subcommand. The
    # SUPPRESS defaults are essential: with parents=, a subparser re-declares
    # these, and an ordinary None default would clobber a value parsed before
    # the subcommand. SUPPRESS leaves the attribute unset when not given, so the
    # value from either position survives.
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--data-dir", default=argparse.SUPPRESS,
                        help="raw data tree (default: sibling data/)")
    common.add_argument("--results-dir", default=argparse.SUPPRESS,
                        help="results tree (default: sibling results/)")

    parser = argparse.ArgumentParser(prog="mvtpy", parents=[common], description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="stage", required=True)

    for name, help_text in (
        ("gps", "assemble the control-vehicle GPS JSON"),
        ("samples", "collect distance-to-AV samples (.npz)"),
        ("fields", "accumulate macroscopic fields (.npz)"),
        ("figures", "render field heatmaps and AV fuel curves"),
        ("all", "gps, samples, fields, figures in order"),
    ):
        stage = sub.add_parser(name, parents=[common], help=help_text)
        stage.add_argument("--day", type=int, choices=DAYS, action="append",
                           help="day to process; repeatable (default: all three)")

    args = parser.parse_args(argv)
    workspace = Workspace.resolve(getattr(args, "data_dir", None),
                                  getattr(args, "results_dir", None))
    days = args.day or list(DAYS)

    runners = {"gps": _run_gps, "samples": _run_samples, "fields": _run_fields,
               "figures": _run_figures, "all": _run_all}
    print(f"[mvtpy] data={workspace.data_dir}  results={workspace.results_dir}")
    for day in days:
        runners[args.stage](workspace, day)
    return 0


# ---------------------------------------------------------------------------


def _run_gps(ws: Workspace, day: int) -> None:
    from . import gpsassemble, matjson

    _require(ws.cars_dir(), "cars data")
    _require(ws.motion_dir(day), "raw MOTION segments")
    timings: dict = {}
    with _timer(f"gps 2022-11-{day}"):
        records = gpsassemble.assemble_day(day, ws.data_dir, timings=timings)
    for phase, seconds in timings.items():
        print(f"    {phase:18s} {seconds:6.1f}s")

    out = ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{day}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(matjson.dumps([_arrays_to_lists(r) for r in records]), encoding="utf-8")
    print(f"    wrote {out}")


def _run_samples(ws: Workspace, day: int) -> None:
    from . import samples

    slim = _require(ws.slim_dir(day), "processed slim JSON")
    with _timer(f"samples 2022-11-{day}"):
        data = samples.collect_samples_from_dir(slim, day)
    _save_npz(ws.figures_dir(day) / f"samples_for_distance_analysis_{day}.npz", data)


def _run_fields(ws: Workspace, day: int) -> None:
    from . import fields

    slim = _require(ws.slim_dir(day), "processed slim JSON")
    with _timer(f"fields 2022-11-{day}"):
        result = fields.macroscopic_fields_from_dir(slim, day)
    flat = {"t": result["t"], "x": result["x"],
            "direction": result["direction"], "lane": result["lane"]}
    for name, value in result["field"].items():
        flat[f"field_{name}"] = value
    _save_npz(ws.figures_dir(day) / f"fields_motion_2022-11-{day}.npz", flat)


def _run_figures(ws: Workspace, day: int) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("    matplotlib not installed; skipping figures "
              "(pip install 'mvtpy[full]')")
        return
    from . import avanalysis, fields, plotting
    from .rawio import iter_trajectories

    slim = _require(ws.slim_dir(day), "processed slim JSON")
    figures = ws.figures_dir(day)
    figures.mkdir(parents=True, exist_ok=True)

    with _timer(f"figures 2022-11-{day}"):
        field_data = fields.macroscopic_fields_from_dir(slim, day)
        gps_file = ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{day}.json"
        gps = list(iter_trajectories(gps_file)) if gps_file.is_file() else None
        for name in ("Rho", "Q", "F", "U", "Phi", "Psi"):
            fig, ax = plt.subplots(figsize=(15, 5))
            plotting.plot_field(field_data, name, gps_records=gps, ax=ax)
            fig.tight_layout()
            fig.savefig(figures / f"fig_field_2022-11-{day}_{name}.png", dpi=120)
            plt.close(fig)

        samples_npz = figures / f"samples_for_distance_analysis_{day}.npz"
        if samples_npz.is_file():
            stats = avanalysis.bin_samples(dict(np.load(samples_npz)))
            fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
            plotting.plot_av_fuel({day: stats}, axes=axes)
            fig.tight_layout()
            fig.savefig(figures / f"fig_av_fuel_2022-11-{day}.png", dpi=120)
            plt.close(fig)
    print(f"    wrote figures to {figures}")


def _run_all(ws: Workspace, day: int) -> None:
    _run_gps(ws, day)
    _run_samples(ws, day)
    _run_fields(ws, day)
    _run_figures(ws, day)


# ---------------------------------------------------------------------------


def _save_npz(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **data)
    print(f"    wrote {path}")


def _arrays_to_lists(record: dict) -> dict:
    return {key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in record.items()}


def _require(path: Path, what: str) -> Path:
    if not path.exists():
        raise SystemExit(f"[mvtpy] required {what} not found: {path}\n"
                         f"        set --data-dir/--results-dir or mount the data.")
    return path


class _timer:
    def __init__(self, label: str):
        self.label = label

    def __enter__(self):
        print(f"[mvtpy] {self.label} ...")
        self.start = time.time()
        return self

    def __exit__(self, *exc):
        print(f"[mvtpy] {self.label} done in {time.time() - self.start:.0f}s")


if __name__ == "__main__":
    sys.exit(main())
