"""Command-line entry point for the Python MVT pipeline.

Runs the ported stages against a data tree and writes outputs into a results
tree, so the whole analysis can run from a shell or inside a container:

    python -m mvtpy gps      --day 16
    python -m mvtpy slim     --day 16
    python -m mvtpy samples  --day 16
    python -m mvtpy fields   --day 16
    python -m mvtpy figures  --day 16
    python -m mvtpy micro    --day 16
    python -m mvtpy all      --day 16

Or let it work out what is stale and run it in parallel, like make:

    python -m mvtpy build -j 8              # everything, all three days
    python -m mvtpy build --day 18 -n       # what is stale, and why
    python -m mvtpy build --target slim -j 6

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
* ``micro`` renders the microscopic trajectory time-space figures (stage 5b).
  It is the last step of ``all``, and the slowest on a cold cache: it reads the
  whole slim tree for the day (~15 GB). Later runs reuse
  ``results/.mvt/cache/*_micro.npz``. Both ``all`` and ``build`` include it.
* ``slim`` produces the released segments from raw MOTION (stage 2). It needs
  the assembled GPS first, takes ~97 s and ~5 GB per segment, and writes
  ~15 GB/day. Split it with ``--segment`` / ``--shard k/N``, or let ``build``
  schedule one process per segment.
* ``build`` is the make: it plans the whole graph, rebuilds only what is stale,
  and runs up to ``-j N`` units at once. Because ``slim`` fans out per segment,
  ``-j`` spreads work across days *and* segments without any shard arithmetic.
  Cap ``-j`` by memory, not cores - ``slim`` peaks near 5 GB per process.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import build
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

    stages = {}
    for name, help_text in (
        ("gps", "assemble the control-vehicle GPS JSON"),
        ("slim", "build the released slim segments from raw MOTION"),
        ("samples", "collect distance-to-AV samples (.npz)"),
        ("fields", "accumulate macroscopic fields (.npz)"),
        ("figures", "render field heatmaps and AV fuel curves"),
        ("micro", "render the microscopic trajectory time-space figures"),
        ("all", "gps, slim, samples, fields, figures, micro in order"),
        ("build", "make-style: rebuild only what is stale, with -j"),
        ("verify", "check generated outputs against expected checksums"),
    ):
        stage = sub.add_parser(name, parents=[common], help=help_text)
        stage.add_argument("--day", type=int, choices=DAYS, action="append",
                           help="day to process; repeatable (default: all three)")
        stages[name] = stage

    stages["slim"].add_argument(
        "--segment", type=int, action="append", metavar="N",
        help="segment index 0..23 to build; repeatable (default: all 24)")
    stages["slim"].add_argument(
        "--shard", metavar="K/N",
        help="build every Nth segment starting at K (1-based), as MATLAB's Shard [k N]")

    stages["build"].add_argument(
        "-j", "--jobs", type=int, default=1, metavar="N",
        help="run up to N units at once, across days and segments (default 1)")
    stages["build"].add_argument(
        "--target", action="append", choices=build.STAGES, metavar="STAGE",
        help="stage to build; repeatable (default: all of "
             + ", ".join(build.DEFAULT_TARGETS) + ")")
    stages["build"].add_argument("-n", "--dry-run", action="store_true",
                                 help="report what is stale and why, build nothing")
    stages["build"].add_argument("-B", "--force", action="store_true",
                                 help="rebuild regardless of timestamps")
    stages["build"].add_argument("-k", "--keep-going", action="store_true",
                                 help="carry on with unaffected units after a failure")
    stages["verify"].add_argument(
        "--reference", metavar="DIR",
        help="compare against another results tree (e.g. the MATLAB one) "
             "instead of the stored checksums; explains how files differ")
    stages["verify"].add_argument(
        "--update", action="store_true",
        help="record the current outputs as the expected checksums")
    stages["verify"].add_argument("-v", "--verbose", action="store_true",
                                  help="list every file, not just failures")

    stages["build"].add_argument(
        "--progress", choices=("auto", "bar", "plain"), default="auto",
        help="'bar' draws a live progress block, 'plain' logs one durable line "
             "per event; 'auto' (default) picks bar on a terminal and plain when "
             "redirected. MVT_NO_PROGRESS=1 also forces plain.")

    stages["micro"].add_argument(
        "--skip-t", type=int, default=None, metavar="N",
        help="sub-sample every Nth trajectory point (default 50; 5 is denser/slower)")
    stages["micro"].add_argument(
        "--no-cache", action="store_true",
        help="re-read the slim JSON instead of reusing results/.mvt/cache/*_micro.npz")

    args = parser.parse_args(argv)
    workspace = Workspace.resolve(getattr(args, "data_dir", None),
                                  getattr(args, "results_dir", None))
    days = args.day or list(DAYS)

    print(f"[mvtpy] data={workspace.data_dir}  results={workspace.results_dir}")

    if args.stage == "build":
        from .progress import Progress

        units = build.plan(workspace, days, args.target or build.DEFAULT_TARGETS,
                           force=args.force, log=lambda *_: None)
        reporter = None
        if not args.dry_run:
            interactive = {"auto": None, "bar": True, "plain": False}[args.progress]
            stale = sum(1 for unit in units if getattr(unit, "stale", True))
            reporter = Progress(total=stale, jobs=args.jobs, interactive=interactive,
                                fresh=len(units) - stale)
        return 1 if build.execute(units, jobs=args.jobs, dry_run=args.dry_run,
                                  keep_going=args.keep_going,
                                  progress=reporter) else 0

    if args.stage == "verify":
        return _run_verify(workspace, days, reference=args.reference,
                           update=args.update, verbose=args.verbose)

    runners = {"gps": _run_gps, "samples": _run_samples, "fields": _run_fields,
               "figures": _run_figures, "micro": _run_micro, "all": _run_all}
    for day in days:
        if args.stage == "micro":
            _run_micro(workspace, day, skip_t=args.skip_t, use_cache=not args.no_cache)
        elif args.stage == "slim":
            _run_slim(workspace, day, segments=args.segment, shard=args.shard)
        else:
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
    _write_dataset_info(ws.gps_dir(), "gps (control-vehicle 10 Hz GPS)", day,
                        info_dir=ws.gps_dir())


def _run_slim(ws: Workspace, day: int, segments=None, shard: str = None) -> None:
    """Build the released slim segments for one day (stage 2).

    Selection is by segment index, so a run can be split across processes:
    ``--segment`` picks specific ones, ``--shard k/N`` takes every Nth starting
    at k (1-based), matching the MATLAB ``Shard [k N]`` option. ``mvtpy build``
    schedules one process per segment instead, which needs no sharding at all.
    """
    from . import segments as segment_map, slim

    _require(ws.motion_dir(day), "raw MOTION segments")
    gps_file = ws.gps_dir() / f"CIRCLES_GPS_10Hz_2022-11-{day}.json"
    _require(gps_file, "assembled GPS JSON (run the gps stage first)")
    grade_csv = _require(ws.models_dir() / "Eastbound_grade_fit.csv", "road-grade fit")

    chosen = segment_map.manifest(ws.data_dir, ws.results_dir, day)
    if segments:
        chosen = [segment for segment in chosen if segment.seq in set(segments)]
    if shard:
        index, count = (int(part) for part in shard.split("/"))
        chosen = [segment for segment in chosen if segment.seq % count == (index - 1) % count]
    if not chosen:
        print("    no segments selected")
        return

    out_dir = ws.slim_dir(day)
    out_dir.mkdir(parents=True, exist_ok=True)
    for segment in chosen:
        with _timer(f"slim 2022-11-{day} #{segment.seq:02d} -> {segment.output_name}"):
            records = slim.build_segment(segment.raw_path, gps_file, grade_csv)
            written = slim.write_segment(records, segment.output_path(out_dir))
        print(f"    wrote {written} ({len(records)} records)")
    _write_dataset_info(out_dir, "slim (I-24 MOTION trajectories)", day)


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


def _run_micro(ws: Workspace, day: int, skip_t: int = None,
               use_cache: bool = True) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("    matplotlib not installed; skipping figures "
              "(pip install 'mvtpy[full]')")
        return
    from . import microplot

    slim = _require(ws.slim_dir(day), "processed slim JSON")
    options = microplot.MicroOptions()
    if skip_t:
        options = replace(options, skip_t_plot=skip_t)

    with _timer(f"micro 2022-11-{day}"):
        written = microplot.plot_microscopic_trajectories(
            slim, day, ws.figures_dir(day), options,
            cache_dir=ws.cache_dir(day) if use_cache else None,
            log=print)
    for path in written:
        print(f"    wrote {path}")


def _run_verify(ws: Workspace, days, reference=None, update: bool = False,
                verbose: bool = False) -> int:
    """Check outputs against stored checksums, or against a reference tree."""
    from . import verify

    failures = 0
    for day in days:
        if update:
            written = verify.write_manifest(ws, day, log=print if verbose else None)
            print(f"[mvtpy] wrote {written}")
            continue

        if reference:
            print(f"[mvtpy] verify 2022-11-{day} against {reference}")
            results = verify.verify_against_reference(ws, day, reference)
        else:
            print(f"[mvtpy] verify 2022-11-{day} against stored checksums")
            results = verify.verify_against_manifest(ws, day)
        failures += verify.summarize(results, verbose=verbose)

    if update:
        return 0
    print("[mvtpy] VERIFIED" if failures == 0
          else f"[mvtpy] {failures} file(s) did not match")
    return 1 if failures else 0


def _run_all(ws: Workspace, day: int) -> None:
    _run_gps(ws, day)
    _run_slim(ws, day)
    _run_samples(ws, day)
    _run_fields(ws, day)
    _run_figures(ws, day)
    _run_micro(ws, day)


# ---------------------------------------------------------------------------


def _write_dataset_info(data_dir, product: str, day: int, info_dir=None) -> None:
    """Record what this product folder holds. Never fatal: a stage that has
    already written its outputs must not fail in the bookkeeping, which is
    exactly how the MATLAB `av` stage once died."""
    from . import datasetinfo

    try:
        written = datasetinfo.write(data_dir, product, day, info_dir=info_dir)
        if written:
            print(f"    wrote {written}")
    except Exception as error:                            # noqa: BLE001
        print(f"    could not write dataset_info.json ({error})")


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
