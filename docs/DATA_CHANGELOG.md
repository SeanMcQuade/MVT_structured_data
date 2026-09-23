# Data set changelog

Versions of the **derived** MVT data products — `gps`, `slim`, `full` and the
figures built from them. This is not the version of the upstream I-24 MOTION
data, which the observatory versions independently and which nothing in this
repository changes.

The version is declared once, in
[`Scripts/+mvt/dataVersion.m`](../Scripts/+mvt/dataVersion.m), and mirrored by
`mvtpy.DATA_VERSION`; a test fails if the two drift apart.

## Versioning scheme

`MAJOR.MINOR.PATCH`, read from the point of view of someone holding the data:

| part | meaning |
| --- | --- |
| MAJOR | the data means something different; re-read the documentation |
| MINOR | fields added, removed, renamed or re-specified; parsing code may need changing |
| PATCH | same fields, same meaning, same format; only the last decimal of a small number of values can differ |

A PATCH bump still changes checksums. Two copies of the data with different
PATCH versions are *not* expected to be byte-identical, but any analysis built
on one holds for the other.

## Telling which version a folder holds

Each product carries a `dataset_info.json` sidecar naming the data-set version,
the product, the copyright and licence, and how it was produced — MATLAB release
or Python version, platform, host, user, and the code commit. `mvt.build` writes
it after each stage; `mvtpy.datasetinfo.write` is the Python equivalent.

It sits at the **product root**, not beside the data:

```
results/slim/dataset_info.json
results/slim/2022-11-18/I-24MOTION_*.json    <- data only
```

A stray `.json` among trajectory JSON is swallowed by any consumer globbing
`*.json`, which is how the `micro` stage broke once. `gps` is the exception —
its three files have no day folder — and that is safe because every consumer
opens them by exact filename.

It is deliberately **excluded from the checksum manifests**: it mixes
reproducible facts with run provenance, so it is machine-specific by design and
would fail a byte comparison for reasons unrelated to the data.

## 2.2

**New derived products, and a folder move.** Three analyses that were run by
hand became pipeline stages, and their outputs joined the data set:

| Product | Where | Size | Written by |
| --- | --- | --- | --- |
| Lane origin/destination sidecars | `results/analysis/2022-11-DD/I-24MOTION_*_orig_dist_lane.mat` | 24/day, 3.6 MB total | `generate_orig_dist_lanes` |
| Lane-change events | `results/analysis/2022-11-DD/LC_data_DD.mat` | 1/day, 40 MB | `extract_lane_changes_v_dist_to_av` |
| Pooled relative speeds | `results/analysis/2022-11-DD/relspeed_data_DD.mat` | 1/day, 690 MB | `relative_speed_histogram` |

At the same time the two existing analysis `.mat` products moved:

```
results/figures/2022-11-DD/samples_for_distance_analysis_DD.mat
results/figures/2022-11-DD/fields_motion_2022-11-DD.mat
                    ↓
results/analysis/2022-11-DD/...
```

`results/figures/` now holds only rendered figures.

*Why MINOR and not PATCH.* No field of `gps`, `slim` or `full` changes, and
their bytes are unchanged — `make verify` passes against the 2.1.1 checksum
manifests without modification. But two published files moved, so code that
opened them by path has to be updated. Under the scheme above that is a
re-specification of the data's layout, which is a MINOR change, not a PATCH.

*Why the move.* The `.mat` files are *inputs* to the figures; the `.png` and
`.pdf` are the outputs. Keeping them in one folder meant that downloading the
inputs needed to rebuild a figure also meant downloading the figure. Separated,
the figure inputs are about 5 GB against the 51 GB of trajectories, which is
what makes it practical to reproduce most of the paper's figures without the
trajectory data at all. See the download routes in the README.

*What to do if you hold 2.1.1.* Move the two files rather than regenerating
them; the contents are identical.

```bash
mkdir -p results/analysis/2022-11-DD
mv results/figures/2022-11-DD/*.mat results/analysis/2022-11-DD/
```

*Not yet in the Python port.* The five new stages are MATLAB-only, and the port
still writes its `.npz` to `results/figures/`. See `docs/PYTHON_PORT.md`.

## 2.1.1

**Deterministic fuel quadrature.** `mvt.neumaierDot` (compensated summation)
replaces MATLAB's `dot` in the trapezoidal fuel integral of
`generate_data_mvt_slim.m` and `generate_data_mvt_full.m`, controlled by
`flag_deterministic_quadrature` (default `true`).

*Why.* `dot` dispatches to whichever BLAS the platform ships — Accelerate on
Apple silicon, MKL on Windows — and those libraries sum in different orders.
Floating-point addition is not associative, so the same script on the same
inputs produced different `total_fuel_consumed_grams` on a Mac and on a PC. This
was not hypothetical: comparing four independent builds of the identical code,
**14 of 72 `slim` files differed between a Mac and a PC**, and the two archived
PC-era trees disagreed with each other on 8. Every one of those differences was
a single trajectory in a single field, e.g. `0.5729` vs `0.573`. Meanwhile all
three `gps` files were byte-identical across every tree — the divergence was
confined to the one stage that calls `dot`.

*What changed in the data.* Nothing structural: no field added, removed,
renamed or re-specified, and no change of units or meaning. The fuel integral
shifts by at most ~4 × 10⁻¹⁶ relative, which after rounding to 4 decimals moves
roughly **one trajectory in ten thousand** by 0.0001 g (0.1 mg). Two full
segments regenerated byte-for-byte identically under the new code. No aggregate
— fleet fuel consumption, energy savings, per-class averages — changes at any
reported precision.

*What it buys.* The pipeline now returns the same bytes on every machine, which
is what makes the released code and data independently reproducible. The new
method is also *more accurate*, not merely more consistent: on 60 real
trajectories it matched the exactly-rounded sum on all 60, where `dot` matched
41.

Full measurements and method: [REPRODUCIBLE_QUADRATURE.md](REPRODUCIBLE_QUADRATURE.md).

*Regenerating.* Outputs built before this change are specific to the machine
that produced them. Rebuild `slim`/`full` and re-derive the checksum manifests
with `python -m mvtpy verify --update`, so the recorded md5s describe a
platform-independent artifact.

## 2.1

The data set used for the Nature submission.
