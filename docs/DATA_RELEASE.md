# Releasing the data: what to publish, and how to lay it out

The team is releasing three things, each one derived from the one before it:
the raw inputs in `data/`, the processed westbound trajectories in
`results/slim/`, and the analysis intermediates in `results/analysis/`. The
middle two are awkward to classify — each is the *output* of one step and the
*input* to the next. This note recommends how to handle that.

Measured sizes of the current tree:

| | size | files | role |
|---|---|---|---|
| `data/i24motion/` | 55 GB | 72 JSON (24 × 10 min × 3 days) | raw |
| `data/cars/` | 4.5 GB | GPS/CAN CSV, VIN map | raw |
| `results/slim/` | 51 GB | 72 JSON | derived trajectories |
| `results/gps/` | 821 MB | 3 JSON | derived trajectories |
| `results/analysis/` | 4.6 GB | 84 `.mat` | figure inputs |
| `results/figures/` | 208 MB | 51 `.png`/`.pdf`/`.fig` | rendered output |

The layering is the point. `analysis/` is reproducible from `slim/` + `gps/`,
which are reproducible from `data/`, and `figures/` is reproducible from
`analysis/`. Each layer is roughly an order of magnitude smaller than the one
it came from, which is what makes publishing more than one of them worthwhile:
a reader who only wants to redraw a figure needs 4.6 GB, not 60.

## Recommendation

**Publish them as three separate, independently cited datasets, each citing
the one it derives from. Do not move `slim` into `data/`, and do not rename
it.**

Concretely:

* **Dataset 1 — MVT raw inputs.** `data/` as it stands: I-24 MOTION segments,
  on-vehicle GPS/CAN, the VIN map, `data/README.md`. ~60 GB.
* **Dataset 2 — MVT derived trajectories, v2.2.** `slim/` plus `gps/`, each
  with its `dataset_info.json`, plus the checksum manifests. ~52 GB. Its
  metadata cites Dataset 1 as its source and names the pipeline commit that
  produced it.
* **Dataset 3 — MVT analysis intermediates, v2.2.** `results/analysis/`: the
  distance samples, macroscopic fields, lane origin/destination sidecars,
  lane-change events and pooled relative speeds, with its `dataset_info.json`.
  ~4.6 GB. Cites Dataset 2 as its source.

Each gets its own DOI, and each landing page states in one line that it is
derived from the previous one by the code at a specific commit and is
reproducible from it.

### Why Dataset 3 is worth separating rather than folding into Dataset 2

It is a different download for a different reader. Someone reproducing the
paper's figures needs 4.6 GB and no trajectory data at all — `make
figures-from-mat` rebuilds every figure but the microscopic trajectory plots
from `analysis/` and `gps/` alone. Someone doing new trajectory work needs the
51 GB of `slim/` and may never open an `.mat`. Publishing them together forces
the first reader to take twelve times more data than they need, and gives them
no identifier for the thing they actually used.

The general rule this follows: **separate a derived layer when reproducing it is
expensive and its audience is distinct.** `analysis/` qualifies on both counts —
rebuilding it means decoding the whole slim tree, and its audience is the
paper's readers. `figures/` qualifies on neither: it is 208 MB, it rebuilds from
`analysis/` in minutes, and the paper is already the citable artifact for those
images. Publish it inside Dataset 3 for convenience if you like, but do not give
it a DOI of its own.

### The derived data is being released first

The raw inputs are not out yet; the DOI being minted now is for the derived
output. That reverses the natural order and creates one problem worth handling
deliberately rather than discovering later.

**Reserve the raw dataset's DOI now, before publishing the derived one.** Every
major repository (Zenodo, Dryad, Figshare, anything on DataCite) will issue a
draft or reserved DOI that resolves only once the record is published. Reserve
it, put it in Dataset 2's metadata as its source, and publish Dataset 1 against
that same identifier whenever it is ready. This costs one step now and avoids
the alternative, which is either a derived dataset that permanently cites
nothing, or a metadata amendment after the fact.

If a reserved DOI is genuinely unavailable, the fallback is to cite the raw
data by name, extent and checksum manifest in Dataset 2's documentation, state
plainly that it is not yet public and will be released separately, and plan to
add the DOI in a metadata update once it exists.

**Say what is verifiable today, and do not overclaim.** Until the raw inputs
are public, an outside reader cannot regenerate `slim` — they have no inputs to
run the pipeline on. What they *can* do is confirm their download is intact
against the published checksums, and read exactly how it was produced. Anyone
who already has the raw data — the consortium, reviewers under embargo, the
I-24 MOTION observatory — can do the full rerun and get a byte-identical tree.
That is a strong claim and it is true; "fully reproducible by anyone" is not,
yet. Word the landing page for the first and let the second become true when
Dataset 1 publishes.

Note also that the I-24 MOTION segments are the observatory's data, versioned
by them and unmodified by this pipeline (the sidecar records this). Whether the
consortium or the observatory publishes them, and under what terms, is worth
confirming before promising a raw release date.

## Why not fold `slim` into `data/`

The instinct — "the analysis reads it, so it is input data" — conflates two
different roles, and the cost shows up later.

*Provenance becomes invisible.* The single most useful fact about `slim` is that
it is derived, and from what, by which code version. A file sitting in `data/`
asserts the opposite: that it is a primary observation. Six months on, nobody
reading the tree can tell that `slim` was computed while `i24motion` was
recorded — and that is precisely the distinction a reader needs in order to
know what can be re-derived and what cannot.

*The build system stops working.* Staleness in this pipeline is decided by
comparing outputs against their inputs and against the code that produced them
(`Scripts/+mvt/isStale.m`). Something in `data/` is, by construction, a leaf
with no inputs. Moving `slim` there means either it is never rebuilt when the
raw data or the fuel models change, or `data/` is no longer a directory that is
safe to treat as read-only. Both are worse than the current arrangement.

*It hides the reproducibility claim you have already earned.* The interesting
property of this release is that 51 GB of it can be regenerated byte for byte
from the other 59 GB, on two platforms and in two languages. Presenting `slim`
as raw input throws that away. Presenting it as a derived product with a pinned
version, a checksum manifest and a `verify` command makes it the strongest part
of the release.

The role `slim` plays — output of one stage, input of the next — is not unusual
and does not need special handling. It is what an *intermediate published
product* is, and the established practice is to publish it as such rather than
to relabel it by whichever role is convenient.

The same argument applies one layer down: `results/analysis/` is derived from
`slim/`, and putting it anywhere that implies otherwise would lose the same
fact. Publishing it as its own dataset that cites Dataset 2 keeps the chain
explicit at every link.

## Why publish it at all, rather than "just publish the code"

Because regenerating it costs 59 GB of download, a MATLAB or Python
installation, and hours of compute. Publishing `slim` means a reader who only
wants to redo the fuel-savings analysis needs 51 GB and no pipeline at all,
while a reader who wants to check the pipeline can still do so. Both audiences
are served, and the second one can *verify* the first one's copy rather than
taking it on trust.

## What makes the derived release credible

Most of this already exists in the repository; it is worth naming so it is not
lost in the release process.

* **A version on the data, distinct from the code version.** `mvt.dataVersion()`
  → `2.2`, with the MAJOR/MINOR/PATCH meaning spelled out in
  `DATA_CHANGELOG.md`. Two copies of `slim` differing only in the last decimal
  of a few fuel totals are otherwise indistinguishable by inspection; the
  version is how a recipient tells them apart.
* **A sidecar per product.** `dataset_info.json` at each product root records
  the version, the file count, the generating code commit, and the run
  environment — including the `host` and `user` that produced it. Those two
  fields are published intentionally: for a dataset whose selling point is that
  two machines produce identical bytes, knowing which machine produced *this*
  copy is provenance, not incidental metadata.
* **Checksums, published with the data.** `python/expected/checksums-*.json`
  lets anyone confirm their download, and lets anyone who reruns the pipeline
  confirm they got the same answer. This is the difference between claiming
  reproducibility and demonstrating it.
* **The code, tagged at the commit that built the release**, so the sidecar's
  `code_commit` resolves to something permanent.
* **A layered chain that can be checked one link at a time.** A reader with
  Dataset 3 can rebuild the figures and compare them to the paper; a reader with
  Dataset 2 can rebuild Dataset 3 and compare content hashes; a reader with
  Dataset 1 can rebuild Dataset 2 and compare bytes. Each claim is testable on
  its own, without the layers below it.

## Naming

Keep `slim`, but do not let it travel alone. Inside the repository it is
established and the documentation uses it consistently. On the release landing
page it means nothing to an outsider, so title each dataset descriptively and
note the directory name in the description:

| Dataset | Suggested title | Directory |
| --- | --- | --- |
| 1 | MVT raw inputs: I-24 MOTION segments and control-vehicle GPS/CAN | `data/` |
| 2 | MVT processed westbound trajectories with fuel estimates (derived, v2.2) | `results/slim/`, `results/gps/` |
| 3 | MVT analysis intermediates: fuel-distance samples, macroscopic fields, lane changes and relative speeds (derived, v2.2) | `results/analysis/` |

For Dataset 2, note that a `full` variant with eastbound and reference
trajectories exists but is not released. Renaming the directories now would
invalidate every path in the code, the docs and the manifests, to fix a problem
a sentence of description solves.

Dataset 3's title is long because its contents are heterogeneous — five
products with different shapes. Resist compressing it to something like "MVT
analysis data", which tells a reader nothing about whether what they want is
inside.

## Suggested release checklist

Ordered for the derived-data-first release actually planned.

1. Tag the code at the commit that builds the release.
2. Build `slim` + `gps` from a clean tree; run `verify` to confirm 75/75.
3. Build `analysis/` from that tree, then the figures from `analysis/`, so the
   published intermediates are demonstrably the ones the figures came from.
4. Confirm each `dataset_info.json` shows `2.2` and the tagged commit — there is
   one per product folder, including `results/analysis/`.
5. **Reserve** the DOI for Dataset 1 (raw), without publishing it.
6. Publish Dataset 2 (derived trajectories) with the checksum manifests
   alongside the data, citing the code tag and the reserved raw DOI as its
   source. State that the raw inputs are to be released separately.
7. Publish Dataset 3 (analysis intermediates), citing Dataset 2's now-real DOI
   and the same code tag. Its landing page should say which figures it
   reproduces and with which command.
8. Have the paper cite Dataset 3 and the code tag — that is the dataset a reader
   redrawing a figure actually needs — and Dataset 2 as the trajectories behind
   it.
9. When the raw inputs are ready, publish Dataset 1 against the reserved DOI.
   The citation in Dataset 2 begins resolving with no edit required.

Step 5 is the one that is easy to skip and expensive to skip. Once Dataset 2 is
published with no source identifier, adding one later means amending a record
that has already been cited.

Two things to settle before step 6, both institutional rather than technical:

* Whoever mints your DOIs may have a view on granularity — some repositories
  discourage splitting a logical release across records. Worth a short
  conversation before committing to three.
* Express the derivation in the metadata rather than only in prose. DataCite's
  `RelatedIdentifier` carries `IsDerivedFrom` and `IsSourceOf`, which is what
  makes the chain machine-readable; the data-set version (`2.2`) belongs in the
  `Version` field, and tracks something different from the DOI's own versioning.
  Expect the two to diverge, and do not try to keep them in step.
