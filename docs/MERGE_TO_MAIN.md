# Merging `tooling` into `main`

Written to answer the question "what would a PR from `tooling` into `main`
actually involve?". Every number here was measured against the repository, not
estimated; the commands are given so they can be re-run when the branches move.

## Where the branches stand

```
merge base            0e23c97
main                  14cf0dc     21 commits ahead of the base
tooling               90ccc5c     53 commits ahead of the base
combined diff         121 files changed, +20,488 / -2,371
```

`rc` is **not** a separate merge target. It was already merged into `main` as
pull request #5 (`13d8411`), so `main` contains it. The only `rc` left anywhere
is a stale local branch in the `MVT_structured_data/` working copy; it has never
been pushed. Any documentation still describing `rc` as "the eventual main" is
out of date.

The 21 commits `main` has that `tooling` lacks are small: README edits, comment
spelling fixes, deletion of superseded `plot_AV_stats*.m` variants, licence-file
renames, and two files added back (`Scripts/mvt_v21_plot_motion_trajectories.m`,
`Scripts/mvt_v21_reduce_data.m`). No pipeline logic among them.

## The merge itself is small

A dry-run merge conflicts in exactly two files:

```
$ git merge-tree --write-tree --name-only main tooling
README.md
Scripts/plot_AV_analysis.m
```

* **`Scripts/plot_AV_analysis.m`** — trivial. `main`'s only change since the
  base is commit `14cf0dc`, "Fixed spelling of several comments": `intrest` →
  `interest`, `analyaze` → `analyze`, `partision` → `partition`, and similar.
  Not one line of executable code differs. Resolve by taking `tooling`'s version
  and re-applying the spelling fixes.
* **`README.md`** — both sides rewrote it. Expect to write the merged version by
  hand rather than resolve hunk by hunk; `tooling`'s describes a workflow that
  no longer resembles `main`'s.

Nothing else conflicts, and `tooling` deletes nothing that `main` has. The two
`mvt_v21_*.m` files that `main` added after the base survive the merge untouched
— whether to keep them is a cleanup decision, not a merge problem, since
`plot_microscopic_trajectories.m` and the reduce cache supersede them.

So the mechanical merge is an afternoon. The reason a PR is still a significant
event is the 20,000 lines it brings, and one change to the numbers.

## What the PR actually introduces

Grouped by what a reviewer would need to think about, roughly in descending
order of consequence.

**1. A change to the output data — the only item that alters results.**
`total_fuel_consumed_grams` is now accumulated by compensated (Neumaier)
summation instead of `dot`, behind `flag_deterministic_quadrature` in
`generate_data_mvt_{slim,full}.m`. `dot` dispatches to the platform BLAS, so its
summation order — and therefore its last bits — differed between Accelerate on
macOS and MKL on Windows. Measured effect: of 72 slim files, 14 differed between
a Mac-built and a PC-built tree; every difference was in this one field, in one
trajectory per file, of order 0.1 mg. It is far below any reported figure, and
it is what now makes the two platforms agree exactly. Data version accordingly
goes to **2.1.1** — a patch, because every field, meaning and format is
unchanged. See `REPRODUCIBLE_QUADRATURE.md` and `DATA_CHANGELOG.md`.

**2. A Python implementation of the whole pipeline** (`python/`, ~50 files).
Independent of MATLAB, with its own make-equivalent, sharding and CLI. It
reproduces the MATLAB JSON byte for byte: 75 of 75 comparable files match across
MATLAB/macOS, MATLAB/Windows and Python/macOS. Figures are not pixel-exact and
are not claimed to be, `.mat` outputs are written as `.npz`, and `full` has no
checksum manifest yet.

**3. Verifiable reproducibility.** `python/expected/checksums-2022-11-{16,17,18}.json`
pin the expected output of every stage; `mvt.verify` (MATLAB) and `mvt verify`
(Python) check a built tree against them. This is what turns "it should
reproduce" into a command anyone can run — and it is what caught the
platform-dependence above.

**4. A build system** (`Scripts/+mvt/`, `Makefile`, `Scripts/make.m`, 27 new
package functions). Make-style staleness — rebuild when outputs are missing,
older than inputs, or older than the code that produced them — replacing
skip-if-exists. Plus process-level parallelism across days and shards within a
day, progress reporting, and timestamped run logs. `make.m` gives Windows users
the same commands without needing GNU make.

**5. Robustness fixes** that stand on their own: the memory blow-up in
`plot_microscopic_trajectories` (see `MEMORY_FIX_microscopic_trajectories.md`),
and the removal of eight hardcoded paths across four scripts that silently
pointed at the wrong tree when the workspace was relocated.

**6. Documentation** (`docs/`, 7 files): algorithms, data dictionary, the
MATLAB JSON number format, the port, the quadrature change, the data changelog.

Crucially for anyone who does not want to adopt any of this: **every original
entry point still exists and still works**. `run_all_scripts`, the six stage
scripts, and their `(processingDay)` signatures are unchanged. The build system
and Python are additive.

## How to make it reviewable

A single 20,000-line PR will not get a real review. The history supports
splitting it into stacked PRs that each stand alone:

1. Path handling and the memory fix — pure bug fixes, uncontroversial, could
   merge today.
2. The build system (`+mvt`, `Makefile`, `make.m`) — additive; no output changes.
3. The quadrature change and version bump — small diff, but the one that needs
   actual scientific sign-off. Give it its own PR so the discussion is not
   buried.
4. Checksum manifests and `verify`.
5. The Python port — largest by line count, lowest risk, since nothing else
   depends on it.
6. Documentation and README reconciliation.

If the group would rather do one merge, item 3 should still be reviewed on its
own terms before the PR opens.

## Questions to expect, and the short answers

* *Does this change our published numbers?* One field, by ~0.1 mg, in 14 of 72
  files, and only where the platforms previously disagreed with each other. No
  figure or reported statistic moves.
* *Do I have to learn a new workflow?* No. The old scripts are untouched.
* *Do I need Python?* No. It is a second implementation used as a cross-check.
* *Why version 2.1.1 and not 2.2?* The data change is smaller than a minor bump
  would imply — no field added, removed, or redefined.
* *Is `main` or `rc` the target?* `main`. `rc` is already in it.

## Open decisions to settle before merging

* ~~**`dataset_info.json` records `host` and `user`.**~~ **Settled: they stay,
  and are published.** For a dataset whose central claim is that two machines
  produce identical bytes, recording which machine produced a given copy is
  provenance rather than incidental metadata.
* **Do the `mvt_v21_*.m` scripts stay?** They are superseded but harmless.
* **Does the Python port live in this repository or its own?** It roughly
  doubles the repository's file count and has a different review audience.
