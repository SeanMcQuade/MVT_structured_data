# Python port status

The MATLAB pipeline under `Scripts/` is the reference implementation. `python/`
reproduces it stage by stage under one hard requirement: **the JSON it writes
must be byte-identical to MATLAB's** for the same inputs.

Everything below is verified against data already on disk — the released
`results/` tree is the oracle — so the port can be developed and checked
without a MATLAB license.

## What is ported and verified

| Piece | Module | Verification | Status |
| --- | --- | --- | --- |
| `jsonencode` output format | `mvtpy.matjson` | ~10M numeric literals from `results/{slim,full,gps}`, all 3 days; whole documents re-encoded | byte-identical |
| MATLAB `round(x, 4)` | `mvtpy.matround` | halves away from zero (NumPy rounds half to even) | done |
| Position, speed, acceleration | `mvtpy.kinematics` | recomputed from raw MOTION samples, compared to released values | exact |
| Road grade from the fit map | `mvtpy.kinematics.GradeMap` | ditto | exact |
| Trapezoidal fuel totals, gallons, mpg | `mvtpy.kinematics` | ditto | exact |
| Fuel models (6 vehicles, 2 families) | `mvtpy.fuel` | constants re-parsed from `Models/*.m`; rates compared to released values | exact |
| Streaming reader for raw segments | `mvtpy.rawio` | full 2 GB segment traversed in ~6 s, 272 MB peak RSS | done |
| Lane identification and clipping | `mvtpy.lanes` | released segmentation reproduced from raw data | exact |
| Sample collection (`generate_data_samples`) | `mvtpy.samples` | all 8 arrays, all 88.8M samples identical to the `.mat` | **exact** |
| Macroscopic fields (`generate_macroscopic_fields`) | `mvtpy.fields` | Rho bit-exact; Q/F/U/Phi/Psi to ~1e-11, NaN layout identical | float-exact |
| AV-effect binning (`plot_AV_analysis` core) | `mvtpy.avanalysis` | median and count bit-exact; effective/mean to ~1e-16 | float-exact |
| Figures (field heatmaps, AV fuel curves) | `mvtpy.plotting` | matplotlib, same colors/colormap/layout, not pixel-exact | visual |
| Microscopic trajectories (`plot_microscopic_trajectories`) | `mvtpy.microplot` | same trajectory count and ribbon geometry per segment; one artist per file | visual |

`mvtpy.plotting` renders the two main figure families in matplotlib using
MATLAB's `parula` colormap, the same per-field color-scale limits, the same
engaged-red / disengaged-white AV overlay, the same per-day blue/red/green, and
a matching layout. It is deliberately visually close rather than pixel-for-pixel
(matplotlib cannot reproduce MATLAB's renderer exactly). matplotlib is imported
lazily, so the rest of `mvtpy` does not depend on it.

`mvtpy.microplot` covers stage 5b, the trajectory time-space figures, and
carries the same memory fix as the MATLAB version — one artist per segment file
rather than one per trajectory. See
[MEMORY_FIX_microscopic_trajectories.md](MEMORY_FIX_microscopic_trajectories.md#the-python-port)
for how the two implementations correspond.
| Distance to upstream/downstream AVs | `mvtpy.avdist` | all four distance fields and their vehicle ids, including empty/null handling | exact |
| **Whole segment: assemble, round, encode, write** | `mvtpy.slim` | **full 409 MB released segment reproduced from raw data, md5 identical** | **byte-identical** |

"Exact" means: run the Python code on the raw I-24 MOTION samples, round to
four decimals the way the pipeline does, and every value equals what MATLAB
wrote — with zero differences in any field.

For lane work specifically, the port reproduces **every** released segment of
the trajectories it processes: the same `-N` segment identifiers, the same
sample windows and first/last timestamps, the same corrected lateral position
(`y_position_corrected_meters`), and the same `lane_number`. That includes the
driving-line estimate, which is a whole-file statistic — hence the streaming
reader, which makes two bounded-memory passes instead of decoding 2 GB at once.

`generate_data_mvt_slim` is fully reproduced. Segment 00 of 2022-11-16 was
rebuilt from the raw MOTION file and matched the released output exactly:

```
409,331,812 bytes   md5 3a7bfc02d454c195f2966caa41d6cdff   (MATLAB and mvtpy)
```

That run took 60 s to build 17,058 trajectory records and 51 s to encode and
write, at 5.8 GB peak RSS.

### The rule that decided it

Getting there required one non-obvious discovery about `round(x, 4)`. MATLAB
snaps near-ties: a scaled value within **one ULP** below the midpoint is
treated as a tie and rounded away from zero, compensating for binary
representation error. It is the same behavior that makes `round(2.675, 2)`
return `2.68` in MATLAB.

For the timestamp 1668600000.1863499 the exact scaled value is
16686000001863.4986877 — below the midpoint, so exact decimal rounding gives
`.1863`, while MATLAB writes `.1864`. Before this was understood, a whole
segment differed from the reference by exactly 73 bytes out of 409 MB.

The tolerance was fitted against 80 values evaluated by MATLAB itself
(`matlab_probes/probe_round.m` → `round_probe_2025b.json`) and is bracketed:
0.5 ULP mismatches 40 of them, 2 ULP mismatches 14, 1 ULP matches all 80.
`tests/test_matround.py` asserts both the rule and its bracketing.

### GPS assembly (`assemble_data_GPS.m`)

Ported in `mvtpy.gpsruns` (run splitting) and `mvtpy.gpsassemble` (everything
past it), and verified against the released GPS file for 2022-11-16:

| Piece | Result |
| --- | --- |
| Run splitting (`parse_gps_data`) | structural: 795/795 runs, per-vehicle counts, containment |
| Control-signal reconstruction (`get_control_car_status`) | **bit-exact**: `control_car` and `control_last30` match on every aligned run (566/566) |
| Lane assignment, direction, `controller_engaged` | exact |
| Connection logic (`get_connection_status` + assembly loop) | faithful: bit-identical to MATLAB on identical inputs (198k samples), 99.96% vs the released file |
| 10 Hz resampling of y / lat / long (`preproc_gps`, `sample_10hz`) | > 99% of samples bit-exact; residual ≤ ~1e-2 m |

### Making the GPS floats portable

GPS data rounds to **6 decimals** and is read from **CSV**, 100x more sensitive
than the JSON-sourced MOTION data. Two MATLAB-specific float behaviors had to be
matched exactly; both were run down to their bit-level cause and fixed, so no
change to the MATLAB pipeline or the released data was needed.

1. **CSV parsing.** pandas' default `read_csv` parser is *not* correctly
   rounded — it disagrees with MATLAB's `readtable` by 1 ULP on ~0.035% of
   values (e.g. the string `30.681164000000006`). MATLAB's `readtable` is the
   correctly-rounded one. pandas' `float_precision='round_trip'` mode is also
   correctly rounded and matches `readtable` **bit-for-bit (0 of 200,000 values
   differ)** across every column. The port uses that mode everywhere it reads a
   CSV.

2. **The colon operator.** MATLAB's `a:step:b` does not compute `a + k·step`; it
   builds the vector from both ends (`a + k·step` for the first half,
   `b − (n−k)·step` for the second) to stay accurate at each end. For a POSIX
   timestamp base the two forms differ in the last bit on ~20% of the 10 Hz
   grid points, which shifts the interp1 queries.

   There is a third case the both-ends description omits: when the number of
   steps is **even** there is an exact middle element, and MATLAB sets it to
   **`(a+b)/2`**, not to either one-sided form. This is not a tie-break that can
   be inferred — the correct midpoint equals the forward form on some grids and
   the backward form on others, and the average on all of them (206/206 of the
   discriminating cases). An earlier note here claimed "0 of 2207 grid points
   differ"; that sample happened to contain no discriminating even-step grid.
   Measured over the 772 real 2022-11-18 run grids, the missing midpoint rule
   was worth **104 wrong grid points per day**.

   `mvtpy.gpsassemble._colon` now reproduces all three cases, verified
   bit-for-bit against `lo:0.1:hi` evaluated in MATLAB: **0 of 3,739,194 points
   differ.**

3. **`interp1`.** MATLAB uses the weighted-blend form `A·(1−w) + B·w`, not
   `A + slope·(x−A)`; they differ at ~1e-14, invisible at 4 decimals but
   decisive at 6. Confirmed against `interp1` on the real inputs (2207/2207).

With all three, **latitude and longitude are bit-exact**. `y_position` (which
also passes through the ft2m multiply) and the unrounded `speed` retain a
sub-ULP residual on a small fraction of samples — a remaining float-order
sensitivity in the resample, not a logic or parse difference. Every field's
*logic* is faithful, verified by feeding identical inputs to the MATLAB rule and
getting bit-identical output.

A useful cross-check along the way: **the current MATLAB code reproduces the
released GPS file byte-for-byte** (all 15 fields, 3.6M samples), validating the
refactor for the GPS stage as well.

### End-to-end timing

`mvtpy.gpsassemble.assemble_day` runs the whole GPS stage. For 2022-11-16:

| Phase | Time |
| --- | --- |
| parse runs (CSV, correctly-rounded) | 40 s |
| preprocess (10 Hz resample) | 0.6 s |
| connection status | 0.4 s |
| **matching bias** | **649 s** |
| assemble records | 13 s |
| encode + write 304 MB JSON | 36 s |
| **total** | **~12 min** |

The matching pass is 92% of the run — which is why it was the optimization
target — and MATLAB's whole `assemble_data_GPS` is a comparable ~14 min. The
assembled fields match MATLAB as above: control/lane/direction exact, lat/long
~99.99%, with the timestamp/x/speed residual being the documented CSV-parse and
grid float sensitivity.

### MOTION-matching bias (`x_position`)

Ported in `mvtpy.gpsmatch`. Each run's `x_position` is corrected by a per-run
median offset (`median_xd`) measured against the MOTION trajectories the vehicle
was observed as — median 1.5 m, up to 5.5 m, so it is not negligible.

The two primitives that had to be reverse-engineered from MATLAB were pinned
against MATLAB output:

* **`smoothdata` gaussian** — MATLAB uses standard deviation `window/5`,
  truncates at `window/2`, and normalizes the weights per point. Matches a
  MATLAB delta probe to ~1e-11.
* **Both-direction lane assignment** — unlike the westbound-only `assign_lanes`
  in `mvtpy.lanes`, this builds a driving line per direction. Matches MATLAB's
  `assign_lanes` on a real segment to ~1e-14.

End-to-end over the full day (all 22 relevant MOTION segments), the per-run
`median_xd` reproduces the value recovered from MATLAB's own output on
**564 of 566 runs to within 1e-6 m** (median difference 2.4e-7 m). The two
outliers differ by 0.045 m and 0.006 m — a run-level constant offset whose
median shifts slightly because a few boundary samples flip in or out of a
matched stretch (the sub-ULP resample residual feeding the match threshold),
not an algorithm difference. `x_position` byte-parity still waits on that
residual, but the bias itself — the piece that makes `x_position` correct rather
than 1.5 m off — is faithful.

The hot paths are vectorized — the lane median filter (sliding window), the
`smoothdata` gaussian (searchsorted-banded, fully vectorized), the match-stretch
walk (run-length on a boolean mask), and per-segment caching of each
trajectory's arrays. This took a full day from **~65 min to ~11 min (6.1×)**
with **bit-identical output** (all 795 runs unchanged); each vectorization was
verified against the original scalar form first. The full-day check is opt-in
(`MVT_RUN_SLOW=1`).

## GPS parity: byte-identical as of 2026-07-25

The Python GPS stage reproduces the released 2022-11-18 file **exactly** —
269,121,567 bytes, md5 `04f301cc5009bebaabd4666c563d040f`. It started 9.3 MB
short.

Five distinct causes, every one identified by dumping MATLAB ground truth and
diffing bit-for-bit, not by reading code and reasoning:

| # | Cause | Cost before the fix |
| --- | --- | --- |
| 1 | `controller_engaged` encoded 0/1, not `true`/`false` | 9.8 MB |
| 2 | Run slicing dropped the first sample outside the testbed | 214 of 772 records short |
| 3 | `a:d:b` midpoint for an even step count is `(a+b)/2` | 104 grid points/day |
| 4 | `interp1` returns the endpoint exactly on a flat segment | 48,547 speed values |
| 5 | `smoothdata` omits NaN; the port propagated it | `median_xd` on 5 runs |

**The method matters more than any single fix.** Write the real inputs to a
file, evaluate the MATLAB expression under `-batch`, dump raw doubles, compare
in Python. Reasoning about what MATLAB "should" do produced two confident and
wrong hypotheses for cause 5 alone — each was killed by a controlled experiment
(re-run the matching, count how many of the 772 runs move) before any code
changed. Both moved zero runs.

Causes 1 and 2, in the order they were found:

1. **`controller_engaged` was encoded as 0/1 instead of `true`/`false`.**
   MATLAB carries it as a logical, so `jsonencode` writes booleans. At 3M
   samples this alone was **9.8 MB** of difference. *Fixed* in
   `gpsassemble.assemble_record`; pinned by
   `tests/test_gpsassemble.py::test_controller_engaged_encodes_as_json_booleans`.

2. **Every run was one 10 Hz sample short — a slicing off-by-one.** MATLAB's
   `runEnd` is the first sample *outside* the testbed and `vehTable(1:runEnd,:)`
   **keeps** it, resuming at `runEnd+1`. The port sliced `rows.iloc[:end]`,
   dropping it, and resumed at `end`. That shortened **214 of 772** records
   (211 by one sample, 3 by two) and caused the 466-sample shortfall in
   `samples_for_distance_analysis_18`. *Fixed* in `gpsruns._split_runs`; pinned
   by `test_run_keeps_the_first_sample_outside_the_testbed`.

   The pre-existing run tests could not have caught this: they check run
   *counts* (795/795), not sample counts.

3. **Two floating-point evaluation differences in the resample.**

   * **The colon operator's midpoint.** For an even number of steps MATLAB sets
     the middle grid element to `(a+b)/2`, not to either one-sided form. This is
     not inferable by tie-breaking: the correct midpoint equals the forward form
     on 102 grids and the backward form on 104, and the average on all 206.
     Verified over the 772 real 2022-11-18 run grids: 0 of 3,739,194 points
     differ.
   * **`interp1` on a flat segment.** Where the bracketing values are equal
     (`A == B`), the blend `A*(1-w) + B*w` rounds its two products
     independently and lands a ULP off; MATLAB returns the value itself. Against
     `interp1` on a real 9,449-point run the plain blend differs on 88 values,
     every one a flat segment; with the guard, 0. Three call sites had the bug,
     and `lanes.interp1_linear_extrap` used the slope form `A + slope*(x-A)`
     entirely, which disagrees with MATLAB on ~25% of points.

   Together these took the file from 484 KB off to **30 bytes**, with nine of
   the ten array fields bit-identical and 767 of 772 records exact.

4. **`smoothdata` omits NaN.** In MATLAB a smoothed point is NaN only when its
   *entire* window is NaN; the port propagated NaN from any window member.
   `dist_to_av` is NaN wherever a MOTION trajectory runs past the AV's own time
   range, so the port poisoned a half-window (1.5 s, ~37 samples at 25 Hz) ahead
   of every NaN. That ended matched stretches early and shifted `median_xd` on 5
   of 772 runs by 2.6-21 mm.

   Two plausible causes were eliminated first, each by re-running the matching
   and counting how many runs moved:

   * *`smoothdata` precision.* MATLAB's `exp` genuinely differs from numpy's in
     the last bit - delta probes recover MATLAB's weight matrix, and the kernel
     *ratios* differ by 4.4e-16 while `np.exp` and `math.exp` agree exactly, so
     no rearrangement can reproduce it. Re-running with a differently-rounded
     equivalent kernel moved **0 of 772** runs. The difference is real and
     irrelevant; no MATLAB-side change is needed.
   * *Filter bounds.* MATLAB's `preproc_gps` overwrites `timestamp` with the
     10 Hz grid but leaves `starting_time`/`ending_time` at the raw values, and
     the matching filters on those. The port used grid bounds, which are
     floor/ceil'd outward by up to 0.1 s. A genuine mismatch, now fixed - and it
     also moved **0 of 772** runs.

   The actual cause was found by instrumenting MATLAB: a copy of
   `assemble_data_GPS.m` restricted to the single segment feeding the smallest
   affected run, dumping its matched pool. Python's pool was a strict subset,
   36 values short, and the per-segment breakdown showed one stretch at `n=90`
   against MATLAB's `n=127` - a 37-sample gap, exactly the half-window.

Ruled out along the way, each against MATLAB ground truth: `round(x,6,'decimals')`
(0 of 18,000 differ, including exact `.5` ties), `median` (0 of 3,000, including
even-length averaging), and FMA in the interpolation (every fused variant ~25x
worse). `mean` does differ from MATLAB's at ~1e-15, but it only gates a
`<= 200 m` test.

## Slim parity: the remaining 0.01%

With the GPS stage byte-identical, a three-day `gps`+`slim` rebuild verifies at
**40 of 75 files byte-identical** (all 3 GPS files, 37 of 72 slim segments).

The residual is one field and one cause. Measured on
`I-24MOTION_2022-11-18_06-39-59.json`:

```
trajectories compared: 20,466
with any difference:        2   (0.0098%)
differing field:            total_fuel_consumed_grams   (both records)
example:                    0.2455  vs  0.2456          (0.1 mg)
```

The fuel total is `integrate = @(t,v) dot(t(2:end)-t(1:end-1), (v(1:end-1)+v(2:end))/2)`.
The port evaluates the identical expression, but MATLAB's `dot` is a BLAS call
whose accumulation order differs from numpy's, by ~1e-16. That only matters when
the integral lands within a ULP of a 4-decimal tie, where MATLAB's near-tie snap
(see "The rule that decided it") rounds one way and the port the other.

Reproducing MATLAB's accumulation was attempted and abandoned. Against `dot`
evaluated in MATLAB on 60 real trajectories:

| accumulation | exact |
| --- | --- |
| `np.dot` (OpenBLAS) | 23/60 |
| `np.sum(a*b)` | 35/60 |
| `math.fsum` / sequential Python loop | 41/60 |
| 8 accumulators, FMA, tree combine | **45/60** (best) |

Blocked (64/128/256), unrolled (2/4/8/16), FMA and non-FMA, sequential and tree
combination were all scanned. Nothing reaches 60/60: MATLAB's BLAS is Accelerate,
numpy's is OpenBLAS, and their `ddot` blocking differs. This is not reachable
from pure Python.

Closing it would require changing the MATLAB side — writing `integrate` as an
explicitly ordered accumulation both languages can reproduce. That is a real
option but not a free one: it changes `generate_data_mvt_slim.m`, so the released
slim data (15 GB/day) would no longer match the code that produced it and would
have to be regenerated along with the checksum manifests.

## What is not ported yet

Every stage of the JSON-producing pipeline is now ported. What remains is
optional polish, not new stages:

| Piece | Notes |
| --- | --- |
| Byte-exact GPS output | Three causes remain, measured below. Not "sub-ULP" as previously recorded — one is a whole missing sample per run. |
| `gpsmatch` performance | Correct but slow; needs vectorization for routine use. |
| `.mat` writers and figures (stages 3–6) | Not started; the analysis/plotting half of the pipeline. |
| Samples, macroscopic fields, figures | stages 3-6 | Later; `.mat` writers and matplotlib equivalents. |

## Running the checks

```bash
cd MVT_structured_data/python
python -m pytest tests -q          # all of the above; skips cleanly without data/ and results/

# JSON parity utilities
python tools/json_parity.py ../../results/slim/2022-11-17/I-24MOTION_2022-11-17_07-59-59.json
python tools/json_parity.py --mode file --records 200 ../../results/gps/CIRCLES_GPS_10Hz_2022-11-16.json
```

Two `jsonencode` rules could not be pinned from the released data alone
(magnitudes between 1e6 and 1e9, and values below 1e-4 or non-finite). Run
`python/matlab_probes/probe_jsonencode.m` once in MATLAB and the skipped probe
tests will pick up the result automatically. See `docs/MATLAB_JSON_FORMAT.md`.

## Design rules for the port

1. **Mirror the operation order, not just the formula.** Every function names
   the MATLAB lines it reproduces. Reassociating arithmetic changes the last
   bits and can change a rounded digit.
2. **Round only where MATLAB rounds.** The pipeline rounds to four decimals
   once, immediately before encoding; intermediate values stay full precision.
3. **No new parameters.** Constants come from the MATLAB sources, and
   `tests/test_fuel.py` re-parses `Models/*.m` to prove they still agree.
4. **Compare against released data, not against expectations.** Any new stage
   lands with a test that reproduces real MATLAB output.
