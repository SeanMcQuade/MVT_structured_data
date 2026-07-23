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

## What is not ported yet

| Piece | Where it lives in MATLAB | Notes |
| --- | --- | --- |
| Segment assembly and file write | body of `generate_data_mvt_slim.m` | Field order, rounding, `jsonencode`; all the parts exist, they need wiring together. |
| GPS assembly | `assemble_data_GPS.m` | Vehicle CSVs, VIN mapping, ping merge, MOTION matching. |
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
