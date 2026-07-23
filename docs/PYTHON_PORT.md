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

"Exact" means: run the Python code on the raw I-24 MOTION samples, round to
four decimals the way the pipeline does, and every value equals what MATLAB
wrote — with zero differences in any field.

For lane work specifically, the port reproduces **every** released segment of
the trajectories it processes: the same `-N` segment identifiers, the same
sample windows and first/last timestamps, the same corrected lateral position
(`y_position_corrected_meters`), and the same `lane_number`. That includes the
driving-line estimate, which is a whole-file statistic — hence the streaming
reader, which makes two bounded-memory passes instead of decoding 2 GB at once.

## What is not ported yet

| Piece | Where it lives in MATLAB | Notes |
| --- | --- | --- |
| Distance to upstream/downstream AVs | `calculate_distance_to_avs` | Needs the assembled GPS file; the last piece before whole-file byte parity. |
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
