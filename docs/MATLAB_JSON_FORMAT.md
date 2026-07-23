# How MATLAB writes the released JSON

The released data files are produced by MATLAB's `jsonencode` followed by a raw
`fwrite` (`Scripts/generate_data_mvt_slim.m`, `..._full.m`,
`Scripts/assemble_data_GPS.m`). For the Python port to produce **byte-identical**
files, it has to reproduce that encoder exactly — `json.dumps` does not, and
neither does any "shortest round-trip" printer.

These rules were derived by measuring the released corpus, and are implemented
in `python/mvtpy/matjson.py` and asserted by `python/tests/test_matjson.py`.

## Structure

| Aspect | Behavior |
| --- | --- |
| Whitespace | None at all: `[{"a":1,"b":[1,2]}]` |
| Key order | MATLAB struct field order, preserved verbatim |
| Booleans | `true` / `false` (the GPS file is full of them) |
| `NaN` | `null` (197 in one slim segment) |
| Empty array `[]` | `[]` (15 232 in one slim segment) |
| Strings | Plain; the released corpus contains no escape sequences at all |

## Numbers

**Digit selection — the rule that matters.** MATLAB formats with 15 significant
digits; if that does not round-trip, it uses 17. **It never uses 16.**

This is why Python's `repr` (shortest round-trip) is not a substitute:

| Value | shortest (`repr`) | MATLAB |
| --- | --- | --- |
| 1668599999.9836 | `1668599999.9836` (14 digits) | `1.6685999999836E+9` |
| 1668691841.800001 | `1668691841.800001` (16 digits) | `1.6686918418000009E+9` (17) |

Using `repr` produced 147 763 mismatches in a single GPS file; the 15-then-17
ladder produces zero across every file measured.

**Notation.** Fixed notation, switching to scientific when the decimal exponent
of the leading digit is `< -4` or `>= 9`. Scientific form is
`<mantissa>E<sign><exponent>` with no zero padding: `1.6685999999836E+9`.

**Integral values** print with no decimal point (`840`, `-1`) — but only in
fixed notation. An integral value above the scientific threshold still goes
scientific: 1668600000 prints as `1.6686E+9`, which appears in the data.

**Negative zero** keeps its sign: `-0`. The `road_grade_radians` arrays contain
long runs of it. Note that `json.loads` decodes `-0` to the *integer* `0` and
loses the sign, so `mvtpy.matjson.loads` decodes every number as a float.

## Evidence and its limits

Measured over ~10 million numeric literals from `results/slim`, `results/full`,
and `results/gps` across all three days: **zero mismatches**, and full
documents (200 trajectory records, 3.9 MB) re-encode byte-for-byte.

Two rules are calibrated rather than directly observed, because the released
data contains no such values:

1. **The upper notation threshold.** The corpus pins it between exponent 5
   (`268983.3831`, fixed) and exponent 9 (`1.6686E+9`, scientific); nothing in
   the data falls between 1e6 and 1e9. The implementation uses 9.
2. **The lower notation threshold and non-finite values.** Four-decimal
   rounding means nothing smaller than `0.0001` (exponent −4, fixed) survives,
   and no infinities occur.

`python/matlab_probes/probe_jsonencode.m` probes exactly these cases. Run it
once in MATLAB:

```matlab
cd MVT_structured_data/python/matlab_probes
probe_jsonencode
```

It writes `jsonencode_probe_<release>.json` next to itself;
`python/tests/test_matjson_probe.py` picks the file up automatically and
asserts that the Python encoder agrees. Until then that test module skips.

## Checking parity yourself

```bash
cd MVT_structured_data/python

# every numeric literal in a released file, re-encoded and compared
python tools/json_parity.py ../../results/slim/2022-11-17/I-24MOTION_2022-11-17_07-59-59.json

# whole-document byte comparison on the first 200 trajectory records
python tools/json_parity.py --mode file --records 200 \
    ../../results/slim/2022-11-17/I-24MOTION_2022-11-17_07-59-59.json

# the test suite (skips the parity tests if results/ is not present)
python -m pytest tests -q
```
