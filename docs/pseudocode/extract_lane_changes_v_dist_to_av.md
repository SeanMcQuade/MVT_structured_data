# `extract_lane_changes_v_dist_to_av.m` — stage `lc`

Collects every lane change in a day, and how far the nearest control vehicle
was when it happened.

**Reads** `results/slim/2022-11-DD/*.json` and the day's 24 lane sidecars
**Writes** `results/analysis/2022-11-DD/LC_data_DD.mat`
**Runs** once per day. Deliberately **not** sharded: it decodes all 24 segments
and concatenates them, so the whole day is one unit of work.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)

slim_files  ← the day's 24 released segments
lane_files  ← the day's 24 lane sidecars
    # Both from mvt.expectedOutputs, not from the raw manifest: this stage
    # reads neither raw file, and must work on a download that has none.
ERROR unless there are 24 of each
STOP with an explanation if any is missing

FOR each of the 24 segments:
    decode the slim segment
    load its sidecar

    CHECK numel(sidecar) == numel(trajectories), else STOP
        # They are paired by index; a mismatch would silently attach the wrong
        # lane history to the wrong vehicle, which is worse than stopping.

    FOR i in 1..numel(trajectories):
        attach sidecar(i).origin_lane and .destination_lane to trajectory i

    # A trajectory is interesting if it entered its lane from somewhere else,
    # or left it for somewhere else, or both
    changed ← trajectories WHERE origin_lane ≠ lane_number
                              OR destination_lane ≠ lane_number

    FOR each trajectory in changed:

        IF origin_lane ≠ lane_number:                    # it merged IN
            # Measured at the trajectory's FIRST timestamp, which is where the
            # clipping placed the change
            IF a downstream AV was present:
                RECORD a merge-in row:
                    lane_number
                    lane_change_at_start  ← lane_number − origin_lane
                    dist_to_av_at_start   ← distance to the downstream AV
                    av_at_start           ← which AV that was
                    dist_to_eng_av_at_start, eng_av_at_start
                                          ← the same for the nearest ENGAGED AV,
                                            left empty if there was none
                    x_position_at_start, t_at_start
            IF an upstream AV was present:
                RECORD a second merge-in row, identically, against that AV

        IF destination_lane ≠ lane_number:               # it merged OUT
            # Measured at the trajectory's LAST timestamp
            IF a downstream AV was present:
                RECORD a merge-out row:
                    lane_number
                    lane_change_at_end ← destination_lane − lane_number
                    … the same fields, with _at_end
            IF an upstream AV was present:
                RECORD a second merge-out row

    accumulate this segment's rows

all_lane_changes_start ← every merge-in row for the day
all_lane_changes_end   ← every merge-out row for the day
WRITE both through a temporary name
```

## The sign convention

Both `lane_change_*` fields are signed, and the sign says which way the vehicle
moved between lane numbers (1 is leftmost, 4 rightmost):

| Field | Definition | Positive means |
|---|---|---|
| `lane_change_at_start` | `lane_number − origin_lane` | it came from a lane to its left |
| `lane_change_at_end` | `destination_lane − lane_number` | it left for a lane to its right |

Zero cannot occur: a row is only emitted when the two differ.

## Notes for review

* **Rows are not vehicles.** One lane change produces one row *per AV it could
  be measured against* — up to two, for the nearest upstream and the nearest
  downstream. A trajectory that merged in *and* out, with both AVs present at
  both ends, contributes four rows across the two arrays. The counts in
  `LC_data_DD.mat` (roughly 160 000–220 000 per array per day) are row counts.

* **Two distances per row.** `dist_to_av_*` is to the nearest control vehicle
  whether or not its controller was engaged; `dist_to_eng_av_*` is to the
  nearest engaged one and is empty when none was present. The figures use the
  engaged distance.

* **Start and end are not symmetric.** A merge-in is recorded at the first
  timestamp of the clipped piece and a merge-out at its last, because that is
  where `clip_lane_changes` put the transition. For a vehicle that changed lane
  twice, the merge-out of one piece and the merge-in of the next describe the
  *same* physical manoeuvre from either side.

* **Not sharded on purpose.** Sharding would mean either 24 partial files and a
  reduce step, or repeated concatenation; the whole day takes about five
  minutes, which is not worth the machinery.
