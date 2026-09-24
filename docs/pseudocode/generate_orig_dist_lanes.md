# `generate_orig_dist_lanes.m` — stage `lanes`

Records, for every trajectory in the released data, which lane it came **from**
and which lane it went **to**.

**Reads** `data/i24motion/2022-11-DD/*_0_*.json` (the raw segments)
**Writes** `results/analysis/2022-11-DD/I-24MOTION_<date>_<time>_orig_dist_lane.mat`
(24 per day, 3.6 MB for all three days)
**Runs** once per day, shardable by segment.

## Why this stage exists

The `slim` stage cuts a lane-changing vehicle into one trajectory per lane and
records, for each piece, only the lane it was driven in. The change itself —
that this piece came from lane 3 and the next went to lane 2 — is discarded.
The lane-change analysis needs exactly that, so this stage re-runs the same
clipping on the raw data and keeps the origin and destination.

It reads the **raw** segments, not `slim`, because the information is not in
`slim` to recover.

## Tunables

Identical to the lane parameters of `slim`, and they must stay identical — see
the fork note at the end.

## Main routine

```
GIVEN a day (16, 17 or 18)

STANDARD PREAMBLE (see assemble_data_GPS.md)
    IF the raw folder is absent:
        STOP with an explanation naming what is missing and how to get it
        # This stage is unavailable on any download without the raw data, which
        # is a normal state, not a fault.

segments ← the manifest: raw file → the sidecar name it produces
outputs  ← the 24 declared sidecar names for this day

FOR each segment owned by this shard:        # worker k of N takes k, k+N, ...
    IF the sidecar exists and is newer than the raw file and the code,
       AND Force is not set:
        report the skip and CONTINUE
    IF Clean: delete it.   IF DryRun: CONTINUE

    decode the raw segment
    CHECK the name derived from the decoded data matches the manifest's

    KEEP only westbound trajectories
    DROP the same 9 unused fields slim drops

    dataLanes ← ASSIGN LANES                   → assign_lanes      (as in slim)
    clipped   ← CLIP LANE CHANGES, tracking
                origin and destination         → clip_lane_changes (extended)

    # One entry per clipped trajectory, in the same order slim produces them,
    # freshly allocated for this segment
    sidecar ← empty array of length(clipped)
    FOR i in 1..length(clipped):
        sidecar(i).origin_lane      ← clipped(i).origin_lane
        sidecar(i).destination_lane ← clipped(i).destination_lane

    WRITE sidecar through a temporary name
```

## How origin and destination are tracked

`clip_lane_changes` is the `slim` routine with two values carried through it:

```
WALKING one trajectory, emitting pieces:

    origin_lane ← round(lane) at the first point where the vehicle is
                  holding its lane
                  # For the first piece this is where the vehicle entered the
                  # observed stretch; for every later piece it is the lane of
                  # the piece before it.

    FOR each piece:
        change_at ← the first point where |lane − current_lane| exceeds
                    LaneChangeThresh
        IF there is no such point:
            destination_lane ← current_lane
            # the trajectory ended without changing lane
        ELSE:
            destination_lane ← current_lane + sign(lane(change_at) − current_lane)
            # ±1: which way it moved, not where it eventually ended up

        EMIT the piece with origin_lane and destination_lane
        origin_lane ← current_lane          # for the next piece
```

So a trajectory that never changes lane has
`origin_lane = destination_lane = lane_number`, and one that did differs in one
or both.

## Notes for review

* **The sidecar is paired with the slim JSON by index.** Entry *i* describes
  released trajectory *i*; the file carries no trajectory id of its own. That
  holds only while both were built from the same raw segment with the same
  clipping parameters. The consumer checks the counts match and refuses to run
  otherwise — which catches a coarse mismatch but not a subtle one.

* **This code is a fork.** `assign_lanes` and `clip_lane_changes` are copies of
  the `slim` versions, extended in place. `mvt.sources` cannot see the coupling,
  so editing one copy does **not** invalidate the other's outputs. Keeping them
  in step is currently a matter of discipline; factoring them into one shared
  file is the obvious follow-up.

* **Fixed bug, worth knowing when reading old output.** The accumulator was
  originally not cleared between segments, so each saved sidecar's length was
  the running maximum of the segment trajectory counts — 17 of 24 segments on
  day 17 carried a tail of stale entries. Nothing read past
  `numel(slim)` entries, so results were unaffected, but the files were wrong.
  Sidecars built before 2026-09-19 should be rebuilt rather than trusted.
