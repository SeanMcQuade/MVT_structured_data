#!/usr/bin/env bash
# Check that the MVT data are laid out where the pipeline expects them.
#
# The pipeline assumes this repository (MVT_structured_data) is a sibling of the
# data/ and results/ folders. This script reports what it finds and tells you
# which workflows can run:
#
#   - PLOT the figures                     needs results/analysis and results/gps
#   - REBUILD the analysis files           needs results/slim and results/gps
#   - BOOTSTRAP from raw data              needs data/cars and data/i24motion
#
# Usage:
#   ./check_data.sh
#   MVT_DATA_DIR=/abs/data MVT_RESULTS_DIR=/abs/results ./check_data.sh
#
# Exits 0 if at least one workflow is runnable, 1 otherwise.
#
# (C) 2026 CIRCLES Consortium. BSD-3-Clause.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE="$(cd "$HERE/.." && pwd)"
DATA="${MVT_DATA_DIR:-$WORKSPACE/data}"
RESULTS="${MVT_RESULTS_DIR:-$WORKSPACE/results}"
DAYS=(16 17 18)

green() { printf '  \033[32m✓\033[0m %s\n' "$1"; }
red()   { printf '  \033[31m✗\033[0m %s\n' "$1"; }
info()  { printf '    %s\n' "$1"; }

# count_glob PATTERN -> number of matching files (0 if none)
count_glob() {
  local n
  n=$(compgen -G "$1" 2>/dev/null | wc -l | tr -d ' ')
  echo "${n:-0}"
}

echo "MVT data check"
echo "  repository : $HERE"
echo "  data       : $DATA"
echo "  results    : $RESULTS"
echo

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
echo "Layout"
layout_ok=1
if [ -d "$DATA" ]; then green "data/ exists"; else red "data/ not found (raw inputs live here)"; layout_ok=0; fi
if [ -d "$RESULTS" ]; then green "results/ exists"; else red "results/ not found (outputs and processed inputs live here)"; fi
echo

# ---------------------------------------------------------------------------
# Analysis files (route 1: plot the figures) : results/analysis + results/gps
#
# This route needs no trajectory data at all, which is why it is checked first
# and separately: a download that has it is ready to make most of the figures
# even though results/slim is empty.
# ---------------------------------------------------------------------------
echo "Analysis files (plot the figures)"
figs_ok=1
gps_n=$(count_glob "$RESULTS/gps/CIRCLES_GPS_10Hz_2022-11-*.json")
if [ "$gps_n" -ge 1 ]; then green "results/gps: $gps_n assembled GPS file(s)"; else red "results/gps: no CIRCLES_GPS_10Hz_*.json"; figs_ok=0; fi
for d in "${DAYS[@]}"; do
  have=""; miss=""
  for f in "fields_motion_2022-11-$d.mat" "LC_data_$d.mat" "relspeed_data_$d.mat" \
           "samples_for_distance_analysis_$d.mat"; do
    if [ -f "$RESULTS/analysis/2022-11-$d/$f" ]; then have="$have ${f%%_*}"; else miss="$miss ${f%%_*}"; fi
  done
  if [ -z "$miss" ]; then green "results/analysis/2022-11-$d: all 4 analysis files"
  elif [ -n "$have" ]; then info "results/analysis/2022-11-$d: have$have; missing$miss"
  else red "results/analysis/2022-11-$d: none"; figs_ok=0; fi
done
echo

# ---------------------------------------------------------------------------
# Processed trajectories (route 2: rebuild the analysis files) : results/slim
# ---------------------------------------------------------------------------
echo "Processed trajectories (rebuild the analysis files)"
plot_ok=1
for d in "${DAYS[@]}"; do
  n=$(count_glob "$RESULTS/slim/2022-11-$d/I-24MOTION_*.json")
  if [ "$n" -ge 24 ]; then green "results/slim/2022-11-$d: $n segments"
  elif [ "$n" -ge 1 ]; then red "results/slim/2022-11-$d: $n segments (need 24)"; plot_ok=0
  else red "results/slim/2022-11-$d: none"; plot_ok=0; fi
done
echo

# ---------------------------------------------------------------------------
# Raw data (to bootstrap from scratch) : data/cars + data/i24motion
# ---------------------------------------------------------------------------
echo "Raw data (bootstrap from scratch)"
raw_ok=1
cars_n=$(count_glob "$DATA/cars/cars_gps/circles_v2_1_car*.csv")
if [ "$cars_n" -ge 1 ]; then green "data/cars/cars_gps: $cars_n vehicle CSV(s)"; else red "data/cars/cars_gps: no circles_v2_1_car*.csv"; raw_ok=0; fi
if [ -f "$DATA/cars/cars_vins.csv" ]; then green "data/cars/cars_vins.csv"; else red "data/cars/cars_vins.csv missing"; raw_ok=0; fi
ping_n=$(count_glob "$DATA/cars/veh_ping_202211*.csv")
if [ "$ping_n" -ge 1 ]; then green "data/cars: $ping_n veh_ping file(s)"; else red "data/cars: no veh_ping_202211*.csv"; raw_ok=0; fi
for d in "${DAYS[@]}"; do
  n=$(count_glob "$DATA/i24motion/2022-11-$d/*_0_*.json")
  if [ "$n" -ge 24 ]; then green "data/i24motion/2022-11-$d: $n raw segments"
  elif [ "$n" -ge 1 ]; then red "data/i24motion/2022-11-$d: $n segments (need 24)"; raw_ok=0
  else red "data/i24motion/2022-11-$d: none"; raw_ok=0; fi
done
echo

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
echo "Verdict"
runnable=0
if [ "$figs_ok" -eq 1 ]; then green "Ready to PLOT THE FIGURES (make figures)."; runnable=1; else info "Not ready to plot: analysis files or GPS incomplete (see above)."; fi
if [ "$plot_ok" -eq 1 ]; then green "Ready to REBUILD THE ANALYSIS FILES from the trajectories (make)."; runnable=1; else info "Not ready to rebuild from trajectories: results/slim incomplete (see above)."; fi
if [ "$raw_ok" -eq 1 ]; then green "Ready to REBUILD EVERYTHING from raw data (make rebuild)."; runnable=1; else info "Not ready to rebuild from raw: raw data incomplete (see above)."; fi

# The lane sidecars come from the raw recordings, so a trajectories-only
# download cannot regenerate LC_data and must keep one or the other.
if [ "$plot_ok" -eq 1 ] && [ "$raw_ok" -eq 0 ]; then
  for d in "${DAYS[@]}"; do
    sc=$(count_glob "$RESULTS/analysis/2022-11-$d/I-24MOTION_*_orig_dist_lane.mat")
    if [ ! -f "$RESULTS/analysis/2022-11-$d/LC_data_$d.mat" ] && [ "$sc" -lt 24 ]; then
      info "2022-11-$d: no LC_data and no lane sidecars; the lane-change figures need one of them (they derive from the raw data)."
    fi
  done
fi

if [ "$runnable" -eq 0 ]; then
  echo
  red "No workflow is runnable yet. Download the data into the sibling folders shown above."
  info "Expected layout:"
  info "  $WORKSPACE/"
  info "    MVT_structured_data/   (this repository)"
  info "    data/                  (raw inputs: cars/, i24motion/)"
  info "    results/               (gps/, analysis/, slim/, figures/)"
  exit 1
fi
exit 0
