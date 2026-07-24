#!/usr/bin/env bash
# Check that the MVT data are laid out where the pipeline expects them.
#
# The pipeline assumes this repository (MVT_structured_data) is a sibling of the
# data/ and results/ folders. This script reports what it finds and tells you
# which workflows can run:
#
#   - PLOT / ANALYZE from processed data   needs results/slim and results/gps
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
# Processed data (for plotting / analysis) : results/slim + results/gps
# ---------------------------------------------------------------------------
echo "Processed data (plot / analyze)"
plot_ok=1
gps_n=$(count_glob "$RESULTS/gps/CIRCLES_GPS_10Hz_2022-11-*.json")
if [ "$gps_n" -ge 1 ]; then green "results/gps: $gps_n assembled GPS file(s)"; else red "results/gps: no CIRCLES_GPS_10Hz_*.json"; plot_ok=0; fi
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
if [ "$plot_ok" -eq 1 ]; then green "Ready to PLOT / ANALYZE (make figures / make av-figs, or mvt fields/figures)."; runnable=1; else info "Not ready to plot: processed data incomplete (see above)."; fi
if [ "$raw_ok" -eq 1 ]; then green "Ready to BOOTSTRAP from raw data (make all, or run_all_scripts)."; runnable=1; else info "Not ready to bootstrap: raw data incomplete (see above)."; fi

if [ "$runnable" -eq 0 ]; then
  echo
  red "No workflow is runnable yet. Download the data into the sibling folders shown above."
  info "Expected layout:"
  info "  $WORKSPACE/"
  info "    MVT_structured_data/   (this repository)"
  info "    data/                  (raw inputs: cars/, i24motion/)"
  info "    results/               (processed inputs and outputs: slim/, gps/, figures/)"
  exit 1
fi
exit 0
