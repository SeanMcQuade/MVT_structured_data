#!/usr/bin/env bash
# Create a virtual environment and install the Python MVT pipeline into it.
#
#   ./setup_venv.sh            # core install (numpy, pandas, scipy, ...)
#   ./setup_venv.sh --full     # also matplotlib + h5py (figures, .mat reading)
#   ./setup_venv.sh --dev      # also pytest (to run the test suite)
#
# Afterward:
#   source .venv/bin/activate
#   python -m mvtpy --help
#
# (C) 2026 CIRCLES Consortium. BSD-3-Clause.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PYTHON="${PYTHON:-python3}"
VENV="${VENV:-.venv}"

extras=""
for arg in "$@"; do
  case "$arg" in
    --full) extras="${extras:+$extras,}full" ;;
    --dev)  extras="${extras:+$extras,}dev" ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

echo "==> creating virtual environment in $VENV (using $PYTHON)"
"$PYTHON" -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

echo "==> upgrading pip"
python -m pip install --upgrade pip >/dev/null

if [ -n "$extras" ]; then
  echo "==> installing mvtpy[$extras] (editable)"
  python -m pip install -e ".[$extras]"
else
  echo "==> installing mvtpy (editable, core dependencies)"
  python -m pip install -e .
fi

echo
echo "Done. Activate with:  source $VENV/bin/activate"
echo "Then run:            python -m mvtpy --help"
