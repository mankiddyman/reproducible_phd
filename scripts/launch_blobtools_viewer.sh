#!/bin/bash
# Launch the BlobToolKit interactive viewer locally.
# Pass the blobdir path as arg, or run with no args to be prompted.
set -euo pipefail

BIN_DIR="$HOME/bin/blobtoolkit"
API_BIN="$BIN_DIR/blobtoolkit-api"
VIEWER_BIN="$BIN_DIR/blobtoolkit-viewer"

if [ ! -x "$API_BIN" ] || [ ! -x "$VIEWER_BIN" ]; then
  echo "Binaries not found at $BIN_DIR — install them per:"
  echo "  https://github.com/genomehubs/blobtoolkit/wiki/Installation"
  exit 1
fi

# The BTK_FILE_PATH must be a directory CONTAINING the blobdir(s) — viewer scans children
if [ "${1:-}" ]; then
  TARGET="$1"
else
  TARGET="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Drosera_paradoxa/blobtoolkit/initial/hap1/output/blobtoolkit"
fi

# BTK_FILE_PATH should be the parent dir that contains all blobdirs we want browsable
# Best: use the repo root and let viewer find all blobdirs
TARGET_DIR="$(dirname "$TARGET")"

echo "Starting API on port 8880, viewer on port 8881"
echo "Serving blobdirs from: $TARGET_DIR"
echo
echo "Open in browser: http://localhost:8881/view/all"
echo "Press Ctrl-C in this terminal to stop."

# Run API in background, viewer in foreground
BTK_API_PORT=8880 BTK_PORT=8881 BTK_FILE_PATH="$TARGET_DIR" "$API_BIN" &
API_PID=$!
trap "kill $API_PID 2>/dev/null" EXIT

sleep 2
BTK_API_PORT=8880 BTK_PORT=8881 "$VIEWER_BIN"
