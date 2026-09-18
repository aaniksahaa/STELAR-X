#!/usr/bin/env bash
# Remove all generated SimPhy data and inferred results. The simulation tools
# recreate this directory automatically when it is needed again.

set -euo pipefail

DRY_RUN=false
ASSUME_YES=false

show_usage() {
  cat <<'EOF'
Usage: ./clear-bulk-simulated.sh [options]

Permanently remove the complete SimPhy data directory:

  $PHYLOGENY_DATA_DIR/simphy/data

This includes simulated datasets, inferred trees, logs, statistics, checkpoints,
and every other file stored beneath that directory. Future simulation commands
will recreate the directory when needed.

Options:
  --dry-run       Show the exact target without deleting it
  --yes, -y       Delete without interactive confirmation
  --help, -h      Show this help

Examples:
  ./clear-bulk-simulated.sh --dry-run
  ./clear-bulk-simulated.sh
  ./clear-bulk-simulated.sh --yes
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --yes|-y) ASSUME_YES=true; shift ;;
    --help|-h) show_usage; exit 0 ;;
    *) echo "Error: unknown option: $1" >&2; show_usage >&2; exit 2 ;;
  esac
done

if [[ -z "${PHYLOGENY_DATA_DIR:-}" ]]; then
  echo "Error: PHYLOGENY_DATA_DIR is not set." >&2
  exit 2
fi

case "$PHYLOGENY_DATA_DIR" in
  /*) ;;
  *)
    echo "Error: PHYLOGENY_DATA_DIR must be an absolute path: $PHYLOGENY_DATA_DIR" >&2
    exit 2
    ;;
esac

BASE_DIR="$(realpath -m -- "$PHYLOGENY_DATA_DIR")"
case "$BASE_DIR" in
  /|"${HOME:-/nonexistent}")
    echo "Error: refusing unsafe PHYLOGENY_DATA_DIR: $BASE_DIR" >&2
    exit 2
    ;;
esac

SIMPHY_DIR="${BASE_DIR}/simphy"
TARGET_DIR="${SIMPHY_DIR}/data"

# Refuse intermediate or target symlinks: following one during recursive
# deletion could remove data outside the explicitly displayed directory tree.
if [[ -L "$SIMPHY_DIR" || -L "$TARGET_DIR" ]]; then
  echo "Error: refusing to delete through a symbolic link: $TARGET_DIR" >&2
  exit 3
fi

if [[ -e "$TARGET_DIR" && ! -d "$TARGET_DIR" ]]; then
  echo "Error: expected a directory but found another file type: $TARGET_DIR" >&2
  exit 3
fi

echo "PHYLOGENY_DATA_DIR: $BASE_DIR"
echo "Delete completely:   $TARGET_DIR"

if [[ ! -e "$TARGET_DIR" ]]; then
  echo "SimPhy data directory does not exist; nothing to remove."
  exit 0
fi

if [[ "$DRY_RUN" == true ]]; then
  echo "Dry run only; nothing was removed."
  exit 0
fi

if [[ "$ASSUME_YES" != true ]]; then
  if [[ ! -t 0 ]]; then
    echo "Refusing non-interactive deletion without --yes. Use --dry-run to preview." >&2
    exit 4
  fi
  echo
  echo "WARNING: this permanently removes all simulated data and inferred results."
  read -r -p "Delete '$TARGET_DIR'? [y/N] " reply
  case "${reply,,}" in
    y|yes) ;;
    *) echo "Cancelled; nothing was removed."; exit 0 ;;
  esac
fi

# TARGET_DIR is constructed above from the validated base and fixed suffix;
# keep this final invariant next to the destructive operation.
if [[ "$TARGET_DIR" != "${BASE_DIR}/simphy/data" ]]; then
  echo "Error: internal target validation failed: $TARGET_DIR" >&2
  exit 3
fi

rm -rf -- "$TARGET_DIR"

if [[ -e "$TARGET_DIR" ]]; then
  echo "Error: deletion did not fully remove: $TARGET_DIR" >&2
  exit 1
fi

echo "Removed all bulk simulated data: $TARGET_DIR"

