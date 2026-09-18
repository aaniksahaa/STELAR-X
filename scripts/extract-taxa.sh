#!/usr/bin/env bash
# Extract union/intersection taxon names using STELAR-X's own Newick tokenizer.

set -euo pipefail

STELARX_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE=""
OUTPUT_FILE=""
SET_MODE="union"
BUILD_FIRST=true

usage() {
  cat <<EOF
Usage: $0 --input <trees.nwk> [--output <taxa.txt>] [--union|--intersection]

Options:
  -i, --input FILE       One or more Newick trees (one tree per non-empty line)
  -o, --output FILE      Output list; stdout when omitted
  --union                Taxa present in any tree (default)
  --intersection         Taxa present in every tree
  --taxa-set MODE        Explicit union | intersection selection
  --no-build             Reuse the existing build directory
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_FILE="$2"; shift 2 ;;
    -o|--output) OUTPUT_FILE="$2"; shift 2 ;;
    --union) SET_MODE="union"; shift ;;
    --intersection) SET_MODE="intersection"; shift ;;
    --taxa-set)
      SET_MODE="$2"
      if [[ "$SET_MODE" != "union" && "$SET_MODE" != "intersection" ]]; then
        echo "Error: --taxa-set expects union or intersection." >&2
        exit 2
      fi
      shift 2
      ;;
    --no-build) BUILD_FIRST=false; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Error: unknown option '$1'." >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: --input is required." >&2
  usage >&2
  exit 2
fi
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "Error: input tree file does not exist: $INPUT_FILE" >&2
  exit 2
fi

INPUT_FILE="$(realpath "$INPUT_FILE")"
if [[ -n "$OUTPUT_FILE" ]]; then
  mkdir -p "$(dirname "$OUTPUT_FILE")"
  OUTPUT_FILE="$(realpath "$OUTPUT_FILE")"
  if [[ "$OUTPUT_FILE" == "$INPUT_FILE" ]]; then
    echo "Error: output taxa file must differ from the input tree file." >&2
    exit 2
  fi
fi

if [[ "$BUILD_FIRST" == true ]]; then
  "$STELARX_ROOT/build.sh"
fi
if [[ ! -f "$STELARX_ROOT/build/stelarx/Main.class" ]]; then
  echo "Error: compiled STELAR-X classes are missing; run ./build.sh first." >&2
  exit 2
fi

ARGS=(--extract-taxa --input "$INPUT_FILE" --taxa-set "$SET_MODE")
if [[ -n "$OUTPUT_FILE" ]]; then
  ARGS+=(--output "$OUTPUT_FILE")
fi

exec java -cp "$STELARX_ROOT/build" stelarx.Main "${ARGS[@]}"
