#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  process_unrooted.sh -i INPUT -o OUTPUT [root_by_outgroups args...] [-- clean.py args...]

This script runs:
  1) root_by_outgroups.py on INPUT -> temp file
  2) clean.py on temp file -> OUTPUT

Examples:
  process_unrooted.sh -i trees.tre -o trees.rooted.cleaned.tre -ogf outgroups.txt
  process_unrooted.sh -i trees.tre -o out.tre -ogf outgroups.txt -- --deterministic

Notes:
- Arguments before '--' are passed to root_by_outgroups.py
- Arguments after '--' are passed to clean.py
- clean.py always receives -i <tmp> -o <OUTPUT>
USAGE
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

root_args=()
clean_args=()
num_workers=""

# Split args on --
mode="root"
for arg in "$@"; do
  if [[ "$arg" == "--" ]]; then
    mode="clean"
    continue
  fi
  if [[ "$arg" == "--num-workers" ]]; then
    mode_tag="$mode"
    read_next="num_workers"
    continue
  fi
  if [[ "${read_next:-}" == "num_workers" ]]; then
    num_workers="$arg"
    read_next=""
    continue
  fi
  if [[ "$mode" == "root" ]]; then
    root_args+=("$arg")
  else
    clean_args+=("$arg")
  fi
done

# Extract -i/--input and -o/--output from root_args (required)
input=""
output=""
for ((i=0; i<${#root_args[@]}; i++)); do
  case "${root_args[$i]}" in
    -i|--input)
      if (( i + 1 < ${#root_args[@]} )); then
        input="${root_args[$((i+1))]}"
      fi
      ;;
    -o|--output)
      if (( i + 1 < ${#root_args[@]} )); then
        output="${root_args[$((i+1))]}"
      fi
      ;;
  esac
  done

if [[ -z "$input" || -z "$output" ]]; then
  echo "Error: -i/--input and -o/--output are required." >&2
  usage
  exit 1
fi

# Make a temp file in /tmp
if tmpfile=$(mktemp /tmp/rooted.XXXXXX.tre 2>/dev/null); then
  :
else
  tmpfile="/tmp/rooted.$$.$RANDOM.tre"
  : > "$tmpfile"
fi

cleanup() {
  rm -f "$tmpfile"
}
trap cleanup EXIT

# Run rooting
if [[ -n "$num_workers" ]]; then
  python root_by_outgroups.py "${root_args[@]}" -o "$tmpfile" --num-workers "$num_workers"
else
  python root_by_outgroups.py "${root_args[@]}" -o "$tmpfile"
fi

# Run cleaning
if [[ -n "$num_workers" ]]; then
  python clean.py -i "$tmpfile" -o "$output" --num-workers "$num_workers" "${clean_args[@]}"
else
  python clean.py -i "$tmpfile" -o "$output" "${clean_args[@]}"
fi
