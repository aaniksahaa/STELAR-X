#!/usr/bin/env bash

# Locate a Python interpreter that can run the Hugging Face transfer tools.
#
# The upload helper (~/utils/hf-data-transfer/hf_upload.py) needs the
# huggingface_hub package. That package usually lives in the conda base
# environment, while an activated project .venv shadows it under "python3".
#
#   stelarx_find_hf_python [REQUESTED]
#
# With REQUESTED set, that interpreter is validated and returned. Otherwise the
# usual candidates are probed in order and the first one that imports
# huggingface_hub is printed. Fails with an actionable message when none works.
stelarx_find_hf_python() {
  local requested="${1:-}"
  local -a candidates=()
  local candidate resolved tried=""

  if [[ -n "$requested" ]]; then
    if [[ "$requested" == */* ]]; then
      resolved="$requested"
    else
      resolved="$(command -v "$requested" 2>/dev/null || true)"
    fi
    if [[ -z "$resolved" || ! -x "$resolved" ]]; then
      echo "Error: Python interpreter is not executable: $requested" >&2
      return 2
    fi
    if ! "$resolved" -c 'import huggingface_hub' >/dev/null 2>&1; then
      echo "Error: $resolved cannot import huggingface_hub." >&2
      echo "  Install it there (pip install huggingface_hub hf_xet) or pass --python with an interpreter that has it." >&2
      return 2
    fi
    printf '%s\n' "$resolved"
    return 0
  fi

  candidates+=(python3 python)
  [[ -n "${CONDA_EXE:-}" ]] && candidates+=("${CONDA_EXE%/bin/conda}/bin/python")
  [[ -n "${CONDA_PREFIX:-}" ]] && candidates+=("${CONDA_PREFIX}/bin/python")
  candidates+=("${HOME}/miniconda3/bin/python" "${HOME}/anaconda3/bin/python" "${HOME}/mambaforge/bin/python" /usr/bin/python3)

  for candidate in "${candidates[@]}"; do
    if [[ "$candidate" == */* ]]; then
      resolved="$candidate"
    else
      resolved="$(command -v "$candidate" 2>/dev/null || true)"
    fi
    [[ -n "$resolved" && -x "$resolved" ]] || continue
    tried+="${tried:+, }${resolved}"
    if "$resolved" -c 'import huggingface_hub' >/dev/null 2>&1; then
      printf '%s\n' "$resolved"
      return 0
    fi
  done

  echo "Error: no Python interpreter with huggingface_hub was found." >&2
  echo "  Tried: ${tried:-none}" >&2
  echo "  Install it (pip install huggingface_hub hf_xet), deactivate the .venv, or pass --python /path/to/python." >&2
  return 2
}
