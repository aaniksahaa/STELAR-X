#!/usr/bin/env bash
# Shared setting-name encoder for the experiment runners.

sanitize_setting_part() {
  local value="$1"
  value="${value// /-}"
  value="${value//\//-}"
  value="${value//:/-}"
  value="${value//=/-}"
  value="${value//,/.-}"
  printf '%s' "$value"
}

canonical_search_space_name() {
  local value="${1,,}"
  value="${value//_/-}"
  case "$value" in
    1|s1|incomplete-local) printf 'S1' ;;
    2|s2|complete-full)    printf 'S2' ;;
    3|s3|exhaustive)       printf 'S3' ;;
    *) sanitize_setting_part "$1" ;;
  esac
}

canonical_intersection_method_name() {
  local value="${1,,}"
  value="${value//_/-}"
  case "$value" in
    1|i1|smaller-side-traversal|smaller-side|smallerside|legacy) printf 'I1' ;;
    2|i2|prefix-sum|prefixsum|prefix)                            printf 'I2' ;;
    3|i3|simple-tree-walk|tree-walk|treewalk|simple)             printf 'I3' ;;
    4|i4|bitset|bitsets|bit-set)                                 printf 'I4' ;;
    *) sanitize_setting_part "$1" ;;
  esac
}

# Encode meaningful options as option_value groups separated by '__'.
# Display-only verbosity flags do not describe an experimental setting and are
# omitted. Friendly and legacy STELAR-X aliases use the compact preset names.
build_setting_name_from_opts() {
  local raw="$1"
  local -a tokens=()
  local -a parts=()
  local i=0 token key value next

  if [[ -z "${raw// }" ]]; then
    printf 'default'
    return
  fi

  read -r -a tokens <<< "$raw"
  while (( i < ${#tokens[@]} )); do
    token="${tokens[$i]}"
    case "$token" in
      -v|-vv|-vvv|-q|--quiet|--verbose)
        ((i+=1))
        continue
        ;;
      --search-space|--intersection-method|--im|--weight-intersection-method)
        if (( i + 1 >= ${#tokens[@]} )); then
          ((i+=1))
          continue
        fi
        value="${tokens[$((i + 1))]}"
        if [[ "$token" == "--search-space" ]]; then
          parts+=("search-space_$(canonical_search_space_name "$value")")
        else
          parts+=("intersection-method_$(canonical_intersection_method_name "$value")")
        fi
        ((i+=2))
        ;;
      --search-space=*|--intersection-method=*|--im=*|--weight-intersection-method=*)
        value="${token#*=}"
        key="${token%%=*}"
        if [[ "$key" == "--search-space" ]]; then
          parts+=("search-space_$(canonical_search_space_name "$value")")
        else
          parts+=("intersection-method_$(canonical_intersection_method_name "$value")")
        fi
        ((i+=1))
        ;;
      -t|-T|--threads|--num-threads)
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("threads_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      -m|--seeds)
        if (( i + 1 < ${#tokens[@]} )); then
          parts+=("seeds_$(sanitize_setting_part "${tokens[$((i + 1))]}")")
          ((i+=2))
        else
          ((i+=1))
        fi
        ;;
      --*=*)
        key="${token%%=*}"
        value="${token#*=}"
        parts+=("$(sanitize_setting_part "${key#--}")_$(sanitize_setting_part "$value")")
        ((i+=1))
        ;;
      --*)
        key="${token#--}"
        if (( i + 1 < ${#tokens[@]} )); then
          next="${tokens[$((i + 1))]}"
        else
          next=""
        fi
        if [[ -n "$next" && ( "$next" != -* || "$next" =~ ^-[0-9] ) ]]; then
          parts+=("$(sanitize_setting_part "$key")_$(sanitize_setting_part "$next")")
          ((i+=2))
        else
          parts+=("$(sanitize_setting_part "$key")_true")
          ((i+=1))
        fi
        ;;
      *)
        # Unknown positional and short-form arguments cannot be interpreted
        # safely without knowing whether they consume another token.
        ((i+=1))
        ;;
    esac
  done

  if [[ ${#parts[@]} -eq 0 ]]; then
    printf 'default'
    return
  fi

  local result="" part
  for part in "${parts[@]}"; do
    if [[ -z "$result" ]]; then
      result="$part"
    else
      result="${result}__${part}"
    fi
  done
  printf '%s' "$result"
}
