#!/usr/bin/env bash
# Verifies the reproducibility mirror of A10K (10k-astral-dataset) run outputs:
#   * default outputs-directory derivation and containment guards,
#   * layout parsing of <data>/10k-simphy/<R>/<method>_outputs/<tree-type>/<setting>,
#   * back-fill through sync-a10k-outputs.sh (dataset record and input
#     fingerprints written, gene trees / species trees never copied),
#   * mirror leaves equal to the source leaves, stale files removed,
#   * upload-a10k-outputs.sh planning, remote paths, filters, and refusals.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${ROOT}/scripts/a10k-outputs-dir.sh"

TMP="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-a10k-mirror-test.XXXXXX")"
trap 'status=$?; rm -rf -- "$TMP"; exit "$status"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# ---------------------------------------------- outputs-directory resolution ---
BASE="${TMP}/phylogeny root"
DATA="${BASE}/10k-astral-dataset"
OUTPUTS="${BASE}/outputs/10k-astral-dataset"
mkdir -p "${DATA}/10k-simphy"

[[ "$(stelarx_default_a10k_outputs_dir "$DATA")" == "$OUTPUTS" ]] || \
  fail "the dataset should mirror into its 'outputs' sibling"
[[ "$(stelarx_default_a10k_outputs_dir "${TMP}/x/outputs/ds")" == "${TMP}/x/outputs/ds_outputs" ]] || \
  fail "a dataset already under an outputs tree should mirror beside itself"
RESOLVED="$(stelarx_prepare_a10k_outputs_dir "" "$DATA")"
[[ "$RESOLVED" == "$OUTPUTS" && -d "$OUTPUTS" ]] || fail "default outputs dir was not created: $RESOLVED"
if stelarx_prepare_a10k_outputs_dir "${DATA}/outputs" "$DATA" >/dev/null 2>"${TMP}/inside.err"; then
  fail "an outputs dir inside the data dir was accepted"
fi
grep -q "inside it" "${TMP}/inside.err" || fail "containment error was unclear"
if stelarx_prepare_a10k_outputs_dir "$DATA" "$DATA" >/dev/null 2>/dev/null; then
  fail "the data dir itself was accepted as the mirror root"
fi
if [[ "$(PHYLOGENY_DATA_DIR="$BASE" stelarx_prepare_a10k_outputs_dir_standalone "")" != "$OUTPUTS" ]]; then
  fail "standalone resolver did not use \$PHYLOGENY_DATA_DIR/outputs/10k-astral-dataset"
fi

# ------------------------------------------------------- synthetic dataset ---
make_replicate_inputs() {
  local replicate="$1" dir="${DATA}/10k-simphy/$1"
  mkdir -p "${dir}/estimatedgenetrees"
  printf '((a,b),(c,d));\n' > "${dir}/s_tree.trees"
  printf '((a,b),(c,d));\n' > "${dir}/truegenetrees"
  printf '((a,b),(c,d));\n' > "${dir}/estimatedgenetrees/estimatedgenetrees.tre"
  printf '((a,b),(c,d));\n' > "${dir}/estimatedgenetrees/estimatedgenetrees.rooted.tre"
}
make_results() {
  local results="$1" tag="$2"
  mkdir -p "$results"
  printf '((a,b),(c,d));\n' > "${results}/out-stelarx.tre"
  printf '# STELAR-X run command\ncd /x && ./run.sh --input i --output o\n' > "${results}/out-stelarx.command"
  printf 'algorithm\nstelar-x\n' > "${results}/out-stelarx_stats.csv"
  printf 'alg,setting\nstelarx,%s\n' "$tag" > "${results}/stat-stelarx.csv"
  printf 'log %s\n' "$tag" > "${results}/.stelarx_run.log"
}

for r in R1 R2; do
  make_replicate_inputs "$r"
  make_results "${DATA}/10k-simphy/${r}/stelarx_outputs/estimated/search-mode_full" "${r}-est-full"
done
make_results "${DATA}/10k-simphy/R1/stelarx_outputs/true/search-mode_full" "R1-true-full"
make_results "${DATA}/10k-simphy/R1/stelarx_outputs/estimated/search-mode_local" "R1-est-local"
make_results "${DATA}/10k-simphy/R1/aster_outputs/estimated/default" "R1-aster"

# ------------------------------------------------------ layout components ---
mapfile -t PARTS < <(stelarx_a10k_results_components "$DATA" "${DATA}/10k-simphy/R1/stelarx_outputs/estimated/search-mode_full")
[[ "${PARTS[0]}" == "R1" && "${PARTS[1]}" == "stelarx_outputs" && "${PARTS[2]}" == "estimated" && "${PARTS[3]}" == "search-mode_full" ]] || \
  fail "results components were parsed wrong: ${PARTS[*]}"
if stelarx_a10k_results_components "$DATA" "${DATA}/10k-simphy/R1/stelarx_outputs/bogus/s" >/dev/null 2>&1; then
  fail "an invalid tree type was accepted"
fi
if stelarx_a10k_results_components "$DATA" "${DATA}/10k-simphy/R1/stelarx_outputs/estimated" >/dev/null 2>&1; then
  fail "a results dir without a setting level was accepted"
fi
[[ "$(stelarx_list_a10k_results_dirs "$DATA" | tr '\0' '\n' | grep -c .)" == "5" ]] || \
  fail "expected 5 results directories in the data tree"

# ------------------------------------------------------------- back-fill ---
SYNC="${ROOT}/sync-a10k-outputs.sh"
"$SYNC" --data-dir "$DATA" --dry-run >"${TMP}/sync-dry.out" 2>&1 || fail "dry run failed: $(cat "${TMP}/sync-dry.out")"
grep -q "would mirror=5" "${TMP}/sync-dry.out" || fail "dry run did not plan 5 leaves: $(cat "${TMP}/sync-dry.out")"
[[ ! -e "${OUTPUTS}/stelarx_outputs/R1" ]] || fail "dry run wrote to the mirror"

"$SYNC" --data-dir "$DATA" >"${TMP}/sync.out" 2>&1 || fail "sync failed: $(cat "${TMP}/sync.out")"
grep -q "mirrored=5 " "${TMP}/sync.out" || fail "sync did not mirror 5 leaves: $(cat "${TMP}/sync.out")"

LEAF="${OUTPUTS}/stelarx_outputs/R1/estimated/search-mode_full"
[[ -s "${LEAF}/out-stelarx.tre" && -s "${LEAF}/stat-stelarx.csv" && -s "${LEAF}/out-stelarx.command" && -s "${LEAF}/.stelarx_run.log" ]] || \
  fail "the mirror leaf is missing run artifacts"
diff -r "${DATA}/10k-simphy/R1/stelarx_outputs/estimated/search-mode_full" "$LEAF" >/dev/null || \
  fail "the mirror leaf does not equal the results leaf"
[[ -s "${OUTPUTS}/stelarx_outputs/R1/true/search-mode_full/out-stelarx.tre" ]] || fail "the true tree type was not mirrored"
[[ -s "${OUTPUTS}/aster_outputs/R1/estimated/default/out-stelarx.tre" ]] || fail "a second method was not mirrored"

# The dataset record and input fingerprints make the mirror reproducible.
RECORD="${OUTPUTS}/stelarx_outputs/${STELARX_A10K_DATASET_RECORD}"
[[ -s "$RECORD" ]] || fail "the dataset record was not written"
grep -q "replicates:   2 (R1 R2)" "$RECORD" || fail "the dataset record does not list the replicates: $(cat "$RECORD")"
grep -q "data_dir:     ${DATA}" "$RECORD" || fail "the dataset record does not name the data directory"
INPUTS="${OUTPUTS}/stelarx_outputs/R1/${STELARX_A10K_INPUTS_RECORD}"
[[ -s "$INPUTS" ]] || fail "the input fingerprints were not written"
grep -qP "10k-simphy/R1/truegenetrees\t15\t" "$INPUTS" || fail "input sizes were not recorded: $(cat "$INPUTS")"
grep -q "estimatedgenetrees.rooted.tre" "$INPUTS" || fail "the rooted gene trees were not recorded"

# Never the inputs themselves.
[[ -z "$(stelarx_a10k_find_forbidden_in_mirror "$OUTPUTS")" ]] || \
  fail "the mirror contains A10K input data: $(stelarx_a10k_find_forbidden_in_mirror "$OUTPUTS")"
[[ -f "${DATA}/10k-simphy/R1/truegenetrees" ]] || fail "the data tree changed"

# Records are stable: a second sync rewrites neither of them.
RECORD_STAMP="$(stat -c '%Y %s' "$RECORD")"
INPUTS_STAMP="$(stat -c '%Y %s' "$INPUTS")"
sleep 1
"$SYNC" --data-dir "$DATA" --quiet >/dev/null 2>&1 || fail "second sync failed"
[[ "$(stat -c '%Y %s' "$RECORD")" == "$RECORD_STAMP" ]] || fail "the dataset record was rewritten without a change"
[[ "$(stat -c '%Y %s' "$INPUTS")" == "$INPUTS_STAMP" ]] || fail "the input fingerprints were rewritten without a change"

# A stale file in the mirror leaf does not survive a refresh.
printf 'stale\n' > "${LEAF}/leftover.txt"
"$SYNC" --data-dir "$DATA" --quiet >/dev/null 2>&1 || fail "refresh sync failed"
[[ ! -e "${LEAF}/leftover.txt" ]] || fail "a stale file survived in the mirror leaf"

# --methods restricts the mirror.
rm -rf "$OUTPUTS"
"$SYNC" --data-dir "$DATA" --methods stelarx --quiet >"${TMP}/sync-filter.out" 2>&1 || fail "filtered sync failed"
grep -q "mirrored=4 filtered-out=1" "${TMP}/sync-filter.out" || fail "--methods did not filter: $(cat "${TMP}/sync-filter.out")"
[[ ! -e "${OUTPUTS}/aster_outputs" ]] || fail "--methods still mirrored another method"
"$SYNC" --data-dir "$DATA" --quiet >/dev/null 2>&1 || fail "restoring sync failed"

# Input data inside a results directory is refused rather than copied.
BAD="${DATA}/10k-simphy/R2/stelarx_outputs/estimated/search-mode_full"
cp "${DATA}/10k-simphy/R2/s_tree.trees" "${BAD}/s_tree.trees"
if stelarx_mirror_a10k_results "$DATA" "$OUTPUTS" "$BAD" >/dev/null 2>"${TMP}/forbidden.err"; then
  fail "a results directory containing input data was mirrored"
fi
grep -q "A10K input data" "${TMP}/forbidden.err" || fail "the refusal message was unclear"
rm -f "${BAD}/s_tree.trees"

# ------------------------------------------------------ uploader dry run ---
UP="${ROOT}/upload-a10k-outputs.sh"
"$SYNC" --data-dir "$DATA" --quiet >/dev/null 2>&1 || fail "pre-upload sync failed"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true >"${TMP}/up.out" 2>&1 || \
  fail "uploader dry run failed: $(cat "${TMP}/up.out")"
grep -q "Plan: upload 2 replicate" "${TMP}/up.out" || fail "uploader did not plan 2 replicate uploads: $(cat "${TMP}/up.out")"
grep -q "plus 1 dataset record" "${TMP}/up.out" || fail "uploader did not plan the dataset record: $(cat "${TMP}/up.out")"
grep -q "ph/d/a10k/outputs/stelarx_outputs/R1/" "${TMP}/up.out" || fail "uploader used an unexpected remote path"
grep -q "aster" "${TMP}/up.out" && fail "the default method filter did not exclude aster"
grep -q "R1-R2 / estimated / search-mode_full" "${TMP}/up.out" || fail "uploader did not collapse replicate ranges: $(cat "${TMP}/up.out")"

# --method all widens the selection; --replicates narrows it.
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true >"${TMP}/up-all.out" 2>&1 || \
  fail "uploader --method all failed"
grep -q "Plan: upload 3 replicate" "${TMP}/up-all.out" || fail "--method all did not add the aster replicate: $(cat "${TMP}/up-all.out")"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --replicates 2 --uploader /bin/true --python /bin/true >"${TMP}/up-r2.out" 2>&1 || \
  fail "uploader --replicates failed"
grep -q "Plan: upload 1 replicate" "${TMP}/up-r2.out" || fail "--replicates did not narrow the plan: $(cat "${TMP}/up-r2.out")"

# A mirror without its dataset record is refused unless explicitly allowed.
mv "$RECORD" "${TMP}/record.saved"
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true >"${TMP}/up-norec.out" 2>&1; then
  fail "a mirror without the dataset record was uploaded"
fi
grep -q "BLOCKED: missing ${STELARX_A10K_DATASET_RECORD}" "${TMP}/up-norec.out" || \
  fail "the missing-record refusal was unclear: $(cat "${TMP}/up-norec.out")"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --allow-missing-command --uploader /bin/true --python /bin/true >"${TMP}/up-allow.out" 2>&1 || \
  fail "--allow-missing-command did not permit the upload"
mv "${TMP}/record.saved" "$RECORD"

# Input data that appeared in the mirror blocks the upload.
cp "${DATA}/10k-simphy/R1/truegenetrees" "${OUTPUTS}/stelarx_outputs/R1/estimated/search-mode_full/truegenetrees"
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true >"${TMP}/up-bad.out" 2>&1; then
  fail "a mirror containing input data was uploaded"
fi
grep -q "BLOCKED: contains truegenetrees" "${TMP}/up-bad.out" || fail "the forbidden-file refusal was unclear: $(cat "${TMP}/up-bad.out")"
rm -f "${OUTPUTS}/stelarx_outputs/R1/estimated/search-mode_full/truegenetrees"

# Input data directly under <method>_outputs/ is uploaded as an individual file,
# so it must be caught by the same guard.
cp "${DATA}/10k-simphy/R1/truegenetrees" "${OUTPUTS}/stelarx_outputs/truegenetrees"
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true >"${TMP}/up-top.out" 2>&1; then
  fail "input data at the method top level was uploaded"
fi
grep -q "BLOCKED: contains truegenetrees" "${TMP}/up-top.out" || \
  fail "the top-level forbidden-file refusal was unclear: $(cat "${TMP}/up-top.out")"
rm -f "${OUTPUTS}/stelarx_outputs/truegenetrees"

# ------------------------------------------------- real STELAR-X run mirror ---
RUN_DATA="${TMP}/run/10k-astral-dataset"
RUN_OUTPUTS="${TMP}/run/outputs/10k-astral-dataset"
mkdir -p "${RUN_DATA}/10k-simphy/R1"
printf '((a,b),(c,d));\n' > "${RUN_DATA}/10k-simphy/R1/truegenetrees"
printf '((a,b),(c,d));\n' > "${RUN_DATA}/10k-simphy/R1/s_tree.trees"

"${ROOT}/run-a10k.sh" --data-dir "$RUN_DATA" --tree-type true --replicates R1 \
  --opts '--search-space S1 --cpu -q' --no-time-monitor --no-gpu-monitor --no-notify \
  >"${TMP}/run1.out" 2>&1 || fail "run-a10k.sh failed: $(tail -20 "${TMP}/run1.out")"

RUN_LEAF="${RUN_OUTPUTS}/stelarx_outputs/R1/true/search-space_S1__cpu_true"
grep -q "outputs mirror: ${RUN_OUTPUTS}" "${TMP}/run1.out" ||   fail "run did not report the default outputs mirror: $(grep -i mirror "${TMP}/run1.out")"
grep -q "Mirrored outputs to: ${RUN_LEAF}" "${TMP}/run1.out" || fail "run did not mirror its outputs"
RUN_SRC="${RUN_DATA}/10k-simphy/R1/stelarx_outputs/true/search-space_S1__cpu_true"
[[ -s "${RUN_LEAF}/out-stelarx.tre" && -s "${RUN_LEAF}/stat-stelarx.csv" && -s "${RUN_LEAF}/out-stelarx_stats.csv" ]] ||   fail "mirrored run is missing tree or CSVs"
[[ -s "${RUN_LEAF}/.stelarx_run.log" ]] || fail "mirrored run lacks the run log"
diff -r "$RUN_SRC" "$RUN_LEAF" >/dev/null || fail "mirror leaf differs from the results dir"
[[ -z "$(stelarx_a10k_find_forbidden_in_mirror "$RUN_OUTPUTS")" ]] || fail "run mirror contains A10K input data"
[[ -s "${RUN_OUTPUTS}/stelarx_outputs/${STELARX_A10K_DATASET_RECORD}" ]] || fail "run did not write the dataset record"
[[ -s "${RUN_OUTPUTS}/stelarx_outputs/R1/${STELARX_A10K_INPUTS_RECORD}" ]] || fail "run did not write the input fingerprints"

# The command record carries the exact run.sh line plus the A10K context.
RUN_CMD_FILE="${RUN_LEAF}/out-stelarx.command"
[[ -s "$RUN_CMD_FILE" ]] || fail "run command record was not mirrored"
grep -q "&& ./run.sh --input .*truegenetrees --output .*out-stelarx.tre --search-space S1 --cpu -q\$" "$RUN_CMD_FILE" ||   fail "command record lacks the exact run.sh invocation: $(cat "$RUN_CMD_FILE")"
grep -q "^# git_commit: " "$RUN_CMD_FILE" || fail "command record lacks the git commit"
grep -q "^# exit_code:    0$" "$RUN_CMD_FILE" || fail "command record lacks the exit code"
grep -q "^# tree type:    true$" "$RUN_CMD_FILE" || fail "command record lacks the tree type"
grep -q "^# setting:      search-space_S1__cpu_true$" "$RUN_CMD_FILE" || fail "command record lacks the setting"
grep -q "^# invoked as: .*run-a10k.sh" "$RUN_CMD_FILE" || fail "command record lacks the outer invocation"

# The skip path rebuilds a deleted mirror for free.
rm -rf "$RUN_OUTPUTS"
"${ROOT}/run-a10k.sh" --data-dir "$RUN_DATA" --tree-type true --replicates R1 \
  --opts '--search-space S1 --cpu -q' --no-time-monitor --no-gpu-monitor --no-notify \
  >"${TMP}/run2.out" 2>&1 || fail "second run-a10k.sh failed"
grep -q "SKIPPING:" "${TMP}/run2.out" || fail "second run did not skip"
[[ -s "${RUN_LEAF}/out-stelarx.tre" ]] || fail "the skip path did not rebuild the mirror"

# --no-outputs-mirror leaves the mirror untouched; an explicit dir is honored.
rm -rf "$RUN_OUTPUTS"
"${ROOT}/run-a10k.sh" --data-dir "$RUN_DATA" --tree-type true --replicates R1 --no-outputs-mirror \
  --opts '--search-space S1 --cpu -q' --no-time-monitor --no-gpu-monitor --no-notify \
  >"${TMP}/run3.out" 2>&1 || fail "run-a10k.sh --no-outputs-mirror failed"
[[ ! -e "$RUN_OUTPUTS" ]] || fail "--no-outputs-mirror still wrote a mirror"
"${ROOT}/run-a10k.sh" --data-dir "$RUN_DATA" --tree-type true --replicates R1 \
  --a10k-outputs-dir "${TMP}/explicit outputs" \
  --opts '--search-space S1 --cpu -q' --no-time-monitor --no-gpu-monitor --no-notify \
  >"${TMP}/run4.out" 2>&1 || fail "run-a10k.sh --a10k-outputs-dir failed"
[[ -s "${TMP}/explicit outputs/stelarx_outputs/R1/true/search-space_S1__cpu_true/out-stelarx.tre" ]] || \
  fail "explicit --a10k-outputs-dir was not used"
if "${ROOT}/run-a10k.sh" --data-dir "$RUN_DATA" --tree-type true --replicates R1 \
  --a10k-outputs-dir "${RUN_DATA}/outputs" --no-time-monitor --no-gpu-monitor --no-notify \
  >"${TMP}/run5.out" 2>&1; then
  fail "an outputs dir inside the data dir was accepted by the runner"
fi

# collect-scores-a10k.sh puts the merged summary in the mirror too.
"${ROOT}/collect-scores-a10k.sh" --data-dir "$RUN_DATA" --start-rep 1 --end-rep 1 \
  >"${TMP}/collect.out" 2>&1 || fail "collect-scores-a10k.sh failed: $(cat "${TMP}/collect.out")"
[[ -s "${RUN_OUTPUTS}/stelarx_outputs/a10k_stelarx_scores_merged.csv" ]] || \
  fail "the merged summary was not mirrored"
cmp -s "${RUN_DATA}/a10k_stelarx_scores_merged.csv" "${RUN_OUTPUTS}/stelarx_outputs/a10k_stelarx_scores_merged.csv" || \
  fail "the mirrored merged summary differs from the data-tree one"

# ------------------------------------------------------ runner integration ---
grep -q -- '--a10k-outputs-dir' "${ROOT}/run-a10k.sh" || fail "run-a10k.sh lacks --a10k-outputs-dir"
grep -q -- '--no-outputs-mirror' "${ROOT}/run-a10k.sh" || fail "run-a10k.sh lacks --no-outputs-mirror"
grep -q 'mirror_results_dir' "${ROOT}/run-a10k.sh" || fail "run-a10k.sh does not mirror its results"
grep -q -- '--no-outputs-mirror' "${ROOT}/collect-scores-a10k.sh" || fail "collect-scores-a10k.sh lacks --no-outputs-mirror"

echo "PASS: A10K outputs mirror"
