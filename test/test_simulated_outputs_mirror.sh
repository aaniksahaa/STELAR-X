#!/usr/bin/env bash
# Verifies the reproducibility mirror of simulated-run outputs:
#   * default outputs-directory derivation and containment guards,
#   * back-fill through sync-simulated-outputs.sh (command files copied, gene
#     trees / species trees / databases never copied, incomplete datasets fall
#     back to the base dataset's command),
#   * automatic mirroring by test-stelarx-simulated.sh on real runs and on the
#     "already completed" skip path, and its --no-outputs-mirror switch,
#   * upload-bulk-simulated-outputs.sh planning, remote paths, and refusals.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${ROOT}/scripts/simphy-outputs-dir.sh"

TMP="$(mktemp -d "${TMPDIR:-/tmp}/stelarx-outputs-mirror-test.XXXXXX")"
trap 'status=$?; rm -rf -- "$TMP"; exit "$status"' EXIT

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# ---------------------------------------------------------------- helpers ---
BASE="${TMP}/phylogeny root"
DATA="${BASE}/simphy/data"
OUTPUTS="${BASE}/outputs/simphy"
mkdir -p "$DATA"

[[ "$(stelarx_default_simphy_outputs_dir "$DATA")" == "$OUTPUTS" ]] || \
  fail "standard simphy/data tree should mirror into outputs/simphy"
[[ "$(stelarx_default_simphy_outputs_dir "${TMP}/other/data")" == "${TMP}/other/outputs" ]] || \
  fail "a non-simphy data dir should mirror into its 'outputs' sibling"
[[ "$(stelarx_default_simphy_outputs_dir "${TMP}/custom")" == "${TMP}/custom_outputs" ]] || \
  fail "custom data dir should mirror into a '<name>_outputs' sibling"
RESOLVED="$(stelarx_prepare_simphy_outputs_dir "" "$DATA")"
[[ "$RESOLVED" == "$OUTPUTS" && -d "$OUTPUTS" ]] || fail "default outputs dir was not created: $RESOLVED"
if stelarx_prepare_simphy_outputs_dir "${DATA}/outputs" "$DATA" >/dev/null 2>"${TMP}/inside.err"; then
  fail "an outputs dir inside the data dir was accepted"
fi
grep -q "inside it" "${TMP}/inside.err" || fail "containment error was unclear"
if stelarx_prepare_simphy_outputs_dir "${BASE}" "$DATA" >/dev/null 2>/dev/null; then
  fail "an outputs dir containing the data dir was accepted"
fi
if PHYLOGENY_DATA_DIR="$BASE" true && [[ "$(PHYLOGENY_DATA_DIR="$BASE" stelarx_prepare_simphy_outputs_dir_standalone "")" != "$OUTPUTS" ]]; then
  fail "standalone resolver did not use \$PHYLOGENY_DATA_DIR/outputs/simphy"
fi

# ------------------------------------------------------- synthetic dataset ---
DS="t_4_g_1_sb_0.000001_spmin_100000_spmax_200000"
DS_INC="${DS}_incomplete"
make_dataset_files() {
  local dir="$1" name="$2"
  mkdir -p "$dir"
  printf 'simphy -sl f:4 -o %s\n' "$name" > "${dir}/${name}.command"
  printf 'params\n' > "${dir}/${name}.params"
  printf 'sqlite\n' > "${dir}/${name}.db"
}
make_results() {
  local results="$1" tag="$2"
  mkdir -p "$results"
  printf '((a,b),(c,d));\n' > "${results}/out-stelarx.tre"
  printf 'alg,setting\nstelarx,%s\n' "$tag" > "${results}/stat-stelarx.csv"
  printf 'algorithm\nstelar-x\n' > "${results}/out-stelarx_stats.csv"
  : > "${results}/.stelarx.lock"
  printf 'exit_code=0\n' > "${results}/.stelarx.success"
  printf 'log %s\n' "$tag" > "${results}/.stelarx_run.log"
}

make_dataset_files "${DATA}/${DS}" "$DS"
for r in R1 R2; do
  mkdir -p "${DATA}/${DS}/${r}"
  printf '((a,b),(c,d));\n' > "${DATA}/${DS}/${r}/all_gt.tre"
  printf '((a,b),(c,d));\n' > "${DATA}/${DS}/${r}/s_tree.trees"
  printf 'num-taxa\n4\n' > "${DATA}/${DS}/${r}/stat-sim.csv"
  make_results "${DATA}/${DS}/${r}/stelarx_outputs/search-space_S1" "${r}-S1"
done
make_results "${DATA}/${DS}/R1/stelarx_outputs/search-space_S2" "R1-S2"
make_results "${DATA}/${DS}/R1/aster_outputs/default" "R1-aster"
# Incomplete variant: no own .command (older sim_incomplete.sh), results present.
mkdir -p "${DATA}/${DS_INC}/R1"
printf '((a,b),c);\n' > "${DATA}/${DS_INC}/R1/all_gt.tre"
printf '((a,b),(c,d));\n' > "${DATA}/${DS_INC}/R1/s_tree.trees"
make_results "${DATA}/${DS_INC}/R1/stelarx_outputs/search-space_S1" "inc-R1-S1"
# Raw SimPhy replicate directory that must be ignored.
mkdir -p "${DATA}/${DS}/1"
printf 'raw\n' > "${DATA}/${DS}/1/g_trees1.trees"

# ------------------------------------------------------------ sync dry run ---
PHYLOGENY_DATA_DIR="$BASE" "${ROOT}/sync-simulated-outputs.sh" --dry-run >"${TMP}/sync-dry.out" 2>&1
grep -q "would mirror=5" "${TMP}/sync-dry.out" || fail "dry run did not count 5 results dirs: $(cat "${TMP}/sync-dry.out")"
[[ -z "$(find "$OUTPUTS" -mindepth 1 -print -quit)" ]] || fail "dry run wrote into the outputs dir"

# ---------------------------------------------------------------- sync run ---
PHYLOGENY_DATA_DIR="$BASE" "${ROOT}/sync-simulated-outputs.sh" >"${TMP}/sync.out" 2>&1
grep -q "mirrored=5 filtered-out=0 failed=0" "${TMP}/sync.out" || fail "sync summary unexpected: $(cat "${TMP}/sync.out")"

M="${OUTPUTS}/stelarx_outputs/${DS}"
[[ -f "${M}/${DS}.command" ]] || fail "dataset .command was not mirrored"
[[ -f "${M}/${DS}.params" ]] || fail "dataset .params was not mirrored"
[[ ! -e "${M}/${DS}.db" ]] || fail "SimPhy database must not be mirrored"
cmp -s "${DATA}/${DS}/${DS}.command" "${M}/${DS}.command" || fail "mirrored .command differs"
for leaf in R1/search-space_S1 R1/search-space_S2 R2/search-space_S1; do
  [[ -s "${M}/${leaf}/out-stelarx.tre" ]] || fail "missing mirrored tree: $leaf"
  [[ -s "${M}/${leaf}/stat-stelarx.csv" ]] || fail "missing mirrored stat csv: $leaf"
  [[ -f "${M}/${leaf}/.stelarx.lock" && -f "${M}/${leaf}/.stelarx.success" ]] || fail "run markers not mirrored: $leaf"
done
[[ -s "${OUTPUTS}/aster_outputs/${DS}/R1/default/out-stelarx.tre" ]] || fail "second method was not mirrored"
[[ -f "${OUTPUTS}/aster_outputs/${DS}/${DS}.command" ]] || fail "second method lacks the .command copy"
[[ -f "${OUTPUTS}/stelarx_outputs/${DS_INC}/${DS}.command" ]] || fail "incomplete dataset lacks base .command"
[[ -s "${OUTPUTS}/stelarx_outputs/${DS_INC}/R1/search-space_S1/out-stelarx.tre" ]] || fail "incomplete results not mirrored"
[[ -z "$(stelarx_simphy_find_forbidden_in_mirror "$OUTPUTS")" ]] || \
  fail "forbidden files leaked into the mirror: $(stelarx_simphy_find_forbidden_in_mirror "$OUTPUTS")"
[[ ! -e "${OUTPUTS}/stelarx_outputs/${DS}/1" ]] || fail "raw SimPhy replicate dir leaked into the mirror"
[[ -z "$(find "$OUTPUTS" -name 'stat-sim.csv' -print -quit)" ]] || fail "stat-sim.csv leaked into the mirror"
# The data tree keeps its exact layout.
[[ -f "${DATA}/${DS}/R1/stelarx_outputs/search-space_S1/out-stelarx.tre" && -f "${DATA}/${DS}/R1/all_gt.tre" ]] || \
  fail "sync altered the data tree"
[[ ! -e "${DATA}/outputs" ]] || fail "sync created an outputs dir inside the data dir"

# Re-sync replaces a stale mirror leaf exactly (stale extra file disappears).
printf 'stale\n' > "${M}/R1/search-space_S1/stale.txt"
printf '((a,c),(b,d));\n' > "${DATA}/${DS}/R1/stelarx_outputs/search-space_S1/out-stelarx.tre"
PHYLOGENY_DATA_DIR="$BASE" "${ROOT}/sync-simulated-outputs.sh" --methods stelarx --quiet >"${TMP}/resync.out" 2>&1
grep -q "mirrored=4 filtered-out=1 failed=0" "${TMP}/resync.out" || fail "method filter summary unexpected: $(cat "${TMP}/resync.out")"
[[ ! -e "${M}/R1/search-space_S1/stale.txt" ]] || fail "stale mirror file survived a re-sync"
cmp -s "${DATA}/${DS}/R1/stelarx_outputs/search-space_S1/out-stelarx.tre" "${M}/R1/search-space_S1/out-stelarx.tre" || \
  fail "re-sync did not refresh the tree"
[[ -z "$(find "${M}/R1" -maxdepth 1 -name '.*mirror*' -print -quit)" ]] || fail "temporary mirror dir left behind"

# Results containing simulated input data are refused.
mkdir -p "${DATA}/${DS}/R2/stelarx_outputs/bad"
printf 'x\n' > "${DATA}/${DS}/R2/stelarx_outputs/bad/all_gt.tre"
if PHYLOGENY_DATA_DIR="$BASE" "${ROOT}/sync-simulated-outputs.sh" --quiet >"${TMP}/bad.out" 2>&1; then
  fail "sync succeeded although a results dir contained all_gt.tre"
fi
grep -q "simulated input data" "${TMP}/bad.out" || fail "refusal reason unclear: $(cat "${TMP}/bad.out")"
[[ ! -e "${M}/R2/bad" ]] || fail "forbidden results dir was mirrored anyway"
rm -rf "${DATA}/${DS}/R2/stelarx_outputs/bad"

# ------------------------------------------------ real STELAR-X run mirror ---
RUN_DATA="${TMP}/run/data"
RUN_OUTPUTS="${TMP}/run/outputs"
RUN_DS="t_4_g_1_sb_0.000001_spmin_100000_spmax_200000"
mkdir -p "${RUN_DATA}/${RUN_DS}/R1"
printf 'simphy -sl f:4\n' > "${RUN_DATA}/${RUN_DS}/${RUN_DS}.command"
printf '((a,b),(c,d));\n' > "${RUN_DATA}/${RUN_DS}/R1/all_gt.tre"
printf '((a,b),(c,d));\n' > "${RUN_DATA}/${RUN_DS}/R1/s_tree.trees"
COMMON=(--simphy-data-dir "$RUN_DATA" -t 4 -g 1 -r R1
  --sb 0.000001 --spmin 100000 --spmax 200000
  --opts '--search-space S1 --cpu -q'
  --no-time-monitor --no-gpu-monitor --no-notify)

env -u PHYLOGENY_DATA_DIR "${ROOT}/test-stelarx-simulated.sh" "${COMMON[@]}" >"${TMP}/run1.out" 2>&1
grep -q "outputs mirror: ${RUN_OUTPUTS}" "${TMP}/run1.out" || fail "run did not report the default outputs mirror: $(grep -i mirror "${TMP}/run1.out")"
RUN_LEAF="${RUN_OUTPUTS}/stelarx_outputs/${RUN_DS}/R1/search-space_S1__cpu_true"
grep -q "Mirrored outputs to: ${RUN_LEAF}" "${TMP}/run1.out" || fail "run did not mirror its outputs"
SRC_LEAF="${RUN_DATA}/${RUN_DS}/R1/stelarx_outputs/search-space_S1__cpu_true"
[[ -s "${RUN_LEAF}/out-stelarx.tre" && -s "${RUN_LEAF}/stat-stelarx.csv" && -s "${RUN_LEAF}/out-stelarx_stats.csv" ]] || \
  fail "mirrored run is missing tree or CSVs"
[[ -f "${RUN_LEAF}/.stelarx.success" && -f "${RUN_LEAF}/.stelarx.lock" ]] || fail "mirrored run lacks markers"
# The exact STELAR-X command is recorded beside the tree and mirrored with it.
[[ -s "${RUN_LEAF}/out-stelarx.command" ]] || fail "run command record was not mirrored"
grep -q "&& ./run.sh --input .*all_gt.tre --output .*out-stelarx.tre --search-space S1 --cpu -q\$" "${RUN_LEAF}/out-stelarx.command" || \
  fail "command record lacks the exact run.sh invocation: $(cat "${RUN_LEAF}/out-stelarx.command")"
grep -q "^# git_commit: " "${RUN_LEAF}/out-stelarx.command" || fail "command record lacks the git commit"
grep -q "^# exit_code:    0$" "${RUN_LEAF}/out-stelarx.command" || fail "command record lacks the exit code"
grep -q "^# invoked as: .*test-stelarx-simulated.sh" "${RUN_LEAF}/out-stelarx.command" || fail "command record lacks the outer invocation"
grep -q "^# setting:      search-space_S1__cpu_true$" "${RUN_LEAF}/out-stelarx.command" || fail "command record lacks the setting"
diff -r "$SRC_LEAF" "$RUN_LEAF" >/dev/null || fail "mirror leaf differs from the results dir"
[[ -f "${RUN_OUTPUTS}/stelarx_outputs/${RUN_DS}/${RUN_DS}.command" ]] || fail "run did not mirror the .command file"
[[ ! -e "${RUN_OUTPUTS}/stelarx_outputs/${RUN_DS}/R1/all_gt.tre" && -z "$(stelarx_simphy_find_forbidden_in_mirror "$RUN_OUTPUTS")" ]] || \
  fail "run mirror contains simulated input data"
[[ -f "${SRC_LEAF}/out-stelarx.tre" && -f "${RUN_DATA}/${RUN_DS}/R1/all_gt.tre" ]] || fail "data tree layout changed"

# The skip path (already completed) rebuilds a deleted mirror.
rm -rf "$RUN_OUTPUTS"
env -u PHYLOGENY_DATA_DIR "${ROOT}/test-stelarx-simulated.sh" "${COMMON[@]}" >"${TMP}/run2.out" 2>&1
grep -q "SKIPPING: successful output already exists" "${TMP}/run2.out" || fail "second run did not skip"
[[ -s "${RUN_LEAF}/out-stelarx.tre" ]] || fail "skip path did not rebuild the mirror"

# --no-outputs-mirror leaves the mirror untouched; explicit dir is honored.
rm -rf "$RUN_OUTPUTS"
env -u PHYLOGENY_DATA_DIR "${ROOT}/test-stelarx-simulated.sh" "${COMMON[@]}" --no-outputs-mirror >"${TMP}/run3.out" 2>&1
[[ ! -e "$RUN_OUTPUTS" ]] || fail "--no-outputs-mirror still wrote a mirror"
grep -q "SKIPPING: successful output already exists" "${TMP}/run3.out" || fail "third run did not take the skip path"
env -u PHYLOGENY_DATA_DIR "${ROOT}/test-stelarx-simulated.sh" "${COMMON[@]}" --simphy-outputs-dir "${TMP}/explicit outputs" >"${TMP}/run4.out" 2>&1
[[ -s "${TMP}/explicit outputs/stelarx_outputs/${RUN_DS}/R1/search-space_S1__cpu_true/out-stelarx.tre" ]] || \
  fail "explicit --simphy-outputs-dir was not used"
if env -u PHYLOGENY_DATA_DIR "${ROOT}/test-stelarx-simulated.sh" "${COMMON[@]}" --simphy-outputs-dir "${RUN_DATA}/outputs" >"${TMP}/run5.out" 2>&1; then
  fail "an outputs dir inside the data dir was accepted by the run script"
fi

# run-bulk-simulated.sh forwards the mirror options.
grep -q -- '--simphy-outputs-dir' "${ROOT}/run-bulk-simulated.sh" || fail "bulk runner lacks --simphy-outputs-dir"
grep -q -- '--no-outputs-mirror' "${ROOT}/run-bulk-simulated.sh" || fail "bulk runner lacks --no-outputs-mirror"

# ------------------------------------------------------- uploader dry run ---
UP="${ROOT}/upload-bulk-simulated-outputs.sh"

# The default publishes only this repository's own results.
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true >"${TMP}/up-default.out" 2>&1 || \
  fail "uploader dry run failed: $(cat "${TMP}/up-default.out")"
grep -q "Plan: upload 2 dataset" "${TMP}/up-default.out" || \
  fail "the default method filter did not plan the 2 stelarx datasets: $(cat "${TMP}/up-default.out")"
grep -q "aster" "${TMP}/up-default.out" && fail "the default method filter did not exclude aster"

# Everything below exercises the full mirror, including other methods.
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true >"${TMP}/up.out" 2>&1 || \
  fail "uploader dry run failed: $(cat "${TMP}/up.out")"
grep -q "Plan: upload 3 dataset" "${TMP}/up.out" || fail "uploader did not plan 3 dataset uploads: $(cat "${TMP}/up.out")"
grep -Fq -- "--path-in-repo ph/d/simulated/outputs/stelarx_outputs/${DS} " "${TMP}/up.out" || \
  fail "remote path for stelarx outputs is wrong"
grep -Fq -- "--path-in-repo ph/d/simulated/outputs/aster_outputs/${DS} " "${TMP}/up.out" || \
  fail "remote path for aster outputs is wrong"
grep -Fq -- "--path-in-repo ph/d/simulated/outputs/stelarx_outputs/${DS_INC} " "${TMP}/up.out" || \
  fail "incomplete dataset was not planned"
grep -Fq -- "--local-path $(printf '%q' "${OUTPUTS}/stelarx_outputs/${DS}") " "${TMP}/up.out" || fail "local path is wrong"
grep -q "nothing was uploaded" "${TMP}/up.out" || fail "dry run did not state that nothing was uploaded"

PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --uploader /bin/true --python /bin/true --methods aster >"${TMP}/up-m.out" 2>&1
grep -q "Plan: upload 1 dataset" "${TMP}/up-m.out" || fail "--methods filter did not narrow the plan"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true --exclude-incomplete >"${TMP}/up-i.out" 2>&1
grep -q "Plan: upload 2 dataset" "${TMP}/up-i.out" || fail "--exclude-incomplete did not drop the incomplete dataset"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true --min-taxa 5 >"${TMP}/up-t.out" 2>&1
grep -q "No <method>_outputs/<dataset> directories matched" "${TMP}/up-t.out" || fail "--min-taxa did not filter everything out"

# --sync refreshes the mirror before planning (a new result appears).
make_results "${DATA}/${DS}/R2/stelarx_outputs/search-space_S3" "R2-S3"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --sync --method all --uploader /bin/true --python /bin/true >"${TMP}/up-s.out" 2>&1 || \
  fail "--sync dry run failed: $(cat "${TMP}/up-s.out")"
grep -q "Refreshing the outputs mirror" "${TMP}/up-s.out" || fail "--sync did not run the mirror sync"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --sync --yes --method all --uploader /bin/true --python /bin/true >"${TMP}/up-s2.out" 2>&1 || \
  fail "--sync upload with a stub uploader failed: $(cat "${TMP}/up-s2.out")"
[[ -s "${M}/R2/search-space_S3/out-stelarx.tre" ]] || fail "--sync did not back-fill the new result"
grep -q "uploaded=3 failed=0" "${TMP}/up-s2.out" || fail "stub upload summary unexpected: $(cat "${TMP}/up-s2.out")"

# Simulated input data inside the mirror blocks the upload.
printf 'leak\n' > "${M}/R1/all_gt.tre"
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true >"${TMP}/up-leak.out" 2>&1; then
  fail "uploader accepted a mirror containing all_gt.tre"
fi
grep -q "BLOCKED: contains all_gt.tre" "${TMP}/up-leak.out" || fail "leak was not reported: $(cat "${TMP}/up-leak.out")"
grep -q "Nothing was uploaded" "${TMP}/up-leak.out" || fail "leak refusal did not state nothing was uploaded"
rm -f "${M}/R1/all_gt.tre"

# A missing .command blocks unless explicitly allowed.
rm -f "${OUTPUTS}/aster_outputs/${DS}/${DS}.command"
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true >"${TMP}/up-cmd.out" 2>&1; then
  fail "uploader accepted a dataset without its .command file"
fi
grep -q "BLOCKED: missing .command" "${TMP}/up-cmd.out" || fail "missing command was not reported"
PHYLOGENY_DATA_DIR="$BASE" "$UP" --dry-run --method all --uploader /bin/true --python /bin/true --allow-missing-command >"${TMP}/up-cmd2.out" 2>&1 || \
  fail "--allow-missing-command did not unblock the plan"
grep -q "upload (NO .command)" "${TMP}/up-cmd2.out" || fail "missing command was not flagged in the plan"

# Non-interactive confirmation is refused without --yes (stdin closed).
if PHYLOGENY_DATA_DIR="$BASE" "$UP" --method all --uploader /bin/true --python /bin/true --allow-missing-command </dev/null >"${TMP}/up-noyes.out" 2>&1; then
  grep -q "Cancelled; nothing was uploaded" "${TMP}/up-noyes.out" || fail "uploader proceeded without confirmation"
fi

echo "Simulated outputs mirror: PASS"
