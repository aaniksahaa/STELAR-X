# Scaling coalescent-based species tree inference to 100,000 taxa with STELAR-X

**STELAR-X** is a scalable, statistically consistent summary method for coalescent-based species tree inference from large collections of gene trees. It combines compact bipartition encodings, fast weight precomputation, GPU-accelerated parallelism, and optimized dynamic programming to analyze datasets as large as **100,000 taxa × 1,000 genes in just 8.5 hours using 86 GB RAM**.

This repository is the **Supplemental Code** for the paper. It contains the complete STELAR-X source code (Java + CUDA), pre-built native GPU libraries, the launchers, and every experiment script used for the simulated, A10K, and biological analyses reported in the paper, together with the tooling that records and publishes the exact command behind each result so that every experiment can be reproduced.

> **Note:** Although STELAR-X provides a CPU fallback option, GPU-enabled execution is recommended to realize its full computational benefits.

## Citation

If you use **STELAR-X**, its source code, or results produced by it in your research, please cite our paper:

> Anik Saha and Md. Shamsuzzoha Bayzid.  
> **Scaling coalescent-based species tree inference to 100,000 taxa with STELAR-X.**  
> Accepted at RECOMB 2026.  
> Genome Research, 2026. https://doi.org/10.1101/gr.282257.126

```bibtex
@article{saha2026stelarx,
  title   = {Scaling coalescent-based species tree inference to
             100,000 taxa with {STELAR-X}},
  author  = {Saha, Anik and Bayzid, Md. Shamsuzzoha},
  journal = {Genome Research},
  year    = {2026},
  doi     = {10.1101/gr.282257.126},
  note    = {Accepted at RECOMB 2026}
}
```

> **Platform:** Developed and tested on **Ubuntu Linux**.

---

## Quick Start

**Prerequisite:** JDK 21 or newer (tested with OpenJDK 21). No Maven, Gradle, or other build tool is needed; the Java sources are compiled directly with `javac`.

```bash
# Install Java if not already installed (Ubuntu/Debian)
sudo apt update && sudo apt install -y openjdk-21-jdk

# Verify
java -version
```

### Build and run

```bash
# 1. Clone the repository
git clone https://github.com/aaniksahaa/STELAR-X-2.git
cd STELAR-X-2

# 2. Compile the Java sources (a few seconds)
./build.sh

# 3. Run on the included example (37-taxon dataset, 200 rooted gene trees)
./stelarx -i example/all_gt_37.tre -o example/out_37.tre
```

That's it. An example gene trees file is included so you can verify it works immediately; you can run inference on any gene tree file given by a relative or absolute path. Pre-built CUDA libraries are shipped in `native/`, so GPU mode is used automatically when an NVIDIA GPU with a working driver is detected; otherwise STELAR-X falls back to the multi-threaded CPU implementation with an explanation.

```bash
# Force CPU mode
./stelarx -i example/all_gt_37.tre -o example/out_37.tre --cpu

# Force GPU mode (falls back to CPU with a warning if CUDA is unusable)
./stelarx -i example/all_gt_37.tre -o example/out_37.tre --gpu

# Run on your own data
./stelarx -i gene_trees.tre -o output.tre

# Score a known rooted species tree against gene trees
./stelarx -i gene_trees.tre -c species_tree.tre

# Print installation and hardware diagnostics (JDK, CPU threads, GPU, CUDA libraries)
./stelarx --diagnose

# See all options
./stelarx --help
```

Each run prints a **Run Summary** with the triplet score, running time, peak CPU RAM, and peak GPU VRAM.

> By default the launcher recompiles the Java sources before every run so that source edits are always picked up. Add `--no-build` to skip this step (for example inside experiment sweeps, or when timing runs).

**Optional** — add to PATH for system-wide access:

```bash
sudo ln -sf $(pwd)/stelarx /usr/local/bin/stelarx
stelarx -i gene_trees.tre -o output.tre    # works from anywhere
```

---

## Input Requirements

- The input file contains **one rooted Newick tree per non-empty line**.
- STELAR-X uses the supplied top-level root exactly as given. It never reroots or invents a root. A top-level node with anything other than two children is rejected, because ordinary Newick has no independent rootedness flag.
- Gene trees may be **incomplete** (missing taxa) and may contain **polytomies**. Polytomies are resolved deterministically during inference by default; the final triplet score always respects the input topology. Use `--keep-polytomy-during-inference` to keep them during inference as well.
- Branch lengths and support values are accepted and ignored. Quoted taxon labels are supported.

If your gene trees are unrooted, `process_unrooted.sh` roots them by a clustered outgroup (with fallbacks) and strips branch lengths and internal labels:

```bash
./process_unrooted.sh -i unrooted_gene_trees.tre -o rooted_gene_trees.tre -ogf outgroup_taxa.txt
```

Run `./process_unrooted.sh --help` and `python3 root_by_outgroups.py --help` for every rooting option.

---

## Building from Source

### Prerequisites

| Dependency | Required | Version | Notes |
|------------|----------|---------|-------|
| **Java (JDK)** | Yes | 21+ | Tested with OpenJDK 21 |
| **NVIDIA CUDA Toolkit** | No | 11.0+ | Only to *rebuild* the GPU libraries; pre-built `.so` files are included |
| **NVIDIA GPU + driver** | No | Compute capability ≥ 5.0 | Required only for GPU mode; the shipped libraries cover `sm_50`–`sm_90` plus PTX; CPU fallback is automatic |
| **Python 3 + DendroPy** | No | DendroPy 5.x | Only for RF-distance evaluation in the experiment scripts |

### Java

```bash
./build.sh
```

`build.sh` compiles every file under `src/stelarx/` into `build/` with plain `javac`. Re-run it after editing any Java source (the `./stelarx` launcher also does this automatically unless `--no-build` is given).

### CUDA libraries

The four native libraries are already built and committed under `native/`:

| Library | Purpose |
|---------|---------|
| `libstelarx_weight.so` | Rooted-triplet weight (intersection) kernels |
| `libstelarx_dp.so` | GPU cross-tree DP transition search |
| `libstelarx_dist.so` | GPU distance matrix (Euler tour + RMQ) for tree completion |
| `libstelarx_sim.so` | GPU similarity matrix for tree completion |

Rebuild them only if you change `src/native/*.cu` or need a different architecture set:

```bash
./build_native.sh                  # all major GPU generations of the installed toolkit (portable)
CUDA_ARCH=native ./build_native.sh # only the GPU in this machine (fast developer build)
CUDA_ARCH=sm_86 ./build_native.sh  # one explicit architecture
```

The library records the oldest compute capability it covers. At start-up, STELAR-X probes the GPU, driver, and library; if any of them is unusable it explains why and selects the CPU implementation. Use `--gpu-strict` to turn that fallback into an immediate failure.

### Installing the CUDA Toolkit (if you want to rebuild the kernels)

```bash
# Check if nvcc is available
nvcc --version

# If missing (Ubuntu/Debian):
sudo apt update
sudo apt install -y nvidia-cuda-toolkit

# Verify
nvcc --version
```

### Development environment (optional)

```bash
./setup_dev.sh              # creates .venv, installs DendroPy, builds Java (+CUDA if nvcc exists), runs CPU tests
./setup_dev.sh --cpu-only   # skip CUDA compilation
./setup_dev.sh --check      # only verify the environment; change nothing
```

### Self-contained release archive

```bash
./build_portable.sh                 # Linux image with bundled Java runtime and CUDA libraries, CPU fallback included
./build_portable.sh --without-cuda  # CPU-only image
```

Artifacts, SHA-256 checksums, and a JSON manifest (version, platform, capability, minimum glibc) are written under `dist/<version>/`. Target machines need neither Java nor CUDA. Run the script on each target platform (`build_portable.ps1` on Windows); macOS builds are CPU-only because CUDA is unavailable there.

---

## Usage

All commands below use the `./stelarx` launcher (`./run.sh` is the same script; `./run-stelarx.sh` is an alias).

### Inference Mode (default)

Infer a rooted species tree from a set of rooted gene trees:

```bash
./stelarx -i <rooted_gene_trees.tre> -o <species_tree.tre> [options]
```

### Score Mode

Calculate the triplet score of a known rooted species tree against the gene trees:

```bash
./stelarx -i <rooted_gene_trees.tre> -c <rooted_species_tree.tre> [options]
```

The machine-readable result is printed as `TRIPLET_SCORE: N`.

### Taxon Extraction

```bash
./stelarx -i <gene_trees.tre> --extract-taxa -o taxa.txt                      # union of all taxa
./stelarx -i <gene_trees.tre> --extract-taxa --taxa-set intersection -o taxa.txt
```

### With Performance Monitoring

Records running time, peak CPU RAM, peak GPU VRAM, the triplet score, and (optionally) the RF rate against a reference tree, and writes a `<output>_stats.csv` and a `<output>.command` file beside the output tree:

```bash
./run-stelarx-with-monitor.sh -i <gene_trees.tre> -o <output.tre> [--reference-species-tree <true_tree.tre>] [--opts "<stelarx options>"]
```

### Examples

```bash
# Automatic mode (GPU if available, otherwise CPU)
./stelarx -i example/all_gt_37.tre -o out-37.tre

# Explicit CPU, 8 threads
./stelarx -i example/all_gt_37.tre -o out-37.tre --cpu --threads 8

# GPU, search space S2 (tree completion + cross-tree recombination)
./stelarx -i example/all_gt_37.tre -o out-37.tre --gpu --search-space S2

# Search space S3 with intersection method I4
./stelarx -i example/all_gt_37.tre -o out-37.tre --search-space S3 --intersection-method I4

# Score-only mode
./stelarx -i example/all_gt_37.tre -c example/true_37.tre

# Restrict inference to a subset of taxa
./stelarx -i example/all_gt_37.tre -o out-sub.tre --taxa-file my_taxa.txt

# Large dataset with a custom Java heap and a log file
./stelarx -i large_dataset.tre -o out.tre --xms 8g --xmx 256g --log-file run.log

# Monitored run with RF rate against the true tree
./run-stelarx-with-monitor.sh -i example/all_gt_37.tre -o out-37.tre --reference-species-tree example/true_37.tre
```

Three larger rooted biological/benchmark inputs used during development are also included at the repository root together with their reference trees: `all_gt_bs_rooted_37.tre`, `all_gt_bs_rooted_48.tre`, and `all_gt_bs_rooted_200.tre` (`true_37.tre`, `true_48.tre`, `true_200.tre`). `bash test_rf.sh 37|48|200|all` runs both DP search modes on them and reports RF distances.

---

## Command-Line Parameters

### Required Parameters

| Flag | Long Form | Description |
|------|-----------|-------------|
| `-i` | `--input <file>` | Input rooted gene trees file in Newick format (one tree per line) |
| `-o` | `--output <file>` | Output species tree file (inference mode; stdout when omitted) |

> In **score mode**, use `-c` instead of `-o`.

### Computation Mode

| Flag | Long Form | Description |
|------|-----------|-------------|
| | `--auto` | Use CUDA when usable, otherwise CPU (**default**) |
| | `--gpu` | Prefer CUDA; warn and fall back to CPU if unavailable |
| | `--gpu-strict` | Require CUDA; fail before reading input if unavailable |
| | `--cpu` | Force multi-threaded CPU execution |
| `-t` | `--threads <n>` | CPU worker threads (default: all available cores) |

### Search Space and Scoring

| Flag | Long Form | Description | Default |
|------|-----------|-------------|---------|
| | `--search-space <S1..S3>` | Candidate search-space preset (`S1` is the baseline; `S2`/`S3` progressively enlarge the search space) | `S1` |
| `--im` | `--intersection-method <I1..I4>` | Weight-calculation (intersection) method | `I2` |
| `-m` | `--seeds <n>` | Number of cluster-hash seeds | `2` |
| | `--keep-polytomy-during-inference` | Keep input polytomies during inference (final scoring always keeps them) | resolve |
| | `--no-prune-search-space` | Disable the DP-reachability weight prune | prune on |
| | `--large-n-score-type <int128\|double>` | Accumulator for very large scores | exact |

### Optional Parameters

| Flag | Long Form | Description | Default |
|------|-----------|-------------|---------|
| `-c` | `--score-species-tree <file>` | Score-only mode: triplet score of the given rooted species tree | — |
| | `--taxa-file <file>` | Restrict inference or scoring to the listed taxa (one per line) | — |
| | `--extract-taxa` | Write the input taxa (one per line) and exit | — |
| | `--taxa-set <union\|intersection>` | Taxon extraction operation | `union` |
| | `--log-file <file>` | Save run messages to a file (progress bars stay terminal-only) | — |
| `-q`, `-v`, `-vv`, `-vvv` | | Quiet / info / debug / trace logging | info |
| | `--xms <size>` | Java minimum heap size | `256m` |
| | `--xmx <size>` | Java maximum heap size | `128g` |
| | `--no-build` | Skip the automatic `build.sh` before running | build |
| `-nn` | `--no-notify` | Disable the optional push notification after score-only runs | — |
| | `--diagnose` | Print runtime/backend diagnostics and exit | — |
| | `--version` | Print the STELAR-X version and exit | — |
| `-h` | `--help` | Show the help message | — |

GPU batching and VRAM controls (`--gpu-batch-size`, `--gpu-vram-occupancy-factor`, `--gpu-sim-vram-cap-mb`, …), tree-completion options, and verification dumps (`--verify-*`, `--dump-clusters`) are documented in `./stelarx --help` and `java -cp build stelarx.Main --help`. Their defaults are the ones used for the paper.

### Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `STELARX_XMS` | Default minimum Java heap size | `STELARX_XMS=8g ./stelarx ...` |
| `STELARX_XMX` | Default maximum Java heap size | `STELARX_XMX=256g ./stelarx ...` |
| `STELARX_CRASH_DIR` | Directory for Java/JVM fatal-error reports | default `crash_logs/` in the checkout |
| `STELARX_PYTHON` | Python interpreter with DendroPy used for RF rates by the monitor script | default `.venv/bin/python`, then `python3` |
| `PHYLOGENY_DATA_DIR` | Root directory for simulated datasets and experiment outputs (see below) | `export PHYLOGENY_DATA_DIR=$HOME/phylogeny-data` |
| `NTFY_CHANNEL_NAME` | Channel for the optional completion notification (disable with `--no-notify`) | — |

---

## Reproducing the Experiments

Every experiment in the paper was run through the scripts in this repository. Two conventions make the runs reproducible:

1. **A shared data root.** Set `PHYLOGENY_DATA_DIR` once; every experiment script places simulated datasets under `$PHYLOGENY_DATA_DIR/simphy/data` and the small run outputs under `$PHYLOGENY_DATA_DIR/outputs/`. An explicit `--simphy-data-dir` / `--data-dir` option always overrides the default.

   ```bash
   echo 'export PHYLOGENY_DATA_DIR="$HOME/phylogeny-data"' >> ~/.bashrc
   source ~/.bashrc
   ```

2. **A command record beside every result.** Each inference writes `out-stelarx.command` next to `out-stelarx.tre`: the exact launcher invocation with every flag and absolute path, the git commit, the wrapper command, the exit code, and the running time. Datasets carry the exact SimPhy command that generated them. Re-running a result is therefore a matter of replaying its command file.

### Simulated datasets (SimPhy)

We use [SimPhy](https://github.com/adamallo/SimPhy) (binary included under `simphy/`) to generate the large-scale simulated datasets. Dataset names encode their parameters: `t_<taxa>_g_<genes>_sb_<speciation rate>_spmin_<min pop. size>_spmax_<max pop. size>`, with replicates `R1`, `R2`, ….

```bash
# Generate a dataset (100 taxa, 200 gene trees, 1 replicate)
./sim.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1

# Run STELAR-X on replicate R1 of that dataset (generates it first if missing)
./test-stelarx-simulated.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 \
    --opts "--search-space S1 --intersection-method I2 -vv"
```

Results are written to `<data>/<dataset>/<replicate>/stelarx_outputs/<setting>/` where `<setting>` is derived from the options (for example `search-space_S1__intersection-method_I2`). Each result directory contains `out-stelarx.tre`, `out-stelarx.command`, `stat-stelarx.csv` (RF rate, triplet score, time, peak RAM/VRAM), `out-stelarx_stats.csv`, and the run log.

**Bulk sweeps.** `run-bulk-simulated.sh` runs the Cartesian product of parameter lists × replicates × settings, prints the complete plan (one line per `<dataset> / <replicates> / <setting>`), asks for confirmation, and skips already-completed runs unless `--fresh` is given:

```bash
./run-bulk-simulated.sh --taxa-list "1000,5000" --genes-list "1000" --num-replicates 5 \
    --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv" --dry-run   # show the plan only
./run-bulk-simulated.sh --taxa-list "1000,5000" --genes-list "1000" --num-replicates 5 \
    --opts-list "--search-space S1 -vv;--search-space S2 -vv;--search-space S3 -vv" --yes
```

The parameter lists used in the paper are kept at the top of the script for reference. Generating a 100,000-taxon dataset takes a long time; the exact datasets used in the paper are therefore published (as ZIPs) in the Hugging Face dataset repository and can be fetched with `./download-bulk-simulated.sh` (see its `--help`).

**Incomplete gene trees.** `sim_incomplete.sh` derives an `<dataset>_incomplete` variant by randomly pruning taxa from a complete dataset; run it with `test-stelarx-simulated.sh --incomplete`.

**Collecting statistics.**

```bash
./collect-stats-simulated.sh --out perf-combined.csv
```

merges every `stat-stelarx.csv` under the data root into one CSV (`alg, setting, num-taxa, gene-trees, replicate, sb, spmin, spmax, rf-rate, optimal-triplet-score, running-time-s, max-cpu-mb, max-gpu-mb, …`).

### A10K dataset (10,000-taxon SimPhy dataset with true and estimated gene trees)

```bash
./run-a10k.sh --data-dir $PHYLOGENY_DATA_DIR/10k-astral-dataset --tree-type "true;estimated" --replicates 1-20 \
    --opts "--search-space S1 --intersection-method I2 -vv"
./collect-scores-a10k.sh --data-dir $PHYLOGENY_DATA_DIR/10k-astral-dataset --start-rep 1 --end-rep 20
```

`--data-dir` must contain the dataset's `10k-simphy/R*/` directories. Results go to `10k-simphy/<R>/stelarx_outputs/<tree-type>/<setting>/`, and the collector writes `a10k_stelarx_scores_merged.csv`.

### Biological datasets

Biological gene trees are analysed with the same launcher, typically through the monitored wrapper so that time, memory, score and the command record are captured:

```bash
./run-stelarx-with-monitor.sh -i all_gt_bs_rooted_48.tre -o out-48.tre --reference-species-tree true_48.tre \
    --opts "--search-space S1 --intersection-method I2 -vv"
```

`run-bulk-standard.sh` and `collect-stats-standard.sh` drive the same runs (and the baseline methods, when their binaries are placed under `baselines/`) across a directory of standard datasets; see `--help` on each.

### Reproducibility mirror of run outputs

Gene trees dominate the size of the data tree while the results are tiny. Every simulated and A10K run therefore also mirrors its results (inferred tree, CSVs, run markers, log, command record) together with the dataset's provenance record (SimPhy `.command`/`.params`, or the A10K dataset record and input fingerprints) into a shareable outputs tree:

```
$PHYLOGENY_DATA_DIR/outputs/simphy/stelarx_outputs/<dataset>/<dataset>.command
$PHYLOGENY_DATA_DIR/outputs/simphy/stelarx_outputs/<dataset>/R1/<setting>/out-stelarx.tre, out-stelarx.command, stat-stelarx.csv, ...
$PHYLOGENY_DATA_DIR/outputs/10k-astral-dataset/stelarx_outputs/a10k-dataset.command
$PHYLOGENY_DATA_DIR/outputs/10k-astral-dataset/stelarx_outputs/R1/<tree-type>/<setting>/...
```

Simulated inputs are never copied. The mirrors can be back-filled from existing results and published to the Hugging Face dataset repository (outputs and command records only):

```bash
./sync-simulated-outputs.sh --dry-run && ./sync-simulated-outputs.sh
./upload-bulk-simulated-outputs.sh --dry-run
./upload-bulk-simulated-outputs.sh --sync

./sync-a10k-outputs.sh --dry-run && ./sync-a10k-outputs.sh
./upload-a10k-outputs.sh --dry-run
```

A dataset whose mirror lacks its provenance record, or that contains simulated input data, is refused before anything is uploaded. Both mirrors share their primitives in `scripts/outputs-mirror-common.sh`; the regression tests `test/test_simulated_outputs_mirror.sh` and `test/test_a10k_outputs_mirror.sh` pin down the layout and the refusal rules.

### Cleaning up

```bash
./clear-bulk-simulated.sh --dry-run     # preview: removes $PHYLOGENY_DATA_DIR/simphy/data completely
./clear-bulk-simulated.sh --yes
./clear-a10k.sh --data-dir $PHYLOGENY_DATA_DIR/10k-astral-dataset --dry-run   # removes only STELAR-X results, keeps inputs
```

### Optional tools

These are not required for STELAR-X itself, but are used by the evaluation scripts:

```bash
pip install dendropy        # RF distance (rf.py, analyze-dataset.py, RF rates in the monitor script)
sudo apt install -y time    # GNU time for peak-RSS monitoring (run-stelarx-with-monitor.sh)
```

`rf.py` computes the normalized Robinson-Foulds distance between two trees; `analyze-dataset.py` reports average gene-tree/gene-tree and gene-tree/species-tree RF distances of a dataset.

---

## Testing

```bash
# Focused regression suite (~1 min): triplet arithmetic across I1–I4 on binary, polytomous and
# incomplete inputs, S1–S3 inference, CLI identity, diagnostics, and the experiment-script tests
test/run_stelarx_tests.sh

# Complete layered suite: Java unit tests, independent randomized Python oracles, malformed inputs,
# end-to-end inference, a CPU-only portable package, and the strict CUDA layer when a GPU is usable
test/run_stelarx_comprehensive_tests.sh              # add --require-gpu, --cpu-only, --quick, --skip-packaging

# Strict GPU layer only (fails if CUDA is unavailable)
test/run_stelarx_gpu_tests.sh

# Accuracy plus wall-time / peak-RSS scaling across the complete S1–S3 × I1–I4 matrix
python3 test/test_stelarx_scalability.py --require-gpu
```

---

## Implementation Notes

- STELAR-X maximizes agreement with the rooted triplets displayed by the input gene trees. Clusters are the descendant clades of the gene-tree nodes; the input root is authoritative.
- Clusters are represented by constant-size, multi-seed hashes with 2^64 wrap-around arithmetic (addition and XOR of per-taxon mixed hashes), so cluster identity and set operations are O(1) regardless of the number of taxa.
- Triplet weights are computed with prefix-sum range intersections over hashed clusters; alternative exact intersection methods are available. Scores are accumulated exactly in `long`/128-bit arithmetic, with an optional `double` mode for astronomically large scores.
- The dynamic programme is hash-binned with reachability pruning. The larger search spaces add hash-subtraction cross-tree recombination, similarity/UPGMA guidance, root-preserving completion of incomplete trees, and consensus enrichment.
- The CUDA kernels (`src/native/*.cu`) cover weight calculation, DP transition search, and distance/similarity matrices, with automatic VRAM-aware batching. The shipped libraries are built with CUDA 12 for every major GPU generation from Maxwell (`sm_50`) to Hopper (`sm_90`) plus forward-compatible PTX.
- Correctness is pinned by the test suites: every intersection method is checked against independent Python re-implementations of the triplet score on binary, polytomous, and incomplete inputs, and the GPU kernels are checked bit-for-bit against the CPU paths.

---

## Troubleshooting

| Problem | Possible Cause | Solution |
|---------|---------------|----------|
| `build/` missing or `ClassNotFoundException` | Project not built | Run `./build.sh` (or drop `--no-build`) |
| `javac: command not found` / unsupported class version | JDK missing or older than 21 | `sudo apt install -y openjdk-21-jdk` |
| Run says CPU mode although a GPU is present | Driver, CUDA runtime, or `native/libstelarx_*.so` unusable | Run `./stelarx --diagnose`; see below |
| `OutOfMemoryError` | Dataset too large for the default heap | Increase heap: `--xms 8g --xmx 256g` |
| GPU out of memory | Batches too large for the VRAM | Lower `--gpu-vram-occupancy-factor`, set `--gpu-batch-size`, or use `--cpu` |
| Input rejected as unrooted | Top-level node has ≠ 2 children | Root the trees first (`process_unrooted.sh`) |
| Unexpected Java failure | — | See the report written under `crash_logs/` (or `$STELARX_CRASH_DIR`) |

### `nvidia-smi` works, but STELAR-X selects CPU mode

`nvidia-smi` only confirms that the driver can see the GPU. STELAR-X additionally probes its native libraries in `native/` and the GPU's compute capability at start-up and prints the reason for any fallback. Check each layer:

```bash
nvidia-smi
./stelarx --diagnose --no-build                     # shows "CUDA usable" or the exact reason it is not
ls -l native/libstelarx_*.so
git status --short -- native/                       # a deleted library shows as "D native/..."
```

If a library is missing, restore the committed one with `git restore native/` or rebuild with `./build_native.sh` (needs `nvcc`). If the GPU is older than the minimum compute capability embedded in the library, rebuild for that GPU, for example `CUDA_ARCH=sm_75 ./build_native.sh`. Use `--gpu-strict` when a silent CPU fallback is unacceptable.

---

## Project Structure

```
STELAR-X-2/
├── stelarx, run.sh, run-stelarx.sh   # Launcher (auto-builds, selects GPU/CPU, sets heap and library path)
├── run-stelarx-with-monitor.sh       # Launcher with time / RAM / VRAM monitoring, RF rate, stats CSV and command record
├── build.sh                          # Compile Java sources into build/
├── build_native.sh (.ps1)            # Rebuild the CUDA JNI libraries into native/
├── build_portable.sh (.ps1)          # Self-contained release image with bundled Java runtime (+ CUDA)
├── setup_dev.sh                      # Development environment: .venv, DendroPy, build, CPU tests
├── src/
│   ├── stelarx/                      # Java sources (package stelarx)
│   │   ├── Main.java, Config.java, CliPresets.java  # CLI, configuration, S1–S3 / I1–I4 presets
│   │   ├── tree/, taxon/             # Newick parsing, rooted trees, taxon registry
│   │   ├── hash/, cluster/, partition/  # Multi-seed hashing, cluster and bipartition tables
│   │   ├── weight/                   # Triplet weight calculation (I1–I4)
│   │   ├── dp/                       # Hash-binned dynamic programming and inference
│   │   ├── completion/, greedy/      # Distance/similarity matrices, UPGMA, tree completion, consensus enrichment
│   │   ├── gpu/                      # Java side of the CUDA bindings and batching
│   │   └── util/                     # Threading, 128-bit arithmetic, progress bars
│   └── native/                       # CUDA kernels: stelarx_weight.cu, stelarx_dp.cu, stelarx_dist.cu, stelarx_similarity.cu
├── native/                           # Pre-built CUDA libraries (libstelarx_weight/dp/dist/sim.so)
├── example/                          # 37-taxon example: all_gt_37.tre (200 rooted gene trees), true_37.tre
├── all_gt_bs_rooted_{37,48,200}.tre, true_{37,48,200}.tre   # Additional rooted inputs with reference trees
├── sim.sh, sim_incomplete.sh, simphy/            # SimPhy simulation (binary and helpers)
├── test-stelarx-simulated.sh, run-bulk-simulated.sh          # Single and bulk simulated experiments
├── run-a10k.sh, collect-scores-a10k.sh                       # A10K experiments
├── run-bulk-standard.sh, collect-stats-standard.sh           # Standard/biological datasets (and baselines)
├── collect-stats-simulated.sh                                # Merge simulated statistics
├── sync-*-outputs.sh, upload-*-outputs.sh                    # Reproducibility mirrors and Hugging Face publication
├── download-bulk-simulated.sh, upload-bulk-simulated.sh      # Fetch / publish the simulated datasets themselves
├── clear-bulk-simulated.sh, clear-a10k.sh, clear-bulk-standard.sh   # Cleanup (always with --dry-run)
├── scripts/                          # Shared helpers (data-dir resolution, mirror primitives, setting names)
├── rf.py, analyze-dataset.py, clean.py, root_by_outgroups.py, process_unrooted.sh, extract-taxa.sh   # Tree utilities
├── test/                             # Regression, oracle, GPU, scalability and script tests (see Testing)
└── packaging/                        # Portable launcher used by build_portable.sh
```
