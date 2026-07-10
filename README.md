# STELAR-X: Scaling Coalescent-Based Species Tree Inference to 100,000 Species and Beyond

**STELAR-X** is a scalable, statistically consistent summary method for coalescent-based species tree inference from large collections of gene trees. It combines compact bipartition encodings, fast weight precomputation, GPU-accelerated parallelism, and optimized dynamic programming to analyze datasets as large as **100,000 taxa × 1,000 genes in just 8.5 hours using 86 GB RAM**.

**The paper  STELAR-X has been accepted for presentation at [RECOMB 2026](https://recomb.org/recomb2026/) in Greece.**

## Citation

If you use **STELAR-X**, its source code, or results produced by it in your research, please cite our paper:

> Anik Saha and Md. Shamsuzzoha Bayzid.  
> **STELAR-X: Scaling Coalescent-Based Species Tree Inference to 100,000 Species and Beyond.**  
> Accepted at RECOMB 2026.  
> bioRxiv, 2025. https://doi.org/10.1101/2025.11.22.689894

```bibtex
@article{saha2025stelarx,
  title   = {{STELAR-X}: Scaling Coalescent-Based Species Tree
             Inference to 100,000 Species and Beyond},
  author  = {Saha, Anik and Bayzid, Md. Shamsuzzoha},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.11.22.689894},
  note    = {Accepted at RECOMB 2026}
}
```

> **Platform:** Developed and tested on **Ubuntu Linux**.

---

## Quick Start (Precompiled Release)

The fastest way to use STELAR-X. **No build tools needed — just Java.**

**Prerequisite:** Java 11+ (tested with OpenJDK 17 and 21)

```bash
# Install Java if not already installed (Ubuntu/Debian)
sudo apt update && sudo apt install -y openjdk-17-jdk

# Verify
java -version
```

### Download

You can download the lightweight, portable `.tar` archive from either source:

* **GitHub Releases:** [https://github.com/aaniksahaa/STELAR-X/releases](https://github.com/aaniksahaa/STELAR-X/releases)
  
  Example: v1.0.0 — [https://github.com/aaniksahaa/STELAR-X/releases/tag/v1.0.0](https://github.com/aaniksahaa/STELAR-X/releases/tag/v1.0.0)
  
* **Google Drive mirror:** [https://drive.google.com/drive/folders/1iyvmd__u_sCLZG1Z5TmzgOOuB1pVXlec?usp=sharing](https://drive.google.com/drive/folders/1iyvmd__u_sCLZG1Z5TmzgOOuB1pVXlec?usp=sharing)

```bash
# 1. Extract the prebuilt release archive
tar xzf stelar-x-1.0.0.tar.gz
cd stelar-x-1.0.0

# 2. Run on the included example (37-taxon dataset, 200 gene trees)
./stelar-x -i examples/all_gt_bs_rooted_37.tre -o examples/out_37.tre
```

That's it. An example gene trees file is included so you can verify it works immediately. Note that, you can run inference for any gene trees with relative or absolute path. GPU mode is used automatically if an NVIDIA GPU is detected; otherwise it falls back to CPU parallel.

```bash
# Force CPU parallel mode
./stelar-x -i examples/all_gt_bs_rooted_37.tre -o examples/out_37.tre --cpu-parallel

# Run on your own data
./stelar-x -i gene_trees.tre -o output.tre

# Score a known species tree against gene trees
./stelar-x -i gene_trees.tre -c species_tree.tre

# See all options
./stelar-x --help
```

**Optional** — add to PATH for system-wide access:

```bash
sudo ln -sf $(pwd)/stelar-x /usr/local/bin/stelar-x
stelar-x -i gene_trees.tre -o output.tre    # works from anywhere
```

---

## Building from Source

For developers who want to edit code, rebuild, or experiment.

### Prerequisites

| Dependency | Required | Version | Notes |
|------------|----------|---------|-------|
| **Java (JDK)** | Yes | 11+ | Tested with OpenJDK 17 and 21 |
| **NVIDIA CUDA Toolkit** | No | 11.0+ | For GPU acceleration; CPU fallback is automatic |
| **NVIDIA GPU** | No | Compute capability ≥ 7.0 | Required only for GPU mode |

> **Maven is NOT required.** The bundled Maven Wrapper (`mvnw`) handles it automatically.

### Install Java (if not installed)

```bash
sudo apt update && sudo apt install -y openjdk-17-jdk

# Verify
java -version
```

### Build

```bash
# 1. Clone the repository
git clone https://github.com/aaniksahaa/STELAR-X.git
cd STELAR-X

# 2. Build everything (one command)
# this install may take time in the first run
./install.sh
```

`install.sh` will:
- Check that Java 11+ (JDK) is available
- Compile the Java project into a fat JAR via the Maven Wrapper (no Maven installation needed)
- Compile the CUDA kernels if `nvcc` is available (otherwise, uses the pre-built library or skips gracefully)

```bash
# Build options
./install.sh                  # Auto-detects CUDA
./install.sh --no-cuda        # Skip CUDA compilation (CPU-only)
./install.sh --force-cuda     # Fail if CUDA compilation fails
./install.sh --clean          # Clean rebuild from scratch
```

### Run

```bash
./run.sh -i all_gt_bs_rooted_37.tre -o out_37.tre
```

After editing source code, rebuild with:

```bash
./build.sh
```

### Build Precompiled Release Archive

To package a standalone distribution (for sharing with others):

```bash
./dist.sh                     # Produces dist/stelar-x-1.0.0.tar.gz (~3.3 MB)
./dist.sh --skip-build        # Package existing build without rebuilding
```

---

## Usage

### Inference Mode (default)

Infer a species tree from a set of gene trees:

```bash
./run.sh -i <gene_trees.tre> -o <output.tre> [options]
```

### Score Mode

Calculate the triplet score of a known species tree against gene trees:

```bash
./run.sh -i <gene_trees.tre> -c <species_tree.tre> [options]
```

### With Performance Monitoring

Records running time, peak CPU RAM, and peak GPU VRAM:

```bash
./run-with-monitor.sh -i <gene_trees.tre> -o <output.tre>
```

### Examples

```bash
# GPU mode (auto-detected)
./run.sh -i all_gt_bs_rooted_37.tre -o out-37.tre

# Explicit CPU parallel
./run.sh -i all_gt_bs_rooted_37.tre -o out-37.tre --cpu-parallel

# GPU with cross-tree recombination (expansion)
./run.sh -i all_gt_bs_rooted_37.tre -o out-37.tre --gpu --expansion

# Score-only mode
./run.sh -i all_gt_bs_rooted_37.tre -c known_species.tre

# Custom Java heap size
./run.sh -i large_dataset.tre -o out.tre --xms 8g --xmx 256g

# Monitored run with stats output
./run-with-monitor.sh -i all_gt_bs_rooted_48.tre -o out-48.tre
```

---

## Command-Line Parameters

### Required Parameters

| Flag | Long Form | Description |
|------|-----------|-------------|
| `-i` | `--input <file>` | Input gene trees file in Newick format |
| `-o` | `--output <file>` | Output species tree file (inference mode) |

> In **score mode**, use `-c` instead of `-o`.

### Computation Mode

| Flag | Long Form | Description |
|------|-----------|-------------|
| | `--gpu` | GPU parallel mode (**default** if NVIDIA GPU is detected) |
| | `--cpu-parallel` | CPU multi-threaded mode (default if no GPU) |
| | `--cpu` | CPU single-threaded mode |
| `-m` | `--mode <mode>` | Explicit mode: `GPU_PARALLEL`, `CPU_PARALLEL`, `CPU_SINGLE` |

### Optional Parameters

| Flag | Long Form | Description | Default |
|------|-----------|-------------|---------|
| `-c` | `--score <file>` | Score-only mode: calculate triplet score for a given species tree | — |
| `-e` | `--expansion` | Enable cross-tree recombination (mixed bipartitions) | Off |
| `-v` | `--verbose` | Verbose expansion output | Off |
| | `--xms <size>` | Java minimum heap size | `4g` |
| | `--xmx <size>` | Java maximum heap size | `128g` |
| `-h` | `--help` | Show help message | — |

### Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `STELAR_XMS` | Default minimum Java heap size | `STELAR_XMS=8g ./run.sh ...` |
| `STELAR_XMX` | Default maximum Java heap size | `STELAR_XMX=256g ./run.sh ...` |

---

## Testing with Simulated Datasets

We use [SimPhy](https://github.com/adamallo/SimPhy) to generate simulated datasets for large-scale benchmarking.

### Generate and test simulated data

```bash
# Generate simulated dataset (100 taxa, 200 gene trees)
./sim.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh

# Run STELAR-X on the simulated data
./test-stelar-simulated.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

### Specifying a base directory

```bash
./sim.sh -b $HOME/research -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -b $HOME/research -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

Here, `-b` should point to the directory *containing* the STELAR-X folder.

### Custom SimPhy data directory

```bash
./sim.sh -b $HOME/research --simphy-data-dir /dev/shm/data -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -b $HOME/research --simphy-data-dir /dev/shm/data -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

### Optional tools

These are not required for STELAR-X itself, but are useful for evaluation:

```bash
pip install dendropy        # RF distance calculation
sudo apt install -y time    # Detailed time/memory monitoring (for run-with-monitor.sh)
```

---

## Building Upon the Codebase

### Modifying Java Code

After editing any Java source files in `src/`, rebuild:

```bash
./build.sh
```

This recompiles the project and regenerates the fat JAR with all dependencies bundled. (`build.sh` and `install.sh` are equivalent — use whichever feels right.)

### Modifying the CUDA Kernel

After editing `cuda/weight_calc.cu`, rebuild:

```bash
cd cuda
make clean && make
cd ..
```

Or simply re-run `./build.sh`, which handles both Java and CUDA compilation.

### Installing the CUDA Toolkit (if not already installed)

```bash
# Check if nvcc is available
nvcc --version

# If missing (Ubuntu/Debian):
sudo apt update
sudo apt install -y nvidia-cuda-toolkit

# Verify
nvcc --version
```

---

## Implementation Notes

- We use the hash functions **addition** and **bitwise XOR**
- We use 2^64 as the modulus since it arises naturally from the wrap-around behavior of the corresponding datatype, providing greater efficiency
- For the single-element hash function, we use a mixing function comprising multiplication and bit-shift operations
- The CUDA kernel is compiled for GPU architectures `sm_70` through `sm_90`, covering NVIDIA Volta, Turing, Ampere, Ada Lovelace, and Hopper generations

---

## Troubleshooting

| Problem | Possible Cause | Solution |
|---------|---------------|----------|
| `JAR not found` | Project not built | Run `./install.sh` |
| CUDA library not loaded | Missing GPU or CUDA runtime | Runs in CPU mode automatically; or install CUDA toolkit |
| `OutOfMemoryError` | Dataset too large for default heap | Increase heap: `--xms 8g --xmx 256g` |
| GPU out of memory | Dataset too large for GPU VRAM | Use `--cpu-parallel` instead |
| `nvcc` not found | CUDA toolkit not installed | `sudo apt install -y nvidia-cuda-toolkit` |
| Slow performance | Running in single-threaded mode | Use `--cpu-parallel` or `--gpu` |

---

## Project Structure

```
STELAR-X/
├── install.sh              # One-command build (Java + CUDA)
├── build.sh                # Same as install.sh (convenient alias for rebuilds)
├── run.sh                  # Main launcher
├── run-with-monitor.sh     # Launcher with performance monitoring
├── dist.sh                 # Build standalone release archive
├── mvnw                    # Maven Wrapper (no Maven installation needed)
├── pom.xml                 # Maven project configuration
├── src/                    # Java source code
│   ├── Main.java
│   ├── core/               # Inference engine, weight calculators, DP
│   ├── preprocessing/      # Gene tree parsing
│   ├── tree/               # Tree and bipartition data structures
│   ├── taxon/              # Taxon representation
│   └── utils/              # Configuration, threading, hashing utilities
├── cuda/
│   ├── weight_calc.cu      # CUDA GPU kernels
│   ├── libweight_calc.so   # Pre-built CUDA library
│   └── Makefile
└── target/
    └── stelar-x-*.jar      # Built fat JAR (after install.sh)
```
