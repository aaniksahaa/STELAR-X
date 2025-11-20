# STELAR-X: Scaling Coalescent-Based Species Tree Inference to 100,000 Species and Beyond

**STELAR-X** is a highly scalable, statistically consistent summary method for species tree inference that reconstructs species trees from large collections of gene trees under the multispecies coalescent model. It achieves near-input-size O(nk) memory usage through a redesigned computational framework built on compact integer-tuple bipartition encodings, fast precomputation of bipartition weights, and GPU-accelerized parallelism, all integrated into an optimized dynamic programming pipeline. With this combination of algorithmic engineering and parallel computation, STELAR-X delivers unprecedented scalability—analyzing 100,000 taxa × 1,000 genes in ~8.5 hours using 86 GB RAM, and 1,000 taxa × 100,000 genes in ~4 minutes using 106 GB RAM. With modest multi-day runtimes and additional memory, it is expected to scale to even larger phylogenomic datasets.

# Prerequisites

STELAR-X requires the following software:

* **Java 11 or higher** (tested with OpenJDK 17 and 21)
* **NVIDIA CUDA Toolkit 11.0 or higher**
* **NVIDIA GPU** with compute capability ≥ 3.5
* **Maven 3.6 or higher**

Before proceeding, please ensure that these dependencies are available on your system.

### 1. Install Java (tested with JDK 17, 21)

```bash
sudo apt update
sudo apt install -y openjdk-21-jdk
```

Verify the installation:

```bash
java -version
```

Sample Output:

```
openjdk version "21.0.8" 2025-07-15
OpenJDK Runtime Environment (build 21.0.8+9-Ubuntu-0ubuntu122.04.1)
OpenJDK 64-Bit Server VM (build 21.0.8+9-Ubuntu-0ubuntu122.04.1, mixed mode, sharing)
```

### 2. Install Maven

```bash
sudo apt install -y maven
```

Check the version:

```bash
mvn -version
```

Sample Output:
```
Apache Maven 3.8.7
Maven home: /usr/share/maven
Java version: 21.0.8, vendor: Ubuntu, runtime: /usr/lib/jvm/java-21-openjdk-amd64
Default locale: en_US, platform encoding: UTF-8
OS name: "linux", version: "6.14.0-35-generic", arch: "amd64", family: "unix"
```

### 3. Install Dendropy

```bash
pip install dendropy
```

---

# Setting Up the Project

### 1. Clone the repository

```bash
git clone https://github.com/aaniksahaa/STELAR-X.git
cd STELAR-X
```

### 2. Grant permission to helper scripts

```bash
chmod +x build.sh run.sh run-with-monitor.sh sim.sh test-stelar-simulated.sh
```

### 3. Build the Java project

```bash
./build.sh
```

This recompiles the Java sources and regenerates the executable JAR.

---

# Running STELAR-X

The recommended way to run the program is through the `run.sh` wrapper:

```bash
./run.sh <input_file> <output_file>
```

Where:

* `<input_file>` — Newick-formatted gene tree file
* `<output_file>` — output path for the inferred species tree

Examples:

```bash
./run.sh all_gt_bs_rooted_37.tre out-37.tre
```
```bash
./run.sh avian-48-gt.tre out-avian-48.tre
```

To record running time and memory usage, use:

```bash
./run-with-monitor.sh <input_file> <output_file>
```

If `time` is not installed:

```bash
sudo apt install -y time
```

Example:

```bash
./run-with-monitor.sh all_gt_bs_rooted_37.tre out-37.tre
```
```bash
./run-with-monitor.sh avian-48-gt.tre out-avian-48.tre
```

---

# Testing With Simulated Datasets

We use **SimPhy** to generate simulated datasets for large-scale benchmarking:
[https://github.com/adamallo/SimPhy](https://github.com/adamallo/SimPhy)

Examples:

```bash
./sim.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

`-t` specifies the number of taxa, and `-g` the number of gene trees.

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

---

# Building Upon the Codebase

## Modifying Java Code

After editing any Java source files, rebuild using:

```bash
./build.sh
```

This recompiles the project and regenerates the JAR.

## Modifying the CUDA Kernel

### 1. Check whether `nvcc` is installed

```bash
nvcc --version
```

If missing:

```bash
sudo apt update
sudo apt-get install nvidia-cuda-toolkit
```

Confirm installation:

```bash
nvcc --version
```

### 2. Rebuild the CUDA kernels

After editing `cuda/*.cu`, rebuild:

```
# Build CUDA code
echo -e "\n${YELLOW}Building CUDA code...${NC}"
cd cuda
make clean
make
if [ $? -ne 0 ]; then
    echo -e "${RED}CUDA compilation failed!${NC}"
    exit 1
fi
cd ..
echo -e "${GREEN}CUDA compilation successful${NC}"
```

---

# Troubleshooting

1. **CUDA errors** may indicate:

* missing or incompatible NVIDIA GPU
* incorrect or incomplete CUDA Toolkit installation
* insufficient GPU memory for the dataset

