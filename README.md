# STELAR-X: Scaling Coalescent-Based Species Tree Inference to 100,000 Species and Beyond

An extended implementation of the STELAR (Species Tree Estimation by maximizing tripLet AgReement) integrating permutation-invariant double hashing and massive parallelism. STELAR-X is capable of analysing ultra-large phylogenetic datasets.

## Prerequisites

- Java 11 or higher (Tested with 17, 21)
- Nvidia CUDA Toolkit 11.0 or higher
- NVIDIA GPU with compute capability 3.5 or higher
- Maven 3.6 or higher

We recommend checking whether the following is supported in your machine, and installing in case they do not exist currently.

1. Install Java. We tested with JDK 21.

```bash
# Update package index
sudo apt update
# Check if OpenJDK 21 is available in default repositories
sudo apt install openjdk-21-jdk
```

```bash
java -version
```

Sample output:
```
openjdk version "21.0.8" 2025-07-15
OpenJDK Runtime Environment (build 21.0.8+9-Ubuntu-0ubuntu122.04.1)
OpenJDK 64-Bit Server VM (build 21.0.8+9-Ubuntu-0ubuntu122.04.1, mixed mode, sharing)
```

2. Install Maven

```bash
sudo apt install maven
```

```bash
mvn -version
```

Sample output:
```
Apache Maven 3.8.7
Maven home: /usr/share/maven
Java version: 21.0.8, vendor: Ubuntu, runtime: /usr/lib/jvm/java-21-openjdk-amd64
Default locale: en_US, platform encoding: UTF-8
OS name: "linux", version: "6.14.0-35-generic", arch: "amd64", family: "unix"
```


## Setting up the Project

1. Clone the repository:
```bash
git clone https://github.com/aaniksahaa/STELAR-X.git
cd STELAR-X
```

2. Grant permissions to the scripts.

```
chmod +x build.sh run.sh run-with-monitor.sh sim.sh test-stelar-simulated.sh
```

3. Build the Java project.

```bash
./build.sh
```

This will re-compile and re-generate the jar file.

## Running the Program

Use the provided `run.sh` script to run the program:

```bash
./run.sh <input_file> <output_file> 
```

Where:
- `<input_file>`: Path to the input gene tree file in Newick format
- `<output_file>`: Path where the output species tree will be written


Example:
```bash
./run.sh all_gt_bs_rooted_37.tre out-37.tre
```
```bash
./run.sh avian-48-gt.tre out-avian-48.tre
```

To also monitor the running time and memory usage, you may use `run-with-monitor.sh`. This may require installing `time`.

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

## Testing with Simulated Datasets

To test STELAR-X with large simulated datasets, we use Simphy (https://github.com/adamallo/SimPhy).

Example Commands for Testing STELAR-X with simulated datasets.

```bash
./sim.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

The options `-t` and `-g` indicate the number of taxa and gene trees respectively.

You may also specify a base-dir as follows.

```bash
./sim.sh -b $HOME/research -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -b $HOME/research -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

Here, the option `-b` expects the base directory where STELAR-X is set up. For instance, if you have set up it at `$HOME/research/STELAR-X`, then "-b" should be set as `$HOME/research`. 

To specify a directory where to store the simulated data (instead of the default location), you may use the option `--simphy-data-dir`

```bash
./sim.sh -b $HOME/research --simphy-data-dir /dev/shm/data -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -rs 1 --fresh
./test-stelar-simulated.sh -b $HOME/research --simphy-data-dir /dev/shm/data -t 100 -g 200 --sb 0.000001 --spmin 100000 --spmax 200000 -r R1 --fresh
```

## Building Upon the Codebase

### Changing the Java code

If the Java codebase is modified, it is sufficient to recompile the sources files and build the Java classes again as follows.

```
./build.sh
```

### Changing the CUDA kernel

1. Check whether nvcc is there.

```bash
nvcc --version
```

If not, install with,

```bash
sudo apt update
sudo apt-get install nvidia-cuda-toolkit
```

And then check,

```bash
nvcc --version
```

2. Make your changes to the Kernel code. Then build the CUDA kernel again.

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

## Troubleshooting

1. If you get a "CUDA error" message:
   - Make sure you have a compatible NVIDIA GPU
   - Verify that CUDA Toolkit is properly installed
   - Check that the GPU has enough memory for your dataset