# STELAR-MP: Massively Parallel Phylogenetic Tree Inference

This project implements a massively parallel version of the STELAR algorithm for phylogenetic tree inference, with support for CPU and GPU computation.

## Prerequisites

- Java 11 or higher (Tested with 17, 21)
- Nvidia CUDA Toolkit 11.0 or higher
- NVIDIA GPU with compute capability 3.5 or higher
- Maven 3.6 or higher

We recommend checking whether the following is supported in your machine, and installing in case they do not exist currently.

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


This may require almost 7GB of disk space. Note that `nvcc` is only needed if you want to test with GPU parallelism.

2. Check whether java is there.

```bash
java -version
```

If not, install Java. We tested with JDK 21.

```bash
# Update package index
sudo apt update
# Check if OpenJDK 21 is available in default repositories
sudo apt install openjdk-21-jdk
```

Then check again,

```bash
java -version
```

Expected output:
openjdk 21.x.x ...


Finally set up JAVA_HOME,

```bash
# Set JAVA_HOME (adjust path if different)
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
# Make it persistent
echo 'export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64' >> ~/.bashrc
source ~/.bashrc
```

3. Install Maven

```bash
# Install Maven
sudo apt install maven
```

```bash
# Verify Maven installation and Java version
mvn -version
```

Expected output:
Apache Maven 3.x.x ...
Java version: 21.x.x ...



## Building the Project

1. Clone the repository:
```bash
git clone https://github.com/yourusername/stelar-mp.git
cd stelar-mp
```

2. Compile the CUDA code:
```bash
cd cuda
make
cd ..
```

3. Build the Java project:
```bash
mvn clean package
```

## Running the Program

Please note that, if you have made any change in the Java codebase, you must run the command,

```bash
chmod +x build.sh
./build.sh
```

This will re-compile and re-generate the jar file.

Use the provided `run.sh` script to run the program:

```bash
./run.sh <input_file> <output_file> [computation_mode]
```

Where:
- `<input_file>`: Path to the input gene tree file in Newick format
- `<output_file>`: Path where the output species tree will be written
- `[computation_mode]`: Optional computation mode (default: GPU_PARALLEL)
  - `CPU_SINGLE`: Single-threaded CPU computation
  - `CPU_PARALLEL`: Multi-threaded CPU computation
  - `GPU_PARALLEL`: GPU-accelerated computation

Example:
```bash
./run.sh input.tre output.tre GPU_PARALLEL
```

## Running in bulk for benchmarking

Please use `bulk_runner.sh` to run in bulk and to benchmark.

## Computation Modes

1. **CPU_SINGLE**: Uses a single CPU thread for weight calculation. This is the slowest but most memory-efficient mode.

2. **CPU_PARALLEL**: Uses multiple CPU threads for weight calculation, with the number of threads equal to the number of available CPU cores. This provides good performance on systems without a GPU.

3. **GPU_PARALLEL**: Uses NVIDIA GPU for weight calculation, providing massive parallelism for large datasets. This is the fastest mode when a compatible GPU is available.


## Generation of Simulated Datasets

To test STELAR-MP with large simulated datasets, we use Simphy (https://github.com/adamallo/SimPhy). Please look at the `simphy` directory to generate simulated datasets. Please read `cmd.txt` there.

## Calculation of RF rate

Please use `RF/getFpFn.py`

## Performance Considerations

- The GPU mode requires a CUDA-capable NVIDIA GPU with sufficient memory to hold the bipartition data.
- For very small datasets, CPU_PARALLEL might be faster than GPU_PARALLEL due to the overhead of data transfer to/from the GPU.
- The program will automatically fall back to CPU_PARALLEL mode if GPU computation fails.

## Troubleshooting

1. If you get a "CUDA error" message:
   - Make sure you have a compatible NVIDIA GPU
   - Verify that CUDA Toolkit is properly installed
   - Check that the GPU has enough memory for your dataset

2. If you get a "JNA error" message:
   - Make sure the CUDA library is in the system library path
   - Try running with `-Djava.library.path=cuda` to specify the library location

3. If the program is slow:
   - Try different computation modes to find the best for your system
   - For large datasets, ensure you have enough system memory
   - Consider using a more powerful GPU for better performance

## License

This project is licensed under the MIT License - see the LICENSE file for details. 