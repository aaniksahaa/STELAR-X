#!/usr/bin/env python3

import subprocess
import time
import os
import sys
from datetime import datetime

# Configuration
INPUT_FILE = "all_gt_bs_rooted_100.tre"  # Change this to your input file
OUTPUT_DIR = "benchmark_results"
PARALLEL_OUTPUT = os.path.join(OUTPUT_DIR, "parallel_out.tre")
ORIGINAL_OUTPUT = os.path.join(OUTPUT_DIR, "original_out.tre")
ORIGINAL_JAR = "/mnt/H/Research/STELAR-extension/STELAR/STELAR.jar"

# Computation modes to test
COMPUTATION_MODES = ["CPU_SINGLE", "CPU_PARALLEL", "GPU_PARALLEL"]

def format_time(seconds):
    """Format time in MM:SS.mmm format"""
    minutes = int(seconds // 60)
    seconds = seconds % 60
    milliseconds = int((seconds - int(seconds)) * 1000)
    return f"{minutes:02d}:{int(seconds):02d}.{milliseconds:03d}"

def run_command(cmd):
    """Run a command and return its output"""
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)

def ensure_output_dir():
    """Ensure the output directory exists"""
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

def run_benchmark(mode):
    """Run benchmark for a specific computation mode"""
    output_file = os.path.join(OUTPUT_DIR, f"{mode.lower()}_out.tre")
    
    print(f"\nRunning {mode} version...")
    start_time = time.time()
    result = run_command(f"./run.sh {INPUT_FILE} {output_file} {mode}")
    execution_time = time.time() - start_time
    
    if result.returncode != 0:
        print(f"{mode} version failed!")
        print(result.stderr)
        return None
    
    print(f"{mode} version completed in {format_time(execution_time)}")
    return execution_time

def main():
    print("=== STELAR-MP Benchmark Script ===")
    print(f"Input file: {INPUT_FILE}")
    print()

    # Ensure output directory exists
    ensure_output_dir()

    # Check if input file exists
    if not os.path.exists(INPUT_FILE):
        print(f"Error: Input file not found at {INPUT_FILE}")
        sys.exit(1)

    # Run benchmarks for each mode
    results = {}
    for mode in COMPUTATION_MODES:
        time_taken = run_benchmark(mode)
        if time_taken is not None:
            results[mode] = time_taken

    # Print performance comparison
    if results:
        print("\n=== Performance Comparison ===")
        baseline = results.get("CPU_SINGLE")
        if baseline:
            print(f"CPU_SINGLE (baseline): {format_time(baseline)}")
            for mode, time_taken in results.items():
                if mode != "CPU_SINGLE":
                    speedup = baseline / time_taken
                    print(f"{mode}: {format_time(time_taken)} (Speedup: {speedup:.2f}x)")

    # Output files
    print("\n=== Output Files ===")
    for mode in COMPUTATION_MODES:
        output_file = os.path.join(OUTPUT_DIR, f"{mode.lower()}_out.tre")
        if os.path.exists(output_file):
            print(f"{mode} output: {output_file}")

if __name__ == "__main__":
    main() 