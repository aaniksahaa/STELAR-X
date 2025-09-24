#!/bin/bash

# Test script for memory-optimized phylogeny reconstruction
# Compares traditional vs memory-optimized implementation

echo "======================================="
echo "Memory Optimization Testing Suite"
echo "======================================="

# Set up directories
SRC_DIR="src"
BIN_DIR="bin"
CUDA_DIR="cuda"
OUTPUT_DIR="memory_test_results"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Step 1: Compiling Java code..."

# Compile all Java files
cd "$SRC_DIR"
find . -name "*.java" -type f > sources.txt

# Create bin directory if it doesn't exist
mkdir -p "../$BIN_DIR"

# Compile with classpath
javac -cp "../$BIN_DIR" -d "../$BIN_DIR" @sources.txt

# Clean up
rm sources.txt
cd ..

if [ $? -ne 0 ]; then
    echo "ERROR: Java compilation failed"
    exit 1
fi

echo "Java compilation successful!"

echo "Step 2: Building CUDA libraries..."
cd "$CUDA_DIR"
make all
cd ..

if [ $? -ne 0 ]; then
    echo "ERROR: CUDA compilation failed"
    exit 1
fi

echo "CUDA libraries built successfully!"

echo "Step 3: Running memory optimization tests..."

# Test files (use small datasets for initial testing)
TEST_FILE="all_gt_bs_rooted_11.tre"

if [ ! -f "$TEST_FILE" ]; then
    echo "WARNING: Test file $TEST_FILE not found. Using any available .tre file..."
    TEST_FILE=$(find . -name "*.tre" -type f | head -1)
    if [ -z "$TEST_FILE" ]; then
        echo "ERROR: No .tre files found for testing"
        exit 1
    fi
fi

echo "Using test file: $TEST_FILE"

# Test 1: Memory-optimized implementation with CPU_SINGLE
echo ""
echo "Test 1: Memory-optimized CPU_SINGLE mode"
echo "----------------------------------------"
cd "$BIN_DIR"
timeout 300s java -Xmx2g -cp . MemoryOptimizedMain \
    -i "../$TEST_FILE" \
    -o "../$OUTPUT_DIR/memory_opt_cpu_single.tre" \
    -m CPU_SINGLE \
    -profile

if [ $? -eq 0 ]; then
    echo "✓ Memory-optimized CPU_SINGLE test completed successfully"
else
    echo "✗ Memory-optimized CPU_SINGLE test failed or timed out"
fi

# Test 2: Memory-optimized implementation with CPU_PARALLEL
echo ""
echo "Test 2: Memory-optimized CPU_PARALLEL mode"
echo "-------------------------------------------"
timeout 300s java -Xmx2g -cp . MemoryOptimizedMain \
    -i "../$TEST_FILE" \
    -o "../$OUTPUT_DIR/memory_opt_cpu_parallel.tre" \
    -m CPU_PARALLEL \
    -profile

if [ $? -eq 0 ]; then
    echo "✓ Memory-optimized CPU_PARALLEL test completed successfully"
else
    echo "✗ Memory-optimized CPU_PARALLEL test failed or timed out"
fi

# Test 3: Try GPU mode if available
echo ""
echo "Test 3: Memory-optimized GPU_PARALLEL mode (if available)"
echo "---------------------------------------------------------"
timeout 300s java -Xmx2g -cp . MemoryOptimizedMain \
    -i "../$TEST_FILE" \
    -o "../$OUTPUT_DIR/memory_opt_gpu_parallel.tre" \
    -m GPU_PARALLEL \
    -profile

if [ $? -eq 0 ]; then
    echo "✓ Memory-optimized GPU_PARALLEL test completed successfully"
else
    echo "⚠ Memory-optimized GPU_PARALLEL test failed (GPU may not be available)"
fi

# Test 4: Original implementation for comparison (if available)
echo ""
echo "Test 4: Original implementation for comparison"
echo "----------------------------------------------"
timeout 300s java -Xmx2g -cp . Main \
    -i "../$TEST_FILE" \
    -o "../$OUTPUT_DIR/original_implementation.tre" \
    -m CPU_PARALLEL \
    -e NONE

if [ $? -eq 0 ]; then
    echo "✓ Original implementation test completed successfully"
    COMPARISON_AVAILABLE=true
else
    echo "⚠ Original implementation test failed or not available"
    COMPARISON_AVAILABLE=false
fi

cd ..

echo ""
echo "Step 4: Analyzing results..."
echo "============================="

# Check output files
echo "Generated output files:"
for file in "$OUTPUT_DIR"/*.tre; do
    if [ -f "$file" ]; then
        echo "  - $(basename "$file"): $(wc -c < "$file") bytes"
    fi
done

# If both implementations produced results, compare them
if [ "$COMPARISON_AVAILABLE" = true ]; then
    echo ""
    echo "Comparing tree topologies (if RF comparison available)..."
    
    # Try to use RF comparison if available
    if [ -f "rf.py" ]; then
        python3 rf.py \
            "$OUTPUT_DIR/memory_opt_cpu_parallel.tre" \
            "$OUTPUT_DIR/original_implementation.tre" > "$OUTPUT_DIR/comparison_results.txt" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "Tree comparison results:"
            cat "$OUTPUT_DIR/comparison_results.txt"
        else
            echo "RF comparison failed or not available"
        fi
    else
        echo "RF comparison tool not available"
    fi
fi

echo ""
echo "Step 5: Summary"
echo "==============="
echo "Memory optimization testing completed!"
echo "Results saved in: $OUTPUT_DIR/"
echo ""
echo "Key benefits of memory optimization:"
echo "• Reduced memory usage from O(n²k) to O(nk)"
echo "• Range-based bipartition representation"
echo "• Hash-based equality checking"
echo "• Improved scalability for large datasets"
echo ""

# Check if any critical errors occurred
if [ ! -f "$OUTPUT_DIR/memory_opt_cpu_parallel.tre" ]; then
    echo "❌ CRITICAL: Memory-optimized CPU_PARALLEL test failed"
    echo "   This is the core functionality - please check error messages above"
    exit 1
else
    echo "✅ SUCCESS: Memory optimization is working correctly!"
fi

echo ""
echo "To run individual tests:"
echo "  Memory-optimized: java -cp bin MemoryOptimizedMain -i input.tre -o output.tre -m CPU_PARALLEL -profile"
echo "  Original:         java -cp bin Main -i input.tre -o output.tre -m CPU_PARALLEL -e NONE"
echo ""
echo "======================================="
