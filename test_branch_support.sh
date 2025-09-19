#!/bin/bash

# Test script for branch support calculation
# This script demonstrates the new branch support functionality

echo "=== Testing Branch Support Calculation ==="

# Compile the project
echo "Compiling project..."
cd /home/aaniksahaa/research/STELAR-MP
./build.sh

if [ $? -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi

echo "Compilation successful!"

# Test with different branch support options
INPUT_FILE="all_gt_bs_rooted_11.tre"
OUTPUT_DIR="branch_support_test_results"

# Create output directory
mkdir -p $OUTPUT_DIR

echo ""
echo "=== Testing different branch support options ==="

# Test 1: No branch support (default)
echo "1. Testing without branch support..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/no_support.tre -m CPU_SINGLE
echo "   Output saved to: $OUTPUT_DIR/no_support.tre"

# Test 2: Posterior probabilities only
echo "2. Testing with posterior probabilities..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/posterior.tre -m CPU_SINGLE -s POSTERIOR
echo "   Output saved to: $OUTPUT_DIR/posterior.tre"

# Test 3: Branch lengths only
echo "3. Testing with branch lengths..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/lengths.tre -m CPU_SINGLE -s LENGTH
echo "   Output saved to: $OUTPUT_DIR/lengths.tre"

# Test 4: Both posterior and branch lengths
echo "4. Testing with both posterior and branch lengths..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/both.tre -m CPU_SINGLE -s BOTH
echo "   Output saved to: $OUTPUT_DIR/both.tre"

# Test 5: Detailed annotation
echo "5. Testing with detailed annotation..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/detailed.tre -m CPU_SINGLE -s DETAILED
echo "   Output saved to: $OUTPUT_DIR/detailed.tre"

# Test 6: P-values for polytomy test
echo "6. Testing with p-values..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/pvalues.tre -m CPU_SINGLE -s PVALUE
echo "   Output saved to: $OUTPUT_DIR/pvalues.tre"

# Test 7: All annotations
echo "7. Testing with all annotations..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/all.tre -m CPU_SINGLE -s ALL
echo "   Output saved to: $OUTPUT_DIR/all.tre"

# Test 8: Custom lambda parameter
echo "8. Testing with custom lambda parameter..."
java -cp bin Main -i $INPUT_FILE -o $OUTPUT_DIR/lambda_1.0.tre -m CPU_SINGLE -s POSTERIOR --lambda 1.0
echo "   Output saved to: $OUTPUT_DIR/lambda_1.0.tre"

echo ""
echo "=== Comparing Results ==="

echo "File sizes:"
ls -la $OUTPUT_DIR/*.tre

echo ""
echo "Sample outputs:"
for file in $OUTPUT_DIR/*.tre; do
    echo "--- $(basename $file) ---"
    head -n 1 $file
    echo ""
done

echo ""
echo "=== Testing Complete ==="
echo "All test results saved in: $OUTPUT_DIR/"
echo ""
echo "Usage examples:"
echo "  Basic usage:                java Main -i input.tre -o output.tre"
echo "  With posterior support:     java Main -i input.tre -o output.tre -s POSTERIOR"
echo "  With branch lengths:        java Main -i input.tre -o output.tre -s LENGTH"
echo "  With both:                  java Main -i input.tre -o output.tre -s BOTH"
echo "  Detailed information:       java Main -i input.tre -o output.tre -s DETAILED"
echo "  All annotations:            java Main -i input.tre -o output.tre -s ALL"
echo "  Custom lambda:              java Main -i input.tre -o output.tre -s POSTERIOR --lambda 1.0"
echo ""
echo "Branch support types:"
echo "  NONE      - No branch support (default)"
echo "  POSTERIOR - Posterior probabilities only"
echo "  LENGTH    - Branch lengths only"
echo "  BOTH      - Posterior probabilities and branch lengths"
echo "  DETAILED  - Detailed quartet frequencies and effective N"
echo "  PVALUE    - P-values for polytomy testing"
echo "  ALL       - All available information"
