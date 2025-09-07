#!/bin/bash

# Test script for bipartition expansion functionality

echo "Testing STELAR-MP Bipartition Expansion"
echo "======================================="

# Compile the project
echo "Compiling project..."
cd src
javac -cp ".:../cuda" *.java */*.java */*/*.java 2>/dev/null
if [ $? -eq 0 ]; then
    echo "✓ Compilation successful"
else
    echo "✗ Compilation failed"
    exit 1
fi

# Test with different expansion methods
echo ""
echo "Testing expansion methods..."

# Test 1: No expansion
echo "Test 1: No expansion"
java -cp ".:../cuda" Main -i ../in.tre -o test_no_expansion.tre --no-expansion -v

# Test 2: Distance-only expansion
echo ""
echo "Test 2: Distance-only expansion"
java -cp ".:../cuda" Main -i ../in.tre -o test_distance_only.tre -e DISTANCE_ONLY -v

# Test 3: Consensus-only expansion
echo ""
echo "Test 3: Consensus-only expansion"
java -cp ".:../cuda" Main -i ../in.tre -o test_consensus_only.tre -e CONSENSUS_ONLY -v

# Test 4: Both distance and consensus
echo ""
echo "Test 4: Distance and consensus expansion"
java -cp ".:../cuda" Main -i ../in.tre -o test_distance_consensus.tre -e DISTANCE_CONSENSUS -v

echo ""
echo "Testing completed. Check output files for results."
echo "Output files:"
echo "  - test_no_expansion.tre"
echo "  - test_distance_only.tre"
echo "  - test_consensus_only.tre"
echo "  - test_distance_consensus.tre"
