#!/bin/bash

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <gene_trees_file> <species_tree_file>"
    echo "Example: $0 all_gt.tre s_tree.trees"
    exit 1
fi

# Run the analyze-dataset.py script and capture output
output=$(python analyze-dataset.py "$1" "$2" 2>&1)

# Display the full output first
echo "$output"

# Extract gt_gt (Average pairwise RF distance gene trees)
gt_gt=$(echo "$output" | grep "Average pairwise discordance (gene trees):" | awk '{print $NF}')

# Extract gt_st (Average RF distance gene trees vs species tree)
gt_st=$(echo "$output" | grep "Average discordance (gene trees vs species tree):" | awk '{print $NF}')

# Display the extracted results
echo ""
echo "============================================================"
echo "EXTRACTED VALUES:"
echo "============================================================"
echo "gt_gt: $gt_gt"
echo "gt_st: $gt_st"
