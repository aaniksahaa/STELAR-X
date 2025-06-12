#!/bin/bash

# Compile all Java files in src directory to bin directory
javac -sourcepath src -d bin $(find src -name "*.java")

input="all_gt_bs_rooted_100.tre"
# input="in.tre"

# Check if compilation was successful
if [ $? -eq 0 ]; then
    # echo "Compilation successful!"
    # echo ""

    # Run the program with in.tre and out.tre
    java -cp bin Main -i ${input} -o out.tre
    
    # Check if execution was successful
    if [ $? -eq 1 ]; then
        echo "Error: Program execution failed!"
        exit 1
    fi
else
    echo "Error: Compilation failed!"
    exit 1
fi 