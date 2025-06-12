#!/bin/bash

# Compile the Java code
javac -sourcepath src -d bin $(find src -name "*.java")

# Define the input tree file paths
TREE_FILE_1="out.tre"
TREE_FILE_2="out-ori.tre"

# Run the CompareTrees program
java -cp bin CompareTrees "$TREE_FILE_1" "$TREE_FILE_2" 