#!/bin/bash

# Compile Java files to bin directory
# echo "Compiling Java files..."
javac -sourcepath src -d bin src/Main.java src/utils/Config.java

# Check if compilation was successful
if [ $? -eq 0 ]; then
    # echo ""
    echo "Compilation successful!"
    echo ""

    # Run the program with in.tre and out.tre
    # echo "Running phylogeny analysis..."
    java -cp bin Main -i in.tre -o out.tre
    
    # Check if execution was successful
    if [ $? -eq 1 ]; then
        echo "Error: Program execution failed!"
        exit 1
    fi
else
    echo "Error: Compilation failed!"
    exit 1
fi 