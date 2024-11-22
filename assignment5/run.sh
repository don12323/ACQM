#!/bin/bash

# Function to compile and run main
run_main() {
    make clean
    make main
    ./main
    mv V* data/
    mv O* data/
}

# Function to compile and run main2
run_main2() {
    make clean
    make main2
    ./main2
    mv k*.txt data/
}

# Check for input parameter
if [ "$1" == "main" ]; then
    run_main
elif [ "$1" == "main2" ]; then
    run_main2
else
    echo "Usage: $0 {main|main2}"
    exit 1
fi
