#!/bin/bash

# Main2 function
run_main() {
    make clean
    make main
    ./main
    mv V* data/
    mv O* data/
}
# Main function
run_main2() {
    make clean
    make main2
    ./main2
    mv k*.txt data/
}

# Check inptu
if [ "$1" == "main" ]; then
    run_main
elif [ "$1" == "main2" ]; then
    run_main2
else
    echo "Usage: $0 {main|main2}"
    exit 1
fi
