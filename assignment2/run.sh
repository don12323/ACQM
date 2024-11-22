#!/bin/bash

# Check if an argument is provided
if [ -z "$1" ]; then
  echo "Usage: $0 {main|main2}"
  exit 1
fi

# Clean and build the project
make clean
make $1

# Check if make was successful
if [ $? -eq 0 ]; then
  # Run the specified executable
  if [ "$1" == "main" ]; then
          ./main
  elif [ "$1" == "main2" ]; then
	  ./main2
	  

  else
    echo "Invalid argument: $1. Use 'main' or 'main2'."
    exit 1
  fi
else
  echo "error"
  exit 1
fi

