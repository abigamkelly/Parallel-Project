#!/bin/bash

# Run the first Python script
python parallel_distance.py
if [ $? -ne 0 ]; then
    echo "Error: parallel_distance.py failed to run"
    exit 1
fi

# Compile the C++ code
g++ -O3 -fopenmp -shared -o c_functions.so -fPIC parallel_c_functions.cpp
if [ $? -ne 0 ]; then
    echo "Error: Failed to compile parallel_c_functions.cpp"
    exit 1
fi

# Run the second Python script
python parallel_colocation.py
if [ $? -ne 0 ]; then
    echo "Error: parallel_colocation.py failed to run"
    exit 1
fi

echo "All tasks completed successfully."
