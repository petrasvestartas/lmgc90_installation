#!/bin/bash
# Build script for OCCT example project

set -e # Exit on any error

# Path to build directory
BUILD_DIR="./build"

# Create build directory if it doesn't exist
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo "=== Configuring OCCT example project ==="
# Configure with CMake
cmake ..

echo "=== Building OCCT example project ==="
# Determine the number of CPU cores for parallel build
NUM_CORES=$(nproc 2>/dev/null || echo 2)
echo "Using $NUM_CORES CPU cores for build"

# Build the project
cmake --build . --config Release --parallel $NUM_CORES

echo "=== Build completed successfully ==="
echo "You can run the example with: ./build/write_step"
