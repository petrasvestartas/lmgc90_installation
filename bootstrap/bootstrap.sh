#!/bin/bash
# LMGC90 Bootstrap Script
# Automatically installs all dependencies locally without admin privileges

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
BOOTSTRAP_DIR="$SCRIPT_DIR"
BUILD_DIR="$BOOTSTRAP_DIR/build"
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "4")

echo -e "${BLUE}🚀 LMGC90 Dependency Bootstrap${NC}"
echo -e "${BLUE}================================${NC}"
echo ""
echo "This script will automatically install all LMGC90 dependencies locally."
echo "No admin privileges required - everything installs to the project directory."
echo ""
echo -e "${YELLOW}System Information:${NC}"
echo "  OS: $(uname -s)"
echo "  Architecture: $(uname -m)"
echo "  CPU Cores: $JOBS"
echo "  Project Root: $PROJECT_ROOT"
echo "  Bootstrap Dir: $BOOTSTRAP_DIR"
echo ""

# Check prerequisites
echo -e "${YELLOW}📋 Checking Prerequisites...${NC}"

# Check for basic build tools
if ! command -v make &> /dev/null; then
    echo -e "${RED}❌ make not found. Please install Xcode Command Line Tools:${NC}"
    echo -e "${BLUE}xcode-select --install${NC}"
    exit 1
else
    echo -e "${GREEN}✅ make found${NC}"
fi

if ! command -v gcc &> /dev/null; then
    echo -e "${RED}❌ gcc not found. Please install Xcode Command Line Tools:${NC}"
    echo -e "${BLUE}xcode-select --install${NC}"
    exit 1
else
    echo -e "${GREEN}✅ gcc found ($(gcc --version | head -n1))${NC}"
fi

# Check CMake - install locally if missing
CMAKE_EXECUTABLE=""
if command -v cmake &> /dev/null; then
    CMAKE_VERSION=$(cmake --version | head -n1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
    CMAKE_MAJOR=$(echo $CMAKE_VERSION | cut -d. -f1)
    CMAKE_MINOR=$(echo $CMAKE_VERSION | cut -d. -f2)
    
    if [ "$CMAKE_MAJOR" -lt 3 ] || ([ "$CMAKE_MAJOR" -eq 3 ] && [ "$CMAKE_MINOR" -lt 20 ]); then
        echo -e "${YELLOW}⚠️  CMake version $CMAKE_VERSION found. Minimum 3.20 required.${NC}"
        echo -e "${BLUE}ℹ️  Will install CMake 3.28 locally.${NC}"
        CMAKE_EXECUTABLE=""
    else
        echo -e "${GREEN}✅ CMake $CMAKE_VERSION found${NC}"
        CMAKE_EXECUTABLE="cmake"
    fi
else
    echo -e "${YELLOW}⚠️  CMake not found. Will install CMake 3.28 locally.${NC}"
    CMAKE_EXECUTABLE=""
fi

# Check Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}❌ Python3 not found. Please install Python 3.8 or later.${NC}"
    exit 1
else
    PYTHON_VERSION=$(python3 --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
    echo -e "${GREEN}✅ Python $PYTHON_VERSION found${NC}"
fi

# Check internet connection
if ! ping -c 1 github.com &> /dev/null; then
    echo -e "${RED}❌ No internet connection. Bootstrap requires internet to download dependencies.${NC}"
    exit 1
else
    echo -e "${GREEN}✅ Internet connection available${NC}"
fi

echo ""

# Function to install CMake locally
install_cmake_locally() {
    echo -e "${YELLOW}📥 Installing CMake 3.28 locally...${NC}"
    
    mkdir -p "$BUILD_DIR/cmake_bootstrap"
    cd "$BUILD_DIR/cmake_bootstrap"
    
    # Download CMake source
    if [ ! -f "cmake-3.28.1.tar.gz" ]; then
        echo -e "${BLUE}Downloading CMake source...${NC}"
        curl -L -o cmake-3.28.1.tar.gz https://github.com/Kitware/CMake/releases/download/v3.28.1/cmake-3.28.1.tar.gz
    fi
    
    # Extract
    if [ ! -d "cmake-3.28.1" ]; then
        echo -e "${BLUE}Extracting CMake...${NC}"
        tar -xzf cmake-3.28.1.tar.gz
    fi
    
    # Configure and build
    cd cmake-3.28.1
    if [ ! -f "Makefile" ]; then
        echo -e "${BLUE}Configuring CMake...${NC}"
        ./configure --prefix="$BUILD_DIR/cmake_local" --parallel=$JOBS
    fi
    
    echo -e "${BLUE}Building CMake (this may take 10-15 minutes)...${NC}"
    make -j$JOBS
    make install
    
    CMAKE_EXECUTABLE="$BUILD_DIR/cmake_local/bin/cmake"
    export PATH="$BUILD_DIR/cmake_local/bin:$PATH"
    
    echo -e "${GREEN}✅ CMake installed locally: $CMAKE_EXECUTABLE${NC}"
    cd "$BOOTSTRAP_DIR"
}

# Install CMake if needed
if [ -z "$CMAKE_EXECUTABLE" ]; then
    install_cmake_locally
fi

# Confirm with user
echo -e "${YELLOW}📦 Dependencies to be installed:${NC}"
echo "  • CMake 3.28 (if not available)"
echo "  • GCC 14.2.0 (C/C++/Fortran compilers)"
echo "  • HDF5 1.14.3 (with Fortran support)"
echo "  • OpenBLAS 0.3.26 (BLAS/LAPACK)"
echo "  • SWIG 3.0.12 (Python bindings)"
echo "  • Python packages (numpy, scipy, matplotlib, vtk, compas)"
echo ""
echo -e "${YELLOW}💾 Estimated disk space: ~3-4 GB${NC}"
echo -e "${YELLOW}⏱️  Estimated time: 45-75 minutes (depending on system)${NC}"
echo ""

read -p "Continue with bootstrap? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${YELLOW}Bootstrap cancelled.${NC}"
    exit 0
fi

echo ""
echo -e "${BLUE}🔧 Starting Bootstrap Process...${NC}"

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure CMake
echo -e "${YELLOW}📝 Configuring CMake...${NC}"
$CMAKE_EXECUTABLE "$BOOTSTRAP_DIR" \
    -DBOOTSTRAP_BUILD_JOBS="$JOBS" \
    -DCMAKE_BUILD_TYPE=Release

if [ $? -ne 0 ]; then
    echo -e "${RED}❌ CMake configuration failed!${NC}"
    exit 1
fi

# Build all dependencies
echo ""
echo -e "${YELLOW}🔨 Building Dependencies...${NC}"
echo -e "${BLUE}This will take a while. Logs are saved to build directory.${NC}"
echo ""

$CMAKE_EXECUTABLE --build . --target bootstrap_complete --parallel "$JOBS"

if [ $? -ne 0 ]; then
    echo ""
    echo -e "${RED}❌ Bootstrap build failed!${NC}"
    echo -e "${YELLOW}Check the build logs in: $BUILD_DIR${NC}"
    exit 1
fi

echo ""
echo -e "${GREEN}🎉 Bootstrap Complete!${NC}"
echo ""
echo -e "${YELLOW}📋 Next Steps:${NC}"
echo "1. Activate the environment:"
echo -e "   ${BLUE}source $BUILD_DIR/activate_bootstrap.sh${NC}"
echo ""
echo "2. Build LMGC90:"
echo -e "   ${BLUE}cd $PROJECT_ROOT${NC}"
echo -e "   ${BLUE}mkdir build_bootstrap && cd build_bootstrap${NC}"
echo -e "   ${BLUE}cmake -DCMAKE_TOOLCHAIN_FILE=$BUILD_DIR/lmgc90_toolchain.cmake ..${NC}"
echo -e "   ${BLUE}make -j$JOBS${NC}"
echo ""
echo "3. Test LMGC90:"
echo -e "   ${BLUE}PYTHONPATH=\$PWD python3 -c \"import pylmgc90; print('✅ SUCCESS!')\"${NC}"
echo ""
echo -e "${GREEN}All dependencies are now installed locally in:${NC}"
echo -e "${BLUE}$BUILD_DIR/deps${NC}"
echo ""
echo -e "${YELLOW}💡 Tip: Add the activation command to your shell profile for convenience!${NC}"
