#!/bin/bash
# LMGC90 Simple Installation Script
# One-command installation for macOS ARM64

set -e

echo "🚀 LMGC90 Simple Installation"
echo "=============================="
echo "This will install all dependencies in ~/spack"
echo ""

# Configuration
SPACK_ROOT="$HOME/spack"
LMGC90_ENV="lmgc90_env"
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "4")

echo "📋 System Info:"
echo "  OS: $(uname -s)"
echo "  Arch: $(uname -m)"
echo "  Cores: $JOBS"
echo "  Spack: $SPACK_ROOT"
echo ""

# Install Spack if not present
if [ ! -d "$SPACK_ROOT" ]; then
    echo "📥 Installing Spack..."
    git clone -c feature.manyFiles=true https://github.com/spack/spack.git "$SPACK_ROOT"
fi

# Source Spack
. "$SPACK_ROOT/share/spack/setup-env.sh"

# Create environment
echo "🔧 Setting up environment..."
if ! spack env list | grep -q "$LMGC90_ENV"; then
    spack env create "$LMGC90_ENV"
fi

# Activate environment
spack env activate "$LMGC90_ENV"

# Add packages
echo "📦 Adding packages..."
spack add gcc@14.2.0
spack add hdf5+fortran
spack add swig@3.0.12
spack add python@3.11
spack add cmake ninja

# Install everything
echo "🔨 Installing packages (this may take 30-60 minutes)..."
spack install --jobs=$JOBS

# Install Python packages
echo "🐍 Installing Python packages..."
spack env activate "$LMGC90_ENV"
python -m ensurepip --upgrade
python -m pip install numpy scipy matplotlib

# Create activation script
echo "📝 Creating activation script..."
cat > activate_lmgc90.sh << 'EOF'
#!/bin/bash
# LMGC90 Environment Activation
. ~/spack/share/spack/setup-env.sh
spack env activate lmgc90_env

# Set compilers
export CC=$(spack location -i gcc)/bin/gcc
export CXX=$(spack location -i gcc)/bin/g++
export FC=$(spack location -i gcc)/bin/gfortran
export F77=$(spack location -i gcc)/bin/gfortran

# Set HDF5
export HDF5_ROOT=$(spack location -i hdf5)
export HDF5_INCLUDE_DIRS="$HDF5_ROOT/include"
export HDF5_LIBRARIES="$HDF5_ROOT/lib"

# Set SWIG
export SWIG_EXECUTABLE=$(spack location -i swig)/bin/swig

# Add to PATH
export PATH=$(spack location -i gcc)/bin:$(spack location -i swig)/bin:$PATH

# Add to library path
export DYLD_LIBRARY_PATH="$HDF5_ROOT/lib:$DYLD_LIBRARY_PATH"

echo "✅ LMGC90 environment activated!"
echo "  GCC: $(gcc --version | head -1)"
echo "  Python: $(python --version)"
echo "  HDF5: $HDF5_ROOT"
echo "  SWIG: $(swig -version | head -1)"
EOF

chmod +x activate_lmgc90.sh

# Create CMake toolchain
echo "🔧 Creating CMake toolchain..."
cat > lmgc90_toolchain.cmake << 'EOF'
# LMGC90 CMake Toolchain
set(CMAKE_C_COMPILER "$(spack location -i gcc)/bin/gcc")
set(CMAKE_CXX_COMPILER "$(spack location -i gcc)/bin/g++")
set(CMAKE_Fortran_COMPILER "$(spack location -i gcc)/bin/gfortran")

set(HDF5_ROOT "$(spack location -i hdf5)")
set(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include")
set(HDF5_LIBRARIES "${HDF5_ROOT}/lib")

set(SWIG_EXECUTABLE "$(spack location -i swig)/bin/swig")
set(PYTHON_EXECUTABLE "$(spack location -i python)/bin/python")
EOF

echo ""
echo "🎉 Installation Complete!"
echo "========================"
echo ""
echo "To activate the environment:"
echo "  source activate_lmgc90.sh"
echo ""
echo "To build LMGC90:"
echo "  mkdir build && cd build"
echo "  cmake .. -DCMAKE_TOOLCHAIN_FILE=../lmgc90_toolchain.cmake"
echo "  make -j$JOBS"
echo ""
echo "📁 Files created:"
echo "  - activate_lmgc90.sh (environment activation)"
echo "  - lmgc90_toolchain.cmake (CMake toolchain)"
echo ""
echo "💾 Installation size: $(du -sh ~/spack | cut -f1)" 