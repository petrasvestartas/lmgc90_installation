#!/bin/bash
# LMGC90 Installation using Spack Package Manager
# This script provides a more reliable alternative to the bootstrap approach

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🚀 LMGC90 Installation using Spack${NC}"
echo -e "${BLUE}====================================${NC}"
echo ""

# Configuration
SPACK_ROOT="$HOME/spack"
LMGC90_ENV="lmgc90_env"
JOBS=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo "4")

echo -e "${YELLOW}System Information:${NC}"
echo "  OS: $(uname -s)"
echo "  Architecture: $(uname -m)"
echo "  CPU Cores: $JOBS"
echo "  Spack Root: $SPACK_ROOT"
echo ""

# Function to install Spack if not present
install_spack() {
    if [ ! -d "$SPACK_ROOT" ]; then
        echo -e "${YELLOW}📥 Installing Spack...${NC}"
        git clone -c feature.manyFiles=true https://github.com/spack/spack.git "$SPACK_ROOT"
        
        # Source Spack
        . "$SPACK_ROOT/share/spack/setup-env.sh"
        
        # Initialize Spack (this sets up the builtin repository properly)
        spack init --scope site
        
        echo -e "${GREEN}✅ Spack installed successfully${NC}"
    else
        echo -e "${GREEN}✅ Spack already installed${NC}"
    fi
}

# Function to setup Spack environment
setup_spack_env() {
    echo -e "${YELLOW}🔧 Setting up Spack environment...${NC}"
    
    # Source Spack
    . "$SPACK_ROOT/share/spack/setup-env.sh"
    
    # Initialize Spack if not already done
    if ! spack repo list | grep -q "builtin"; then
        echo -e "${BLUE}Initializing Spack repositories...${NC}"
        spack init --scope site
    fi
    
    # Create environment if it doesn't exist
    if ! spack env list | grep -q "$LMGC90_ENV"; then
        echo -e "${BLUE}Creating Spack environment: $LMGC90_ENV${NC}"
        spack env create "$LMGC90_ENV"
    fi
    
    # Activate environment
    spack env activate "$LMGC90_ENV"
    
    echo -e "${GREEN}✅ Spack environment ready${NC}"
}

# Function to install dependencies
install_dependencies() {
    echo -e "${YELLOW}📦 Installing LMGC90 Dependencies...${NC}"
    
    # Source Spack and activate environment
    . "$SPACK_ROOT/share/spack/setup-env.sh"
    spack env activate "$LMGC90_ENV"
    
    # Install core dependencies
    echo -e "${BLUE}Installing GCC (C/C++/Fortran compilers)...${NC}"
    spack add gcc@14.2.0
    
    echo -e "${BLUE}Installing HDF5 with Fortran support...${NC}"
    spack add hdf5@1.14.3 +fortran +hl +cxx
    
    echo -e "${BLUE}Installing OpenBLAS (BLAS/LAPACK)...${NC}"
    spack add openblas@0.3.26
    
    echo -e "${BLUE}Installing SWIG 3.0.12 (critical version)...${NC}"
    spack add swig@3.0.12
    
    echo -e "${BLUE}Installing Python and packages...${NC}"
    spack add python@3.11
    spack add py-numpy py-scipy py-matplotlib
    spack add vtk@9.2.6
    
    # Install all dependencies
    echo -e "${BLUE}Building all dependencies (this may take 30-60 minutes)...${NC}"
    spack install --jobs="$JOBS"
    
    echo -e "${GREEN}✅ All dependencies installed successfully${NC}"
}

# Function to create activation script
create_activation_script() {
    echo -e "${YELLOW}📝 Creating activation script...${NC}"
    
    # Source Spack and activate environment
    . "$SPACK_ROOT/share/spack/setup-env.sh"
    spack env activate "$LMGC90_ENV"
    
    # Get environment path
    ENV_PATH=$(spack env path "$LMGC90_ENV")
    
    # Create activation script
    cat > "activate_lmgc90_spack.sh" << EOF
#!/bin/bash
# LMGC90 Spack Environment Activation Script

# Source Spack
. "$SPACK_ROOT/share/spack/setup-env.sh"

# Activate LMGC90 environment
spack env activate "$LMGC90_ENV"

# Set environment variables for LMGC90
export LMGC90_SPACK_ENV="$ENV_PATH"
export CC=\$(spack location -i gcc)/bin/gcc
export CXX=\$(spack location -i gcc)/bin/g++
export FC=\$(spack location -i gcc)/bin/gfortran
export F77=\$(spack location -i gcc)/bin/gfortran

# Add Spack packages to PATH
export PATH=\$(spack location -i swig)/bin:\$PATH
export PATH=\$(spack location -i python)/bin:\$PATH

# Add libraries to LD_LIBRARY_PATH/DYLD_LIBRARY_PATH
if [[ "\$OSTYPE" == "darwin"* ]]; then
    export DYLD_LIBRARY_PATH=\$(spack location -i hdf5)/lib:\$(spack location -i openblas)/lib:\$DYLD_LIBRARY_PATH
else
    export LD_LIBRARY_PATH=\$(spack location -i hdf5)/lib:\$(spack location -i openblas)/lib:\$LD_LIBRARY_PATH
fi

# Add Python packages to PYTHONPATH
export PYTHONPATH=\$(spack location -i py-numpy)/lib/python*/site-packages:\$PYTHONPATH
export PYTHONPATH=\$(spack location -i py-scipy)/lib/python*/site-packages:\$PYTHONPATH
export PYTHONPATH=\$(spack location -i py-matplotlib)/lib/python*/site-packages:\$PYTHONPATH
export PYTHONPATH=\$(spack location -i vtk)/lib/python*/site-packages:\$PYTHONPATH

echo "✅ LMGC90 Spack environment activated"
echo "GCC: \$CC"
echo "Python: \$(which python)"
echo "SWIG: \$(which swig)"
EOF
    
    chmod +x "activate_lmgc90_spack.sh"
    echo -e "${GREEN}✅ Activation script created: activate_lmgc90_spack.sh${NC}"
}

# Function to create CMake toolchain
create_cmake_toolchain() {
    echo -e "${YELLOW}📝 Creating CMake toolchain file...${NC}"
    
    # Source Spack and activate environment
    . "$SPACK_ROOT/share/spack/setup-env.sh"
    spack env activate "$LMGC90_ENV"
    
    cat > "lmgc90_spack_toolchain.cmake" << EOF
# LMGC90 Spack Toolchain File
# Generated automatically by install_lmgc90_spack.sh

# Set compilers
set(CMAKE_C_COMPILER "$(spack location -i gcc)/bin/gcc")
set(CMAKE_CXX_COMPILER "$(spack location -i gcc)/bin/g++")
set(CMAKE_Fortran_COMPILER "$(spack location -i gcc)/bin/gfortran")

# Set HDF5
set(HDF5_ROOT "$(spack location -i hdf5)")
set(HDF5_INCLUDE_DIRS "\${HDF5_ROOT}/include")
set(HDF5_LIBRARIES "\${HDF5_ROOT}/lib")

# Set OpenBLAS
set(BLAS_LIBRARIES "$(spack location -i openblas)/lib")
set(LAPACK_LIBRARIES "$(spack location -i openblas)/lib")

# Set SWIG
set(SWIG_EXECUTABLE "$(spack location -i swig)/bin/swig")

# Set Python
set(PYTHON_EXECUTABLE "$(spack location -i python)/bin/python")
set(PYTHON_INCLUDE_DIRS "$(spack location -i python)/include/python*")
set(PYTHON_LIBRARIES "$(spack location -i python)/lib")

# Set VTK
set(VTK_DIR "$(spack location -i vtk)/lib/cmake/vtk")

# Add to CMAKE_PREFIX_PATH
list(APPEND CMAKE_PREFIX_PATH "$(spack location -i hdf5)")
list(APPEND CMAKE_PREFIX_PATH "$(spack location -i openblas)")

# Set find_package hints
set(HDF5_DIR "$(spack location -i hdf5)/lib/cmake/hdf5")
set(BLAS_DIR "$(spack location -i openblas)/lib/cmake/blas")
set(LAPACK_DIR "$(spack location -i openblas)/lib/cmake/lapack")
EOF
    
    echo -e "${GREEN}✅ CMake toolchain created: lmgc90_spack_toolchain.cmake${NC}"
}

# Function to verify installation
verify_installation() {
    echo -e "${YELLOW}🔍 Verifying installation...${NC}"
    
    # Source Spack and activate environment
    . "$SPACK_ROOT/share/spack/setup-env.sh"
    spack env activate "$LMGC90_ENV"
    
    echo -e "${BLUE}Checking GCC...${NC}"
    $(spack location -i gcc)/bin/gcc --version
    
    echo -e "${BLUE}Checking GFortran...${NC}"
    $(spack location -i gcc)/bin/gfortran --version
    
    echo -e "${BLUE}Checking SWIG version...${NC}"
    $(spack location -i swig)/bin/swig -version
    
    echo -e "${BLUE}Checking Python...${NC}"
    $(spack location -i python)/bin/python --version
    
    echo -e "${BLUE}Checking HDF5...${NC}"
    ls -la $(spack location -i hdf5)/lib/libhdf5*
    
    echo -e "${BLUE}Checking OpenBLAS...${NC}"
    ls -la $(spack location -i openblas)/lib/libopenblas*
    
    echo -e "${GREEN}✅ Installation verification complete${NC}"
}

# Function to show usage instructions
show_instructions() {
    echo ""
    echo -e "${GREEN}🎉 LMGC90 Spack Installation Complete!${NC}"
    echo ""
    echo -e "${YELLOW}📋 Next Steps:${NC}"
    echo "1. Activate the environment:"
    echo -e "   ${BLUE}source activate_lmgc90_spack.sh${NC}"
    echo ""
    echo "2. Build LMGC90:"
    echo -e "   ${BLUE}mkdir build_spack && cd build_spack${NC}"
    echo -e "   ${BLUE}cmake -DCMAKE_TOOLCHAIN_FILE=../lmgc90_spack_toolchain.cmake ..${NC}"
    echo -e "   ${BLUE}make -j$JOBS${NC}"
    echo ""
    echo "3. Test LMGC90:"
    echo -e "   ${BLUE}PYTHONPATH=\$PWD python3 -c \"import pylmgc90; print('✅ SUCCESS!')\"${NC}"
    echo ""
    echo -e "${YELLOW}💡 Useful Spack Commands:${NC}"
    echo "  • List installed packages: spack find"
    echo "  • Show environment: spack env status"
    echo "  • Deactivate environment: spack env deactivate"
    echo "  • Remove environment: spack env remove $LMGC90_ENV"
    echo ""
    echo -e "${GREEN}All dependencies are now managed by Spack in: $SPACK_ROOT${NC}"
}

# Main execution
main() {
    echo -e "${YELLOW}📦 Dependencies to be installed:${NC}"
    echo "  • GCC 14.2.0 (C/C++/Fortran compilers)"
    echo "  • HDF5 1.14.3 (with Fortran support)"
    echo "  • OpenBLAS 0.3.26 (BLAS/LAPACK)"
    echo "  • SWIG 3.0.12 (Python bindings)"
    echo "  • Python 3.11 with scientific packages"
    echo ""
    echo -e "${YELLOW}💾 Estimated disk space: ~4-5 GB${NC}"
    echo -e "${YELLOW}⏱️  Estimated time: 45-90 minutes (depending on system)${NC}"
    echo ""

    read -p "Continue with Spack installation? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo -e "${YELLOW}Installation cancelled.${NC}"
        exit 0
    fi

    install_spack
    setup_spack_env
    install_dependencies
    create_activation_script
    create_cmake_toolchain
    verify_installation
    show_instructions
}

# Run main function
main "$@" 