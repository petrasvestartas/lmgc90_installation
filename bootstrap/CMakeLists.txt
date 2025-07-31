cmake_minimum_required(VERSION 3.20)
project(LMGC90_Bootstrap)

# Bootstrap configuration
set(BOOTSTRAP_PREFIX "${CMAKE_BINARY_DIR}/deps" CACHE PATH "Local dependency installation prefix")
set(BOOTSTRAP_DOWNLOADS "${CMAKE_BINARY_DIR}/downloads" CACHE PATH "Download cache directory")
set(BOOTSTRAP_BUILD_JOBS "4" CACHE STRING "Number of parallel build jobs")

# Detect system architecture and OS
if(APPLE)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(ARCH STREQUAL "arm64")
        set(MACOS_ARCH "ARM64")
        set(MACOS_ARCH_SUFFIX "arm64")
    else()
        set(MACOS_ARCH "Intel")
        set(MACOS_ARCH_SUFFIX "x86_64")
    endif()
    message(STATUS "Detected macOS ${MACOS_ARCH} (${ARCH})")
endif()

include(ExternalProject)
include(FetchContent)

# Create directories
file(MAKE_DIRECTORY ${BOOTSTRAP_PREFIX})
file(MAKE_DIRECTORY ${BOOTSTRAP_DOWNLOADS})
file(MAKE_DIRECTORY ${BOOTSTRAP_PREFIX}/bin)
file(MAKE_DIRECTORY ${BOOTSTRAP_PREFIX}/lib)
file(MAKE_DIRECTORY ${BOOTSTRAP_PREFIX}/include)

# =============================================================================
# 1. GCC/GFortran Compiler Toolchain
# =============================================================================

if(APPLE AND MACOS_ARCH STREQUAL "ARM64")
    # Install prebuilt GCC for macOS ARM64 using custom script
    set(GCC_VERSION "14.2.0")
    
    ExternalProject_Add(gcc_toolchain
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND 
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/install_gcc.sh ${BOOTSTRAP_DOWNLOADS} ${BOOTSTRAP_PREFIX} arm64
        LOG_INSTALL ON
    )
    
    set(BOOTSTRAP_CC "${BOOTSTRAP_PREFIX}/gcc/bin/gcc")
    set(BOOTSTRAP_CXX "${BOOTSTRAP_PREFIX}/gcc/bin/g++")
    set(BOOTSTRAP_FC "${BOOTSTRAP_PREFIX}/gcc/bin/gfortran")
    
elseif(APPLE AND MACOS_ARCH STREQUAL "Intel")
    # Install prebuilt GCC for macOS Intel using custom script
    
    ExternalProject_Add(gcc_toolchain
        DOWNLOAD_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND 
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/install_gcc.sh ${BOOTSTRAP_DOWNLOADS} ${BOOTSTRAP_PREFIX} x86_64
        LOG_INSTALL ON
    )
    
    set(BOOTSTRAP_CC "${BOOTSTRAP_PREFIX}/gcc/bin/gcc")
    set(BOOTSTRAP_CXX "${BOOTSTRAP_PREFIX}/gcc/bin/g++")
    set(BOOTSTRAP_FC "${BOOTSTRAP_PREFIX}/gcc/bin/gfortran")
    
else()
    # Build GCC from source for Linux or other systems
    ExternalProject_Add(gcc_toolchain
        URL "https://gcc.gnu.org/releases/gcc-14.2.0/gcc-14.2.0.tar.xz"
        DOWNLOAD_DIR ${BOOTSTRAP_DOWNLOADS}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure 
            --prefix=${BOOTSTRAP_PREFIX}/gcc
            --enable-languages=c,c++,fortran
            --disable-multilib
            --enable-shared
            --enable-threads=posix
        BUILD_COMMAND make -j${BOOTSTRAP_BUILD_JOBS}
        INSTALL_COMMAND make install
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )
    
    set(BOOTSTRAP_CC "${BOOTSTRAP_PREFIX}/gcc/bin/gcc")
    set(BOOTSTRAP_CXX "${BOOTSTRAP_PREFIX}/gcc/bin/g++")
    set(BOOTSTRAP_FC "${BOOTSTRAP_PREFIX}/gcc/bin/gfortran")
endif()

# =============================================================================
# 2. CMake (if system version is too old)
# =============================================================================

if(CMAKE_VERSION VERSION_LESS "3.20")
    ExternalProject_Add(cmake_bootstrap
        URL "https://github.com/Kitware/CMake/releases/download/v3.28.1/cmake-3.28.1.tar.gz"
        DOWNLOAD_DIR ${BOOTSTRAP_DOWNLOADS}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${BOOTSTRAP_PREFIX}/cmake
        BUILD_COMMAND make -j${BOOTSTRAP_BUILD_JOBS}
        INSTALL_COMMAND make install
        LOG_DOWNLOAD ON
        LOG_CONFIGURE ON
        LOG_BUILD ON
        LOG_INSTALL ON
    )
    set(BOOTSTRAP_CMAKE "${BOOTSTRAP_PREFIX}/cmake/bin/cmake")
else()
    set(BOOTSTRAP_CMAKE "${CMAKE_COMMAND}")
endif()

# =============================================================================
# 3. HDF5 with Fortran Support
# =============================================================================

ExternalProject_Add(hdf5
    DEPENDS gcc_toolchain
    URL "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_14_3.tar.gz"
    DOWNLOAD_DIR ${BOOTSTRAP_DOWNLOADS}
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${BOOTSTRAP_PREFIX}
        -DCMAKE_BUILD_TYPE=Release
        -DHDF5_BUILD_FORTRAN=ON
        -DHDF5_BUILD_CPP_LIB=ON
        -DHDF5_ENABLE_PARALLEL=OFF
        -DBUILD_SHARED_LIBS=ON
        -DHDF5_BUILD_TOOLS=ON
        -DCMAKE_C_COMPILER=${BOOTSTRAP_CC}
        -DCMAKE_CXX_COMPILER=${BOOTSTRAP_CXX}
        -DCMAKE_Fortran_COMPILER=${BOOTSTRAP_FC}
    BUILD_COMMAND ${BOOTSTRAP_CMAKE} --build . --parallel ${BOOTSTRAP_BUILD_JOBS}
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
    LOG_INSTALL ON
)

# =============================================================================
# 4. OpenBLAS (BLAS/LAPACK implementation)
# =============================================================================
# OpenBLAS - Using system package manager for reliable prebuilt binaries
ExternalProject_Add(openblas
    DEPENDS gcc_toolchain
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/install_openblas.sh ${BOOTSTRAP_PREFIX}
    LOG_INSTALL ON
)

# =============================================================================
# 5. SWIG 3.0.12 (Critical version for LMGC90)
# =============================================================================

# SWIG 3.0.12 - Using system package manager with version check
ExternalProject_Add(swig
    DEPENDS gcc_toolchain
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/install_swig.sh ${BOOTSTRAP_PREFIX}
    LOG_INSTALL ON
)

# =============================================================================
# 6. Python Dependencies (using local Python)
# =============================================================================

find_program(PYTHON_EXECUTABLE python3 REQUIRED)

ExternalProject_Add(python_deps
    DOWNLOAD_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND 
        ${PYTHON_EXECUTABLE} -m pip install --user --upgrade pip
        COMMAND ${PYTHON_EXECUTABLE} -m pip install --user numpy scipy matplotlib
        COMMAND ${PYTHON_EXECUTABLE} -m pip install --user vtk==9.2.6
        COMMAND ${PYTHON_EXECUTABLE} -m pip install --user compas compas_dem compas_viewer
    LOG_INSTALL ON
)

# =============================================================================
# 7. Environment Setup Scripts
# =============================================================================

# Create activation script for the bootstrap environment
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/activate_bootstrap.sh.in
    ${CMAKE_BINARY_DIR}/activate_bootstrap.sh
    @ONLY
)

# Create CMake toolchain file for LMGC90 build
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/lmgc90_toolchain.cmake.in
    ${CMAKE_BINARY_DIR}/lmgc90_toolchain.cmake
    @ONLY
)

# =============================================================================
# 8. Verification Target
# =============================================================================

add_custom_target(verify_bootstrap
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/verify_bootstrap.cmake
    DEPENDS gcc_toolchain hdf5 openblas swig python_deps
    COMMENT "Verifying bootstrap installation"
)

# =============================================================================
# 9. Summary and Instructions
# =============================================================================

add_custom_target(bootstrap_complete
    DEPENDS gcc_toolchain hdf5 openblas swig python_deps verify_bootstrap
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "🎉 LMGC90 Bootstrap Complete!"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "Next steps:"
    COMMAND ${CMAKE_COMMAND} -E echo "1. Source the environment: source ${CMAKE_BINARY_DIR}/activate_bootstrap.sh"
    COMMAND ${CMAKE_COMMAND} -E echo "2. Build LMGC90: cmake -DCMAKE_TOOLCHAIN_FILE=${CMAKE_BINARY_DIR}/lmgc90_toolchain.cmake .."
    COMMAND ${CMAKE_COMMAND} -E echo "3. Compile: make -j${BOOTSTRAP_BUILD_JOBS}"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "All dependencies installed in: ${BOOTSTRAP_PREFIX}"
    COMMENT "Bootstrap installation complete"
)
