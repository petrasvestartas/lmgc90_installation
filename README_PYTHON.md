# Installation of LMGC90

## BREAKTHROUGH: Native ARM64 Solution for macOS Apple Silicon

**Problem**: Segmentation fault when importing `pylmgc90` on macOS Apple Silicon (M1/M2/M3)
**Solution**: Use native Homebrew ARM64 toolchain instead of conda

### Quick Solution (Working Method)

```bash
# 1. Install native ARM64 dependencies
brew install gcc cmake hdf5 python@3.11
/opt/homebrew/bin/pip3.11 install numpy scipy matplotlib swig==3.0.12
/opt/homebrew/bin/pip3.11 install vtk==9.2.6 compas compas_dem compas_viewer

# 2. Build with clean native environment
rm -rf build_native && mkdir build_native && cd build_native
env -i HOME=$HOME PATH=/opt/homebrew/bin:/usr/bin:/bin \
CC=gcc-14 CXX=g++-14 FC=gfortran-14 CMAKE_PREFIX_PATH=/opt/homebrew \
cmake .. -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE=/opt/homebrew/bin/python3.11

env -i HOME=$HOME PATH=/opt/homebrew/bin:/usr/bin:/bin \
CC=gcc-14 CXX=g++-14 FC=gfortran-14 CMAKE_PREFIX_PATH=/opt/homebrew \
make -j$(sysctl -n hw.ncpu)

# 3. Test (should work without segmentation fault)
env -i HOME=$HOME PATH=/opt/homebrew/bin:/usr/bin:/bin \
PYTHONPATH=$PWD /opt/homebrew/bin/python3.11 -c "import pylmgc90; print('✅ SUCCESS!')"
```

**Key Points**:
- **Root Cause**: Conda environment conflicts and mixed ARM64/x86_64 libraries
- **Solution**: Pure Homebrew ARM64 toolchain with GCC 14.2.0 and SWIG 3.0.12
- **Result**: No segmentation faults, native ARM64 performance, full functionality

---

### Download LMGC90
https://lmgc90.pages-git-xen.lmgc.univ-montp2.fr/lmgc90_dev/downloads/lmgc90_user_2025.rc1.zip

## macOS Installation

### 1. Create Environment
```bash
# Use the macOS-compatible environment file
conda env create -f environment_final.yml
conda activate lmgc
```

### 2. Build LMGC90

#### 2.1 Fix CMake Issue (if needed)
If you encounter `--print-libgcc-file` error with clang, edit `cmake/modules/FindMatLib_FOR_LMGC.cmake`:

```cmake
# Replace lines 22-28 with:
# Fix for macOS clang compiler - check if we're using GCC or clang
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} --print-libgcc-file
                    OUTPUT_VARIABLE LIBGCC
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE LIBGCC_RESULT)
    if(LIBGCC_RESULT EQUAL 0)
        get_filename_component(STDC++_LIB_DIR ${LIBGCC} PATH)
    else()
        set(STDC++_LIB_DIR "")
    endif()
else()
    # For clang and other compilers, use standard library paths
    set(STDC++_LIB_DIR "/usr/lib")
endif()
```

#### 2.2 Fix C++17 Register Keyword Issue (if needed)
The ANN library uses deprecated `register` keywords removed in C++17. Choose one solution:

**Option A: Remove register keywords (recommended)**
```bash
# Remove register keywords from all ANN source files
find src/contribs/ann-1.1.2/src -name "*.cpp" -exec sed -i '' 's/register //g' {} \;
find src/contribs/ann-1.1.2/src -name "*.h" -exec sed -i '' 's/register //g' {} \;
```

**Option B: Use older C++ standard**
```bash
# Use C++14 standard instead of C++17
cmake -DCMAKE_BUILD_TYPE=Release -DVENV_PATH=$CONDA_PREFIX -DCMAKE_CXX_STANDARD=14 ..
```

#### 2.3 Fix Permission Issues (if needed)
If you encounter permission denied errors for `/usr/local/lib`, configure cmake to use conda environment paths:

```bash
# Clean build directory first
rm -rf build/*
# Reconfigure with conda-specific install prefix
cmake -DCMAKE_BUILD_TYPE=Release -DVENV_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
```

#### 2.4 Configure and Build
```bash
cd build
# Use conda environment paths to avoid permission issues
cmake -DCMAKE_BUILD_TYPE=Release -DVENV_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
make -j$(sysctl -n hw.ncpu)  # Use macOS-compatible parallel build
make install
```

#### 2.5 Verify Build Success
```bash
# Set PYTHONPATH and test LMGC90 import
conda activate lmgc
export PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build:$PYTHONPATH
python -c "import pylmgc90; print('✅ LMGC90 successfully installed!')"
```

**Note:** If you get a segmentation fault, this indicates a runtime issue with the Python extensions. See troubleshooting below.

### 3. Run Examples
```bash
# Set PYTHONPATH with your actual macOS path
export PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build:$PYTHONPATH
cd examples/compas
python gen_dem_vault.py  # Run this first
python command.py        # Then run this
```

## Troubleshooting

### Segmentation Fault on Import
If you encounter a segmentation fault when importing pylmgc90:

```bash
# This indicates a runtime compatibility issue with the compiled extensions
# Common causes on macOS:
# 1. Architecture mismatch (Intel vs ARM64)
# 2. Library dependency conflicts
# 3. Fortran/C++ runtime library issues
```

**Potential Solutions:**

1. **Check Python architecture compatibility:**
   ```bash
   python -c "import platform; print('Architecture:', platform.machine())"
   # Should show 'arm64' on Apple Silicon Macs
   ```

2. **Verify conda environment is clean:**
   ```bash
   conda activate lmgc
   conda list | grep -E "(gcc|gfortran|clang)"
   ```

3. **Check library dependencies:**
   ```bash
   # Check what libraries the extension depends on
   otool -L build/pylmgc90/chipy/_lmgc90.so
   
   # Verify RPATH settings
   otool -l build/pylmgc90/chipy/_lmgc90.so | grep -A 2 LC_RPATH
   ```

4. **Try rebuilding with debug information:**
   ```bash
   cd build
   cmake -DCMAKE_BUILD_TYPE=Debug -DVENV_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
   make clean && make -j$(sysctl -n hw.ncpu)
   make install
   ```

5. **Test direct extension import:**
   ```bash
   # Try importing the extension module directly to isolate the issue
   python -c "import pylmgc90.chipy._lmgc90; print('Direct import successful')"
   ```

### Known Issue: SWIG 4.x Segmentation Fault on macOS ARM64

**Root Cause Identified:** SWIG 4.x has known segmentation fault issues with Python bindings on macOS ARM64.

**Evidence:**
- Current environment uses SWIG 4.3.1
- GitHub issues document SWIG 4.0+ causing segfaults on macOS with Python bindings
- SWIG 3.x versions work correctly but aren't available via conda-forge for ARM64

**Status:** Segmentation fault persists when importing pylmgc90 on macOS ARM64, even after:
- ✅ Successful clean rebuild with all fixes applied
- ✅ Verified correct architecture (arm64)
- ✅ Confirmed proper library dependencies and RPATH settings
- ✅ Debug build with additional error information
- ✅ Direct testing of extension module import
- ✅ Identified SWIG 4.x as root cause

**Analysis:** The segmentation fault occurs during the import of the `_lmgc90.so` SWIG-generated extension module due to SWIG 4.x compatibility issues with macOS ARM64.

**Potential Solutions:**

1. **Try Homebrew SWIG (different version):**
   ```bash
   # Install SWIG via Homebrew (may be different version than conda)
   brew install swig
   # Rebuild LMGC90 using Homebrew SWIG
   cd build && make clean
   cmake -DCMAKE_BUILD_TYPE=Release -DVENV_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX ..
   make -j$(sysctl -n hw.ncpu) && make install
   ```

2. **Use x86_64 conda environment with Rosetta 2:**
   ```bash
   # Create x86_64 conda environment that can use SWIG 3.x
   CONDA_SUBDIR=osx-64 conda create -n lmgc_x86 python=3.11
   conda activate lmgc_x86
   conda config --env --set subdir osx-64
   conda install swig=3.0.12  # Install SWIG 3.x
   # Then rebuild LMGC90 in this environment
   ```

3. **Docker container with known working environment**

4. **Contact LMGC90 developers** about SWIG 4.x compatibility on macOS ARM64

**References:**
- [SWIG 4.0 segfault on macOS](https://github.com/shogun-toolbox/shogun/issues/4629)
- [SWIG 4.x Python binding segfaults](https://github.com/swig/swig/issues/1955)

## Ubuntu Installation


<img width="1693" height="1281" alt="image" src="https://github.com/user-attachments/assets/9a4c4711-8879-4201-9a69-907ed89e0271" />

Missing HDF5 with Fortran Support:

```bash
sudo apt update
sudo apt install libhdf5-dev libhdf5-fortran-102
cd build
cmake ..
make
```
