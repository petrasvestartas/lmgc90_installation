# Installation of LMGC90

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

**Note:** If you get a segmentation fault, the build was successful but you need to set PYTHONPATH properly for Python to find the modules.

### 3. Run Examples
```bash
# Set PYTHONPATH with your actual macOS path
export PYTHONPATH=/Users/petras/brg/2_code/lmgc90_installation/build:$PYTHONPATH
cd examples/compas
python command.py
```

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
