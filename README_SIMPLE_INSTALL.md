# LMGC90 Simple Installation

## 🚀 One-Command Installation

This script installs all LMGC90 dependencies using Spack package manager.

### Quick Start

```bash
# Download and run the installation script
curl -O https://raw.githubusercontent.com/lmgc90/lmgc90/main/install_lmgc90_simple.sh
chmod +x install_lmgc90_simple.sh
./install_lmgc90_simple.sh
```

### What It Does

1. **Installs Spack** (if not present)
2. **Creates environment** with all dependencies
3. **Installs packages**:
   - GCC 14.2.0 (C/C++/Fortran compilers)
   - HDF5 1.14.3 (with Fortran support)
   - SWIG 3.0.12 (Python bindings)
   - Python 3.11.11
   - numpy, scipy, matplotlib
4. **Creates activation script** and CMake toolchain

### Usage After Installation

```bash
# Activate environment
source activate_lmgc90.sh

# Build LMGC90
mkdir build && cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=../lmgc90_toolchain.cmake
make -j$(nproc)
```

### System Requirements

- macOS ARM64 (Apple Silicon) or x86_64
- Git
- Internet connection
- ~4GB disk space

### Installation Time

- **First run**: 30-60 minutes (compiles from source)
- **Subsequent runs**: Instant (uses cached packages)

### Files Created

- `~/spack/` - All installed packages (3.9GB)
- `activate_lmgc90.sh` - Environment activation script
- `lmgc90_toolchain.cmake` - CMake toolchain file

### Troubleshooting

If installation fails:
1. Check internet connection
2. Ensure you have sufficient disk space
3. Try running with fewer jobs: `JOBS=4 ./install_lmgc90_simple.sh`

### Replicability

This installation is **completely local** and can be replicated on any macOS system with the same script. 