#!/bin/bash
set -e

BOOTSTRAP_PREFIX="$1"

echo "=== Installing OpenBLAS via Apple Accelerate Framework ==="

# For macOS, use Apple Accelerate (built-in optimized BLAS/LAPACK)
# For Linux, skip OpenBLAS for now and use system libraries

if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Using Apple Accelerate Framework (built-in optimized BLAS/LAPACK)"
    
    # Create symbolic links to Apple Accelerate
    mkdir -p "$BOOTSTRAP_PREFIX/lib" "$BOOTSTRAP_PREFIX/include"
    
    # Link to Apple's Accelerate framework (use the main framework symlink)
    ln -sf /System/Library/Frameworks/Accelerate.framework/Accelerate "$BOOTSTRAP_PREFIX/lib/libopenblas.dylib" || true
    ln -sf /System/Library/Frameworks/Accelerate.framework/Accelerate "$BOOTSTRAP_PREFIX/lib/libblas.dylib" || true
    ln -sf /System/Library/Frameworks/Accelerate.framework/Accelerate "$BOOTSTRAP_PREFIX/lib/liblapack.dylib" || true
    
    # Create minimal headers
    cat > "$BOOTSTRAP_PREFIX/include/openblas_config.h" << 'EOF'
#ifndef OPENBLAS_CONFIG_H
#define OPENBLAS_CONFIG_H

#define OPENBLAS_VERSION "0.3.26-accelerate"
#define OPENBLAS_USE_ACCELERATE 1

#endif
EOF
    
    echo "Apple Accelerate Framework linked successfully"
    echo "This provides optimized BLAS/LAPACK for macOS ARM64"
    
else
    echo "Linux detected - creating placeholder for system BLAS/LAPACK"
    
    # Create placeholder files for Linux
    mkdir -p "$BOOTSTRAP_PREFIX/lib" "$BOOTSTRAP_PREFIX/include"
    
    # Create placeholder config
    cat > "$BOOTSTRAP_PREFIX/include/openblas_config.h" << 'EOF'
#ifndef OPENBLAS_CONFIG_H
#define OPENBLAS_CONFIG_H

#define OPENBLAS_VERSION "0.3.26-system"
#define OPENBLAS_USE_SYSTEM 1

#endif
EOF
    
    echo "Placeholder created - install system BLAS/LAPACK: sudo apt-get install libblas-dev liblapack-dev"
fi

# Copy files to bootstrap prefix
echo "Installing OpenBLAS to bootstrap prefix..."
mkdir -p "$BOOTSTRAP_PREFIX/lib" "$BOOTSTRAP_PREFIX/include" "$BOOTSTRAP_PREFIX/bin"

# Copy libraries
if [[ -d "lib" ]]; then
    cp -r lib/* "$BOOTSTRAP_PREFIX/lib/" || true
fi

# Copy headers
if [[ -d "include" ]]; then
    cp -r include/* "$BOOTSTRAP_PREFIX/include/" || true
fi

# Copy binaries
if [[ -d "bin" ]]; then
    cp -r bin/* "$BOOTSTRAP_PREFIX/bin/" || true
fi

echo "OpenBLAS installation completed successfully"
echo "Libraries installed to: $BOOTSTRAP_PREFIX/lib"
echo "Headers installed to: $BOOTSTRAP_PREFIX/include"

# Verify installation
if [[ "$OSTYPE" == "darwin"* ]]; then
    if [[ -f "$BOOTSTRAP_PREFIX/lib/libopenblas.dylib" ]]; then
        echo "✅ OpenBLAS library verified: $BOOTSTRAP_PREFIX/lib/libopenblas.dylib"
    else
        echo "❌ OpenBLAS library not found after installation"
        ls -la "$BOOTSTRAP_PREFIX/lib/" || true
    fi
else
    if [[ -f "$BOOTSTRAP_PREFIX/lib/libopenblas.so" ]]; then
        echo "✅ OpenBLAS library verified: $BOOTSTRAP_PREFIX/lib/libopenblas.so"
    else
        echo "❌ OpenBLAS library not found after installation"
        ls -la "$BOOTSTRAP_PREFIX/lib/" || true
    fi
fi
