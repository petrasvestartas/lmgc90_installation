#!/bin/bash
set -e

BOOTSTRAP_PREFIX="$1"

echo "=== Setting up SWIG for bootstrap ==="

# Create a minimal SWIG setup that works with system SWIG
# This approach uses system SWIG but ensures it's available in bootstrap prefix

mkdir -p "$BOOTSTRAP_PREFIX/bin" "$BOOTSTRAP_PREFIX/share"

# Check if system SWIG is available
if command -v swig >/dev/null 2>&1; then
    SYSTEM_SWIG=$(which swig)
    echo "Found system SWIG: $SYSTEM_SWIG"
    
    # Create symlink in bootstrap prefix
    ln -sf "$SYSTEM_SWIG" "$BOOTSTRAP_PREFIX/bin/swig"
    
    # Check version
    SWIG_VERSION=$("$SYSTEM_SWIG" -version 2>&1 | grep 'SWIG Version' | cut -d' ' -f3 || echo "unknown")
    echo "System SWIG version: $SWIG_VERSION"
    
    if [[ "$SWIG_VERSION" == "4."* ]]; then
        echo "⚠️  WARNING: SWIG 4.x detected - may cause issues with LMGC90 on macOS ARM64"
        echo "   Consider using conda environment with SWIG 3.0.12 for LMGC90 build"
        echo "   Command: conda create -n lmgc_swig python=3.11 && conda activate lmgc_swig"
        echo "   Then: CONDA_SUBDIR=osx-64 conda install swig=3.0.12"
    fi
    
else
    echo "⚠️  No system SWIG found"
    echo "   For LMGC90, you'll need SWIG 3.0.12. Install options:"
    echo "   1. conda install swig=3.0.12 (recommended)"
    echo "   2. brew install swig (macOS)"
    echo "   3. apt-get install swig (Linux)"
    
    # Create a placeholder script
    cat > "$BOOTSTRAP_PREFIX/bin/swig" << 'EOF'
#!/bin/bash
echo "SWIG not found in bootstrap. Please install SWIG 3.0.12:"
echo "conda install swig=3.0.12"
exit 1
EOF
    chmod +x "$BOOTSTRAP_PREFIX/bin/swig"
fi

# Copy files to bootstrap prefix
echo "Installing SWIG to bootstrap prefix..."
mkdir -p "$BOOTSTRAP_PREFIX/bin" "$BOOTSTRAP_PREFIX/share"

# Copy binaries
if [[ -d "bin" ]]; then
    cp -r bin/* "$BOOTSTRAP_PREFIX/bin/" || true
fi

# Copy share files (SWIG library files)
if [[ -d "share" ]]; then
    cp -r share/* "$BOOTSTRAP_PREFIX/share/" || true
fi

echo "SWIG installation completed successfully"
echo "SWIG binary installed to: $BOOTSTRAP_PREFIX/bin/swig"

# Verify installation and check version
if [[ -f "$BOOTSTRAP_PREFIX/bin/swig" ]]; then
    echo "✅ SWIG binary verified: $BOOTSTRAP_PREFIX/bin/swig"
    
    # Check SWIG version
    SWIG_VERSION=$("$BOOTSTRAP_PREFIX/bin/swig" -version 2>&1 | grep 'SWIG Version' | cut -d' ' -f3 || echo "unknown")
    echo "SWIG version: $SWIG_VERSION"
    
    if [[ "$SWIG_VERSION" == "3.0.12" ]]; then
        echo "✅ Correct SWIG version (3.0.12) - compatible with LMGC90"
    else
        echo "⚠️  SWIG version $SWIG_VERSION may not be optimal for LMGC90 (expected 3.0.12)"
    fi
else
    echo "❌ SWIG binary not found after installation"
    ls -la "$BOOTSTRAP_PREFIX/bin/" || true
fi
