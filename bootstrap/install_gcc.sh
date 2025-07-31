#!/bin/bash
# GCC Installation Script for macOS
set -e

DOWNLOADS_DIR="$(cd "$1" && pwd)"
INSTALL_PREFIX="$(cd "$2" && pwd)"
ARCH="$3"

echo "Installing GCC for macOS $ARCH..."
echo "Downloads: $DOWNLOADS_DIR"
echo "Install prefix: $INSTALL_PREFIX"

# Create directories
mkdir -p "$DOWNLOADS_DIR"
mkdir -p "$INSTALL_PREFIX"

# Determine download URL based on architecture
if [ "$ARCH" = "arm64" ]; then
    DMG_URL="https://github.com/fxcoudert/gfortran-for-macOS/releases/download/14.2-sonoma/gfortran-ARM-14.2-Sonoma.dmg"
    DMG_FILE="$DOWNLOADS_DIR/gfortran-ARM-14.2-Sonoma.dmg"
else
    DMG_URL="https://github.com/fxcoudert/gfortran-for-macOS/releases/download/14.2-sonoma/gfortran-Intel-14.2-Sonoma.dmg"
    DMG_FILE="$DOWNLOADS_DIR/gfortran-Intel-14.2-Sonoma.dmg"
fi

# Download DMG if not already present
if [ ! -f "$DMG_FILE" ]; then
    echo "Downloading GCC DMG..."
    curl -L -o "$DMG_FILE" "$DMG_URL"
else
    echo "DMG already downloaded: $DMG_FILE"
fi

# Mount DMG
echo "Mounting DMG..."
MOUNT_POINT="/tmp/gfortran_mount_$$"
hdiutil attach "$DMG_FILE" -mountpoint "$MOUNT_POINT" -nobrowse -quiet

# Extract PKG to temporary location
echo "Extracting PKG installer..."
TEMP_EXTRACT="/tmp/gfortran_extract_$$"
mkdir -p "$TEMP_EXTRACT"
xar -xf "$MOUNT_POINT/gfortran.pkg" -C "$TEMP_EXTRACT"

# Extract payload
echo "Extracting payload..."
cd "$TEMP_EXTRACT"
cat Payload | gunzip -dc | cpio -i

# Copy files to install prefix using absolute paths
echo "Installing GCC to $INSTALL_PREFIX/gcc..."
mkdir -p "$INSTALL_PREFIX/gcc"
if [ -d "$TEMP_EXTRACT/usr/local/gfortran" ]; then
    # Use absolute paths and preserve directory structure
    cp -R "$TEMP_EXTRACT/usr/local/gfortran/"* "$INSTALL_PREFIX/gcc/"
    echo "Debug: Files copied to $INSTALL_PREFIX/gcc"
    echo "Debug: Verifying copy was successful..."
    ls -la "$INSTALL_PREFIX/gcc/bin/" 2>/dev/null || echo "Warning: bin directory not found"
else
    echo "Error: usr/local/gfortran directory not found after extraction"
    echo "Debug: Contents of temp extract:"
    find "$TEMP_EXTRACT" -name "*gfortran*" -o -name "*gcc*" | head -5
    exit 1
fi

# Clean up
echo "Cleaning up..."
rm -rf "$TEMP_EXTRACT"

# Unmount DMG
echo "Unmounting DMG..."
hdiutil detach "$MOUNT_POINT" -quiet

# Verify installation
if [ -f "$INSTALL_PREFIX/gcc/bin/gcc" ] && [ -f "$INSTALL_PREFIX/gcc/bin/gfortran" ]; then
    echo "✅ GCC installation successful!"
    echo "GCC version: $($INSTALL_PREFIX/gcc/bin/gcc --version | head -n1)"
    echo "GFortran version: $($INSTALL_PREFIX/gcc/bin/gfortran --version | head -n1)"
else
    echo "❌ GCC installation failed!"
    echo "Debug: Checking what was actually installed..."
    ls -la "$INSTALL_PREFIX/gcc/" 2>/dev/null || echo "GCC directory not found"
    exit 1
fi
