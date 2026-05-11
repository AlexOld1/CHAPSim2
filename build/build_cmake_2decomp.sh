#!/bin/bash
set -euo pipefail

# =============================================================================
# 2decomp-fft Library Build Script with Architecture Auto-Detection
# - On ARCHER2/Cray: auto-load cray-fftw (if available) and use FFTW_F03
# - Otherwise: if FFTW_ROOT exists and looks valid, use FFTW_F03; else generic
# =============================================================================

echo "========================================================================="
echo "  2decomp-fft Library Build"
echo "========================================================================="
echo ""

# -----------------------------------------------------------------------------
# Detect Architecture and Platform
# -----------------------------------------------------------------------------
detect_architecture() {
    local arch=""
    local platform=""

    if [[ "$OSTYPE" == "darwin"* ]]; then
        platform="macOS"
        if [[ $(uname -m) == "arm64" ]]; then
            arch="arm64"
            echo "Detected: Apple Silicon (M1/M2/M3/M4)"
        elif [[ $(uname -m) == "x86_64" ]]; then
            arch="x86_64"
            echo "Detected: Intel Mac"
        else
            arch=$(uname -m)
            echo "Detected: macOS - $(uname -m)"
        fi
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        platform="Linux"
        arch=$(uname -m)
        echo "🐧 Detected: Linux - $arch"
    else
        platform="Unknown"
        arch=$(uname -m)
        echo "❓ Detected: $OSTYPE - $arch"
    fi

    echo "Platform: $platform"
    echo "Architecture: $arch"
    echo ""
}

# -----------------------------------------------------------------------------
# Set Architecture-Specific Compiler Flags
# -----------------------------------------------------------------------------
set_compiler_flags() {
    local arch="$1"
    local platform="$2"

    # default (may be overridden later on Cray)
    export FC=${FC:-mpif90}

    if [[ "$platform" == "macOS" ]]; then
        if [[ "$arch" == "arm64" ]]; then
            export FFLAGS="-arch arm64 -fallow-argument-mismatch"
            export FCFLAGS="-arch arm64 -fallow-argument-mismatch"
            export CFLAGS="-arch arm64"
            export CXXFLAGS="-arch arm64"
            echo "✅ Using Apple Silicon flags: -arch arm64"
        elif [[ "$arch" == "x86_64" ]]; then
            export FFLAGS="-arch x86_64 -fallow-argument-mismatch"
            export FCFLAGS="-arch x86_64 -fallow-argument-mismatch"
            export CFLAGS="-arch x86_64"
            export CXXFLAGS="-arch x86_64"
            echo "✅ Using Intel Mac flags: -arch x86_64"
        fi
    else
        export FFLAGS="-fallow-argument-mismatch"
        export FCFLAGS="-fallow-argument-mismatch"
        echo "✅ Using standard Unix flags"
    fi

    echo "FC (initial): $FC"
    echo "FFLAGS: ${FFLAGS:-none}"
    echo ""
}

# -----------------------------------------------------------------------------
# ARCHER2 / Cray environment setup:
# - detect Cray PE
# - load cray-fftw module if module system exists
# - use Cray wrapper compilers (ftn/cc/CC)
# - set FFTW_ROOT from CRAY_FFTW_PREFIX when available
# -----------------------------------------------------------------------------
setup_archer2_cray_env() {
    # module is often a shell function; source init if needed
    if ! command -v module >/dev/null 2>&1; then
        if [[ -f /etc/profile.d/modules.sh ]]; then
            # shellcheck disable=SC1091
            source /etc/profile.d/modules.sh
        fi
    fi

    local is_cray=false
    if [[ -n "${CRAYPE_VERSION:-}" ]] || [[ -n "${PE_ENV:-}" ]] || [[ -n "${CRAY_FFTW_PREFIX:-}" ]]; then
        is_cray=true
    fi

    if [[ "$is_cray" == true ]]; then
        echo "🛰️  Detected Cray PE / ARCHER2-like environment"
        echo "   CRAYPE_VERSION=${CRAYPE_VERSION:-unknown}  PE_ENV=${PE_ENV:-unknown}"

        # Prefer Cray wrappers unless user already forced something
        export FC=${FC:-ftn}
        export CC=${CC:-cc}
        export CXX=${CXX:-CC}

        # Load cray-fftw if possible and not already loaded
        if command -v module >/dev/null 2>&1; then
            if ! module -t list 2>&1 | grep -q '^cray-fftw/'; then
                echo "Loading module: cray-fftw"
                module load cray-fftw || echo "⚠️  Warning: module load cray-fftw failed"
            else
                echo "✅ cray-fftw already loaded"
            fi
        fi

        # Prefer CRAY_FFTW_PREFIX as FFTW_ROOT
        if [[ -n "${CRAY_FFTW_PREFIX:-}" ]] && [[ -d "${CRAY_FFTW_PREFIX}" ]]; then
            export FFTW_ROOT="${FFTW_ROOT:-$CRAY_FFTW_PREFIX}"
            echo "✅ FFTW_ROOT set from CRAY_FFTW_PREFIX: $FFTW_ROOT"
        fi

        echo "FC (Cray): $FC"
        echo "CC (Cray): $CC"
        echo ""
    fi
}

# -----------------------------------------------------------------------------
# Choose FFT backend:
# - if FFTW_ROOT exists and looks valid (include + lib), use FFTW_F03
# - else fallback to generic
# -----------------------------------------------------------------------------
choose_fft_backend() {
    # allow user override; else fall back to your mac path
    local fftw_root="${FFTW_ROOT:-/usr/local/}"

    FFT_CHOICE="generic"
    FFTW_CMAKE_ARGS=""

    local inc_ok=false
    local lib_ok=false

    if [[ -d "$fftw_root" ]]; then
        if [[ -d "$fftw_root/include" ]] && \
           ([[ -f "$fftw_root/include/fftw3.f03" ]] || [[ -f "$fftw_root/include/fftw3.h" ]]); then
            inc_ok=true
        fi
        if compgen -G "$fftw_root/lib/libfftw3*" >/dev/null 2>&1 || \
           compgen -G "$fftw_root/lib64/libfftw3*" >/dev/null 2>&1; then
            lib_ok=true
        fi
    fi

    if [[ "$inc_ok" == true && "$lib_ok" == true ]]; then
        FFT_CHOICE="FFTW_F03"
        FFTW_CMAKE_ARGS="-DFFT_Choice=FFTW_F03 -DFFTW_ROOT=$fftw_root"
        echo "✅ FFTW detected at: $fftw_root"
        echo "   -> Using FFT backend: FFTW_F03"
    else
        echo "ℹ️  FFTW not found / incomplete at: $fftw_root"
        echo "   include ok? $inc_ok   lib ok? $lib_ok"
        echo "   -> Using FFT backend: generic"
    fi
    echo ""
}

# -----------------------------------------------------------------------------
# Clean Previous Build
# -----------------------------------------------------------------------------
clean_build() {
    echo "Cleaning previous build artifacts..."
    rm -rf CMakeCache.txt CMakeFiles/ cmake_install.cmake Makefile
    rm -rf opt/ lib/ lib64/ include/ bin/
    rm -rf src/CMakeFiles/ examples/CMakeFiles/
    echo "✅ Clean complete"
    echo ""
}

# -----------------------------------------------------------------------------
# Configure with CMake
# -----------------------------------------------------------------------------
configure_cmake() {
    echo "Configuring with CMake..."

    local install_prefix="$(pwd)/opt"

    cmake -S ../ -B ./ \
        -DCMAKE_INSTALL_PREFIX="$install_prefix" \
        -DCMAKE_Fortran_COMPILER="$FC" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_Fortran_FLAGS="${FFLAGS:-}" \
        -DCMAKE_C_FLAGS="${CFLAGS:-}" \
        -DCMAKE_CXX_FLAGS="${CXXFLAGS:-}" \
        ${FFTW_CMAKE_ARGS} \
        || { echo "❌ CMake configuration failed"; return 1; }

    echo "✅ CMake configuration complete"
    echo ""
}

# -----------------------------------------------------------------------------
# Build the Library
# -----------------------------------------------------------------------------
build_library() {
    echo "Building library..."
    local log_file="build_2decomp.log"
    if cmake --build ./ --parallel >"$log_file" 2>&1; then
        echo "✅ Build complete"
    else
        echo "❌ Build failed"
        echo "Last 30 lines from $log_file:"
        tail -30 "$log_file"
        return 1
    fi
    echo "Build log saved to: $(pwd)/$log_file"
    echo ""
}

# -----------------------------------------------------------------------------
# Install the Library
# -----------------------------------------------------------------------------
install_library() {
    echo "Installing library..."
    cmake --install ./ || { echo "❌ Installation failed"; return 1; }
    echo "✅ Installation complete"
    echo ""
}

# -----------------------------------------------------------------------------
# Verify Installation
# -----------------------------------------------------------------------------
verify_installation() {
    echo "Verifying installation..."

    local lib_path=""

    if [ -f "./opt/lib64/libdecomp2d.a" ]; then
        lib_path="./opt/lib64/libdecomp2d.a"
    elif [ -f "./opt/lib/libdecomp2d.a" ]; then
        lib_path="./opt/lib/libdecomp2d.a"
    else
        echo "❌ ERROR: libdecomp2d.a not found in opt/lib or opt/lib64"
        return 1
    fi

    echo "📍 Library location: $lib_path"

    echo ""
    echo "File information:"
    file "$lib_path"

    echo ""
    echo "Archive contents (first 10 entries):"
    if ar -t "$lib_path" > /dev/null 2>&1; then
        ar -t "$lib_path" | head -10
        local obj_count
        obj_count=$(ar -t "$lib_path" | wc -l | tr -d ' ')
        echo "..."
        echo "Total object files: $obj_count"
    else
        echo "❌ ERROR: Cannot read archive contents"
        return 1
    fi

    echo ""
    echo "Checking for suspicious entries..."
    if ar -t "$lib_path" | grep -E "^/$|^//" > /dev/null 2>&1; then
        echo "⚠️  WARNING: Found suspicious '/' entries in archive!"
        ar -t "$lib_path" | grep -E "^/$|^//"
        return 1
    else
        echo "✅ No suspicious entries found"
    fi

    if command -v lipo &> /dev/null; then
        echo ""
        echo "Architecture information:"
        lipo -info "$lib_path" 2>&1 || echo "Note: lipo info not available for static libraries"
    fi

    echo ""
    echo "Rebuilding archive index with ranlib..."
    ranlib "$lib_path" || echo "⚠️  Warning: ranlib failed"

    echo ""
    echo "✅ Verification complete"
    echo ""
}

# =============================================================================
# Main Execution
# =============================================================================

ARCH=$(uname -m)
PLATFORM=""
if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="macOS"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PLATFORM="Linux"
else
    PLATFORM="Unknown"
fi

detect_architecture
set_compiler_flags "$ARCH" "$PLATFORM"

# ARCHER2/Cray optional setup (loads cray-fftw and sets FFTW_ROOT automatically)
setup_archer2_cray_env

# Decide FFT backend (FFTW_F03 if FFTW_ROOT is valid; else generic)
choose_fft_backend

clean_build
configure_cmake || exit 1
build_library || exit 1
install_library || exit 1
verify_installation || exit 1

echo "========================================================================="
echo "✅ 2decomp-fft library build completed successfully!"
echo "========================================================================="
echo ""
echo "Library installed to: $(pwd)/opt"
echo "FFT backend selected: ${FFT_CHOICE}"
echo ""

if [ -f "./opt/lib64/libdecomp2d.a" ]; then
    echo "Use this path in your Makefile: ../lib/2decomp-fft/build/opt/lib64"
elif [ -f "./opt/lib/libdecomp2d.a" ]; then
    echo "Use this path in your Makefile: ../lib/2decomp-fft/build/opt/lib"
fi
echo ""
