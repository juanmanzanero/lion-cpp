#!/usr/bin/env bash
# Offline smoke test for the CMake oneMKL bootstrap; does NOT download Intel MKL.
# Builds tiny fake libraries that export pardiso, plus a fake offline installer.
set -euo pipefail

module_dir="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

mkdir -p "$work/payload/mkl/2026.2/include" \
         "$work/payload/mkl/2026.2/lib/cmake/mkl" \
         "$work/src"
printf '%s\n' 'extern "C" void pardiso(void);' > "$work/payload/mkl/2026.2/include/mkl.h"

printf '%s\n' 'void pardiso(void) {}' > "$work/pardiso.c"
printf '%s\n' 'void dummy(void) {}' > "$work/dummy.c"
gcc -fPIC -shared "$work/pardiso.c" -o "$work/payload/mkl/2026.2/lib/libmkl_core.so"
gcc -fPIC -shared "$work/dummy.c" -o "$work/payload/mkl/2026.2/lib/libmkl_gf_lp64.so"
gcc -fPIC -shared "$work/dummy.c" -o "$work/payload/mkl/2026.2/lib/libmkl_sequential.so"
gcc -fPIC -shared "$work/dummy.c" -o "$work/payload/mkl/2026.2/lib/libmkl_gnu_thread.so"
gcc -c "$work/pardiso.c" -o "$work/pardiso.o"
gcc -c "$work/dummy.c" -o "$work/dummy.o"
ar rcs "$work/payload/mkl/2026.2/lib/libmkl_core.a" "$work/pardiso.o"
ar rcs "$work/payload/mkl/2026.2/lib/libmkl_gf_lp64.a" "$work/dummy.o"
ar rcs "$work/payload/mkl/2026.2/lib/libmkl_sequential.a" "$work/dummy.o"
ar rcs "$work/payload/mkl/2026.2/lib/libmkl_gnu_thread.a" "$work/dummy.o"

cat > "$work/payload/mkl/2026.2/lib/cmake/mkl/MKLConfig.cmake" <<'CMAKE'
get_filename_component(_fake_mkl_root "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
if(NOT TARGET MKL::MKL)
    add_library(MKL::MKL INTERFACE IMPORTED)
    if(MKL_LINK STREQUAL "static")
        set(_fake_suffix "a")
    else()
        set(_fake_suffix "so")
    endif()
    set_target_properties(MKL::MKL PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_fake_mkl_root}/include"
        INTERFACE_LINK_LIBRARIES "${_fake_mkl_root}/lib/libmkl_gf_lp64.${_fake_suffix};${_fake_mkl_root}/lib/libmkl_${MKL_THREADING}.${_fake_suffix};${_fake_mkl_root}/lib/libmkl_core.${_fake_suffix}")
endif()
set(MKL_ROOT "${_fake_mkl_root}")
set(MKL_FOUND TRUE)
CMAKE

cat > "$work/installer.sh" <<'SH'
#!/bin/sh
set -eu
prefix=""
for arg in "$@"; do
    case "$arg" in
        --install-dir=*) prefix="${arg#--install-dir=}" ;;
    esac
done
[ -n "$prefix" ] || exit 2
mkdir -p "$prefix"
printf '%s\n' "called" >> "$prefix/install-count.txt"
cp -R "$FAKE_MKL_PAYLOAD/mkl" "$prefix/"
SH
chmod +x "$work/installer.sh"

cat > "$work/src/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.10)
project(mkl_bootstrap_smoke LANGUAGES C CXX Fortran)
include("$ENV{LION_TEST_MKL_FILE}")
if(NOT blaslapack_FOUND)
    message(FATAL_ERROR "MKL bootstrap did not create the BLAS/LAPACK target")
endif()
CMAKE

# Download disabled: actionable failure rather than a silent fallback.
if env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/missing" \
    -DMKL_THREADING=sequential -DLION_MKL_AUTO_INSTALL=OFF > "$work/missing.log" 2>&1; then
    echo 'FAIL: missing MKL should be rejected with auto-install disabled' >&2
    exit 1
fi
grep -q 'LION_MKL_AUTO_INSTALL=OFF' "$work/missing.log"

# First-time provisioning must not accept the Intel EULA implicitly.
if env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/no-consent" \
    -DMKL_THREADING=sequential -DLION_MKL_INSTALLER="$work/installer.sh" > "$work/no-consent.log" 2>&1; then
    echo 'FAIL: missing explicit EULA acceptance' >&2
    exit 1
fi
grep -q 'LION_MKL_ACCEPT_EULA=ON' "$work/no-consent.log"

# Provision into a build-local directory; both PARDISO link probes must work.
env FAKE_MKL_PAYLOAD="$work/payload" LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/build" \
    -DMKL_THREADING=sequential \
    -DLION_MKL_INSTALLER="$work/installer.sh" -DLION_MKL_ACCEPT_EULA=ON \
    -DLION_MKL_INSTALL_PREFIX="$work/local-oneapi" > "$work/install.log" 2>&1 || {
    cat "$work/install.log" >&2
    exit 1
}
grep -q 'LION_MKL_TARGET_HAS_PARDISO:INTERNAL=1' "$work/build/CMakeCache.txt"
grep -q 'LION_MKL_AUTOTOOLS_HAS_PARDISO:INTERNAL=1' "$work/build/CMakeCache.txt"
[ "$(wc -l < "$work/local-oneapi/install-count.txt")" -eq 1 ]

# Reconfigure and fresh-build reuse must not invoke the installer again.
env FAKE_MKL_PAYLOAD="$work/payload" LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/build" > "$work/reconfigure.log" 2>&1
env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/reuse" \
    -DLION_MKL_INSTALL_PREFIX="$work/local-oneapi" \
    -DLION_MKL_AUTO_INSTALL=OFF -DMKL_THREADING=sequential > "$work/reuse.log" 2>&1
[ "$(wc -l < "$work/local-oneapi/install-count.txt")" -eq 1 ]

# Explicit external root, no installer and no network.
env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/explicit" \
    -DLION_MKL_ROOT="$work/local-oneapi/mkl/2026.2" \
    -DLION_MKL_AUTO_INSTALL=OFF -DMKL_THREADING=sequential > "$work/explicit.log" 2>&1

# Standard CMake discovery via CMAKE_PREFIX_PATH is also supported.
env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/cmake-prefix" \
    -DCMAKE_PREFIX_PATH="$work/local-oneapi/mkl/2026.2" \
    -DLION_MKL_INSTALL_PREFIX="$work/empty-oneapi" \
    -DLION_MKL_AUTO_INSTALL=OFF -DMKL_THREADING=sequential > "$work/cmake-prefix.log" 2>&1 || {
    cat "$work/cmake-prefix.log" >&2
    exit 1
}

# Exercise the GNU OpenMP and static-library link configurations too.
env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/gnu-thread" \
    -DLION_MKL_ROOT="$work/local-oneapi/mkl/2026.2" \
    -DLION_MKL_AUTO_INSTALL=OFF -DMKL_THREADING=gnu_thread > "$work/gnu-thread.log" 2>&1 || {
    cat "$work/gnu-thread.log" >&2
    exit 1
}
env LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/static" \
    -DLION_MKL_ROOT="$work/local-oneapi/mkl/2026.2" \
    -DLION_MKL_AUTO_INSTALL=OFF -DMKL_THREADING=sequential -DMKL_LINK=static > "$work/static.log" 2>&1 || {
    cat "$work/static.log" >&2
    exit 1
}

# Reject a wrong digest before executing any installer.
if env FAKE_MKL_PAYLOAD="$work/payload" LION_TEST_MKL_FILE="$module_dir/mkl.cmake" \
    cmake -S "$work/src" -B "$work/bad-hash" \
    -DMKL_THREADING=sequential \
    -DLION_MKL_INSTALLER="$work/installer.sh" -DLION_MKL_ACCEPT_EULA=ON \
    -DLION_MKL_INSTALLER_SHA384="$(printf '%096d' 0)" \
    -DLION_MKL_INSTALL_PREFIX="$work/bad-oneapi" > "$work/bad-hash.log" 2>&1; then
    echo 'FAIL: wrong installer digest accepted' >&2
    exit 1
fi
grep -q 'SHA384 mismatch' "$work/bad-hash.log"
[ ! -e "$work/bad-oneapi/install-count.txt" ]

echo 'PASS: MKL bootstrap provision/reuse/consent/root/integrity/dynamic/static/threading smoke tests'
