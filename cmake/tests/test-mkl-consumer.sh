#!/usr/bin/env bash
# Verify that an installed MKL-enabled lion can be used by a C/C++ project
# without enabling Fortran, and that its BLAS implementation links and runs.
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <lion install prefix> [oneMKL root]" >&2
    exit 2
fi

lion_prefix="$(cd "$1" && pwd)"
workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

cat > "$workdir/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.13)
project(lion_mkl_consumer LANGUAGES C CXX)

find_package(lion CONFIG REQUIRED PATHS "${LION_INSTALL_PREFIX}" NO_DEFAULT_PATH)
add_executable(lion_mkl_smoke main.cpp)
target_link_libraries(lion_mkl_smoke PRIVATE lion::blaslapack)
CMAKE

cat > "$workdir/main.cpp" <<'CPP'
#include <mkl.h>

int main()
{
    const double a[1] = {2.0};
    const double b[1] = {3.0};
    double c[1] = {0.0};
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                1, 1, 1, 1.0, a, 1, b, 1, 0.0, c, 1);
    return c[0] == 6.0 ? 0 : 1;
}
CPP

cmake_args=("-DLION_INSTALL_PREFIX=$lion_prefix")
if [ "$#" -eq 2 ]; then
    cmake_args+=("-DLION_MKL_ROOT=$2")
fi

cmake -S "$workdir" -B "$workdir/build" "${cmake_args[@]}"
cmake --build "$workdir/build" --parallel
"$workdir/build/lion_mkl_smoke"
echo "PASS: a C/C++-only consumer of lion linked and used oneMKL BLAS."
