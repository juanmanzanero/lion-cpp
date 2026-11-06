# lion: lightweigh interfaces for optimization and numerics
[![Latest commit](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake.yml/badge.svg)](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake.yml)
[![MacOS (latest)](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake_macos.yml/badge.svg)](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake_macos.yml)
[![Windows](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake_windows.yml/badge.svg)](https://github.com/juanmanzanero/lion-cpp/actions/workflows/cmake_windows.yml)
[![codecov](https://codecov.io/gh/juanmanzanero/lion-cpp/branch/main/graph/badge.svg?token=YNPC90CWH0)](https://codecov.io/gh/juanmanzanero/lion-cpp)


## Optional Intel oneMKL + Ipopt PARDISO on Linux

The default build continues to use reference BLAS/LAPACK and MUMPS. On Linux
with GCC/gfortran, `WITH_MKL=ON` uses oneMKL for **all** BLAS/LAPACK calls
(including those from MUMPS, HSL if enabled, and Ipopt) and enables Ipopt's
`pardisomkl` linear solver. No reference BLAS/LAPACK libraries are built or
linked by lion in that configuration.

Install the Intel oneAPI MKL **development** package first, then configure
from a **new build directory**:

```bash
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-mkl -DCMAKE_BUILD_TYPE=Release -DWITH_MKL=ON
cmake --build build-mkl --parallel
cmake --install build-mkl --prefix /your/install/prefix
```

If CMake cannot locate MKL, pass
`-DMKL_DIR=/opt/intel/oneapi/mkl/latest/lib/cmake/mkl`. When consuming an
installed lion package, MKL must remain installed and discoverable through
`MKL_DIR` / `CMAKE_PREFIX_PATH` (and, for shared MKL, available to the runtime
linker, e.g. after sourcing `setvars.sh`). Consumers should enable the
Fortran language before `find_package(lion)` (for example,
`project(my_app LANGUAGES CXX Fortran)`) so Intel's `MKLConfig.cmake` selects
the same GNU Fortran interface.

Supported MKL settings for this integration are:

- `MKL_INTERFACE_FULL=gf_lp64` (required: GNU Fortran, 32-bit BLAS integers).
- `MKL_THREADING=gnu_thread` (default) or `sequential`.
- `MKL_LINK=dynamic` (default) or `static`.

For example, use `-DMKL_THREADING=sequential` to avoid nested MKL OpenMP
threads during parallel CI runs. The build checks that both the CMake MKL
target and the Ipopt/Autotools linker flags resolve PARDISO. A different
BLAS/LAPACK selection or MKL mode requires a separate build directory to
prevent stale Ipopt/MUMPS binaries from being reused.

To use PARDISO from the Ipopt C++ API:

```cpp
app->Options()->SetStringValue("linear_solver", "pardisomkl");
```

MUMPS remains available via `linear_solver=mumps`. Verify a completed Ipopt
solve using PARDISO (not merely that `SetStringValue` succeeds) and compare
against MUMPS. For a dynamically linked executable, inspect dependencies with
`ldd` and make sure no unrelated reference BLAS/LAPACK shared libraries are
loaded. Other application dependencies may still load their own BLAS provider;
keep the final application link graph consistent.
