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

CMake looks for an existing oneMKL development installation first. If missing,
it can download Intel's pinned **2026.2.0.890** standalone binary installer and
install it **without sudo** under `build-mkl/lion/thirdparty/oneapi`. This is a
configure-time step because lion builds Ipopt/MUMPS during CMake configuration.
**At `cmake --install` time**, lion also copies the complete oneMKL development
installation to `<prefix>/oneapi/mkl/<version>` (by default). Installed
`lionConfig.cmake` resolves that copy *relative to its own location* rather
than referring to the original build tree. The installed HSL and MUMPS shared
libraries also receive an `$ORIGIN`-relative MKL RUNPATH. This makes the lion
installation usable by Jenkins and after the original build directory is removed.

First-time local install (review and explicitly accept Intel's license):

```bash
cmake -S . -B build-mkl -DCMAKE_BUILD_TYPE=Release \
  -DWITH_MKL=ON -DLION_MKL_ACCEPT_EULA=ON
cmake --build build-mkl --parallel
sudo cmake --install build-mkl --prefix /opt/lion-cpp-mkl
```

An existing installation is reused automatically on subsequent configurations;
no EULA flag is needed when no installation occurs. To prohibit any download:

```bash
cmake -S . -B build-mkl -DWITH_MKL=ON -DLION_MKL_AUTO_INSTALL=OFF \
  -DLION_MKL_ROOT=/opt/intel/oneapi/mkl/latest
```

You may alternatively pass `-DMKL_DIR=/path/to/mkl/lib/cmake/mkl` or use the
`MKLROOT` environment variable. To share a local installation across multiple
independent build trees, set `-DLION_MKL_INSTALL_PREFIX=/your/cache/oneapi`.
A cached installer can be supplied with `-DLION_MKL_INSTALLER=/path/offline.sh`
to avoid network access. If Intel updates the distribution, set both
`LION_MKL_VERSION` and `LION_MKL_INSTALLER_URL` (the default is deliberately
version-pinned, not a moving `latest` URL). For additional download integrity
verification, set `LION_MKL_INSTALLER_SHA384` to a trusted digest.

**Notes:** The standalone installer is large and the first configuration may
need substantial disk space. The download uses HTTPS, but no authoritative
checksum was available for the standalone 2026.2 installer when this integration
was authored; the SHA384 check is optional rather than a fabricated pin.
The installation copies the *entire* oneMKL development directory, including
headers, shared and static libraries and Intel's CMake package. Allow sufficient
disk space. You need write permissions (or `sudo`) only for `cmake --install`.
The build and Intel download still run as your normal user.

By default, consumers need only `find_package(lion CONFIG REQUIRED)` with
`CMAKE_PREFIX_PATH=/opt/lion-cpp-mkl`; no `MKLROOT`, additional MKL install or
access to the original build tree is required. An explicit consumer-side
`-DLION_MKL_ROOT=/another/mkl/root` can override the bundled copy. To keep MKL
externally managed instead, configure lion with
`-DLION_INSTALL_MKL_WITH_LION=OFF` (the previous behavior).
C/C++-only consumers do **not** need to enable Fortran. The exported package
provides the same `gf_lp64` MKL link interface when Fortran is disabled,
while Fortran-enabled consumers continue to use Intel's `MKLConfig.cmake`.
OpenMP is requested only for the languages enabled in the consuming project.
Verify with `cmake/tests/test-mkl-consumer.sh <lion install prefix>`.
The offline `cmake/tests/test-mkl-bundled-install.sh` also exercises install
relocation and HSL/MUMPS runtime paths using small fake MKL libraries.

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
