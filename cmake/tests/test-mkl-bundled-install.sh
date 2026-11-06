#!/usr/bin/env bash
# Offline integration check: oneMKL bundle, relocatable CMake config, HSL/MUMPS
# RUNPATH changes and a C/C++ consumer without an enabled Fortran language.
set -euo pipefail

cmake_dir="$(cd "$(dirname "$0")/.." && pwd)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
root="$work/build_directory_with_a_long_name/thirdparty/oneapi/mkl/2026.2"
mkdir -p "$root/include" "$root/lib/cmake/mkl" "$work/src" "$work/consumer"
cat > "$root/include/mkl.h" <<'HEADER'
#ifdef __cplusplus
extern "C" {
#endif
int fake_mkl(void);
#ifdef __cplusplus
}
#endif
HEADER
cat > "$work/mkl.c" <<'C'
int fake_mkl(void) { return 17; }
C
cat > "$work/consumer.c" <<'C'
extern int fake_mkl(void);
int hsl(void) { return fake_mkl(); }
int mumps(void) { return fake_mkl(); }
C
for lib in gf_lp64 gnu_thread core; do
  gcc -fPIC -shared "$work/mkl.c" -Wl,-soname,"libmkl_$lib.so.3" \
      -o "$root/lib/libmkl_$lib.so.3"
  ln -s "libmkl_$lib.so.3" "$root/lib/libmkl_$lib.so"
done
cat > "$root/lib/cmake/mkl/MKLConfig.cmake" <<'CMAKE'
get_filename_component(_fake_root "${CMAKE_CURRENT_LIST_DIR}/../../.." ABSOLUTE)
if(NOT TARGET MKL::MKL)
    add_library(MKL::MKL INTERFACE IMPORTED)
    set_target_properties(MKL::MKL PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_fake_root}/include"
        INTERFACE_LINK_LIBRARIES "${_fake_root}/lib/libmkl_gf_lp64.so")
endif()
CMAKE

for lib in coinhsl coinmumps; do
    gcc -fPIC -shared "$work/consumer.c" -Wl,-soname,"lib$lib.so.2" \
        -Wl,--no-as-needed -L"$root/lib" -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core \
        -Wl,-rpath,"$root/lib" -o "$work/src/lib$lib.so.2"
done
cat > "$work/src/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.13)
project(lion VERSION 1.0 LANGUAGES C CXX)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(WITH_MKL ON)
set(MKL_LINK "dynamic")
set(MKL_THREADING "gnu_thread")
set(MKL_INTERFACE_FULL "gf_lp64")
set(MKL_ROOT "@FAKE_ROOT@")
set(MKL_LIBDIR "${MKL_ROOT}/lib")
set(project_name "lion")
set(targets_export_name "lionTargets")
# third-party ELF libraries are installed before RPATH relocation is run
install(FILES "@FAKE_WORK@/src/libcoinhsl.so.2" "@FAKE_WORK@/src/libcoinmumps.so.2"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}")
include("@MODULE_DIR@/install-mkl.cmake")
add_library(blaslapack INTERFACE)
target_link_libraries(blaslapack INTERFACE MKL::MKL)
install(TARGETS blaslapack EXPORT lionTargets)
install(EXPORT lionTargets NAMESPACE lion:: DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lion")
configure_package_config_file("@MODULE_DIR@/Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lionConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lion")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lionConfig.cmake" DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lion")
CMAKE
# Don't expand CMake ${vars} in the template, only our @...@ placeholders.
FAKE_ROOT="$root" FAKE_WORK="$work" MODULE_DIR="$cmake_dir" python - <<'PY'
from pathlib import Path
import os
p=Path(os.environ['FAKE_WORK'])/'src'/'CMakeLists.txt'
s=p.read_text()
for key in ('FAKE_ROOT','FAKE_WORK','MODULE_DIR'):
    s=s.replace('@'+key+'@',os.environ[key])
p.write_text(s)
PY
cmake -S "$work/src" -B "$work/cmake-build" -DCMAKE_INSTALL_PREFIX="$work/original_prefix" >/dev/null
cmake --install "$work/cmake-build" --prefix "$work/install_a" >/dev/null
lib="$work/install_a/lib/libcoinhsl.so.2"
readelf -d "$lib" | grep -F '$ORIGIN' >/dev/null || {
    echo "FAIL: HSL RPATH was not rewritten to an origin-relative location" >&2
    readelf -d "$lib" >&2
    exit 1
}
# Removing the entire original 'build directory' must not break the install.
rm -rf "$work/build_directory_with_a_long_name"
mv "$work/install_a" "$work/relocated_install"
readelf -d "$work/relocated_install/lib/libcoinmumps.so.2" | grep -F '$ORIGIN' >/dev/null
ldd "$work/relocated_install/lib/libcoinhsl.so.2" | grep -F "oneapi/mkl/2026.2/lib/" >/dev/null
# The lionConfig file must not embed the original build directory at all.
if grep -F "$root" "$work/relocated_install/lib/cmake/lion/lionConfig.cmake"; then
    echo "FAIL: installed lionConfig.cmake still embeds the build-time MKL root" >&2
    exit 1
fi
cat > "$work/consumer/CMakeLists.txt" <<'CMAKE'
cmake_minimum_required(VERSION 3.13)
project(mkl_consumer LANGUAGES C CXX)
find_package(lion CONFIG REQUIRED PATHS "${LION_INSTALL_PREFIX}" NO_DEFAULT_PATH)
add_executable(mkltest main.cpp)
target_link_libraries(mkltest PRIVATE lion::blaslapack)
CMAKE
cat > "$work/consumer/main.cpp" <<'CPP'
#include <mkl.h>
int main() { return fake_mkl() == 17 ? 0 : 1; }
CPP
cmake -S "$work/consumer" -B "$work/consumer-build" \
  -DLION_INSTALL_PREFIX="$work/relocated_install" >/dev/null
cmake --build "$work/consumer-build" >/dev/null
"$work/consumer-build/mkltest"
# The Fortran-enabled branch should also resolve MKL from the installed tree.
sed -i 's/LANGUAGES C CXX)/LANGUAGES C CXX Fortran)/' "$work/consumer/CMakeLists.txt"
cmake -S "$work/consumer" -B "$work/consumer-fortran" \
  -DLION_INSTALL_PREFIX="$work/relocated_install" >/dev/null
cmake --build "$work/consumer-fortran" >/dev/null
"$work/consumer-fortran/mkltest"
echo "PASS: bundled MKL survives relocation; HSL/MUMPS RUNPATH and C/C++/Fortran-enabled consumers work."
