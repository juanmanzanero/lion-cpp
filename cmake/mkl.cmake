# Linux + gfortran integration of Intel oneMKL as the sole BLAS/LAPACK provider.
# MKL::MKL is used for CMake consumers; MKL_FLAGS mirrors that selection for
# the Autotools-based Ipopt, MUMPS and (if enabled) HSL builds.
if(NOT CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
    message(FATAL_ERROR "WITH_MKL is currently supported only on Linux")
endif()
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "WITH_MKL currently requires GNU Fortran (gfortran) and the LP64 interface")
endif()

set(MKL_LINK "dynamic" CACHE STRING "MKL link mode: dynamic or static")
set_property(CACHE MKL_LINK PROPERTY STRINGS dynamic static)
set(MKL_THREADING "gnu_thread" CACHE STRING "MKL threading: gnu_thread or sequential")
set_property(CACHE MKL_THREADING PROPERTY STRINGS gnu_thread sequential)
set(MKL_INTERFACE_FULL "gf_lp64" CACHE STRING "MKL Fortran LP64 interface")

if(NOT MKL_LINK STREQUAL "dynamic" AND NOT MKL_LINK STREQUAL "static")
    message(FATAL_ERROR "WITH_MKL requires MKL_LINK=dynamic or static (got '${MKL_LINK}')")
endif()
if(NOT MKL_THREADING STREQUAL "gnu_thread" AND NOT MKL_THREADING STREQUAL "sequential")
    message(FATAL_ERROR "WITH_MKL requires MKL_THREADING=gnu_thread or sequential (got '${MKL_THREADING}')")
endif()
if(NOT MKL_INTERFACE_FULL STREQUAL "gf_lp64")
    message(FATAL_ERROR "WITH_MKL requires MKL_INTERFACE_FULL=gf_lp64 for gfortran (got '${MKL_INTERFACE_FULL}')")
endif()

if(MKL_THREADING STREQUAL "gnu_thread")
    find_package(OpenMP REQUIRED COMPONENTS C CXX Fortran)
endif()
find_package(MKL CONFIG REQUIRED)
if(NOT TARGET MKL::MKL)
    message(FATAL_ERROR "MKLConfig.cmake did not provide target MKL::MKL")
endif()

# Newer oneAPI releases use lib/ instead of lib/intel64/.
if(NOT MKL_ROOT)
    get_filename_component(MKL_ROOT "${MKL_DIR}/../../.." ABSOLUTE)
endif()
set(MKL_INCDIR "${MKL_ROOT}/include")
if(NOT EXISTS "${MKL_INCDIR}/mkl.h")
    message(FATAL_ERROR "Could not find MKL headers at ${MKL_INCDIR}; set MKL_DIR or MKLROOT")
endif()

if(MKL_LINK STREQUAL "dynamic")
    set(_lion_mkl_suffix "so")
else()
    set(_lion_mkl_suffix "a")
endif()
set(MKL_LIBDIR "")
foreach(_libdir "${MKL_ROOT}/lib" "${MKL_ROOT}/lib/intel64" "${MKL_ROOT}/lib/intel64_lin")
    if(EXISTS "${_libdir}/libmkl_gf_lp64.${_lion_mkl_suffix}"
       AND EXISTS "${_libdir}/libmkl_${MKL_THREADING}.${_lion_mkl_suffix}"
       AND EXISTS "${_libdir}/libmkl_core.${_lion_mkl_suffix}")
        set(MKL_LIBDIR "${_libdir}")
        break()
    endif()
endforeach()
if(NOT MKL_LIBDIR)
    message(FATAL_ERROR "Missing ${MKL_LINK} MKL LP64/${MKL_THREADING} libraries below ${MKL_ROOT}/lib; install the oneMKL development libraries")
endif()
if(MKL_LIBDIR MATCHES "[ ;]")
    message(FATAL_ERROR "MKL library paths containing whitespace or semicolons are not supported by the Ipopt Autotools build: ${MKL_LIBDIR}")
endif()

if(MKL_LINK STREQUAL "dynamic")
    set(MKL_FLAGS "-L${MKL_LIBDIR} -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_${MKL_THREADING} -lmkl_core")
else()
    # A single linker option prevents libtool from reordering the cyclic archives.
    set(MKL_FLAGS "-Wl,--start-group,${MKL_LIBDIR}/libmkl_gf_lp64.a,${MKL_LIBDIR}/libmkl_${MKL_THREADING}.a,${MKL_LIBDIR}/libmkl_core.a,--end-group")
endif()
if(MKL_THREADING STREQUAL "gnu_thread")
    string(APPEND MKL_FLAGS " -lgomp")
endif()
string(APPEND MKL_FLAGS " -lpthread -lm -ldl")

# Detect a missing PARDISO symbol before attempting to build Ipopt. Verify
# both the exported CMake target and the flags passed to Autotools.
include(CheckCXXSourceCompiles)
set(_lion_pardiso_probe [=[
#include <mkl.h>
int main() {
    auto* volatile pardiso_symbol = &pardiso;
    return pardiso_symbol == nullptr;
}
]=])
set(_lion_old_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
set(_lion_old_required_includes "${CMAKE_REQUIRED_INCLUDES}")

set(CMAKE_REQUIRED_LIBRARIES MKL::MKL)
unset(LION_MKL_TARGET_HAS_PARDISO CACHE)
check_cxx_source_compiles("${_lion_pardiso_probe}" LION_MKL_TARGET_HAS_PARDISO)
if(NOT LION_MKL_TARGET_HAS_PARDISO)
    message(FATAL_ERROR "MKL::MKL cannot link PARDISO; check MKL_LINK, MKL_THREADING and oneMKL installation")
endif()

separate_arguments(_lion_mkl_link_args NATIVE_COMMAND "${MKL_FLAGS}")
set(CMAKE_REQUIRED_LIBRARIES ${_lion_mkl_link_args})
set(CMAKE_REQUIRED_INCLUDES "${MKL_INCDIR}")
unset(LION_MKL_AUTOTOOLS_HAS_PARDISO CACHE)
check_cxx_source_compiles("${_lion_pardiso_probe}" LION_MKL_AUTOTOOLS_HAS_PARDISO)
set(CMAKE_REQUIRED_LIBRARIES "${_lion_old_required_libraries}")
set(CMAKE_REQUIRED_INCLUDES "${_lion_old_required_includes}")
if(NOT LION_MKL_AUTOTOOLS_HAS_PARDISO)
    message(FATAL_ERROR "Autotools MKL_FLAGS cannot link PARDISO; inspect CMakeFiles/CMakeError.log")
endif()

add_library(blaslapack INTERFACE)
target_link_libraries(blaslapack INTERFACE MKL::MKL)
set(blaslapack_FOUND TRUE)
message(STATUS "lion: using oneMKL (${MKL_LINK}, ${MKL_THREADING}, ${MKL_INTERFACE_FULL}) from ${MKL_ROOT}")
