# Select exactly one BLAS/LAPACK provider for lion and all of its dependencies.
set(BUILD_LAPACK NO)
if (WITH_MKL)
    include(mkl)
else()
    find_package(blaslapack)
    if (NOT blaslapack_FOUND)
        set(BUILD_LAPACK YES)
    endif()
endif()

# A previously configured Ipopt/MUMPS cannot safely change BLAS/LAPACK ABI or
# linking mode in place. Require a separate build directory on any change.
set(_lion_math_backend "reference|WITH_HSL=${WITH_HSL}")
if (WITH_MKL)
    set(_lion_math_backend "MKL|${MKL_ROOT}|${MKL_LINK}|${MKL_THREADING}|${MKL_INTERFACE_FULL}|WITH_HSL=${WITH_HSL}")
endif()
set(_lion_math_backend_file "${CMAKE_BINARY_DIR}/lion-math-backend.txt")
if (EXISTS "${_lion_math_backend_file}")
    file(READ "${_lion_math_backend_file}" _lion_previous_math_backend)
    if (NOT _lion_previous_math_backend STREQUAL _lion_math_backend)
        message(FATAL_ERROR "BLAS/LAPACK configuration has changed in this build directory. Use a fresh build directory for WITH_MKL, MKL settings or WITH_HSL changes.")
    endif()
elseif(WITH_MKL AND EXISTS "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libipopt.a")
    message(FATAL_ERROR "Found an existing Ipopt build without a recorded BLAS/LAPACK provider. Use a fresh build directory for WITH_MKL=ON.")
else()
    file(WRITE "${_lion_math_backend_file}" "${_lion_math_backend}")
endif()


if ( ${ENABLE_TEST} )
    # Google tests 
    find_package(GTest PATHS ${CMAKE_BINARY_DIR} NO_DEFAULT_PATH)
    
    if ( NOT ${GTest_FOUND})
        set(BUILD_GTEST YES)
    endif()
else()
	set(BUILD_GTEST NO)
endif()

# Tinyxml2
find_package(tinyxml2)

if (NOT ${tinyxml2_FOUND})
    set(BUILD_TINYXML YES)
endif()

# nlohmann-json
find_package(nlohmann_json)

if (NOT ${nlohmann_json_FOUND})
    set(BUILD_NLOHMANNJSON YES)
endif()


# Ipopt + dependencies
find_package(ipopt)

if (NOT ${ipopt_FOUND})
    set(BUILD_IPOPT YES)
endif()

# Logger cpp
find_package(loggercpp)

if (NOT ${loggercpp_FOUND})
    set(BUILD_LOGGERCPP YES)
endif()

# Cpp Automatic Differentiation
find_package(cppad)

if (NOT ${cppad_FOUND})
    set(BUILD_CPPAD YES)
endif()


# HDF5
find_package(hdf5)

if (NOT ${hdf5_FOUND})
    set(BUILD_HDF5 YES)
endif()


# cppcodec
find_package(cppcodec)

if (NOT ${cppcodec_FOUND})
    set(BUILD_CPPCODEC YES)
endif()

# ZSTD
find_package(zstd)

if (NOT ${zstd_FOUND})
    set(BUILD_ZSTD YES)
endif()

# CLI11 (header-only)
find_package(cli11)

if (NOT cli11_FOUND)
    set(BUILD_CLI11 YES)
endif()

#######################################################################

##### BUILD ALL REQUIRED THIRD PARTY LIBRARIES ##########
message(STATUS "Compilation of the required third party libraries")
configure_file(cmake/third-party/CMakeLists.txt ${CMAKE_BINARY_DIR}/thirdparty/CMakeLists.txt)

execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/thirdparty"
    RESULT_VARIABLE _lion_thirdparty_configure_status
)
if(NOT "${_lion_thirdparty_configure_status}" STREQUAL "0")
    message(FATAL_ERROR "lion: third-party CMake configuration failed: ${_lion_thirdparty_configure_status}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/thirdparty"
    RESULT_VARIABLE _lion_thirdparty_build_status
)
if(NOT "${_lion_thirdparty_build_status}" STREQUAL "0")
    message(FATAL_ERROR "lion: third-party build failed: ${_lion_thirdparty_build_status}. See the preceding ExternalProject error.")
endif()
message("")
message(STATUS "Configuration of fastest-lap")
#######################################################################

if ( ${ENABLE_TEST} )
    find_package(GTest PATHS ${CMAKE_BINARY_DIR}/thirdparty REQUIRED NO_DEFAULT_PATH)
endif()

if (NOT WITH_MKL)
    find_package(blaslapack REQUIRED)
endif()
find_package(tinyxml2 REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(ipopt REQUIRED)
find_package(loggercpp REQUIRED)
find_package(cppad REQUIRED)
find_package(hdf5 REQUIRED)
find_package(cppcodec REQUIRED)
find_package(zstd REQUIRED)
find_package(cli11 REQUIRED)

if (MSYS)
	file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-tinyxml2-9.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
        file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/libcoinmumps-3.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
	if (HAS_HSL)
            file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-coinhsl-2.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
	endif()
endif()
