cmake_minimum_required(VERSION 3.10)

# set the project name
project(fastest-lap-third-party VERSION 1.0)

SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(THIRD_PARTY_DIR ${CMAKE_BINARY_DIR})
set(PATCH_DIR ${CMAKE_SOURCE_DIR}/cmake/patch)

set(BUILD_LAPACK ${BUILD_LAPACK})
set(BUILD_GTEST ${BUILD_GTEST})
set(BUILD_TINYXML ${BUILD_TINYXML})
set(BUILD_IPOPT ${BUILD_IPOPT})
set(BUILD_LOGGERCPP ${BUILD_LOGGERCPP})
set(BUILD_CPPAD ${BUILD_CPPAD})

set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_FORTRAN_COMPILER ${CMAKE_FORTRAN_COMPILER})


## Configure doxygen
include(third-party/lapack)
include(third-party/gtest)
include(third-party/tinyxml)
include(third-party/ipopt)
include(third-party/loggercpp)
include(third-party/cppad)
