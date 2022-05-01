# Find lapack
set(BUILD_LAPACK NO)
find_package(blaslapack) 

if (NOT ${blaslapack_FOUND})
    set(BUILD_LAPACK YES)
endif()

if ( ${ENABLE_TEST} )
    # Google tests 
    find_package(GTest PATHS ${CMAKE_BINARY_DIR})
    
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


#######################################################################




##### BUILD ALL REQUIRED THIRD PARTY LIBRARIES ##########
message(STATUS "Compilation of the required third party libraries")
configure_file(cmake/third-party/CMakeLists.txt ${CMAKE_BINARY_DIR}/thirdparty/CMakeLists.txt)

execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/thirdparty"
)
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/thirdparty"
)
message("")
message(STATUS "Configuration of fastest-lap")
#######################################################################

if ( ${ENABLE_TEST} )
	find_package(GTest PATHS ${CMAKE_BINARY_DIR}/thirdparty REQUIRED)
endif()

find_package(blaslapack REQUIRED)
find_package(tinyxml2 REQUIRED)
find_package(ipopt REQUIRED)
find_package(loggercpp REQUIRED)
find_package(cppad REQUIRED)

if (MSYS)
	file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-blas-3.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
	file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-lapack-3.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
	file(COPY ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-tinyxml2-9.dll DESTINATION ${CMAKE_BINARY_DIR}/lion/thirdparty/lib)
endif()
