include(python)

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

# Python and matplotlib cpp
if ( ${Python3_FOUND} AND ${Python3_MATPLOTLIB} )
    find_package(matplotlibcpp)
    if (NOT ${matplotlibcpp_FOUND})
        set(BUILD_MATPLOTLIBCPP YES)
    endif()
else()
    configure_file(cmake/third-party/fake_matplotlibcpp.h ${CMAKE_BINARY_DIR}/thirdparty/include/matplotlibcpp.h COPYONLY)
endif()

# Logger cpp
find_package(loggercpp)

if (NOT ${loggercpp_FOUND})
    set(BUILD_LOGGERCPP YES)
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

find_package(tinyxml2)
find_package(ipopt)
if ( ${Python3_FOUND} AND ${Python3_MATPLOTLIB} )
    find_package(matplotlibcpp)
endif()
find_package(loggercpp)
