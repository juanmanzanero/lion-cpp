macro(new_test)
    get_filename_component(binary_name ${CMAKE_CURRENT_LIST_DIR} NAME)
    string(REPLACE " " "_" binary_name ${binary_name})
    
    set(BINARY ${binary_name}_test)
    
    file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp ../main.cpp)
    
    set(SOURCES ${TEST_SOURCES})
    
    add_executable(${BINARY} ${TEST_SOURCES})
    
    add_test(NAME ${BINARY} COMMAND ${BINARY})

    if ( NOT APPLE)
        target_link_options(${BINARY} PUBLIC -Wl,--no-as-needed -ldl)
    endif()

    # Link libraries
    target_link_libraries(${BINARY} LINK_PRIVATE GTest::gtest lion)
    
    # Copy required files
    add_custom_target(${BINARY}_link_data ALL
                      COMMAND ${CMAKE_COMMAND} -E create_symlink 
                      ${CMAKE_CURRENT_SOURCE_DIR}/data 
                      ${CMAKE_BINARY_DIR}/test/${binary_name}/data
                     )
    if (MSYS)
	ADD_CUSTOM_TARGET(${BINARY}_link_tinyxml_windows ALL
	                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/msys-tinyxml2-9.dll ${CMAKE_BINARY_DIR}/test/${binary_name}/msys-tinyxml2-9.dll)
	ADD_CUSTOM_TARGET(${BINARY}_link_mumps_windows ALL
	                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_BINARY_DIR}/lion/thirdparty/bin/libcoinmumps-3.dll ${CMAKE_BINARY_DIR}/test/${binary_name}/libcoinmumps-3.dll)			  
    endif()
endmacro()

add_subdirectory(./foundation)
add_subdirectory(./frame)
add_subdirectory(./io)
add_subdirectory(./math)
add_subdirectory(./propagators)


# Create a test if code coverage
if (CODE_COVERAGE)
	include(code-coverage)

	add_code_coverage_all_targets(EXCLUDE test/* googletest-src/*)

	set(BINARY coverage_test)

	file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false **/*.cpp main.cpp)

	set(SOURCES ${TEST_SOURCES})

	add_executable(${BINARY} ${TEST_SOURCES})

    target_link_libraries(${BINARY} LINK_PRIVATE GTest::gtest lion)

	set_target_properties(${BINARY} PROPERTIES 
        EXCLUDE_FROM_ALL 1 
        EXCLUDE_FROM_DEFAULT_BUILD 1
    )

	target_code_coverage(${BINARY} AUTO ALL EXCLUDE test/* build/googletest-src/*) 

    # Copy required files
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test/data)

    add_custom_command(TARGET ${BINARY} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/test/foundation/data/ ${CMAKE_BINARY_DIR}/test/data)

    add_custom_command(TARGET ${BINARY} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/test/frame/data/ ${CMAKE_BINARY_DIR}/test/data)

    add_custom_command(TARGET ${BINARY} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/test/io/data/ ${CMAKE_BINARY_DIR}/test/data)

    add_custom_command(TARGET ${BINARY} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/test/math/data/ ${CMAKE_BINARY_DIR}/test/data)

    add_custom_command(TARGET ${BINARY} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/test/propagators/data/ ${CMAKE_BINARY_DIR}/test/data)

endif()
