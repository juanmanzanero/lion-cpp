if (NOT cppad_FOUND)
    find_path(CPPAD_INCLUDE_DIR cppad.hpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include/cppad)

    find_library(CPPAD_LIBRARY NAMES cppad_lib PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 

    if (CPPAD_INCLUDE_DIR AND CPPAD_LIBRARY)
        set(cppad_FOUND YES)
	
	    get_filename_component(CPPAD_LIB_NAME ${CPPAD_LIBRARY} NAME)

	    if (APPLE)
		    # Change relative paths by absolute path
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${CPPAD_LIB_NAME}" ${CPPAD_LIBRARY})
    	endif()


	    add_library(cppad INTERFACE)


	    file(GLOB CPPAD_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libcppad_lib*")
	
	    set(CPPAD_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${CPPAD_LIB_NAME})

        target_link_libraries(cppad INTERFACE 
               "$<BUILD_INTERFACE:${CPPAD_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${CPPAD_LIB_NAME}>")


        foreach(t ${CPPAD_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(cppad_FOUND NO)
    endif()
endif()
