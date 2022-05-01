if (NOT blaslapack_FOUND)
    # Try to find it in the thirdparty folder
    find_library(BLAS_LIBRARY NAMES blas PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH)
    find_library(LAPACK_LIBRARY NAMES lapack PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 

    if (BLAS_LIBRARY AND LAPACK_LIBRARY)
        set(blaslapack_FOUND YES)

        get_filename_component(LAPACK_LIB_NAME ${LAPACK_LIBRARY} NAME)
        get_filename_component(BLAS_LIB_NAME ${BLAS_LIBRARY} NAME)


        add_library(blaslapack INTERFACE)
        target_link_libraries(blaslapack INTERFACE ${CMAKE_DL_LIBS})

        set(IPOPT_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${IPOPT_LIB_NAME})

        target_link_libraries(blaslapack INTERFACE 
                              "$<BUILD_INTERFACE:${LAPACK_LIBRARY}>"
                              "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${LAPACK_LIB_NAME}>")

        target_link_libraries(blaslapack INTERFACE 
                              "$<BUILD_INTERFACE:${BLAS_LIBRARY}>"
				  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${BLAS_LIB_NAME}>")

        # Install them
	file(GLOB LAPACK_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/liblapack*")
	file(GLOB BLAS_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libblas*")
	
        foreach(t ${LAPACK_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

        foreach(t ${BLAS_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(blaslapack_FOUND NO)
    endif()
endif()
