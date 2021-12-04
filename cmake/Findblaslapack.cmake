if (NOT blaslapack_FOUND)
    # Look for Lapack in the system
    if ( NOT BUILD_LAPACK )
        find_package(LAPACK)
    endif()

    if ( LAPACK_FOUND )
        set(blaslapack_FOUND YES)

        add_library(blaslapack INTERFACE)
	target_link_libraries(blaslapack INTERFACE LAPACK_LIBRARIES)

    else()

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
            install(FILES "${LAPACK_LIBRARY}" TYPE LIB) 
	    install(FILES "${BLAS_LIBRARY}" TYPE LIB) 

        else()
            set(blaslapack_FOUND NO)
        endif()
    endif()
endif()
