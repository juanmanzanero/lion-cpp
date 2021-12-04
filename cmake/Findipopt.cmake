if (NOT ipopt_FOUND)
    find_path(IPOPT_INCLUDE_DIR IpIpoptApplication.hpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include/coin-or NO_DEFAULT_PATH)
    find_library(IPOPT_LIBRARY NAMES ipopt PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH)
    find_library(MUMPS_LIBRARY NAMES coinmumps PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 


    if (IPOPT_INCLUDE_DIR AND IPOPT_LIBRARY)
        set(ipopt_FOUND YES)
	
	get_filename_component(IPOPT_LIB_NAME ${IPOPT_LIBRARY} NAME)
	get_filename_component(MUMPS_LIB_NAME ${MUMPS_LIBRARY} NAME)

	if (APPLE)
		# Change relative paths by absolute path
		execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${IPOPT_LIB_NAME}" ${IPOPT_LIBRARY})
	endif()
    
	add_library(ipopt INTERFACE)
	target_link_libraries(ipopt INTERFACE ${CMAKE_DL_LIBS})
	target_include_directories(ipopt SYSTEM INTERFACE
		  "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
		  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

	file(GLOB IPOPT_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libipopt*")
	file(GLOB MUMPS_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libcoinmumps*")
	
	set(IPOPT_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${IPOPT_LIB_NAME})

        target_link_libraries(ipopt INTERFACE 
		"$<BUILD_INTERFACE:${IPOPT_LIBRARY}>"
 		"$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${IPOPT_LIB_NAME}>")

        target_link_libraries(ipopt INTERFACE 
               "$<BUILD_INTERFACE:${MUMPS_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${MUMPS_LIB_NAME}>")

        target_link_directories(ipopt INTERFACE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
        target_link_libraries(ipopt INTERFACE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

        foreach(t ${IPOPT_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

        foreach(t ${MUMPS_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(ipopt_FOUND NO)
    endif()
endif()
