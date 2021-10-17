if (NOT tinyxml2_FOUND)
    find_path(TINYXML_INCLUDE_DIR tinyxml2.h PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include NO_DEFAULT_PATH)
    find_library(TINYXML_LIBRARY NAMES tinyxml2 PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH)
    
    if (TINYXML_INCLUDE_DIR AND TINYXML_LIBRARY)
        set(tinyxml2_FOUND YES)
    
        add_library(tinyxml2 INTERFACE)

	target_include_directories(tinyxml2 SYSTEM INTERFACE
		  "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
		  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

	file(GLOB TINYXML_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libtinyxml*")

	get_filename_component(TINYXML_LIB_NAME ${TINYXML_LIBRARY} NAME)
	set(TINYXML_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${TINYXML_LIB_NAME})

        target_link_libraries(tinyxml2 INTERFACE 
		"$<BUILD_INTERFACE:${TINYXML_LIBRARY}>"
 		"$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${TINYXML_LIB_NAME}>")

        foreach(t ${TINYXML_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(tinyxml2_FOUND NO)
    endif()
endif()
