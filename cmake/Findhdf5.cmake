if (NOT hdf5_FOUND)
    find_path(HDF5_INCLUDE_DIR hdf5.h PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include)

    find_library(HDF5_LIBRARY NAMES hdf5 PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 
    find_library(HDF5_HL_LIBRARY NAMES hdf5_hl PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 
    find_library(HDF5_TOOLS_LIBRARY NAMES hdf5_tools PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 
    find_library(HDF5_CPP_LIBRARY NAMES hdf5_cpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 
    find_library(HDF5_CPP_HL_LIBRARY NAMES hdf5_hl_cpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib HINTS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH) 

    find_package(ZLIB)
    if (WITH_SZIP)
    	find_package(SZIP)
    endif()

    if (HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
        set(hdf5_FOUND YES)
	
	    get_filename_component(HDF5_LIB_NAME ${HDF5_LIBRARY} NAME)
	    get_filename_component(HDF5_HL_LIB_NAME ${HDF5_HL_LIBRARY} NAME)
	    get_filename_component(HDF5_TOOLS_LIB_NAME ${HDF5_TOOLS_LIBRARY} NAME)
	    get_filename_component(HDF5_CPP_LIB_NAME ${HDF5_CPP_LIBRARY} NAME)
	    get_filename_component(HDF5_CPP_HL_LIB_NAME ${HDF5_CPP_HL_LIBRARY} NAME)

	    if (APPLE)
		    # Change relative paths by absolute path
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${HDF5_LIB_NAME}" ${HDF5_LIBRARY})
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${HDF5_HL_LIB_NAME}" ${HDF5_HL_LIBRARY})
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${HDF5_TOOLS_LIB_NAME}" ${HDF5_TOOLS_LIBRARY})
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${HDF5_CPP_LIB_NAME}" ${HDF5_CPP_LIBRARY})
		    execute_process(COMMAND ${CMAKE_INSTALL_NAME_TOOL} -id "@rpath/${HDF5_CPP_HL_LIB_NAME}" ${HDF5_CPP_HL_LIBRARY})
    	endif()


	    add_library(hdf5 INTERFACE)


	    file(GLOB HDF5_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libhdf5**")
	
	    set(HDF5_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${HDF5_LIB_NAME})
	    set(HDF5_HL_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${HDF5_HL_LIB_NAME})
	    set(HDF5_TOOLS_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${HDF5_TOOLS_LIB_NAME})
	    set(HDF5_CPP_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${HDF5_CPP_LIB_NAME})
	    set(HDF5_CPP_HL_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${HDF5_CPP_HL_LIB_NAME})

        target_link_libraries(hdf5 INTERFACE 
               "$<BUILD_INTERFACE:${HDF5_CPP_HL_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${HDF5_CPP_HL_LIB_NAME}>")

        target_link_libraries(hdf5 INTERFACE 
               "$<BUILD_INTERFACE:${HDF5_CPP_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${HDF5_CPP_LIB_NAME}>")

        target_link_libraries(hdf5 INTERFACE 
               "$<BUILD_INTERFACE:${HDF5_TOOLS_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${HDF5_TOOLS_LIB_NAME}>")

        target_link_libraries(hdf5 INTERFACE 
               "$<BUILD_INTERFACE:${HDF5_HL_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${HDF5_HL_LIB_NAME}>")

        target_link_libraries(hdf5 INTERFACE 
               "$<BUILD_INTERFACE:${HDF5_LIBRARY}>"
               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${HDF5_LIB_NAME}>")

        target_link_libraries(hdf5 INTERFACE ZLIB::ZLIB)

	if (WITH_SZIP)
        	target_link_libraries(hdf5 INTERFACE SZIP::SZIP)
	endif()


        foreach(t ${HDF5_LIB_ALL_FILES})
    	    install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(hdf5_FOUND NO)
    endif()
endif()
