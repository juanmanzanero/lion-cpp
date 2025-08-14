if (NOT zstd_FOUND)
    find_path(ZSTD_INCLUDE_DIR zstd/zstd.h PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include NO_DEFAULT_PATH)
    find_library(ZSTD_LIBRARY NAMES zstd PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/lib NO_DEFAULT_PATH)

    if (ZSTD_INCLUDE_DIR AND ZSTD_LIBRARY)
        set(zstd_FOUND YES)
        add_library(zstd INTERFACE)

    target_include_directories(zstd SYSTEM INTERFACE
         "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
         "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

    file(GLOB ZSTD_LIB_ALL_FILES "${CMAKE_BINARY_DIR}/lion/thirdparty/lib/libzstd*")

    get_filename_component(ZSTD_LIB_NAME ${ZSTD_LIBRARY} NAME)
    set(ZSTD_LIBRARY_INSTALL ${CMAKE_INSTALL_FULL_LIBDIR}/${ZSTD_LIB_NAME})

        target_link_libraries(zstd INTERFACE 
        "$<BUILD_INTERFACE:${ZSTD_LIBRARY}>"
        "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_LIBDIR}/${ZSTD_LIB_NAME}>")

        foreach(t ${ZSTD_LIB_ALL_FILES})
            install(FILES "${t}" TYPE LIB) 
        endforeach()

    else()
        set(zstd_FOUND NO)
    endif()
endif()
