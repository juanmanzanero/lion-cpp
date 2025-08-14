if (NOT cppcodec_FOUND)
    find_path(CPPCODEC_INCLUDE_DIR cppcodec/base64_rfc4648.hpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include NO_DEFAULT_PATH)

    if (CPPCODEC_INCLUDE_DIR)
        set(cppcodec_FOUND YES)
        add_library(cppcodec INTERFACE)

        target_include_directories(cppcodec SYSTEM INTERFACE
            "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
            "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

    else()
        set(cppcodec_FOUND NO)
    endif()
endif()
