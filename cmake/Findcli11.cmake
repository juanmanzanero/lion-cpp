# CLI11 is provided as a single header in the lion third-party include tree.
find_path(CLI11_INCLUDE_DIR CLI/CLI11.hpp
    PATHS "${CMAKE_BINARY_DIR}/lion/thirdparty/include"
    NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cli11 DEFAULT_MSG CLI11_INCLUDE_DIR)

if (cli11_FOUND AND NOT TARGET cli11)
    add_library(cli11 INTERFACE)

    target_include_directories(cli11 SYSTEM INTERFACE
        "$<BUILD_INTERFACE:${CLI11_INCLUDE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/lion/thirdparty/include>")
endif()
