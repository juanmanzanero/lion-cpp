if (BUILD_CLI11)
    include(ExternalProject)

    # CLI11 v2.7.2 amalgamated header. No configure or compilation needed.
    ExternalProject_Add(cli11
        URL https://github.com/CLIUtils/CLI11/releases/download/v2.7.2/CLI11.hpp
        URL_HASH SHA256=ffa9a30da295c5858fb5f91f9f45771bab09471d7010a34c7c68c857a330dd76
        DOWNLOAD_NO_EXTRACT TRUE
        PREFIX "${THIRD_PARTY_DIR}/cli11"
        INSTALL_DIR "${THIRD_PARTY_DIR}/lion/thirdparty"
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND "${CMAKE_COMMAND}" -E make_directory "${THIRD_PARTY_DIR}/lion/thirdparty/include/CLI"
            COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                <DOWNLOADED_FILE> "${THIRD_PARTY_DIR}/lion/thirdparty/include/CLI/CLI11.hpp"
    )
endif()
