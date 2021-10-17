if ( ${BUILD_MATPLOTLIBCPP} )
    include (ExternalProject)
    ExternalProject_Add(matplotlibcpp
      GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
      GIT_TAG master
      PREFIX "${THIRD_PARTY_DIR}/matplotlibcpp"
      SOURCE_DIR ${THIRD_PARTY_DIR}/matplotlibcpp/source
      BINARY_DIR ${THIRD_PARTY_DIR}/matplotlibcpp/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND cp ${THIRD_PARTY_DIR}/matplotlibcpp/source/matplotlibcpp.h ${THIRD_PARTY_DIR}/lion/thirdparty/include && cd ${THIRD_PARTY_DIR}/lion/thirdparty/include && patch -p1 < ${PATCH_DIR}/matplotlibcpp.patch
    )

endif()
