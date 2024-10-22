if ( ${BUILD_LOGGERCPP} )
    include (ExternalProject)
    ExternalProject_Add(loggercpp
      GIT_REPOSITORY https://github.com/juanmanzanero/logger-cpp.git
      GIT_TAG main
      PREFIX "${THIRD_PARTY_DIR}/loggercpp"
      SOURCE_DIR ${THIRD_PARTY_DIR}/loggercpp/source
      BINARY_DIR ${THIRD_PARTY_DIR}/loggercpp/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND cp ${THIRD_PARTY_DIR}/loggercpp/source/logger.h ${THIRD_PARTY_DIR}/lion/thirdparty/include
    )

endif()
