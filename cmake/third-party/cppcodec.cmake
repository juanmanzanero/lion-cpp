if (${BUILD_CPPCODEC})

include (ExternalProject)
ExternalProject_Add(cppcodec
  GIT_REPOSITORY https://github.com/tplgy/cppcodec.git
  GIT_TAG v0.2
  PREFIX "${THIRD_PARTY_DIR}/cppcodec"
  SOURCE_DIR "${THIRD_PARTY_DIR}/cppcodec/source" 
  BINARY_DIR "${THIRD_PARTY_DIR}/cppcodec/build" 
  INSTALL_DIR "${THIRD_PARTY_DIR}/lion/thirdparty"

  CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -DBUILD_TESTING=Off
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${THIRD_PARTY_DIR}/lion/thirdparty
                -DCMAKE_INSTALL_BINDIR:STRING=lib
                -DCMAKE_BUILD_TYPE:STRING=Release
                -DCMAKE_CXX_STANDARD=11
                -DCMAKE_CXX_FLAGS=${windows_flag}
                -DBUILD_TESTING=OFF
                -DCPPCODEC_BUILD_TOOLS=OFF
                ${THIRD_PARTY_DIR}/cppcodec/source
)

endif()
