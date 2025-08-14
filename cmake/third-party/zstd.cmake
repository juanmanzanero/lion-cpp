if ( ${BUILD_ZSTD} )

include (ExternalProject)
ExternalProject_Add(zstd
  GIT_REPOSITORY https://github.com/facebook/zstd.git
  GIT_TAG v1.5.7
  PREFIX "${THIRD_PARTY_DIR}/zstd"
  SOURCE_DIR "${THIRD_PARTY_DIR}/zstd/source"
  BINARY_DIR "${THIRD_PARTY_DIR}/zstd/build"
  INSTALL_DIR "${THIRD_PARTY_DIR}/lion/thirdparty"

  # we cannot use ${CMAKE_COMMAND} here!! Only a simple "cmake" works
  CONFIGURE_COMMAND cmake -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${THIRD_PARTY_DIR}/lion/thirdparty
                -DCMAKE_INSTALL_INCLUDEDIR:STRING=include/zstd
                -DCMAKE_INSTALL_LIBDIR:STRING=lib
                -DCMAKE_BUILD_TYPE:STRING=Release
                -DCMAKE_CXX_STANDARD=11
                -DCMAKE_CXX_FLAGS=${windows_flag}
                -DZSTD_BUILD_STATIC=ON
                -DZSTD_BUILD_SHARED=OFF
                -DZSTD_BUILD_PROGRAMS=OFF
                ${THIRD_PARTY_DIR}/zstd/source/build/cmake
)

endif()
