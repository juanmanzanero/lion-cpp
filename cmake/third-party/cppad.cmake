if ( ${BUILD_CPPAD} )
    include (ExternalProject)
    ExternalProject_Add(cppad
      GIT_REPOSITORY https://github.com/coin-or/CppAD.git
      GIT_TAG 6856fdde535d2a7634a6fc36d8787a72166dce91
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${THIRD_PARTY_DIR}/lion/thirdparty -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      PREFIX "${THIRD_PARTY_DIR}/cppad"
      SOURCE_DIR ${THIRD_PARTY_DIR}/cppad/source
      BINARY_DIR ${THIRD_PARTY_DIR}/cppad/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
    )
    add_dependencies(cppad ipopt)

endif()
