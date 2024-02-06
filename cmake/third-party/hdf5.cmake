if ( ${BUILD_HDF5} )
    include (ExternalProject)

    ExternalProject_Add(hdf5
      GIT_REPOSITORY https://github.com/HDFGroup/hdf5.git
      GIT_TAG hdf5-1_13_2 # we use this branch because latest ones require cmake 3.18...
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${THIRD_PARTY_DIR}/lion/thirdparty -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DBUILD_TESTING=No -DBUILD_SHARED_LIBS=No -DBUILD_STATIC_LIBS=Yes -DHDF5_BUILD_CPP_LIB=Yes
      PREFIX "${THIRD_PARTY_DIR}/hdf5"
      SOURCE_DIR ${THIRD_PARTY_DIR}/hdf5/source
      BINARY_DIR ${THIRD_PARTY_DIR}/hdf5/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
    )

endif()
