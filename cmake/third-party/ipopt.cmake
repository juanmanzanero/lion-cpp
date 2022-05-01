if (${BUILD_IPOPT})
    set(lapack_flags "--with-lapack-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib/ -llapack -lblas")
    if (${BUILD_LAPACK})
	set(depends "lapack")
    endif()

    if (MSYS)
	set(disable_dependency_tracking "--disable-dependency-tracking")
    endif()

    include (ExternalProject)
    ExternalProject_Add(mumps
      GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
      GIT_TAG 1063386179a14b34bf762cdefca051606bd09fc4
      PREFIX "${THIRD_PARTY_DIR}/mumps"
      SOURCE_DIR ${THIRD_PARTY_DIR}/mumps/source
      BINARY_DIR ${THIRD_PARTY_DIR}/mumps/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source && ./get.Mumps && ./configure --prefix=${THIRD_PARTY_DIR}/lion/thirdparty --enable-static=yes --enable-shared=no ${lapack_flags} ${disable_dependency_tracking} --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib --without-metis
      BUILD_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source && make && make install
      INSTALL_COMMAND ""
      UPDATE_COMMAND ""
      DEPENDS ${depends}
    )

    if (NOT MSYS)
	    ExternalProject_Add(ipopt
	      GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
	      GIT_TAG releases/3.14.3 
	      PREFIX "${THIRD_PARTY_DIR}/ipopt"
	      SOURCE_DIR ${THIRD_PARTY_DIR}/ipopt/source
	      BINARY_DIR ${THIRD_PARTY_DIR}/ipopt/build 
	      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
	      PATCH_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/source && git apply ${PATCH_DIR}/ipopt.patch --reject
	      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/build &&
			../source/configure CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_FORTRAN_COMPILER}
					   --disable-java --with-mumps-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include/coin-or/mumps --enable-static=yes --enable-shared=no --with-lapack ${lapack_flags}
					       --with-mumps-lflags="-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lcoinmumps" --prefix=${THIRD_PARTY_DIR}/lion/thirdparty
                      --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib
	      DEPENDS mumps
	    )
    else()
	    ExternalProject_Add(ipopt
	      GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
	      GIT_TAG releases/3.14.3 
	      PREFIX "${THIRD_PARTY_DIR}/ipopt"
	      SOURCE_DIR ${THIRD_PARTY_DIR}/ipopt/source
	      BINARY_DIR ${THIRD_PARTY_DIR}/ipopt/build 
	      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
	      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/build &&
			../source/configure CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_FORTRAN_COMPILER}
					   --disable-java --with-mumps-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include/coin-or/mumps --enable-static=yes --enable-shared=no ${lapack_flags}
					       --with-mumps-lflags="-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lcoinmumps" --prefix=${THIRD_PARTY_DIR}/lion/thirdparty
	      DEPENDS mumps
	    )    
    endif()
endif()

