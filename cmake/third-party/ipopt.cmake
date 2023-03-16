if (${BUILD_IPOPT})
    set(fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    set(fortran_dirs ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    list(TRANSFORM fortran_libs PREPEND "-l")
    list(TRANSFORM fortran_dirs PREPEND "-L")
    string (REPLACE ";" " " fortran_libs_str "${fortran_libs}")
    string (REPLACE ";" " " fortran_dirs_str "${fortran_dirs}")

    set(lapack_flags "--with-lapack-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib/ -llapack -lblas ${fortran_dirs_str} ${fortran_libs_str}")
    if (${BUILD_LAPACK})
	set(depends "lapack")
    endif()

    if (MSYS)
	set(disable_dependency_tracking "--disable-dependency-tracking")
    endif()

    if (NOT MSYS)
	set(ipopt_patch_command cd ${THIRD_PARTY_DIR}/ipopt/source && git apply ${PATCH_DIR}/ipopt.patch --reject)
    else()
        set(ipopt_patch_command "")
    endif()

    include (ExternalProject)
    ExternalProject_Add(GKLib
      GIT_REPOSITORY https://github.com/KarypisLab/GKlib.git
      PREFIX "${THIRD_PARTY_DIR}/gklib"
      SOURCE_DIR ${THIRD_PARTY_DIR}/gklib/source
      BINARY_DIR ${THIRD_PARTY_DIR}/gklib/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/gklib/source && make config prefix=${THIRD_PARTY_DIR}/lion/thirdparty cc=${CMAKE_C_COMPILER}
      BUILD_COMMAND cd ${THIRD_PARTY_DIR}/gklib/source && make && make install
      INSTALL_COMMAND ""
      UPDATE_COMMAND ""
    )

    ExternalProject_Add(metis
      GIT_REPOSITORY https://github.com/KarypisLab/METIS.git
      PREFIX "${THIRD_PARTY_DIR}/metis"
      SOURCE_DIR ${THIRD_PARTY_DIR}/metis/source
      BINARY_DIR ${THIRD_PARTY_DIR}/metis/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND ""
      BUILD_COMMAND 
        cd ${THIRD_PARTY_DIR}/metis/source 
        && make config prefix=${THIRD_PARTY_DIR}/lion/thirdparty gklib_path=${THIRD_PARTY_DIR}/lion/thirdparty cc=${CMAKE_C_COMPILER}
        && cd build && make && make install
      INSTALL_COMMAND ""
      UPDATE_COMMAND ""
      DEPENDS GKLib
    )

    if (WITH_HSL)
        ExternalProject_Add(hsl
          GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-HSL.git
          GIT_TAG stable/2.2
          PREFIX "${THIRD_PARTY_DIR}/hsl"
          SOURCE_DIR ${THIRD_PARTY_DIR}/hsl/source
          BINARY_DIR ${THIRD_PARTY_DIR}/hsl/build 
          INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
          CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/hsl/source 
    		&& cp ${WITH_HSL} . 
    		&& unzip coinhsl-2022.11.09.zip && mv coinhsl-2022.11.09 coinhsl 
    		&& ./configure 
    		--prefix=${THIRD_PARTY_DIR}/lion/thirdparty 
    		--enable-static=no --enable-shared=yes ${lapack_flags} ${disable_dependency_tracking} --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib 
    		--with-metis-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include
    		"--with-metis-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lmetis -lGKlib -lm"
          BUILD_COMMAND cd ${THIRD_PARTY_DIR}/hsl/source && make && make install
          INSTALL_COMMAND ""
          UPDATE_COMMAND ""
          DEPENDS ${depends} metis
        )

        set(hsl_c_flags "--with-hsl-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include/coin-or/hsl")
        set(hsl_l_flags "--with-hsl-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lcoinhsl")
        set(depends_hsl "hsl")
    else()
        set(hsl_c_flags "")
        set(hsl_l_flags "")
        set(depends_hsl "")

    endif()

    ExternalProject_Add(mumps
      GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
      GIT_TAG releases/3.0.2
      PREFIX "${THIRD_PARTY_DIR}/mumps"
      SOURCE_DIR ${THIRD_PARTY_DIR}/mumps/source
      BINARY_DIR ${THIRD_PARTY_DIR}/mumps/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source 
 	&& ./get.Mumps 
	&& ./configure 
		--prefix=${THIRD_PARTY_DIR}/lion/thirdparty 
		--enable-static=no --enable-shared=yes ${lapack_flags} ${disable_dependency_tracking} --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib
		--with-metis-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include
		"--with-metis-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lmetis -lGKlib -lm"
      BUILD_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source && make && make install
      INSTALL_COMMAND ""
      UPDATE_COMMAND ""
      DEPENDS ${depends} metis
    )

    ExternalProject_Add(ipopt
      GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
      GIT_TAG stable/3.14
      PREFIX "${THIRD_PARTY_DIR}/ipopt"
      SOURCE_DIR ${THIRD_PARTY_DIR}/ipopt/source
      BINARY_DIR ${THIRD_PARTY_DIR}/ipopt/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      PATCH_COMMAND "${ipopt_patch_command}"
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/build &&
	../source/configure CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_FORTRAN_COMPILER}
                 --disable-java --enable-static=yes --enable-shared=no --with-lapack ${lapack_flags}
                 --with-mumps-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include/coin-or/mumps 
                 --with-mumps-lflags="-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lcoinmumps" 
                 ${hsl_c_flags} ${hsl_l_flags}
                 --prefix=${THIRD_PARTY_DIR}/lion/thirdparty
                 --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib
                 --without-asl
                 DEPENDS mumps ${depends_hsl}
    )
endif()
