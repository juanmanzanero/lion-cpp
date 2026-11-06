if (${BUILD_IPOPT})
    include (ExternalProject)

    set(fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    set(fortran_dirs ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    list(TRANSFORM fortran_libs PREPEND "-l")
    list(TRANSFORM fortran_dirs PREPEND "-L")
    string (REPLACE ";" " " fortran_libs_str "${fortran_libs}")
    string (REPLACE ";" " " fortran_dirs_str "${fortran_dirs}")
    
    # Shared HSL/MUMPS need a persistent RUNPATH to the locally installed MKL.
    # Ipopt's standalone HSL link probe needs rpath-link to resolve MKL DSOs
    # that are needed indirectly by libcoinhsl.so (its own LAPACK flags are not
    # included in that probe).
    # This file is included by lion's separate third-party CMake process.
    # Fail immediately if the parent did not forward the MKL settings.
    if (WITH_MKL)
        if (NOT MKL_LINK STREQUAL "dynamic" AND NOT MKL_LINK STREQUAL "static")
            message(FATAL_ERROR "lion: missing or invalid MKL_LINK in nested third-party build")
        endif()
        if (NOT IS_DIRECTORY "${MKL_LIBDIR}")
            message(FATAL_ERROR "lion: missing MKL_LIBDIR in nested third-party build: '${MKL_LIBDIR}'")
        endif()
    endif()
    set(lion_autotools_mkl_ldflags "")
    if (WITH_MKL AND MKL_LINK STREQUAL "dynamic")
        set(lion_autotools_mkl_ldflags "LDFLAGS=$ENV{LDFLAGS} -Wl,-rpath-link,${MKL_LIBDIR} -Wl,-rpath,${MKL_LIBDIR}")
    endif()
    set(lion_mumps_lflags "--with-mumps-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lcoinmumps")

    set(mumps_configure_command ./configure)
    set(pardisomkl_flag "")
    if (WITH_MKL)
      set(lapack_flags "--with-lapack-lflags=${MKL_FLAGS}")
      set(MUMPS_CPPFLAGS "$ENV{CPPFLAGS} -I${MKL_INCDIR}")
      set(mumps_configure_command "${CMAKE_COMMAND}" -E env "CPPFLAGS=${MUMPS_CPPFLAGS}" ./configure)
      set(pardisomkl_flag "--enable-pardisomkl")
      message(STATUS "${lapack_flags}")
	
    else()
        set(lapack_flags "--with-lapack-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib/ -llapack -lblas ${fortran_dirs_str} ${fortran_libs_str}")
        if (${BUILD_LAPACK})
            set(depends "lapack")
        endif()
    endif()


    if (WITH_HSL)
        ExternalProject_Add(GKLib
          GIT_REPOSITORY https://github.com/KarypisLab/GKlib.git
          PREFIX "${THIRD_PARTY_DIR}/gklib"
          SOURCE_DIR ${THIRD_PARTY_DIR}/gklib/source
          BINARY_DIR ${THIRD_PARTY_DIR}/gklib/build 
          INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
	  CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/gklib/source && make config prefix=${THIRD_PARTY_DIR}/lion/thirdparty cc=${CMAKE_C_COMPILER} CFLAGS=-fPIC
          BUILD_COMMAND cd ${THIRD_PARTY_DIR}/gklib/source && make && make install
          INSTALL_COMMAND ""
          UPDATE_COMMAND ""
        )

   
        ExternalProject_Add(metis
          GIT_REPOSITORY https://github.com/KarypisLab/METIS.git
	  GIT_TAG 076b60ea8f90b3e213db22ad7d8f706b62c687b5
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

	set(metis_c_flags "--with-metis-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include")
	set(metis_l_flags "--with-metis-lflags=-L${THIRD_PARTY_DIR}/lion/thirdparty/lib -lmetis -lGKlib -lm")

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
            ${lion_autotools_mkl_ldflags}
            --prefix=${THIRD_PARTY_DIR}/lion/thirdparty 
            --enable-static=no --enable-shared=yes --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib 
	    ${lapack_flags} ${metis_c_flags} ${metis_l_flags}
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
	set(metis_c_flags "--without-metis")
	set(metis_l_flags "")
    endif()

    ExternalProject_Add(mumps
      GIT_REPOSITORY https://github.com/coin-or-tools/ThirdParty-Mumps.git
      GIT_TAG stable/3.0
      PREFIX "${THIRD_PARTY_DIR}/mumps"
      SOURCE_DIR ${THIRD_PARTY_DIR}/mumps/source
      BINARY_DIR ${THIRD_PARTY_DIR}/mumps/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source 
         && ./get.Mumps 
         && ${mumps_configure_command}
           ${lion_autotools_mkl_ldflags}
           --prefix=${THIRD_PARTY_DIR}/lion/thirdparty 
	   --enable-static=no --enable-shared=yes --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib
	   --without-metis
	   ${lapack_flags}
      BUILD_COMMAND cd ${THIRD_PARTY_DIR}/mumps/source && make && make install
      INSTALL_COMMAND ""
      UPDATE_COMMAND ""
      DEPENDS ${depends}
    )


    ExternalProject_Add(ipopt
      GIT_REPOSITORY https://github.com/coin-or/Ipopt.git
      GIT_TAG stable/3.14
      PREFIX "${THIRD_PARTY_DIR}/ipopt"
      SOURCE_DIR ${THIRD_PARTY_DIR}/ipopt/source
      BINARY_DIR ${THIRD_PARTY_DIR}/ipopt/build 
      INSTALL_DIR ${THIRD_PARTY_DIR}/lion/thirdparty
      PATCH_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/source && git apply ${PATCH_DIR}/ipopt.patch --reject
      CONFIGURE_COMMAND cd ${THIRD_PARTY_DIR}/ipopt/build &&
      ../source/configure CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER} F77=${CMAKE_Fortran_COMPILER} CFLAGS=-O3\ -fno-fast-math\ -m64 CXXFLAGS=-O3\ -fno-fast-math\ -m64
          --disable-java --enable-static=yes --enable-shared=no --with-lapack ${lapack_flags}
          ${lion_autotools_mkl_ldflags}
          --with-mumps-cflags=-I${THIRD_PARTY_DIR}/lion/thirdparty/include/coin-or/mumps 
          ${lion_mumps_lflags}
          ${hsl_c_flags} ${hsl_l_flags}
          --prefix=${THIRD_PARTY_DIR}/lion/thirdparty
          --libdir=${THIRD_PARTY_DIR}/lion/thirdparty/lib
          --without-asl
          ${pardisomkl_flag}
          DEPENDS mumps ${depends_hsl}
    )

endif()
