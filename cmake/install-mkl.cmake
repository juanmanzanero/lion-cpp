# Install the full oneMKL development tree next to lion, so an installed lion
# package does not depend on the CMake build directory or the original MKL root.
# This is deliberately an install-time copy: configuring/building lion never
# needs write access to /opt or another protected installation prefix.
option(LION_INSTALL_MKL_WITH_LION
    "Install the oneMKL development tree inside the lion installation prefix" ON)

set(LION_MKL_INSTALL_SUBDIR "")
set(LION_MKL_EXTERNAL_ROOT "")
if(WITH_MKL)
    if(LION_INSTALL_MKL_WITH_LION)
        get_filename_component(_lion_mkl_source_root "${MKL_ROOT}" REALPATH)
        get_filename_component(_lion_mkl_version "${_lion_mkl_source_root}" NAME)
        if(NOT EXISTS "${_lion_mkl_source_root}/include/mkl.h" OR
           NOT EXISTS "${_lion_mkl_source_root}/lib/cmake/mkl/MKLConfig.cmake")
            message(FATAL_ERROR "Cannot bundle incomplete oneMKL development tree: ${_lion_mkl_source_root}")
        endif()
        if(_lion_mkl_version MATCHES "[/\\\\; ]" OR _lion_mkl_version STREQUAL "latest")
            message(FATAL_ERROR "Expected a versioned oneMKL root, got: ${_lion_mkl_source_root}")
        endif()

        set(LION_MKL_INSTALL_SUBDIR "oneapi/mkl/${_lion_mkl_version}")
        install(DIRECTORY "${_lion_mkl_source_root}/"
            DESTINATION "${LION_MKL_INSTALL_SUBDIR}"
            USE_SOURCE_PERMISSIONS)

        if(MKL_LINK STREQUAL "dynamic")
            # Relative path between lion's installed shared libraries and the
            # installed oneMKL libraries, independent of --prefix/DESTDIR.
            file(RELATIVE_PATH _lion_mkl_library_subdir
                "${_lion_mkl_source_root}" "${MKL_LIBDIR}")
            if(_lion_mkl_library_subdir MATCHES "^\\.\\.")
                message(FATAL_ERROR "MKL_LIBDIR is outside MKL_ROOT: ${MKL_LIBDIR}")
            endif()
            file(RELATIVE_PATH _lion_mkl_relative_libdir
                "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
                "${CMAKE_INSTALL_PREFIX}/${LION_MKL_INSTALL_SUBDIR}/${_lion_mkl_library_subdir}")
            set(LION_MKL_NEW_RPATH "\$ORIGIN/${_lion_mkl_relative_libdir}")
            set(LION_MKL_OLD_RPATH "${MKL_LIBDIR}")
            set(LION_MKL_INSTALLED_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
            set(_lion_mkl_rpath_script "${CMAKE_CURRENT_BINARY_DIR}/generated/lion-relocate-mkl-rpath.cmake")
            file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/generated")
            configure_file("${CMAKE_CURRENT_LIST_DIR}/mkl-relocate-rpath.cmake.in"
                "${_lion_mkl_rpath_script}" @ONLY)
            install(SCRIPT "${_lion_mkl_rpath_script}")
        endif()
        message(STATUS "lion: oneMKL will be installed at <prefix>/${LION_MKL_INSTALL_SUBDIR}")
    else()
        # Preserve the existing external-MKL behavior when explicitly chosen.
        set(LION_MKL_EXTERNAL_ROOT "${MKL_ROOT}")
        message(STATUS "lion: oneMKL will remain external at ${MKL_ROOT}")
    endif()
endif()
