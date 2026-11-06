# Locate a complete Intel oneMKL development installation. If requested and no
# installation is discoverable, provision a pinned binary release locally.
#
# This runs at *configure time* by design: lion builds third-party projects
# during configuration, and MKLConfig.cmake must already exist when the
# Autotools arguments for Ipopt/MUMPS are constructed. ExternalProject_Add()
# alone would run too late in lion's current build architecture.
#
# Neither the Intel installer nor its libraries are copied into lion's install
# tree. An installed lion package continues to require oneMKL at consumption.

set(LION_MKL_AUTO_INSTALL ON CACHE BOOL
    "Download/install oneMKL locally when WITH_MKL=ON and no installation is found")
set(LION_MKL_ACCEPT_EULA OFF CACHE BOOL
    "Explicitly accept the Intel oneMKL license for an automatic installation")
set(LION_MKL_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/lion/thirdparty/oneapi" CACHE PATH
    "Private oneAPI installation prefix (may be shared across build directories)")
set(LION_MKL_VERSION "2026.2.0.890" CACHE STRING
    "Pinned Intel oneMKL standalone installer version")
set(LION_MKL_INSTALLER_URL
    "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/97415287-271a-4d14-967f-72e2d9ceeb1e/intel-onemkl-2026.2.0.890_offline.sh"
    CACHE STRING "HTTPS URL of a pinned Intel oneMKL offline installer")
set(LION_MKL_INSTALLER "" CACHE FILEPATH
    "Use this pre-downloaded Intel oneMKL offline installer instead of downloading")
set(LION_MKL_INSTALLER_SHA384 "" CACHE STRING
    "Optional SHA384 hex digest for verifying an Intel oneMKL installer")
set(LION_MKL_DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/thirdparty/downloads" CACHE PATH
    "Directory for cached oneMKL installer downloads")

# Unlike setting MKL_DIR, this accepts a root such as
# /opt/intel/oneapi/mkl/latest or /other/prefix/mkl/2026.2.
set(LION_MKL_ROOT "" CACHE PATH "Use an existing oneMKL development root")

function(_lion_mkl_validate_root root out_root)
    set(_root "")
    if(NOT "${root}" STREQUAL "" AND EXISTS "${root}")
        # Resolve Intel's 'latest' symlink to the actual version, so changing
        # that symlink cannot silently switch an existing build to another MKL.
        get_filename_component(_candidate "${root}" REALPATH)
        if(EXISTS "${_candidate}/include/mkl.h"
           AND EXISTS "${_candidate}/lib/cmake/mkl/MKLConfig.cmake")
            if(MKL_LINK STREQUAL "static")
                set(_suffix "a")
            else()
                set(_suffix "so")
            endif()
            foreach(_libdir "${_candidate}/lib"
                            "${_candidate}/lib/intel64"
                            "${_candidate}/lib/intel64_lin")
                if(EXISTS "${_libdir}/libmkl_gf_lp64.${_suffix}"
                   AND EXISTS "${_libdir}/libmkl_${MKL_THREADING}.${_suffix}"
                   AND EXISTS "${_libdir}/libmkl_core.${_suffix}")
                    set(_root "${_candidate}")
                    break()
                endif()
            endforeach()
        endif()
    endif()
    set(${out_root} "${_root}" PARENT_SCOPE)
endfunction()

function(_lion_mkl_find_version_in_prefix prefix out_root)
    # Ignore 'latest': a previous, unrelated version can own that symlink.
    # Accept e.g. mkl/2026.2 and mkl/2026.2.0 for version 2026.2.0.890.
    set(_root "")
    if(NOT LION_MKL_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.")
        message(FATAL_ERROR "Invalid LION_MKL_VERSION: ${LION_MKL_VERSION}")
    endif()
    set(_release "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}")
    string(REPLACE "." "\\." _release_pattern "${_release}")
    file(GLOB _configs "${prefix}/mkl/*/lib/cmake/mkl/MKLConfig.cmake")
    foreach(_config IN LISTS _configs)
        get_filename_component(_config_dir "${_config}" DIRECTORY)
        get_filename_component(_candidate "${_config_dir}/../../.." ABSOLUTE)
        get_filename_component(_version_name "${_candidate}" NAME)
        if(_version_name MATCHES "^${_release_pattern}($|\\.)")
            _lion_mkl_validate_root("${_candidate}" _valid)
            if(_valid)
                set(_root "${_valid}")
            endif()
        endif()
    endforeach()
    set(${out_root} "${_root}" PARENT_SCOPE)
endfunction()

function(_lion_mkl_select_root root)
    set(MKL_DIR "${root}/lib/cmake/mkl" CACHE PATH
        "Intel oneMKL CMake package directory" FORCE)
    set(MKL_ROOT "${root}" PARENT_SCOPE)
    message(STATUS "lion: found Intel oneMKL development files: ${root}")
endfunction()

function(_lion_mkl_bootstrap)
    if(LION_MKL_ROOT)
        _lion_mkl_validate_root("${LION_MKL_ROOT}" _existing)
        if(NOT _existing)
            message(FATAL_ERROR
                "LION_MKL_ROOT=${LION_MKL_ROOT} is not a complete oneMKL development root "
                "(headers, CMake package or ${MKL_LINK}/${MKL_THREADING} libraries missing)")
        endif()
        _lion_mkl_select_root("${_existing}")
        set(MKL_ROOT "${_existing}" PARENT_SCOPE)
        return()
    endif()

    # Honor an explicitly supplied MKL_DIR, including a prior CMake cache.
    # A stale cache path is discarded, allowing automatic reprovisioning.
    if(MKL_DIR AND EXISTS "${MKL_DIR}/MKLConfig.cmake")
        get_filename_component(_candidate "${MKL_DIR}/../../.." ABSOLUTE)
        _lion_mkl_validate_root("${_candidate}" _existing)
        if(_existing)
            set(MKL_ROOT "${_existing}" PARENT_SCOPE)
            message(STATUS "lion: using Intel oneMKL via MKL_DIR=${MKL_DIR}")
            return()
        endif()
        message(FATAL_ERROR "MKL_DIR=${MKL_DIR} is missing oneMKL development headers or selected libraries")
    elseif(MKL_DIR AND NOT MKL_DIR MATCHES "-NOTFOUND$")
        message(STATUS "lion: invalid/stale MKL_DIR=${MKL_DIR}; searching again")
        unset(MKL_DIR CACHE)
    endif()

    foreach(_candidate "${MKL_ROOT}" "$ENV{MKLROOT}" "/opt/intel/oneapi/mkl/latest")
        _lion_mkl_validate_root("${_candidate}" _existing)
        if(_existing)
            _lion_mkl_select_root("${_existing}")
            set(MKL_ROOT "${_existing}" PARENT_SCOPE)
            return()
        endif()
    endforeach()

    # A locally provisioned release is preferred over an unrelated system
    # installation, and can be reused without accepting the EULA again.
    _lion_mkl_find_version_in_prefix("${LION_MKL_INSTALL_PREFIX}" _existing)
    if(_existing)
        _lion_mkl_select_root("${_existing}")
        set(MKL_ROOT "${_existing}" PARENT_SCOPE)
        return()
    endif()

    # Also honor CMAKE_PREFIX_PATH, toolchains, and other standard CMake hints.
    find_package(MKL CONFIG QUIET)
    if(MKL_FOUND AND TARGET MKL::MKL)
        message(STATUS "lion: using oneMKL found by CMake: ${MKL_DIR}")
        return()
    endif()

    if(NOT LION_MKL_AUTO_INSTALL)
        message(FATAL_ERROR
            "Intel oneMKL development files not found and LION_MKL_AUTO_INSTALL=OFF. "
            "Set LION_MKL_ROOT or MKL_DIR, or enable LION_MKL_AUTO_INSTALL.")
    endif()
    if(NOT LION_MKL_ACCEPT_EULA)
        message(FATAL_ERROR
            "Intel oneMKL was not found. Automatic installation requires explicit "
            "acceptance of Intel's license: set -DLION_MKL_ACCEPT_EULA=ON after "
            "reviewing https://www.intel.com/content/www/us/en/developer/articles/license/end-user-license-agreement.html "
            "or provide LION_MKL_ROOT/MKL_DIR to reuse an existing installation.")
    endif()
    if(NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64)$")
        message(FATAL_ERROR
            "The pinned Intel oneMKL Linux installer is x86-64 only (host: ${CMAKE_HOST_SYSTEM_PROCESSOR})")
    endif()
    if("${LION_MKL_INSTALL_PREFIX}" MATCHES "[ ;]" OR
       "${LION_MKL_DOWNLOAD_DIR}" MATCHES "[;]")
        message(FATAL_ERROR "oneMKL installation paths must not contain spaces or semicolons")
    endif()
    if(NOT LION_MKL_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$")
        message(FATAL_ERROR "LION_MKL_VERSION must be a pinned Intel version like 2026.2.0.890")
    endif()
    if(NOT LION_MKL_INSTALLER AND
       NOT LION_MKL_VERSION STREQUAL "2026.2.0.890" AND
       LION_MKL_INSTALLER_URL MATCHES "intel-onemkl-2026\\.2\\.0\\.890_offline\\.sh$")
        message(FATAL_ERROR
            "LION_MKL_VERSION changed but the default 2026.2 installer URL did not. "
            "Set LION_MKL_INSTALLER_URL to the matching pinned Intel release.")
    endif()

    file(MAKE_DIRECTORY "${LION_MKL_INSTALL_PREFIX}")
    # Safe for multiple independent CMake builds using the same install prefix.
    # GUARD FUNCTION releases the lock as soon as provisioning is complete.
    file(LOCK "${LION_MKL_INSTALL_PREFIX}/.lion-mkl-install.lock"
        GUARD FUNCTION TIMEOUT 3600 RESULT_VARIABLE _lock_result)
    if(NOT _lock_result STREQUAL "0")
        message(FATAL_ERROR "Could not lock oneMKL install prefix: ${_lock_result}")
    endif()

    # Recheck under the lock; a concurrent CMake process may have installed MKL.
    _lion_mkl_find_version_in_prefix("${LION_MKL_INSTALL_PREFIX}" _existing)
    if(_existing)
        _lion_mkl_select_root("${_existing}")
        set(MKL_ROOT "${_existing}" PARENT_SCOPE)
        return()
    endif()

    if(LION_MKL_INSTALLER)
        if(NOT EXISTS "${LION_MKL_INSTALLER}")
            message(FATAL_ERROR "LION_MKL_INSTALLER does not exist: ${LION_MKL_INSTALLER}")
        endif()
        set(_installer "${LION_MKL_INSTALLER}")
    else()
        if(NOT LION_MKL_INSTALLER_URL MATCHES "^https://")
            message(FATAL_ERROR
                "LION_MKL_INSTALLER_URL must use HTTPS; use LION_MKL_INSTALLER for a local file")
        endif()
        file(MAKE_DIRECTORY "${LION_MKL_DOWNLOAD_DIR}")
        string(MD5 _url_id "${LION_MKL_INSTALLER_URL}")
        set(_installer "${LION_MKL_DOWNLOAD_DIR}/intel-onemkl-${LION_MKL_VERSION}-${_url_id}_offline.sh")
        if(NOT EXISTS "${_installer}")
            message(STATUS "lion: downloading pinned Intel oneMKL ${LION_MKL_VERSION} from Intel")
            file(DOWNLOAD "${LION_MKL_INSTALLER_URL}" "${_installer}.part"
                SHOW_PROGRESS TLS_VERIFY ON STATUS _download_status)
            list(GET _download_status 0 _download_code)
            if(NOT _download_code EQUAL 0)
                file(REMOVE "${_installer}.part")
                message(FATAL_ERROR "Intel oneMKL download failed: ${_download_status}")
            endif()
            file(RENAME "${_installer}.part" "${_installer}")
        endif()
    endif()

    # Intel currently does not publish an accessible checksum alongside the
    # standalone installer URL. SHA384 is therefore optional; HTTPS is always
    # enforced for downloads and users may supply the checksum for CI pinning.
    if(LION_MKL_INSTALLER_SHA384)
        string(LENGTH "${LION_MKL_INSTALLER_SHA384}" _hash_length)
        if(NOT _hash_length EQUAL 96 OR
           NOT LION_MKL_INSTALLER_SHA384 MATCHES "^[0123456789abcdefABCDEF]+$")
            message(FATAL_ERROR "LION_MKL_INSTALLER_SHA384 must contain 96 hexadecimal characters")
        endif()
        file(SHA384 "${_installer}" _actual_sha384)
        string(TOLOWER "${LION_MKL_INSTALLER_SHA384}" _expected_sha384)
        if(NOT _actual_sha384 STREQUAL _expected_sha384)
            message(FATAL_ERROR
                "Intel oneMKL installer SHA384 mismatch: ${_installer}. "
                "Delete the cached installer before retrying.")
        endif()
    endif()

    # Compatible with lion's CMake 3.10 baseline (file(SIZE) needs 3.14).
    file(READ "${_installer}" _installer_first_byte LIMIT 1)
    if(_installer_first_byte STREQUAL "")
        message(FATAL_ERROR "Intel oneMKL installer is empty: ${_installer}")
    endif()

    # 'instance' is unique for the target prefix, avoiding interference with
    # other per-user Intel oneAPI installations on the same Linux machine.
    string(MD5 _instance_hash "${LION_MKL_INSTALL_PREFIX}")
    set(_instance "lionmkl${_instance_hash}")
    set(_log "${LION_MKL_INSTALL_PREFIX}/lion-mkl-installer.log")
    message(STATUS "lion: installing Intel oneMKL ${LION_MKL_VERSION} into ${LION_MKL_INSTALL_PREFIX}")
    execute_process(
        COMMAND sh "${_installer}" -a --silent --eula accept
            "--install-dir=${LION_MKL_INSTALL_PREFIX}" "--instance=${_instance}"
        RESULT_VARIABLE _install_result
        OUTPUT_FILE "${_log}"
        ERROR_FILE "${_log}")
    if(NOT _install_result STREQUAL "0")
        message(FATAL_ERROR
            "Intel oneMKL installation failed (${_install_result}). See ${_log}")
    endif()
    _lion_mkl_find_version_in_prefix("${LION_MKL_INSTALL_PREFIX}" _installed)
    if(NOT _installed)
        message(FATAL_ERROR
            "Intel oneMKL installer completed but development headers, "
            "MKLConfig.cmake or selected libraries are missing under ${LION_MKL_INSTALL_PREFIX}/mkl. "
            "See ${_log}")
    endif()
    _lion_mkl_select_root("${_installed}")
    set(MKL_ROOT "${_installed}" PARENT_SCOPE)
endfunction()

_lion_mkl_bootstrap()
