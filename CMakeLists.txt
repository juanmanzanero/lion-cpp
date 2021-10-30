cmake_minimum_required(VERSION 3.10)

# set the project name
project(lioncpp VERSION 1.0 LANGUAGES C CXX Fortran)

include(CMakeDependentOption)
include(GNUInstallDirs)

# Options
option(ENABLE_TEST "Enables tests" ON)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/thirdparty)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/thirdparty/include)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lion/thirdparty)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lion/thirdparty/include)
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Download third party libraries
include(get-third-party)

# Configure doxygen
include(Doxygen)

# Set C++ 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(./)

include(compilerflags)

if ( ${ENABLE_TEST} )
	# Google tests
	enable_testing()
	list(APPEND CMAKE_CTEST_ARGUMENTS "--verbose")

	# Compile the subdirectories
	add_subdirectory(./test)
endif()

############# INSTALLATION #################

# Create a lion interface to be used as front end
add_library(lion INTERFACE)
target_include_directories(lion SYSTEM INTERFACE
	"$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
	"$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/lion/thirdparty/include>"
	"$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>"
	"$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>"
	"$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}/lion/thirdparty/include>"
)

target_link_libraries(lion INTERFACE ipopt tinyxml2)

# Install lion
install(DIRECTORY "${CMAKE_SOURCE_DIR}/lion" TYPE INCLUDE)

# Install third party libraries
install(DIRECTORY "${CMAKE_BINARY_DIR}/lion/thirdparty/include/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/lion/thirdparty/include/")


include(CMakePackageConfigHelpers)

# Set names
set(cmake_package_name lion CACHE INTERNAL "")
set(targets_export_name ${cmake_package_name}Targets CACHE INTERNAL "")
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated" CACHE INTERNAL "")
set(cmake_files_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/${cmake_package_name}")
set(version_file "${generated_dir}/${cmake_package_name}ConfigVersion.cmake")

# Create the package
install(EXPORT ${targets_export_name}
  NAMESPACE ${cmake_package_name}::
  DESTINATION ${cmake_files_install_dir})

# Create the config file (the one to be included for its use)
set(config_file "${generated_dir}/${cmake_package_name}Config.cmake")
configure_package_config_file("${CMAKE_SOURCE_DIR}/cmake/Config.cmake.in"
  "${config_file}" INSTALL_DESTINATION ${cmake_files_install_dir})
write_basic_package_version_file(${version_file} VERSION 1.0 COMPATIBILITY AnyNewerVersion)
install(FILES ${version_file} ${config_file}
  DESTINATION ${cmake_files_install_dir})

# Export the targets
install(TARGETS lion ipopt tinyxml2
  EXPORT ${targets_export_name}
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
