if (${BUILD_TINYXML})

include (ExternalProject)
ExternalProject_Add(tinyxml
  GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
  GIT_TAG master
  PREFIX "${THIRD_PARTY_DIR}/tinyxml"
  SOURCE_DIR "${THIRD_PARTY_DIR}/tinyxml/source" 
  BINARY_DIR "${THIRD_PARTY_DIR}/tinyxml/build" 
  INSTALL_DIR "${THIRD_PARTY_DIR}/lion/thirdparty"
  
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${THIRD_PARTY_DIR}/lion/thirdparty
                -DCMAKE_CXX_FLAGS:STRING=-DTIXML_USE_STL
		-DCMAKE_BUILD_TYPE:STRING=Release
		-DBUILD_SHARED_LIBS:BOOL=On
                ${THIRD_PARTY_DIR}/tinyxml/source
)

endif()
