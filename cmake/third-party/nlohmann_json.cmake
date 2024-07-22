if (${BUILD_NLOHMANNJSON})

include (ExternalProject)
ExternalProject_Add(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git 
  GIT_TAG master
  PREFIX "${THIRD_PARTY_DIR}/nlohmann_json"
  SOURCE_DIR "${THIRD_PARTY_DIR}/nlohmann_json/source" 
  BINARY_DIR "${THIRD_PARTY_DIR}/nlohmann_json/build" 
  INSTALL_DIR "${THIRD_PARTY_DIR}/lion/thirdparty"
  
  CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -DJSON_BuildTests=OFF
                -DBUILD_TESTING=Off
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_INSTALL_PREFIX:PATH=${THIRD_PARTY_DIR}/lion/thirdparty
		-DCMAKE_BUILD_TYPE:STRING=Release
		-DCMAKE_CXX_STANDARD=11
		-DCMAKE_CXX_FLAGS=${windows_flag}
                ${THIRD_PARTY_DIR}/nlohmann_json/source
)

endif()
