if (NOT nlohmann_json_FOUND)
    find_path(NLOHMANNJSON_INCLUDE_DIR nlohmann/json.hpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include NO_DEFAULT_PATH)
    
    if (NLOHMANNJSON_INCLUDE_DIR)
        set(nlohmann_json_FOUND YES)
    
        add_library(nlohmann_json INTERFACE)

	    target_include_directories(nlohmann_json SYSTEM INTERFACE
		  "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>"
		  "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

    else()
        set(nlohmann_json_FOUND NO)
    endif()
endif()
