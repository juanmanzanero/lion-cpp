if (NOT loggercpp_FOUND)
    find_path(LOGGERCPP_INCLUDE_DIR logger.h PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include NO_DEFAULT_PATH)
 
    if (LOGGERCPP_INCLUDE_DIR)
        set(loggercpp_FOUND YES)
    
        add_library(loggercpp INTERFACE)
    	target_include_directories(loggercpp INTERFACE ${CMAKE_BINARY_DIR})
    else()
        set(loggercpp_FOUND NO)
    endif()
endif()
