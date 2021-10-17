if (NOT matplotlibcpp_FOUND)
    find_path(MATPLOTLIBCPP_INCLUDE_DIR matplotlibcpp.h PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include)
 
    if (MATPLOTLIBCPP_INCLUDE_DIR)
        set(matplotlibcpp_FOUND YES)
    
    	target_include_directories(python INTERFACE ${MATPLOTLIBCPP_INCLUDE_DIR})
    	target_compile_definitions(python INTERFACE _MATPLOTLIB_CPP_)
    else()
        set(matplotlibcpp_FOUND NO)
    endif()
endif()
