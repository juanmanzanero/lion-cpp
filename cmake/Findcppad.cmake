if (NOT cppad_FOUND)
    find_path(CPPAD_INCLUDE_DIR cppad.hpp PATHS ${CMAKE_BINARY_DIR}/lion/thirdparty/include/cppad)
 
    if (CPPAD_INCLUDE_DIR)
        set(cppad_FOUND YES)
    else()
        set(cppad_FOUND NO)
    endif()
endif()
