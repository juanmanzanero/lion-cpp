file(GLOB_RECURSE math_src "*.cpp")

add_library(math STATIC ${math_src})
