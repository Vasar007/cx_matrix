cmake_minimum_required(VERSION 3.8)

project(cx_matrix)

# Set compiler options.
if (MSVC)
    add_compile_options(/W4)
else ()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add cx_math library.
add_library(cx_math STATIC
    cx_math/cx_math.h
)
target_include_directories(cx_math PUBLIC
    cx_math
)
set_target_properties(cx_math PROPERTIES LINKER_LANGUAGE CXX)

# Add cx_options library.
add_library(cx_options STATIC
    cx_options/cx_algorthms.hpp
)
target_include_directories(cx_options PUBLIC
    cx_options
)
set_target_properties(cx_options PROPERTIES LINKER_LANGUAGE CXX)

# Include sub-projects.
add_subdirectory(cx_matrix)
