cmake_minimum_required(VERSION 3.8)

project(cx_matrix_impl)

add_library(${PROJECT_NAME}
    cx_matrix.hpp
)
target_link_libraries(${PROJECT_NAME}
    cx_math
    cx_options
)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
