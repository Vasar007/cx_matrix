set(LIBS_NAME, "cx_matrix_impl")

add_library(${LIBS_NAME}
    cx_matrix.hpp
)
target_link_libraries(${LIBS_NAME}
    cx_math
    cx_options
)
set_target_properties(${LIBS_NAME} PROPERTIES LINKER_LANGUAGE CXX)
