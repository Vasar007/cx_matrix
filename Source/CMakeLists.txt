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
