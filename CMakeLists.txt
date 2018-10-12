cmake_minimum_required(VERSION 3.8)

project(cx_matrix)

# Set compiler options.
if (MSVC)
    add_compile_options(/W4)
else ()
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options(-Wall -Wextra -pedantic -fconstexpr-steps=1000000000)
    else()
        add_compile_options(-Wall -Wextra -pedantic -fconstexpr-loop-limit=1000000000)
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if ("${TARGET_CPU}" STREQUAL "amd64")
    set(CMAKE_SIZEOF_VOID_P 8)

    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)
    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB32_PATHS FALSE)
elseif ("${TARGET_CPU}" STREQUAL "x86")
    set(CMAKE_SIZEOF_VOID_P 4)

    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS FALSE)
    set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB32_PATHS TRUE)

    if (GCC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -m32")
    endif()
else ()
    message (FATAL_ERROR "Unsupported CPU: ${TARGET_CPU}")
endif()

# Include sub-projects.
add_subdirectory(Source)
