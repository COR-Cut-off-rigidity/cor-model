cmake_minimum_required(VERSION 3.10)
project(cor VERSION 3.2.0 LANGUAGES C)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

add_subdirectory(lib)
add_subdirectory(src)
