cmake_minimum_required(VERSION 3.10)
project(chaos C)
add_library(${PROJECT_NAME} chaos.c)
set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -O3 -fno-math-errno -fno-trapping-math -ffp-contract=off -march=native -finline-limit=300) # -g

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
