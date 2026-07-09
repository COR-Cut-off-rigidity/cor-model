cmake_minimum_required(VERSION 3.10)
project(ts05 C)

set(SRC
    ts05.c
    ts05.h
    deformed.c
    deformed.h
    birk_tot.c
    birk_tot.h
    full_rc.c
    full_rc.h
)

add_library(${PROJECT_NAME} ${SRC})
target_include_directories(${PROJECT_NAME} PUBLIC ../common)
set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -O3 -fno-math-errno -fno-trapping-math -ffp-contract=off -march=native -finline-limit=300)
target_link_libraries(${PROJECT_NAME} PRIVATE m)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
