cmake_minimum_required(VERSION 3.10)
project(cor C)

add_subdirectory(lib/internal)
add_subdirectory(lib/ts05)
add_subdirectory(lib/t96)
add_subdirectory(lib/chaos)

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()

file(GLOB MMUL_SRCS CONFIGURE_DEPENDS "src/*.h" "src/*.c")
add_executable(${PROJECT_NAME} ${MMUL_SRCS})
set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -O3 -fno-math-errno -fno-trapping-math -ffp-contract=off -march=native -finline-limit=300) #-g
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/lib)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

target_link_libraries(${PROJECT_NAME} PRIVATE m PRIVATE internal)
target_link_libraries(${PROJECT_NAME} PRIVATE m PRIVATE ts05)
target_link_libraries(${PROJECT_NAME} PRIVATE m PRIVATE t96)
target_link_libraries(${PROJECT_NAME} PRIVATE m PRIVATE chaos)
