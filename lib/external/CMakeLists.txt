add_subdirectory(t96)
add_subdirectory(ts05)

set(SRC
    external.c
    external.h
)

add_library(external ${SRC})
target_include_directories(external PUBLIC ../common)
target_include_directories(external PUBLIC common)
set_target_properties(external PROPERTIES C_STANDARD 17)
target_compile_options(external PRIVATE -Wall -Werror -O3 -fno-math-errno -fno-trapping-math -ffp-contract=off -march=native -finline-limit=300)
