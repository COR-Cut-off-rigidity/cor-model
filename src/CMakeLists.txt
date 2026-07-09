set(SRC
    args.c
    args.h
    field.c
    field.h
    io.c
    io.h
    main.c
    models.c
    models.h
    trajectory.c
    trajectory.h
    utils.c
    utils.h
)

add_executable(${PROJECT_NAME} ${SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES C_STANDARD 17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -O3 -fno-math-errno -fno-trapping-math -ffp-contract=off -march=native -finline-limit=300)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/lib/internal)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/lib/external)
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/lib/crustal)

target_link_libraries(${PROJECT_NAME} PRIVATE internal)
target_link_libraries(${PROJECT_NAME} PRIVATE external)
target_link_libraries(${PROJECT_NAME} PRIVATE ts05)
target_link_libraries(${PROJECT_NAME} PRIVATE t96)
target_link_libraries(${PROJECT_NAME} PRIVATE crustal)

if(DEFINED ENV{COR_ENGINE_AUTHOR_EMAIL})
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEF_COR_ENGINE_AUTHOR_EMAIL="$ENV{COR_ENGINE_AUTHOR_EMAIL}")
endif()
