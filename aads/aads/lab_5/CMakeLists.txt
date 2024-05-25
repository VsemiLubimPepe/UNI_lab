cmake_minimum_required(VERSION 3.5)

set(CMAKE_C_STADARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

project(create_matrix)

project(lab_5)

add_library(matrix STATIC matrix.c)

add_executable(bench bench.c)

target_link_libraries(bench matrix)

set(create_matrix_SOURCE create_matrix.c)

add_executable(create_matrix ${create_matrix_SOURCE})

set(ORDER "2")

add_custom_command(
    TARGET create_matrix
    POST_BUILD
    COMMAND $<TARGET_FILE:create_matrix> ${ORDER}    
)