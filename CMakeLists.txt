cmake_minimum_required(VERSION 3.15)
project(triangular_matrix_solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")

include_directories(.)
include_directories(matrix-reading)

find_package(OpenMP REQUIRED)
add_executable (
        triangular_matrix_solver
        matrix-reading/CSC.cpp
        matrix-reading/CSC.h
        matrix-reading/VectorRead.cpp
        matrix-reading/VectorRead.h
        main.cpp
        MatrixSolver.cpp
        MatrixSolver.h
        matrix-reading/Utils.cpp 
        matrix-reading/Utils.h )

#target_link_libraries(triangular_matrix_solver PRIVATE "${OpenMP_CXX_FLAGS}")
#target_compile_options(triangular_matrix_solver PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(triangular_matrix_solver PRIVATE OpenMP::OpenMP_CXX)
