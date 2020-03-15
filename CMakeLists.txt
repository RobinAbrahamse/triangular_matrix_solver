project(triangular_matrix_solver)

cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(.)
include_directories(matrix-reading)

add_executable (
        triangular_matrix_solver
        matrix-reading/CSC.cpp
        matrix-reading/CSC.h
        matrix-reading/VectorRead.cpp
        matrix-reading/VectorRead.h
        main.cpp
        matrix-operations/MatrixSolver.cpp
        matrix-operations/MatrixSolver.h
        matrix-reading/Utils.cpp
        matrix-reading/Utils.h )
