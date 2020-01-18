cmake_minimum_required(VERSION 3.15)
project(triangular_matrix_solver)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
include_directories(matrix-reading)

find_package(OpenMP REQUIRED)
add_executable(triangular_matrix_solver
        .idea/misc.xml
        .idea/modules.xml
        .idea/triangular_matrix_solver.iml
        .idea/vcs.xml
        .idea/workspace.xml
        cmake-build-debug/cmake_install.cmake
        cmake-build-debug/CMakeCache.txt
        cmake-build-debug/Makefile
        cmake-build-debug/triangular_matrix_solver.cbp
        data/dense_vector.mtx
        data/sparse_triangular_matrix.mtx
        data/sparse_vector.mtx
        matrix-reading/CSC.cpp
        matrix-reading/CSC.h
        matrix-reading/VectorRead.cpp
        matrix-reading/VectorRead.h
        .gitignore
        CMakeLists.txt
        main.cpp
        MatrixSolver.cpp
        MatrixSolver.h
        README.md)
target_compile_options(triangular_matrix_solver PRIVATE -Wall ${OpenMp_CXX_FLAGS})
target_link_libraries(triangular_matrix_solver PRIVATE ${OpenMp_CXX_FLAGS})
