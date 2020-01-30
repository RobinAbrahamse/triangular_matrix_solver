#!/bin/bash
FILE=CMakeCache.txt
if [ ! -f "$FILE" ]; then
	cmake CMakeLists.txt
fi
make
./triangular_matrix_solver ./data/sparse_triangular_matrix.mtx ./data/dense_vector.mtx

