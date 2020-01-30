#!/bin/bash
FILE=CMakeCache.txt
if [ ! -f "$FILE" ]; then
	cmake CMakeLists.txt
fi
make
./triangular_matrix_solver ./data/matrix.mtx ./data/vector.mtx

