#!/bin/bash
FILE=CMakeCache.txt
if [ ! -f "$FILE" ]; then
	cmake CMakeLists.txt
	echo
fi
if [ ! -d data ]; then
	mkdir data
fi
if [ ! -f data/matrix.mtx ] && [ ! -f data/vector.mtx ]; then
	wget -O data/af_0_k101.tar.gz https://sparse.tamu.edu/MM/Schenk_AFE/af_0_k101.tar.gz
	tar -xvzf data/af_0_k101.tar.gz -C data
	rm -r data/af_0_k101.tar.gz
	mv data/af_0_k101/af_0_k101.mtx data/matrix.mtx
	mv data/af_0_k101/af_0_k101_b.mtx data/vector.mtx
	rmdir data/af_0_k101
	echo
fi
make
echo
echo
./triangular_matrix_solver ./data/matrix.mtx ./data/vector.mtx

