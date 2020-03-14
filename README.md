triangular_matrix_solver

To run program:
 - Navigate to main folder
 - Execute script named "run" (this will compile the program and run it)
 - Set "-D" flag for compiling in debug mode
 - The script takes 2 input arguments: (1) the (matrix market) matrix file; (2) the (matrix market) vector file
 - If no input arguments are supplied, the script will download an example matrix+vector to 'data' directory and solve for them
   (example can be found here: https://sparse.tamu.edu/MM/Schenk_AFE/af_0_k101.tar.gz)

 example: ./run -D data/Matrix.mtx data/Vector.mtx
