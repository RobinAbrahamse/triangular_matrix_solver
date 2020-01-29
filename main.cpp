#include <omp.h>
#include "matrix-reading/CSC.h"
#include "matrix-reading/VectorRead.h"
#include "MatrixSolver.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) return 1;
    const auto matrix_file_path = argv[1];
    const auto vector_file_path = argv[2];
    int* col = nullptr;
    int* row = nullptr;
    double* val = nullptr;
    auto t1 = omp_get_wtime();
    int n = CSC(matrix_file_path, col, row, val);
    t1 = omp_get_wtime() - t1;
    printf("Time to read matrix: %fs\n", t1);

    int* ind = nullptr;
    double* b = nullptr;
    bool denseVector = isDense(vector_file_path);
    auto t2 = omp_get_wtime();
    readVector(vector_file_path, b);
    t2 = omp_get_wtime() - t2;
    printf("Time to read dense vector: %fs\n", t2);

    auto t3 = omp_get_wtime();
    analyse(n, col, row);
    t3 = omp_get_wtime() - t3;
    printf("Time to analyse dense vector: %fs\n", t3);

    auto t4 = omp_get_wtime();
    solve(n, col, row, val, b);
    t4 = omp_get_wtime() - t4;
    printf("Time to solve for dense vector: %fs\n", t4);

    double* b_copy = nullptr;
    readVector(vector_file_path, b_copy);
    int valid = verify(n, col, row, val, b, b_copy);
    printf("Solution valid: %d\n", valid);
    return 0;
}