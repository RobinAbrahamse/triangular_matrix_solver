#include <omp.h>
#include "matrix-reading/CSC.h"
#include "matrix-reading/VectorRead.h"
#include "matrix-operations/MatrixSolver.h"

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
    auto t = new int[1];
    analyseDependencies(n, col, row, t);
    t3 = omp_get_wtime() - t3;
    printf("Time to analyse dense vector: %fs\n", t3);

    auto t4 = omp_get_wtime();
    createLevelsets(n, col, row, t, t);
    t4 = omp_get_wtime() - t4;
    printf("Time to analyse levels for dense vector: %fs\n", t4);
    delete[] t;

    auto t5 = omp_get_wtime();
    solve(n, col, row, val, b);
    t5 = omp_get_wtime() - t5;
    printf("Time to solve for dense vector: %fs\n", t5);

    int valid = verifyMult(n, col, row, val, b);
    printf("Multiplication valid: %d\n", valid);

    double* b_copy = nullptr;
    readVector(vector_file_path, b_copy);
    valid = verifySolve(n, col, row, val, b, b_copy);
    printf("Solution valid: %d\n", valid);
    return 0;
}