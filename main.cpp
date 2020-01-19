#include <omp.h>
#include "matrix-reading/CSC.h"
#include "matrix-reading/VectorRead.h"
#include "MatrixSolver.h"

using namespace std;

static string MATRIX_FILE_PATH = "..\\data\\sparse_triangular_matrix.mtx";
static string DENSE_VECTOR_FILE_PATH = "..\\data\\dense_vector.mtx";
//static string MATRIX_FILE_PATH = "..\\data\\matrix2.mtx";
//static string DENSE_VECTOR_FILE_PATH = "..\\data\\dense_vector2.mtx";
//static string MATRIX_FILE_PATH = "..\\data\\matrix3.mtx";
//static string DENSE_VECTOR_FILE_PATH = "..\\data\\dense_vector3.mtx";
static string SPARSE_VECTOR_FILE_PATH = "..\\data\\sparse_vector.mtx";

int main() {
    int* col = nullptr;
    int* row = nullptr;
    double* val = nullptr;
    auto t1 = omp_get_wtime();
    int n = CSC(MATRIX_FILE_PATH, col, row, val);
    auto t2 = omp_get_wtime();
    printf("Time to read matrix: %fs\n", t2-t1);
//    printf("First element: %d, %d, %f\n", col[0], row[col[0]], val[col[0]]);
//    printf("Last element: %d, %d, %f\n", col[n-1], row[col[n-1]], val[col[n-1]]);

    int* ind = nullptr;
    double* b = nullptr;
    t1 = omp_get_wtime();
    readVector(DENSE_VECTOR_FILE_PATH, ind, b);
    t2 = omp_get_wtime();
    printf("Time to read dense vector: %fs\n", t2-t1);
//    printf("Dense vector: %f, %f, %f, ... %f, %f\n", b[0], b[1], b[2], b[n-2], b[n-1]);
//    readVector(SPARSE_VECTOR_FILE_PATH, ind, b);
//    printf("Sparse vector: (%i, %f), (%i, %f), (%i, %f), (%i, %f), (%i, %f), (%i, %f)\n", ind[1], b[1], ind[2], b[2], ind[3], b[3], ind[4], b[4], ind[ind[0]-1], b[ind[0]-1], ind[ind[0]], b[ind[0]]);

    t1 = omp_get_wtime();
    solve(n, col, row, val, b);
    t2 = omp_get_wtime();
    printf("Time to solve for dense vector: %fs\n", t2-t1);

    double* b_copy = nullptr;
    readVector(DENSE_VECTOR_FILE_PATH, ind, b_copy);
    int success = verify(n, col, row, val, b, b_copy);
    printf("Solution valid: %d\n", success);
    return 0;
}