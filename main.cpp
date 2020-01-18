#include <omp.h>
#include "matrix-reading/CSC.h"
#include "matrix-reading/VectorRead.h"
#include "MatrixSolver.h"

using namespace std;

static string MATRIX_FILE_PATH = "..\\data\\sparse_triangular_matrix.mtx";
static string DENSE_VECTOR_FILE_PATH = "..\\data\\dense_vector.mtx";
static string SPARSE_VECTOR_FILE_PATH = "..\\data\\sparse_vector.mtx";

int main() {
    auto t1 = omp_get_wtime();
    #pragma omp parallel default(none)
    {
        #pragma omp single
        { printf("Number of threads: %d\n", omp_get_num_threads()); }
//        printf(" Hello ");
//        printf(" World \n");
    }
    auto t2 = omp_get_wtime();
    printf("Time to initialise threads: %fs\n", t2-t1);

    t1 = omp_get_wtime();
    int* col = nullptr;
    int* row = nullptr;
    double* val = nullptr;
    int n = CSC(MATRIX_FILE_PATH, col, row, val);
    t2 = omp_get_wtime();
    printf("Time to read matrix: %fs\n", t2-t1);
    printf("First element: %d, %d, %f\n", col[1], row[col[1]], val[col[1]]);
    printf("Last element: %d, %d, %f\n", col[n], row[col[n]], val[col[n]]);

    int* ind = nullptr;
    double* b = nullptr;
    t1 = omp_get_wtime();
    readVector(DENSE_VECTOR_FILE_PATH, ind, b);
    t2 = omp_get_wtime();
    printf("Time to read dense vector: %fs\n", t2-t1);
    printf("Dense vector: %f, %f, %f, ... %f, %f\n", b[1], b[2], b[3], b[n-1], b[n]);
    readVector(SPARSE_VECTOR_FILE_PATH, ind, b);
    printf("Sparse vector: (%i, %f), (%i, %f), (%i, %f), (%i, %f), (%i, %f), (%i, %f)\n", ind[1], b[1], ind[2], b[2], ind[3], b[3], ind[4], b[4], ind[ind[0]-1], b[ind[0]-1], ind[ind[0]], b[ind[0]]);

    solve(n, col, row, val, b);
    delete[] col;
    delete[] row;
    delete[] val;
    delete[] ind;
    delete[] b;
    return 0;
}