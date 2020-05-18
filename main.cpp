#include <omp.h>
#include <cstring>
#include <thread>
#include "matrix-reading/CSC.h"
#include "matrix-reading/VectorRead.h"
#include "matrix-operations/MatrixSolver.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) return 1;
    printf("Number of cores: %d\n", std::thread::hardware_concurrency());
    const auto matrix_file_path = argv[1];
    const auto vector_file_path = argv[2];
    int* col = nullptr;
    int* row = nullptr;
    double* val = nullptr;
    auto t = omp_get_wtime();
    const int n = CSC(matrix_file_path, col, row, val);
    t = omp_get_wtime() - t;
    printf("Time to read matrix: %fs\n", t);

    // // int* ind = nullptr;
    // double* b = nullptr;
    // // bool denseVector = isDense(vector_file_path);
    // t = omp_get_wtime();
    // readVector(vector_file_path, b);
    // t = omp_get_wtime() - t;
    // printf("Time to read vector: %fs\n", t);
    auto b = new double[n];
    for (int i = 0; i < n; i++) {
      b[i] = 1;
    }

    int v = 3;
    int tn = 5;
    double*** bt = new double**[v];
    for (int i = 0; i < v; i++){
      bt[i] = new double*[tn];
      for(int j = 0; j < tn; ++j) {
        bt[i][j] = new double[n];
        memcpy(bt[i][j], b, n*8);
        memcpy(bt[i][j], b, n*8);
      }
    }

    //    solveParallelLevels(n, col, row, val, b);
    //    solveParallelTasks(n, col, row, val, b);
    //    solveParallelFor(n, col, row, val, b);
    //    solveSerial(n, col, row, val, b);

    t = omp_get_wtime();
    for (int i = 0; i < tn; i++) {
      solveParallelLevels(n, col, row, val, bt[0][i]);
    }
    t = (omp_get_wtime() - t) / tn;
    printf("Average time to solve (parallel level sets) for dense vector: %fs\n", t);

    // t = omp_get_wtime();
    // for (int i = 0; i < tn; i++) {
    //   solveParallelTasks(n, col, row, val, bt[1][i]);
    // }
    // t = (omp_get_wtime() - t) / tn;
    // printf("Average time to solve (parallel sync-free) for dense vector: %fs\n", t);

    t = omp_get_wtime();
    for (int i = 0; i < tn; i++) {
      solveSerial(n, col, row, val, bt[2][i]);
    }
    t = (omp_get_wtime() - t) / tn;
    printf("Average time to solve (serial) for dense vector: %fs\n", t);

    int valid = verifyMult(n, col, row, val, b);
    printf("Multiplication valid: %d\n", valid);

    valid = verifySolve(n, col, row, val, bt[0][0], b);
    printf("Solution valid: %d\n", valid);
    return 0;
}
