//
// Created by Robin on 11-Jan-20.
//

#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include "MatrixSolver.h"

using namespace std;

/*
* Lower triangular solver Lx=b
* L is stored in the compressed column storage format
* Inputs are:
* n : the matrix dimension
* col : the column pointer of L
* row : the row index of L
* val : the values of L
* In/Out:
* b : the right hand-side b at start and the solution x at the end.
*/
int solve(int n, const int *col, const int *row, const double *val, double *b) {
    if (!col || !row || !b) return 0;
    int p, j;
    auto d_in_degree = new int[n];
    auto s_in_degree = new int[n];
    auto d_sum = new double[n];
    auto s_sum = new double[n];
    analyse(n, col, row, d_in_degree);
#pragma omp parallel default(none) shared(n, col, row, val, b, d_in_degree, s_in_degree, d_sum, s_sum, j) private(p)
    {
#pragma omp single
        {
            for (j = 0; j < n; j++) {
#pragma omp task default(none) shared(col, row, val, b, d_in_degree, s_in_degree, d_sum, s_sum) private(p, j)
                {
                    while (s_in_degree[j] + 1 != d_in_degree[j]) {
                        printf("Thread: %d, task: %d, expected: %d, got: %d\n", omp_get_thread_num(), j, d_in_degree[j+1], s_in_degree[j]);
#pragma omp taskyield
                    }      //Lock-wait section
                    b[j] = (b[j] - d_sum[j] - s_sum[j]) / val[col[j]];
                    double t;
                    for (p = col[j] + 1; p < col[j + 1]; p++) {
                        t = val[row[p]] * b[j];
#pragma omp atomic
                        d_sum[j] += t;                                                  //Critical section
#pragma omp atomic
                        d_in_degree[j]--;                                               //Lock-update section
                    }
                };
            }
        };
    }
    delete[] d_in_degree;
    delete[] s_in_degree;
    delete[] d_sum;
    delete[] s_sum;
    return 1;
}

//int solve(int n, const int *col, const int *row, const double *val, double *b) {
//    if (!col || !row || !b) return 0;
//    int p, j;
//    auto d_in_degree = new int[n];
//    auto s_in_degree = new int[n];
//    auto d_sum = new double[n];
//    auto s_sum = new double[n];
//    analyse(n, col, row, d_in_degree);
//#pragma omp parallel default(none) shared(n, col, row, val, b, d_in_degree, s_in_degree, d_sum, s_sum, j) private(p)
//    {
//        #pragma omp for
//        for (j = 0; j < n; j++) {
//            while (s_in_degree[j] + 1 != d_in_degree[j]) {
////                        printf("Thread: %d, task: %d, expected: %d, got: %d\n", omp_get_thread_num(), j, d_in_degree[j+1], s_in_degree[j]);
//            }      //Lock-wait section
//            b[j] = (b[j] - d_sum[j] - s_sum[j]) / val[col[j]];
//            double t;
//            for (p = col[j] + 1; p < col[j + 1]; p++) {
//                t = val[row[p]] * b[j];
//                #pragma omp atomic
//                d_sum[j] += t;                                                  //Critical section
//                #pragma omp atomic
//                d_in_degree[j]--;                                               //Lock-update section
//            }
//        }
//    }
//    delete[] d_in_degree;
//    delete[] s_in_degree;
//    delete[] d_sum;
//    delete[] s_sum;
//    return 1;
//}

//int solve(int n, const int *col, const int *row, const double *val, double *b) {
//    int p, j;
//    if (!col || !row || !b) return 0;
//    for (j = 0; j < n; j++) {
//        if (abs(b[j]) >= DBL_MIN) {
//            b[j] /= val[col[j]];
//            for (p = col[j] + 1; p < col[j + 1]; p++) {
//                b[row[p]] -= val[p] * b[j];
//            }
//        }
//    }
//    return 1;
//}

int analyse(int n, const int *col, const int *row, int *d) {
    int j;
#pragma omp parallel for default(none) shared(n, col, row) private(j) reduction(+:d[:n])
    for (j = 0; j < col[n]; j++) {
        d[row[j]]++;
    }
    return 1;
}

/*
* Sparse matrix-vector multiply: y = A*x
* A is stored in the compressed column storage format
* Inputs:
* col : the column pointer of A
* row : the row index of A
* val : the values of A
* x : is a dense vector
* Output:
* y : is a dense vector that stores the result of multiplication
*/
int mult(int n, int *col, int *row, double *val, double *x, double *y) {
    int p, j;
    double xj;
    if (!col || !x || !y) return (0);
#pragma omp parallel for default(none) shared(n, col, row, val, x, y) private(p, j, xj)
    for (j = 0; j < n; j++) {
        xj = x[j];
        for (p = col[j]; p < col[j + 1]; p++) {
            #pragma omp atomic                  //Atomic instead of reduction because array size (stack overflow)
            y[row[p]] += val[p] * xj;
        }
    }
    return (1);
}

int serial_mult(int n, int *col, int *row, double *val, double *x, double *y) {
    int p, j;
    double xj;
    if (!col || !x || !y) return (0);
    for (j = 0; j < n; j++) {
        xj = x[j];
        for (p = col[j]; p < col[j + 1]; p++) {
            y[row[p]] += val[p] * xj;
        }
    }
    return (1);
}

int verify(int n, int *col, int *row, double *val, double *x, double *b) {
    auto *y = new double[n];
    auto t1 = omp_get_wtime();
    mult(n, col, row, val, x, y);
    t1 = omp_get_wtime() - t1;
    printf("Time to parallel multiply: %f\n", t1);

    y = new double[n];
    auto t2 = omp_get_wtime();
    serial_mult(n, col, row, val, x, y);
    t2 = omp_get_wtime() - t2;
    printf("Time to serial multiply: %f\n", t2);

    for (int i = 0; i < n; i++) {
        if (!nearly_equal(y[i], b[i])) {
            cout << setprecision(15) << "Expected: " << b[i] << ", but got: " << y[i] << ", index: " << i << endl; return 0;
        }
    }
    return 1;
}

bool nearly_equal(double a, double b, double epsilon, double relth) {
    if (a == b) return true;

    auto diff = abs(a-b);
    auto norm = min((abs(a) + abs(b)), numeric_limits<double>::max());
    return diff < max(relth, epsilon * norm);
}