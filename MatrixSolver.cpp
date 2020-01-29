//
// Created by Robin on 11-Jan-20.
//

#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include "MatrixSolver.h"

using namespace std;

int analyse(int n, const int *col, const int *row) {
    int p, j;
    if (!col || !row) return 0;
    auto x = new int[n];
#pragma omp parallel for default(none) shared(n, col, row, x) private(j, p)
    for (j = 0; j < n; j++) {
        for (p = col[j] + 1; p < col[j + 1]; p++) {
            x[row[p]] = j;
        }
    }
    return 1;
}

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
    int p, j;
    if (!col || !row || !b) return 0;
    for (j = 0; j < n; j++) {
        if (abs(b[j]) >= DBL_MIN) {
            b[j] /= val[col[j]];
            for (p = col[j] + 1; p < col[j + 1]; p++) {
                b[row[p]] -= val[p] * b[j];
            }
        }
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