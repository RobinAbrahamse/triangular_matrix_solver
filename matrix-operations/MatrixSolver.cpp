//
// Created by Robin on 11-Jan-20.
//

#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <algorithm>
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
int solveParallelTasks(int n, const int *col, const int *row, const double *val, double *b) {
    if (!col || !row || !b) return 0;
    int p, j;
    auto d_in_degree = new int[n];
    analyseDependencies(n, col, row, d_in_degree);
    #pragma omp parallel default(none) shared(n, col, row, val, b, d_in_degree, j) private(p)
    {
        #pragma omp single nowait
        {
            for (j = 0; j < n; j++) {
                #pragma omp task default(none) shared(col, row, val, b, d_in_degree) private(p, j)
                {
                    while (d_in_degree[j] != 1) {
                        #pragma omp taskyield
                    }      //Lock-wait section
                    b[j] /= val[col[j]];
                    double t;
                    int rid;
                    for (p = col[j] + 1; p < col[j + 1]; p++) {
                        rid = row[p];
                        t = val[p] * b[j];
                        #pragma omp atomic
                        b[rid] -= t;                                                  //Critical section
                        #pragma omp atomic
                        d_in_degree[rid]--;                                           //Lock-update section
                    }
                };
            }
        };
    };
    delete[] d_in_degree;
    return 1;
}

int solveParallelFor(int n, const int *col, const int *row, const double *val, double *b) {
    if (!col || !row || !b) return 0;
    int p, j;
    auto d_in_degree = new int[n]();
    analyseDependencies(n, col, row, d_in_degree);
    #pragma omp parallel for default(none) shared(n, col, row, val, b, d_in_degree) private(p, j)
    for (j = 0; j < n; j++) {
        while (d_in_degree[j] != 1) {/*Spin thread*/}                       //Lock-wait section
        b[j] /= val[col[j]];
        double t;
        int rid;
        for (p = col[j] + 1; p < col[j + 1]; p++) {
            rid = row[p];
            t = val[p] * b[j];
            #pragma omp atomic
            b[rid] -= t;                                                    //Critical section
            #pragma omp atomic
            d_in_degree[rid]--;                                             //Lock-update section
        }
    }
    delete[] d_in_degree;
    return 1;
}

int solveParallelGPU(int n, const int *col, const int *row, const double *val, double *b) {
    if (!col || !row || !b) return 0;
    int p, j;
    auto d_in_degree = new int[n];
    analyseDependencies(n, col, row, d_in_degree);
    #pragma omp target teams distribute default(none) shared(n, col, row, val, b, d_in_degree) private(p, j)
    for (j = 0; j < n; j++) {
        while (d_in_degree[j] != 1) {/*Spin thread*/}                      //Lock-wait section
        b[j] /= val[col[j]];
        double t;
        int rid;
        #pragma omp parallel for default(none) shared(col, row, val, b, d_in_degree, j) private(p, t, rid)
        for (p = col[j] + 1; p < col[j + 1]; p++) {
            rid = row[p];
            t = val[p] * b[j];
            #pragma omp atomic
            b[rid] -= t;                                                  //Critical section
            #pragma omp atomic
            d_in_degree[rid]--;                                           //Lock-update section
        }
    }
    delete[] d_in_degree;
    return 1;
}

int analyseDependencies(int n, const int *col, const int *row, int *&d) {
    int j;
    d = new int[n]();
#pragma omp parallel for default(none) shared(n, col, row) private(j) reduction(+:d[:n])
    for (j = 0; j < col[n]; j++) {
        d[row[j]]++;
    }
    return 1;
}

int solveParallelLevels(int n, const int *col, const int *row, const double *val, double *b) {
    if (!col || !row || !b) return 0;
    int p, j, k, l, nlev;
    int* levels = nullptr;
    int* levelptrs = nullptr;
    nlev = createLevelsets(n, col, row, levels, levelptrs);
    for (k = 0; k < nlev; k++) {
//#pragma omp parallel for default(none) shared(n, col, row, val, b, levels, levelptrs, k) private(p, j, l)
        for (l = levelptrs[k]; l < levelptrs[k+1]; l++) {
            j = levels[l];
            b[j] /= val[col[j]];
            for (p = col[j] + 1; p < col[j+1]; p++) {
//                #pragma omp critical
                b[row[p]] -= val[p] * b[j];
            }
        }
    }
    delete[] levels;
    delete[] levelptrs;
    return 1;
}

int createLevelsets(int n, const int *col, const int *row, int *&levels, int *&levelptrs) {
    levels = new int[n]();
    int nlev = analyseLevels(n, col, row, levels);
    levelptrs = new int[nlev + 1];
    sortLevels(n, nlev, levels, levelptrs);
    return nlev;
}

int analyseLevels(int n, const int *col, const int *row, int *levels) {
    int p, j, t, nlev = 0;
//#pragma omp parallel for default(none) shared(n, col, row) private(j, p, t) reduction(max:levels[:n], nlev)
    for (j = 0; j < n; j++) {
        for (p = col[j]; p < col[j + 1]; p++) {
            t = max(levels[j] + 1, levels[row[p]]);
            nlev = max(nlev, t);
            levels[row[p]] = t;
        }
    }
    return nlev;
}

void sortLevels(int n, int nlev, int *levels, int *levelptrs) {
    int t, j;
    auto perm = new long long[n]();
//#pragma omp parallel for default(none) shared(n, levels, perm) private(j)
    for (j = 0; j < n; j++) {
        perm[j] = ((long long) j << 32) | levels[j];
    }
    sort(perm, perm+n, [](long long a, long long b) {
        return (int) (a & 0xFFFFFFFF) > (int) (b & 0xFFFFFFFF);
    });

    for (j = 0; j < nlev; j++) {
        levelptrs[j] = INT_MAX;
    }
//#pragma omp parallel for default(none) shared(n, nlev, levels, t, perm) private(j) reduction(min:levelptrs[:nlev])
    for (j = 0; j < n; j++) {
        t = perm[j] & 0xFFFFFFFF;
        levelptrs[t] = min(levelptrs[t], j);
        levels[j] = (int) (perm[j] >> 32);
    }
    delete[] perm;
}

int solveSerial(int n, const int *col, const int *row, const double *val, double *b) {
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
* Parallelized sparse matrix-vector multiply: y = A*x
* A is stored in the compressed column storage format
* Inputs:
* n : the matrix dimension
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

/*
* Serial sparse matrix-vector multiply: y = A*x
* A is stored in the compressed column storage format
* Inputs:
* n : the matrix dimension
* col : the column pointer of A
* row : the row index of A
* val : the values of A
* x : is a dense vector
* Output:
* y : is a dense vector that stores the result of multiplication
*/
int serialMult(int n, int *col, int *row, double *val, double *x, double *y) {
    int p, j;
    if (!col || !x || !y) return (0);
    for (j = 0; j < n; j++) {
        for (p = col[j]; p < col[j + 1]; p++) {
            y[row[p]] += val[p] * x[j];
        }
    }
    return (1);
}

/*
* Verification function for the matrix solver
* Multiplies the matrix stored in [col], [row], [val] by the solution vector [x]
* Then compares results of multiplication to the result vector [b]
* Inputs:
* n : the matrix dimension
* col : the column pointer of A
* row : the row index of A
* val : the values of A
* x : the solution vector
* b : the vector to be solved for
* Output:
* Returns the validity of the solution
*/
int verifySolve(int n, int *col, int *row, double *val, double *x, double *b) {
    auto *y = new double[n];
    mult(n, col, row, val, x, y);
    for (int i = 0; i < n; i++) {
        if (!nearlyEqual(y[i], b[i])) {
            cout << setprecision(15) << "Expected: " << b[i] << ", but got: " << y[i] << ", index: " << i << endl;
            return 0;
        }
    }
    return 1;
}

/*
* Verification function for the matrix-vector multiplier
* Compares output of the parallel and pre-verified serial multiplier functions
* Inputs:
* n : the matrix dimension
* col : the column pointer of A
* row : the row index of A
* val : the values of A
* x : the dense vector
* Output:
* Returns the validity of the multiplication
*/
int verifyMult(int n, int *col, int *row, double *val, double *x) {
    auto y1 = new double[n];
    auto t1 = omp_get_wtime();
    mult(n, col, row, val, x, y1);
    t1 = omp_get_wtime() - t1;
    printf("Time to parallel multiply: %f\n", t1);

    auto y2 = new double[n];
    auto t2 = omp_get_wtime();
    serialMult(n, col, row, val, x, y2);
    t2 = omp_get_wtime() - t2;
    printf("Time to serial multiply: %f\n", t2);

    for (int i = 0; i < n; i++) {
        if (!nearlyEqual(y1[i], y2[i])) {
            cout << setprecision(15) << "Expected: " << y2[i] << ", but got: " << y1[i] << ", index: " << i << endl;
            return 0;
        }
    }
    return 1;
}

/*
* Fuzzy equals method for double precision numbers
* Input:
* a : first number
* b : second number
* epsilon : relative precision parameter
* relth : absolute precision parameter
* Output:
* Returns the fuzzy equality of [a] and [b]
*/
int nearlyEqual(double a, double b, double epsilon, double relth) {
    if (a == b) return true;

    auto diff = abs(a-b);
    auto norm = min((abs(a) + abs(b)), numeric_limits<double>::max());
    return diff < max(relth, epsilon * norm);
}