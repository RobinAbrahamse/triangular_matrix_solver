//
// Created by Robin on 17-Jan-20.
//

#ifndef TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
#define TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H

#include <cfloat>

int analyseDependencies(int n, const int *col, const int *row, int *&d);

void createLevelsets(int n, const int *col, const int *row, int *&levels, int *&levelpntrs);

int analyseLevels(int n, const int *col, const int *row, int *levels);

void sortLevels(int n, int nlev, int *levels, int *levelpntrs);

int solve(int n, const int* col, const int* row, const double* val, double *b);

int mult(int n, int *col, int *row, double *val, double *x, double *y);

int serialMult(int n, int *col, int *row, double *val, double *x, double *y);

int verifyMult(int n, int *col, int *row, double *val, double *x);

int verifySolve(int n, int* col, int* row, double* val, double* x, double* b);

int nearlyEqual(double a, double b, double epsilon = 128 * FLT_EPSILON, double relth = FLT_MIN);

#endif //TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
