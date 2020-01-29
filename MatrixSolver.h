//
// Created by Robin on 17-Jan-20.
//

#ifndef TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
#define TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H

#include <cfloat>

int analyse(int n, const int *col, const int *row);

int solve(int n, const int* col, const int* row, const double* val, double *b);

int mult(int n, int *col, int *row, double *val, double *x, double *y);
int serial_mult(int n, int *col, int *row, double *val, double *x, double *y);

int verify(int n, int* col, int* row, double* val, double* x, double* b);

bool nearly_equal(double a, double b, double epsilon = 128 * FLT_EPSILON, double relth = FLT_MIN);

#endif //TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
