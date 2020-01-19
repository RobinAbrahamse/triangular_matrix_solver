//
// Created by Robin on 17-Jan-20.
//

#ifndef TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
#define TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H

int solve (int n, const int* col, const int* row, const double* val, double *b);

int mult(int n, int *col, int *row, double *val, double *x, double *y);

bool verify(int n, int* col, int* row, double* val, double* x, double* b);

#endif //TRIANGULAR_MATRIX_SOLVER_MATRIXSOLVER_H
