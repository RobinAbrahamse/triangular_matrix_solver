//
// Created by Robin on 11-Jan-20.
//

#include <cmath>

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
    if (!col || !row || !b) return (0);
//#pragma omp parallel for default(none) shared(n, b, val, col, row) private(p)
    for (j = 0; j < n; j++)
    {
        b[j] /= val[col[j]];
        for (p = col[j] + 1; p < col[j + 1]; p++)
        {
            b[row[p]] -= val[p] * b[j];
        }
    }
    return (1);
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
    if (!col || !x || !y) return (0);
    for (j = 0; j < n; j++)
    {
        for (p = col[j]; p < col[j + 1]; p++)
        {
            y[row[p]] += val[p] * x[j];
        }
    }
    return (1);
}

int verify(int n, int *col, int *row, double *val, double *x, double *b) {
    auto *y = new double[n+1];
    mult(n, col, row, val, x, y);
    for (int i = 0; i < n; i++)
    {
        if (fabs(y[i] - b[i]) > 0.0000001)
        return 0;
    }
    return 1;
}
