//
// Created by Robin on 17-Jan-20.
//

#include <string>
#include "CSC.h"

using namespace std;

#ifndef TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
#define TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H

bool isDense(char *filepath);

void readVector(char *filepath, double *&val);

void readVector(char *filepath, int* &ind, double* &val);

#endif //TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
