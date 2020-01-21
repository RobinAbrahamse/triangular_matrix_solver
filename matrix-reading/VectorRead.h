//
// Created by Robin on 17-Jan-20.
//

#include <string>
#include "CSC.h"

using namespace std;

#ifndef TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
#define TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H

bool isDense(string &filepath);

void readVector(string &filepath, double *&val);

void readVector(string &filepath, int* &ind, double* &val);

#endif //TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
