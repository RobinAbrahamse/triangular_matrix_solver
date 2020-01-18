//
// Created by Robin on 17-Jan-20.
//

#include <string>
#include "CSC.h"

using namespace std;

#ifndef TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
#define TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H

bool readVector(string &filepath, int* &ind, double* &val);

void readDense(basic_ifstream<char> &file, double* &val);

void readSparse(basic_ifstream<char> &file, int* &ind, double* &val);

extern void ignoreComments(basic_ifstream<char> &file);

#endif //TRIANGULAR_MATRIX_SOLVER_VECTORREAD_H
