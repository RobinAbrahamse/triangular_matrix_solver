//
// Created by Robin on 21-Jan-20.
//

#include <fstream>

using namespace std;

#ifndef TRIANGULAR_MATRIX_SOLVER_UTILS_H
#define TRIANGULAR_MATRIX_SOLVER_UTILS_H

int fast_atoi(const char *str);

double fast_atof(const char *c);

void ignoreComments(basic_ifstream<char> &file);

#endif //TRIANGULAR_MATRIX_SOLVER_UTILS_H
