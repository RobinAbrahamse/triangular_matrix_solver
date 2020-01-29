//
// Created by Robin on 15-Jan-20.
//

#include <string>

using namespace std;

int CSC(char *filepath, int* &col, int* &row, double* &val);

void ignoreComments(basic_ifstream<char> &file);

int matrixDimensions(basic_ifstream<char> &file);

int nonzeroElements(basic_ifstream<char> &file);

void compress(basic_ifstream<char> &file, int lines, int* col, int* row, double* val);

int fast_atoi(const char *str);

double fast_atof (const char *c);