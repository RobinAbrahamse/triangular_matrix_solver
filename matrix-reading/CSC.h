//
// Created by Robin on 15-Jan-20.
//

#include <string>

using namespace std;

int CSC(char *filepath, int* &col, int* &row, double* &val);

void matrixData(basic_ifstream<char> &file, int &n, int &Nz);

void compress(basic_ifstream<char> &file, int lines, int* col, int* row, double* val);