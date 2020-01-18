//
// Created by Robin on 15-Jan-20.
//

#include <fstream>
#include <string>
#include <iostream>
#include "CSC.h"

using namespace std;

int CSC(string &filepath, int* &col, int* &row, double* &val) {
    ifstream file(filepath);
    ignoreComments(file);
    int n = matrixDimensions(file);
    int Nz = nonzeroElements(file);
    col = new int[n+2];
    row = new int[Nz+1];
    val = new double[Nz+1];
    compress(file, Nz, col, row, val);
    file.close();
    return n;
}

void ignoreComments(basic_ifstream<char> &file) {
    while (file.peek() == '%') file.ignore(2048, '\n');
}

int matrixDimensions(basic_ifstream<char> &file) {
    int n = 0;
    file.ignore(2048, ' ');
    file >> n;
    cout << "Matrix dimension: " << n << endl;
    return n;
}

int nonzeroElements(basic_ifstream<char> &file) {
    int Nz = 0;
    file >> Nz;
    cout << "Non-zero elements: " << Nz << endl;
    return Nz;
}

void compress(basic_ifstream<char> &file, int lines, int* col, int* row, double* val) {
    int j = 0, j_prev = 0; //j_prev instead of col[j-1] to prevent caching the col array (reading col is not required)

    row[0] = lines;
    val[0] = lines;
    for (int k = 1; k <= lines; k++) {
        file >> row[k] >> j >> val[k];
        if (j > j_prev) {
            col[j] = k;
            j_prev = j;
        }
    }
    col[j+1] = lines + 1;
    col[0] = j;
}
