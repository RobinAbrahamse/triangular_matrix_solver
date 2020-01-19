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
    col = new int[n+1];
    row = new int[Nz];
    val = new double[Nz];
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
    int i = 0, j = 0, j_prev = 0; //j_prev instead of col[j-1] to prevent caching the col array (reading col is not required)

    for (int k = 0; k < lines; k++) {
        file >> i >> j >> val[k];
        row[k] = i - 1;
        if (j > j_prev) {
//        if (j <= j_prev) {} else {
            col[j-1] = k;
            j_prev = j;
        }
    }
    col[j] = lines;
}
