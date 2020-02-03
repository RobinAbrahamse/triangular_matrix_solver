//
// Created by Robin on 15-Jan-20.
//

#include <fstream>
#include <string>
#include <iostream>
#include "Utils.h"
#include "CSC.h"

using namespace std;

int CSC(char *filepath, int* &col, int* &row, double* &val) {
    int n, Nz;
    ifstream file(filepath);
    ignoreComments(file);
    matrixData(file, n, Nz);
    col = new int[n+1];
    row = new int[Nz];
    val = new double[Nz];
    compress(file, Nz, col, row, val);
    file.close();
    return n;
}

void matrixData(basic_ifstream<char> &file, int &n, int &Nz) {
    string s;
    file.ignore(2048, ' ');
    getline(file, s, ' ');
    n = fast_atoi(s.c_str());
    getline(file, s);
    Nz = fast_atoi(s.c_str());
    cout << "Matrix dimension: " << n << endl;
    cout << "Non-zero elements: " << Nz << endl;
}

void compress(basic_ifstream<char> &file, int lines, int* col, int* row, double* val) {
    int i = 0, j = 0, j_prev = 0; //j_prev instead of col[j-1] to prevent caching the col array (reading col is not required)
    double x = 0.0;
    string s;

    for (int k = 0; k < lines; k++) {
        getline(file, s, ' ');
        i = fast_atoi(s.c_str());
        getline(file, s, ' ');
        j = fast_atoi(s.c_str());
        getline(file, s);
        x = fast_atof(s.c_str());

        row[k] = i - 1;
        val[k] = x;
        if (j > j_prev) {
            col[j-1] = k;
            j_prev = j;
        }
    }
    col[j] = lines;
}
