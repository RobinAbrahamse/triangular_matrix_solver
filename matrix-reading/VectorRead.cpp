//
// Created by Robin on 17-Jan-20.
//

#include <fstream>
#include <string>
#include <iostream>
#include "VectorRead.h"

using namespace std;

bool readVector(string &filepath, int *&ind, double *&val) {
    ifstream file(filepath);
    string header;
    bool isDense;
    file.ignore(22); //ignore irrelevant header part
    file >> header;
    file.ignore(2048, '\n');
    ignoreComments(file);
    if (header == "array") {
        readDense(file, val);
        isDense = true;
    } else {
        readSparse(file, ind, val);
        isDense = false;
    }
    file.close();
    return isDense;
}

void readDense(basic_ifstream<char> &file, double* &val) {
    int n = 0;
    file >> n;
    file.ignore(2048, '\n');
    val = new double[n+1];

    val[0] = n;
    for (int k = 1; k <= n; k++) {
        file >> val[k];
    }
}

void readSparse(basic_ifstream<char> &file, int* &ind, double* &val) {
    int Nz = 0, _ = 0;
    file.ignore(2048, ' ');
    file.ignore(2048, ' ');
    file >> Nz;
    ind = new int[Nz+1];
    val = new double[Nz+1];

    ind[0] = Nz;
    val[0] = Nz;
    for (int k = 1; k <= Nz; k++) {
        file >> ind[k] >> _ >> val[k];
    }
}
