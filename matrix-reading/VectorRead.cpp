//
// Created by Robin on 17-Jan-20.
//

#include <fstream>
#include <string>
#include "Utils.h"
#include "VectorRead.h"

using namespace std;

bool isDense(char *filepath) {
    ifstream file(filepath);
    string header;
    file.ignore(22); //ignore irrelevant header part
    getline(file, header, ' '); //get format type (array or coordinate)
    file.close();
    return header == "array";
}

void readVector(char *filepath, double *&val) {
    ifstream file(filepath);
    ignoreComments(file);
    string s;
    double x;
    getline(file, s, ' ');
    int n = fast_atoi(s.c_str());
    file.ignore(2048, '\n');
    val = new double[n];

    for (int k = 0; k < n; k++) {
        getline(file, s);
        val[k] = fast_atof(s.c_str());
    }
    file.close();
}

void readVector(char *filepath, int *&ind, double *&val) {
    ifstream file(filepath);
    ignoreComments(file);
    string s;
    file.ignore(2048, ' ');
    file.ignore(2048, ' ');
    getline(file, s);
    int Nz = fast_atoi(s.c_str());
    ind = new int[Nz];
    val = new double[Nz];

    for (int k = 0; k < Nz; k++) {
        getline(file, s, ' ');
        ind[k] = fast_atoi(s.c_str());
        file.ignore(2048, ' ');
        getline(file, s);
        val[k] = fast_atof(s.c_str());
    }
    file.close();
}
