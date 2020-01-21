//
// Created by Robin on 21-Jan-20.
//

#include <fstream>

using namespace std;

/**
 * Source: http://www.leapsecond.com/tools/fast_atof.c
 *
 * @param c
 * c-string containing double
 * @return
 * double number contained in c
 */
#define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')
double fast_atof(const char *c) {
    int frac;
    double sign, value, scale;

    // Skip leading white space, if any.
    while (white_space(*c) ) {
        c += 1;
    }

    // Get sign, if any.
    sign = 1.0;
    if (*c == '-') {
        sign = -1.0;
        c += 1;

    } else if (*c == '+') {
        c += 1;
    }

    // Get digits before decimal point or exponent, if any.
    for (value = 0.0; valid_digit(*c); c += 1) {
        value = value * 10.0 + (*c - '0');
    }

    // Get digits after decimal point, if any.
    if (*c == '.') {
        double pow10 = 10.0;
        c += 1;
        while (valid_digit(*c)) {
            value += (*c - '0') / pow10;
            pow10 *= 10.0;
            c += 1;
        }
    }

    // Handle exponent, if any.
    frac = 0;
    scale = 1.0;
    if ((*c == 'e') || (*c == 'E')) {
        unsigned int expon;

        // Get sign of exponent, if any.
        c += 1;
        if (*c == '-') {
            frac = 1;
            c += 1;
        } else if (*c == '+') {
            c += 1;
        }

        // Get digits of exponent, if any.
        for (expon = 0; valid_digit(*c); c += 1) {
            expon = expon * 10 + (*c - '0');
        }
        if (expon > 308) expon = 308;

        // Calculate scaling factor.
        while (expon >= 50) { scale *= 1E50; expon -= 50; }
        while (expon >=  8) { scale *= 1E8;  expon -=  8; }
        while (expon >   0) { scale *= 10.0; expon -=  1; }
    }

    // Return signed and scaled floating point result.
    return sign * (frac ? (value / scale) : (value * scale));
}

/**
 * Source: https://stackoverflow.com/questions/16826422/c-most-efficient-way-to-convert-string-to-int-faster-than-atoi
 *
 * @param c
 * c-string containing integer
 * @return
 * integer number contained in c
 */
int fast_atoi(const char *c) {
    int x = 0;
    while(*c) {
        x = x*10 + (*c++ - '0');
    }
    return x;
}

/**
 * Ignores the comment lines before the MM format matrix values
 *
 * @param file
 * ifstream containing MM format matrix
 */
void ignoreComments(basic_ifstream<char> &file) {
    while (file.peek() == '%') file.ignore(2048, '\n');
}
