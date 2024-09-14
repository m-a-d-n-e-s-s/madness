#ifndef MAD_TESTPC_H
#define MAD_TESTPC_H

#include <iostream>
#include <complex>
#include <numeric>
#include <vector>
#include <cmath>

#include <lapacke.h>

#include <madness/misc/gnuplot.h>

const double pi = 3.14159265358979323846;
const double L = 1.0; // must be integer for cosine test to work, and > 1 for gaussian tests to work
const double xshift = 0.0; // shift from origin of charge distributions to test the periodicity in [-0.5*L,0.5*L]
const double yshift = 0.0; 
const double zshift = 0.0;

// Will be assigned by main program based on test_case selection
extern double (*f)(double, double, double);
extern double (*exact)(double, double, double);

void set_test_case(int test_case);

// print the matrix M[m,n] in row-major order
template <typename T>
void print(size_t m, size_t n, const std::vector<T>& M) {
    for (size_t i=0; i<m; i++) {
        for (size_t j=0; j<n; j++) {
            std::cout << i << " " << M[i*n+j] << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<double> tabulate(double(*f)(double, double, double), std::vector<double> x, std::vector<double> y, std::vector<double> z);

std::vector<double> linspace(double a, double b, size_t n, bool include_right_endpoint = true);

#endif
