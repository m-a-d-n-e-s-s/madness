#include <iostream>
#include <cmath>
#include <vector>
#include <misc/interpolation_1d.h>

using namespace std;

/// A simple program for testing the CubicInterpolationTable class.

double func(double x) {
    return sin(x);
}

int main() {
    cout.precision(12);

    // Uniform mesh for sin(x)
    CubicInterpolationTable<double> fit(-10.0,30.0,100000,func);
    cout << "maxerr " << fit.err(func) << endl;

    return 0;
}
