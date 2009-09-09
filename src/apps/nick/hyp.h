//By: Robert Harrison
#include <iostream>
#include <complex>
#include <algorithm>
#include <cstdio> //NEEDED
#include <cmath>
#include <nick/mpreal.h>


typedef mpfr::mpreal extended_real;
typedef std::complex<extended_real> extended_complex;
typedef std::complex<double> complexd;


/// Computes 1F1(a,b,z) internally using extended precision

/// If result is larger than 1.0, result should be accurate to
/// full double precision, otherwise result is accurate to
/// somewhat better than 1e-17.  However, the termination
/// test is not very smart so there may be failures.
complexd conhyp(const complexd& a_arg,
                const complexd& b_arg,
                const complexd& z_arg);
