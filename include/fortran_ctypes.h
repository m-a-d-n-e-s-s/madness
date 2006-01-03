#ifndef FORTRAN_CTYPES_H
#define FORTRAN_CTYPES_H

/// \file fortran_ctypes.h
/// \brief Corresponding C and Fortran types

#include <complex>

/// Fortran integer
typedef long integer;

/// Fortran double precision
typedef double real8;
typedef double double_precision ;

/// Fortran single precision
typedef float real4;
typedef float single_precision;

/// Fortran double complex
typedef std::complex<double> complex_real8;
typedef std::complex<double> double_precision_complex;

/// Fortran single complex
typedef std::complex<float> complex_real4;
typedef std::complex<float> single_precision_complex;

/// Type of variable appended to argument list for length of fortran character strings
#ifdef _CRAY
    typedef long char_len;
#else
    typedef int char_len;
#endif


#endif
