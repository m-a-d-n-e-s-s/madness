/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov 
  tel:   865-241-3937
  fax:   865-572-0680

  
  $Id$
*/

  
#ifndef FORTRAN_CTYPES_H
#define FORTRAN_CTYPES_H

/// \file fortran_ctypes.h
/// \brief Corresponding C and Fortran types

#include <complex>

/// Fortran integer
#ifdef _CRAY
#  ifndef MADNESS_FORINT
#    define MADNESS_FORINT int
#  endif
#endif

#ifdef MADNESS_FORINT
typedef MADNESS_FORINT integer;
#else
typedef long integer;
#endif


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
