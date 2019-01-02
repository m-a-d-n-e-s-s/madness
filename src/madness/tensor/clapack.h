/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  Copyright (C) 2018 Virginia Tech

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

  
#ifndef MADNESS_LINALG_CLAPACK_H__INCLUDED
#define MADNESS_LINALG_CLAPACK_H__INCLUDED

/// \file clapack.h
/// \brief C++ interface to LAPACK, either directly via Fortran API (see clapack_fortran.h) or via LAPACKE (see clapack_lapacke.h)

#ifdef MADNESS_LINALG_USE_LAPACKE
#  if HAVE_INTEL_MKL
#    if !__has_include(<mkl_lapacke.h>)
#      error "INTEL MKL detected at configure time, and MADNESS_LINALG_USE_LAPACKE defined, but mkl_lapacke.h not found. Provide -I/path/to/mkl/include to the compiler."
#    endif
#    if !__has_include(<mkl_lapack.h>)
#      error "INTEL MKL detected at configure time, and MADNESS_LINALG_USE_LAPACKE defined, but mkl_lapack.h not found. Provide -I/path/to/mkl/include to the compiler."
#    endif
#    include <mkl_lapacke.h>
#    include <mkl_lapack.h>
#  else
#    if !__has_include(<lapacke.h>)
#      error "MADNESS_LINALG_USE_LAPACKE defined, but lapacke.h not found. Provide -I/path/to/lapacke/dot/h to the compiler."
#    endif
#    include <lapacke.h>
#  endif
#else
#  include "clapack_fortran.h"
#endif

#endif // MADNESS_LINALG_CLAPACK_H__INCLUDED
