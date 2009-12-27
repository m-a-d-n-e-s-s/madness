/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
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
#include <madness_config.h>

#include <tensor/aligned.h>
#include <tensor/tensor.h>

namespace madness {

#if (defined(X86_32) || defined(X86_64))
    template <> 
    void aligned_zero<double_complex>(long n, double_complex* a) {
      aligned_zero(2*n, (double *) a);
    }
#endif

    void aligned_add(long n, double_complex* restrict a, const double_complex* restrict b) {
        aligned_add(2*n, (double*) a, (const double*) b);
    }

    void aligned_sub(long n, double_complex* restrict a, const double_complex* restrict b) {
        aligned_sub(2*n, (double*) a, (const double*) b);
    }
}


