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

  $LastChangedDate$
  $Rev$
*/

  
#ifndef MTRAND_H
#define MTRAND_H

namespace madness {

    /// Initialize Mersenne twister random number generator
    extern void init_genrand(unsigned long s);

    /// Initialize Mersenne twister random number generator

    /// Initialize by an array with array-length
    /// init_key is the array for initializing keys
    /// key_length is its length
    extern void init_by_array(unsigned long init_key[], int key_length);

    /// Generates a random number on [0,0xffffffff]-interval
    extern unsigned long genrand_int32(void);

    /// Generates a random number on [0,0x7fffffff]-interval
    extern long genrand_int31(void);

    /// Generates a random number on [0,1]-real-interval
    extern double genrand_real1(void);

    /// Generates a random number on [0,1)-real-interval
    extern double genrand_real2(void);

    /// Generates a random number on (0,1)-real-interval
    extern double genrand_real3(void);

    /// Generates a random number on [0,1) with 53-bit resolution
    double genrand_res53(void) ;

}

#endif
