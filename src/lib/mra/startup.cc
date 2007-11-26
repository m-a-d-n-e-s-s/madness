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

  
/// \file mra/startup.cc

#include <mra/mra.h>
#include <iomanip>

namespace madness {
    void startup(World& world, int argc, char** argv) {
        for (int arg=1; arg<argc; arg++) {
            if (std::strcmp(argv[arg],"-dx")==0) 
                xterm_debug("world", 0);
            else if (std::strcmp(argv[arg],"-dn") ==0 && 
                     std::atoi(argv[arg+1])==world.rank()) 
                xterm_debug("world",0);
            else if (std::strcmp(argv[arg],"-dam")==0) 
                world.am.set_debug(true);
            else if (std::strcmp(argv[arg],"-dmpi")==0) 
                world.mpi.set_debug(true);
            else if (std::strcmp(argv[arg],"-rio")==0)
                redirectio(world);
        }

        world.gop.fence();

        std::cout << std::boolalpha;  // Pretty printing of booleans
        std::cout << std::scientific;
        std::cout << std::showpoint; 
        //std::cout << std::showpos;
        std::cout << std::setprecision(8);

#ifdef FUNCTION_INSTANTIATE_1
        FunctionDefaults<1>::set_defaults(world);
#endif
#ifdef FUNCTION_INSTANTIATE_2
        FunctionDefaults<2>::set_defaults(world);
#endif
#ifdef FUNCTION_INSTANTIATE_3
        FunctionDefaults<3>::set_defaults(world);
#endif
#ifdef FUNCTION_INSTANTIATE_4
        FunctionDefaults<4>::set_defaults(world);
#endif
#ifdef FUNCTION_INSTANTIATE_5
        FunctionDefaults<5>::set_defaults(world);
#endif
#ifdef FUNCTION_INSTANTIATE_6
        FunctionDefaults<6>::set_defaults(world);
#endif


        if (world.rank() == 0) print("loading coeffs, etc.");

        load_coeffs(world);
        load_quadrature(world);

        if (world.rank() == 0) print("testing coeffs, etc.");
        MADNESS_ASSERT(gauss_legendre_test());
        MADNESS_ASSERT(test_two_scale_coefficients());

        if (world.rank() == 0) print("done with startup");

        world.gop.fence();
    }
}
