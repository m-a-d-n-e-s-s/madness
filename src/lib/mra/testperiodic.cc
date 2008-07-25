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

/// \file testperiodic.cc
/// \brief test periodic convolutiosn

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

typedef Vector<double,3> coordT;

double source(const coordT& r) {
    return cos(2*constants::pi*r[0])*cos(2*constants::pi*r[1])*cos(2*constants::pi*r[2]);
}

double potential(const coordT& r) {
    return -source(r)/3.0;
}

void test_periodic(World& world) {
    const long k = 10;
    const double thresh = 1e-8;
    const double L = 0.5;
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_thresh(thresh);
    
    Function<double,3> f = FunctionFactory<double,3>(world).f(source);

    std::vector< SharedPtr< Convolution1D<double> > > ops(1);

    ops[0] = SharedPtr< Convolution1D<double> >(new PeriodicGaussianConvolution1D<double>(k, 12, 100.0, 80000.0));

    SeparatedConvolution<double,3> op(world, k, ops, false, true);

    Function<double,3> opf = apply(op,f);

    coordT r(0.49);

    opf.reconstruct();
    f.reconstruct();

    cout.precision(10);

    print(r, f(r), source(r), opf(r));

    world.gop.fence();

}


int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    try {
        startup(world,argc,argv);
        
        test_periodic(world);

    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    MPI::Finalize();

    return 0;
}

