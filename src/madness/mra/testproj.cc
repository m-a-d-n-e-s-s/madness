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

/// \file testproj.cc
/// \brief test box size dependence of projection etc.

#include <madness/mra/mra.h>
#include <madness/constants.h>
#include <vector>

using namespace madness;

bool smalltest = true;

template <typename T, std::size_t NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;

    Gaussian(const coordT& center, double exponent, T coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {}

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (std::size_t i=0; i<NDIM; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    }

    std::vector<coordT> special_points() const {
    	return std::vector<coordT>(1,center);
    }
};



template <typename T, std::size_t NDIM>
int test_proj(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef std::shared_ptr< FunctionFunctorInterface<T,NDIM> > functorT;

    FunctionDefaults<NDIM>::set_k(6);
    FunctionDefaults<NDIM>::set_thresh(1e-8);
    FunctionDefaults<NDIM>::set_initial_level(5);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_autorefine(true);
    FunctionDefaults<NDIM>::set_truncate_mode(1);
    FunctionDefaults<NDIM>::set_truncate_on_project(false);
    int success=0;

    const double expnt = 100.0;
    const double coeff = pow(expnt/constants::pi,1.5);
    double x = 1.0/11.0; // A non-dyadic point
    //double x = 0.0;
    const coordT origin(x);

    int ilo=2, ihi=2;
    if (!smalltest) {ilo=0; ihi=7;}
    for (int i=ilo; i<=ihi; i+=1) {
        double L = pow(2.0,double(i));
        FunctionDefaults<NDIM>::set_cubic_cell(-L,L);
        print("I think the cell volume is", FunctionDefaults<NDIM>::get_cell_volume());

        functorT gaussfunctor(new Gaussian<T,NDIM>(origin, expnt, coeff));
        Gaussian<T,NDIM> gauss(origin, expnt, coeff);
        Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(gaussfunctor);
        f.truncate();
        f.reconstruct();
        print("L",L,f.trace(),f.norm2(),f.size()/6/6/6,f.max_depth());
        double err=f.err(gauss);
        if (world.rank()==0) print("error in ",NDIM, "dimensions ", err);
        if (err>1.2*FunctionDefaults<3>::get_thresh()) success++;

        //f.print_tree();

    }
    return success;
}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);
    int success=0;

    try {
        startup(world,argc,argv);
        if (getenv("MAD_SMALL_TESTS")) smalltest=true;
        for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
        std::cout << "small test : " << smalltest << std::endl;

        std::cout.precision(8);

        success+=test_proj<double,1>(world);
        success+=test_proj<double,2>(world);
        if (!smalltest) success+=test_proj<double,3>(world);

    }
    catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
//    catch (const char* s) {
//        print(s);
//        error("caught a c-string exception");
//    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    world.gop.fence();
    finalize();

    return success;
}

