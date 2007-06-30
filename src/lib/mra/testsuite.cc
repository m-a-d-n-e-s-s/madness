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

  
  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/

/// \file testsuite.cc
/// \brief The QA/test suite for Function

  
#include <mra/mra.h>

using namespace madness;

template <typename T, int NDIM>
class GaussianFunctor : public FunctionFunctorInterface<T,NDIM> {
private:
    typedef Vector<double,NDIM> coordT;            ///< Type of vector holding coordinates
    const coordT center;
    const double exponent;
    const T coefficient;
    
public:
    GaussianFunctor(const coordT& center, double exponent, T coefficient) 
        : center(center), exponent(exponent), coefficient(coefficient) {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (int i=0; i<NDIM; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    };
};

template <typename T, int NDIM>
void test_basic(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) 
        print("Test compression of a normalized gaussian at origin, type =",archive::get_type_name<T>(),", ndim =",NDIM);

//     for (int i=0; i<NDIM; i++) {
//         FunctionDefaults<NDIM>::cell(i,0) = -10.0;
//         FunctionDefaults<NDIM>::cell(i,1) =  10.0;
//     }
    FunctionDefaults<NDIM>::k = 7;
    FunctionDefaults<NDIM>::thresh = 1e-12;
    FunctionDefaults<NDIM>::compress = false;
    FunctionDefaults<NDIM>::refine = true;
    FunctionDefaults<NDIM>::initial_level = 2;
    
    double used;
    const double PI = 3.1415926535897932384;
    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(1.0/PI,0.25*NDIM);

    functorT functor(new GaussianFunctor<T,NDIM>(origin, expnt, coeff));

    used = -wall_time();
    Function<T,NDIM> f = FunctionFactory<double,3>(world).functor(functor);
    used += wall_time();
    cout << "numerical(0.4,0.5,0.55) " << f(origin) << endl;
    cout << " analytic(0.4,0.5,0.55) " << (*functor)(origin) << endl;

    print(f.err(*functor));

//     //if (f.err(myg) > Function::defaults.thresh) error("failed");
//     cout << " used " << used << "s" << endl;
//     cout << f.norm2() << endl;
//     //f.coeff->summarize();
    
//     cout << "Test compression of a simple Gaussian using vector interface" << endl;
//     used = -wall_time();
//     Function vf = function_factory().vf(vector_myg).nocompress().initial_level(2).thresh(1e-9);
//     used += wall_time();
//     cout << "numerical(0.4,0.5,0.55) " << vf(0.4,0.5,0.55) << endl;
//     cout << " analytic(0.4,0.5,0.55) " << myg(0.4,0.5,0.55) << endl;
//     //if (vf.err(myg) > Function::defaults.thresh) error("failed");
//     cout << " used " << used << "s" << endl;
//     cout << vf.norm2() << endl;
//     //vf.coeff->summarize();
//     vf = Function();
    
//     cout << "Testing compress and reconstruct" << endl;
//     used = -wall_time();
//     f.compress();
//     used += wall_time();
//     cout << " compress used " << used << "s" << endl;
    
//     //cout << "compressed function" << endl;
//     //f.coeff->summarize();
//     used = -wall_time();
//     f.reconstruct();
//     used += wall_time();
//     cout << " reconstruct used " << used << "s" << endl;
//     //cout << "reconstructed function" << endl;
//     //f.coeff->summarize();
//     if (f.err(myg) > Function::defaults.thresh) error("failed");
    
//     cout << "Testing truncate" << endl;
//     used = -wall_time();
//     f.truncate();
//     used += wall_time();
//     cout << " truncate used " << used << "s" << endl;
//     cout << "numerical(0.4,0.5,0.55) " << f(0.4,0.5,0.55) << endl;
//     cout << " analytic(0.4,0.5,0.55) " << myg(0.4,0.5,0.55) << endl;
//     if (f.err(myg) > Function::defaults.thresh*10.0) error("failed");
//     //std::exit(0);
}

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);
        test_basic<double,3>(world);
    } catch (const MPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();

    return 0;
}
