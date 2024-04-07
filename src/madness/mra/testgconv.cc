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

/// \file mra/testgconv.cc
/// \brief Test convolution with Gaussian * polyn

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>

using namespace madness;

static const int k = 10;
static const double thresh = 1e-6;
static const double L = 17;

// exp(-r^2) / sqrt(pi) = normalized gaussian
double g(const coord_1d& r) {
    static const double fac = 1.0/sqrt(constants::pi);
    return exp(-r[0]*r[0]) * fac;
}


// g' = exp(-r^2) / sqrt(pi) = derivative of a normalized gaussian
double gprime(const coord_1d& r) {
    static const double fac = 1.0/sqrt(constants::pi);
    return -2.* r[0] * exp(-r[0]*r[0]) * fac;
}

// conv g()
double convg(const coord_1d& r) {
    static const double fac = 1.0/sqrt(2.0*constants::pi);
    return exp(-r[0]*r[0]*0.5) * fac;
}

// sqrt(8)*x*exp(-x^2)
double h(const coord_1d& r) {
    static const double fac = sqrt(8.0);
    return exp(-r[0]*r[0]) * fac * r[0];
}

// g() conv h() == h() conv g()
double gconvh(const coord_1d& r) {
    return exp(-0.5*r[0]*r[0]) * r[0];
}

// D(conv) (g())
double conv_prime_g(const coord_1d& r) {
    static const double fac = 1.0/sqrt(2.0*constants::pi);
    return -exp(-0.5*r[0]*r[0]) * r[0] *fac;
}


int test_gconv(World& world) {
    coord_1d origin(0.0), lo(-L), hi(L);
    double width = 2.0*L;
    int success=0;

    if (world.rank() == 0) print("Test gconv operation");

    const real_function_1d f = real_factory_1d(world).f(g);
    double error=f.trace()-1.0;
    print("error in integral(g) ", error);
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 0 ", success);

    // convolve with a normalized Gaussian kernel
    std::vector< std::shared_ptr< Convolution1D<double> > > ops(1);
    ops[0].reset(new GaussianConvolution1D<double>(k, width/sqrt(constants::pi),
            width*width, 0, false));
    real_convolution_1d op(world, ops);

    real_function_1d opf = op(f);
    error=opf.trace()-1.0;
    print("error in integral(op(g)) ", error);
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 1 ", success);

    real_function_1d exact = real_factory_1d(world).f(convg);
    print("norm2(g conv g - exact)", (opf-exact).norm2());
    error=(opf-exact).norm2();
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 2 ", success);

    real_function_1d q = real_factory_1d(world).f(h);
    error=q.trace();
    print("error in integral(h) ", error);
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 3 ", success);

    error=q.norm2() - sqrt(sqrt(2.0*constants::pi));
    print("error in norm2(h)", error);
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 4 ", success);

    real_function_1d opq = op(q);
    exact = real_factory_1d(world).f(gconvh);
    error=(opq-exact).norm2();
    print("norm2(g conv h - exact)", error);
    if (error>FunctionDefaults<1>::get_thresh()) success++;
    print("success 5 ", success);



    // test the convolution with a derivative Gaussian:
    // result(y) = \int g'(x-y) f(x) = \int g(x-y) f'(x)
    // where we use
    // f(x)      = exp(-x^2) / sqrt(pi)
    // f'(x)     = -2 x exp(-x^2) / sqrt(pi)
    // result(y) = -y exp(-y/2) / sqrt(2 pi)

    // the derivative Gaussian convolution kernel:
    // note the scaling of the coeffs because the derivative operator brings
    // down the scaling factor of the exponent
    ops[0].reset(new GaussianConvolution1D<double>(k, 1.0/sqrt(constants::pi),
            width*width, 1, false));

    real_convolution_1d oph(world, ops);

    // this is the result hardwired
    const real_function_1d convpg=real_factory_1d(world).f(conv_prime_g);

    // apply the derivative Gaussian on f
    opq = oph(f);

    // apply the Gaussian on the derivative of f
    const real_function_1d fp=real_factory_1d(world).f(gprime);
    real_function_1d opfp=op(fp);

    // the error
    const double error1=(opq-convpg).norm2();
    const double error2=(opfp-convpg).norm2();
    print("norm2(conv' g - exact)", error1);
    print("norm2(conv g' - exact)", error2);
    if (error1>FunctionDefaults<1>::get_thresh()) success++;
    print("success 6a ", success);
    if (error2>FunctionDefaults<1>::get_thresh()) success++;
    print("success 6b ", success);

    plot_line("opf.dat", 1001, lo, hi, q, opq, convpg);

    world.gop.fence();
    return success;
}


int main(int argc, char**argv) {
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);

    int success=0;
    try {
        startup(world,argc,argv);

        FunctionDefaults<1>::set_cubic_cell(-L,L);
        FunctionDefaults<1>::set_k(k);
        FunctionDefaults<1>::set_thresh(thresh);
        FunctionDefaults<1>::set_truncate_mode(1);
        FunctionDefaults<1>::set_initial_level(5);
        if (world.rank()==0) {
        	print(" threshold  ", thresh);
        	print(" polynomial ", k,"\n");
        }
        success+=test_gconv(world);

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

