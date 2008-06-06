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

/// \file testbsh.cc
/// \brief test the bsh operator

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>

using namespace madness;

template <typename T, int NDIM>
class Gaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;
    
    Gaussian(const coordT& center, double exponent, T coefficient) 
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

double q (double r)
{
    const double pi = constants::pi;
    const double fac = pow(2.0/pi,0.25*3)/(4.0*pi);
    if (r < 1e-4) {
        return fac*(2.854819526231167-1.618591848021335*r*r);
    }
    else {
        return fac*pow(pi, 0.3e1 / 0.2e1) * (-exp((double) (2 * r)) + 0.1e1 + erf((double) r - 0.1e1 / 0.2e1) + erf((double) r + 0.1e1 / 0.2e1) * exp((double) (2 * r))) * exp(0.1e1 / 0.4e1 - (double) r) / (double) r / 0.2e1;

    }

}

struct Qfunc : public FunctionFunctorInterface<double,3> {
    double operator()(const Vector<double,3>& x) const {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        return q(r);
    }
};

template <typename T>
void test_bsh(World& world) {
    double mu = 1.0;
    vector<long> npt(3,201);
    typedef Vector<double,3> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,3> > functorT;

    int nn = 1001;
    double lo = -2.0, hi=2.0, range=hi-lo;

    if (world.rank() == 0) 
        print("Test BSH operation, type =",
              archive::get_type_name<T>(),", ndim =",3);

    FunctionDefaults<3>::set_cubic_cell(-20,20);
    FunctionDefaults<3>::set_k(9);
    FunctionDefaults<3>::set_thresh(1e-7);
    FunctionDefaults<3>::set_initial_level(2);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_autorefine(true);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_on_project(false);
    
    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(2.0/constants::pi,0.25*3);

    Function<T,3> f = FunctionFactory<T,3>(world).functor(functorT(new Gaussian<T,3>(origin, expnt, coeff)));
    f.truncate();
    f.reconstruct();
    print("before",f.size());
    f.refine();
    print("after ",f.size());

//     for (int i=0; i<nn; i++) {
//         double z=lo + i*range/double(nn-1);
//         coordT p(z);
//         double  exact = Gaussian<T,3>(origin, expnt, coeff)(p);
//         print(z, f(p), exact, f(p)-exact);
//     }

//     plotdx(f, "f.dx", FunctionDefaults<3>::get_cell(), npt);
//     for (int i=0; i<101; i++) {
//         double z=-4 + 0.08*i;
//         coordT p(z);
//         double exact = Gaussian<T,3>(origin, expnt, coeff)(p);
//         print(z, f(p), f(p)-exact);
//     }

    double norm = f.norm2();
    print("norm of initial function", norm, f.err(Gaussian<T,3>(origin, expnt, coeff)));

    SeparatedConvolution<T,3> op = BSHOperator<T,3>(world, 
                                                    mu, 
                                                    FunctionDefaults<3>::get_k(), 
                                                    1e-5, 
                                                    1e-12);
    //op.doleaves = true;

    print("applying - 1");
    Function<T,3> opf = apply(op,f);
    print("done");
    print("err in opf", opf.err(Qfunc()));

//     Function<double,3> qf = FunctionFactory<T,3>(world).functor(functorT(new Qfunc()));
//     Function<double,3> xerror = qf-opf;
//     //xerror.truncate();
//     xerror.reconstruct();
//     xerror.print_tree();
//     //return;

//     opf.reconstruct();
//     coordT pt(0.1);
//     print("compare");
//     print(q(0.1*sqrt(3.0)), opf(pt), opf(pt)/q(0.1*sqrt(3.0)));
//     for (int i=0; i<nn; i++) {
//         double z=lo + i*range/double(nn-1);
//         double r = fabs(z)*sqrt(3.0);
//         coordT p(z);
//         print(z, opf(p), q(r), q(r)/opf(p), q(r)-opf(p));
//     }

//     plotdx(opf, "opf.dx", FunctionDefaults<3>::get_cell(), npt);

    Function<T,3> opinvopf = opf*(mu*mu);
    for (int axis=0; axis<3; axis++) {
        print("diffing",axis);
        opinvopf = opinvopf - diff(diff(opf,axis),axis);
    }

//     plotdx(opinvopf, "opinvopf.dx", FunctionDefaults<3>::get_cell(), npt);

    print("norm of (-del^2+mu^2)*G*f", opinvopf.norm2());
    Function<T,3> error = (f-opinvopf);
    
//     error.reconstruct();
//     plotdx(error, "err.dx", FunctionDefaults<3>::get_cell(), npt);

    print("error",error.norm2());

    opinvopf.reconstruct();
    f.reconstruct();
    error.reconstruct();

//     for (int i=0; i<101; i++) {
//         double z=-4 + 0.08*i;
//         coordT p(z);
//         print(z, opinvopf(p), f(p), opinvopf(p)/f(p), error(p));
//     }

    opf.clear(); opinvopf.clear();

    Function<T,3> g = (mu*mu)*f;
    for (int axis=0; axis<3; axis++) {
        g = g - diff(diff(f,axis),axis);
    }
    g = apply(op,g);
    print("norm of G*(-del^2+mu^2)*f",g.norm2());
    print("error",(g-f).norm2());

    world.gop.fence();

}

template <typename T, int NDIM>
void test_op(World& world) {

    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("\nTest separated operators - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.5);
    const double expnt = 4096.0; //1.0;
    const double coeff = pow(2.0*expnt/constants::pi,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::set_k(16);
    FunctionDefaults<NDIM>::set_thresh(1e-13);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(10);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    FunctionDefaults<NDIM>::set_cubic_cell(-100.0,100.0);
    
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    f.truncate();
    f.reconstruct();

//     print("F INITIAL PROJECTION");
//     f.print_tree();
//     print("F COMPRESSED");
//     f.compress();
//     f.print_tree();
//     print("F RECONSTRUCTED");
//     f.reconstruct();
//     f.print_tree();
//     print("F NONSTANDARD");
//     f.nonstandard();
//     f.print_tree();

//     f.standard();

    //f.truncate();
//     print("fsize before refine1 ", f.size());
//     f.refine();
//     print("fsize before refine2 ", f.size());
//     f.refine();
//     print("fsize  after refine2 ", f.size());
//     f.refine();
//     print("fsize  after refine3 ", f.size());

    double n2 = f.norm2();
    double t2 = f.trace();
    double e2 = f.err(*functor);
    if (world.rank() == 0) {
        print("         f  norm is", n2);
        print("         f trace is", t2);
        print("      f total error", e2);
    }

    // Convolution exp(-a*x^2) with exp(-b*x^2) is
    // exp(-x^2*a*b/(a+b))* (Pi/(a+b))^(NDIM/2)

    Tensor<double> coeffs(1), exponents(1);
    exponents(0L) = 1.0; //4096.0; //0.015625;
    //    while (exponents(0L) < 1.1e6) {
        coeffs(0L) = pow(exponents(0L)/constants::pi, 0.5*NDIM);

        SeparatedConvolution<T,NDIM> op(world, FunctionDefaults<NDIM>::get_k(), coeffs, exponents);
        //op.doleaves = true;
        //GenericConvolution1D<double,GaussianGenericFunctor<double> > gen(FunctionDefaults<NDIM>::get_k(),
        //                                                                         GaussianGenericFunctor<double>(coeffs[0L],exponents[0L]));


        Function<T,NDIM> r = apply(op,f);
        r.verify_tree();
        f.verify_tree();

        double newexpnt = expnt*exponents(0L)/(expnt+exponents(0L));
        double newcoeff = pow(constants::pi/(expnt+exponents(0L)),0.5*NDIM)*coeff*coeffs(0L);
        functorT fexact(new Gaussian<T,NDIM>(origin, newexpnt, newcoeff));

        double rn = r.norm2();
        double re = r.err(*fexact);
        if (world.rank() == 0) print("exponent",exponents[0L],"norm",rn,"err",re);

        exponents[0L] *= 2.0;
        //    }
}


int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    try {
        startup(world,argc,argv);
        
        //test_op<double,1>(world);
        test_bsh<double>(world);

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

