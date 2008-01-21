#include <unistd.h>
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

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <cstdio>

#include <tensor/random.h>

const double PI = 3.1415926535897932384;

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

template <typename T, int NDIM>
class DerivativeGaussian : public FunctionFunctorInterface<T,NDIM> {
public:
    typedef Vector<double,NDIM> coordT;
    const coordT center;
    const double exponent;
    const T coefficient;
    const int axis;
    
    DerivativeGaussian(const coordT& center, double exponent, T coefficient, int axis) 
        : center(center), exponent(exponent), coefficient(coefficient), axis(axis) 
    {};

    DerivativeGaussian(const Gaussian<T,NDIM>& g, int axis) 
        : center(g.center), exponent(g.exponent), coefficient(g.coefficient), axis(axis)
    {};

    T operator()(const coordT& x) const {
        double sum = 0.0;
        for (int i=0; i<NDIM; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return -2.0*exponent*(x[axis]-center[axis])*coefficient*exp(-exponent*sum);
    };
};


template <typename T, typename L, typename R>
inline T product(L l, R r) {
    return T(l*r);
}

template <typename T, typename L, typename R>
inline T sum(L l, R r) {
    return T(l+r);
}


/// Makes a square-normalized Gaussian with random origin and exponent
template <typename T, int NDIM>
Gaussian<T,NDIM>*
RandomGaussian(const Tensor<double> cell, double expntmax=1e5) {
    typedef Vector<double,NDIM> coordT;
    coordT origin;
    for (int i=0; i<NDIM; i++) {
        origin[i] = RandomNumber<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(0.1);
    double hi = log(expntmax);
    double expnt = exp(RandomNumber<double>()*(hi-lo) + lo);
    T coeff = pow(2.0*expnt/PI,0.25*NDIM);            
    return new Gaussian<T,NDIM>(origin,expnt,coeff);
}

/// Returns a new functor combining two functors via operation op(left,right)
template <typename resultT, typename L, typename R, typename opT, int NDIM>
class BinaryOp : public FunctionFunctorInterface<resultT,NDIM> {
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<L,NDIM> > functorL;
    typedef SharedPtr< FunctionFunctorInterface<R,NDIM> > functorR;

    functorL left;
    functorR right;
    opT op;

public:
    BinaryOp(functorL& left, functorR& right, opT& op) 
        : left(left), right(right), op(op)
    {};

    resultT operator()(const coordT& x) const {
        return op((*left)(x),(*right)(x));
    };
};


#define CHECK(value, threshold, message) \
        do { \
             if (world.rank() == 0) { \
                bool status = abs(value) < threshold; \
                const char* msgs[2] = {"FAILED","OK"}; \
                std::printf("%20.20s :%5d :%30.30s : %10.2e  < %10.2e : %s\n", \
                            __FUNCTION__,__LINE__,message,abs(value),threshold, msgs[status]); \
                if (!status) ok = false; \
             } \
        } while (0) 

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt)


template <typename T, int NDIM>
void test_basic(World& world) {
    bool ok = true;
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) 
        print("Test compression of a normalized gaussian at origin, type =",
              archive::get_type_name<T>(),", ndim =",NDIM);

    for (int i=0; i<NDIM; i++) {
        FunctionDefaults<NDIM>::cell(i,0) = -11.0-2*i;  // Deliberately asymmetric bounding box
        FunctionDefaults<NDIM>::cell(i,1) =  10.0+i;
    }
    FunctionDefaults<NDIM>::k = 7;
    FunctionDefaults<NDIM>::thresh = 1e-5;
    FunctionDefaults<NDIM>::refine = true;
    FunctionDefaults<NDIM>::initial_level = 2;
    
    const coordT origin(0.0);
    coordT point;
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);

    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    for (int i=0; i<NDIM; i++) point[i] = 0.1*i;

    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);
    double norm = f.norm2();
    double err = f.err(*functor);
    T val = f(point);
    CHECK(abs(norm-1.0), 1e-10, "norm");
    CHECK(err, 3e-7, "err");
    CHECK(val-(*functor)(point), 1e-8, "error at a point");

    f.compress();
    double new_norm = f.norm2();
    CHECK(new_norm-norm, 1e-14, "new_norm");
    
    f.reconstruct();
    new_norm = f.norm2();
    double new_err = f.err(*functor);
    CHECK(new_norm-norm, 1e-14, "new_norm");
    CHECK(new_err-err, 1e-14, "new_err");
    
    f.compress();
    new_norm = f.norm2();
    CHECK(new_norm-norm, 1e-14, "new_norm");

    f.truncate();
    new_norm = f.norm2();
    new_err = f.err(*functor);
    CHECK(new_norm-norm, 1e-9, "new_norm");
    CHECK(new_err, 3e-5, "new_err");

    //MADNESS_ASSERT(ok);
    
    world.gop.fence();
    if (world.rank() == 0) print("projection, compression, reconstruction, truncation OK\n\n");
}

template <typename T, int NDIM>
void test_conv(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("Test convergence - log(err)/(n*k) should be roughly const, a least for each value of k");
        print("                 - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    for (int i=0; i<NDIM; i++) {
        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    }

    for (int k=1; k<=15; k+=2) {
	if (world.rank() == 0) printf("k=%d\n", k);
	int ntop = 5;
	if (NDIM > 2 && k>5) ntop = 4;
	for (int n=1; n<=ntop; n++) {
	    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor).norefine().initial_level(n).k(k);
	    double err2 = f.err(*functor);
            std::size_t size = f.size();
            if (world.rank() == 0) 
                printf("   n=%d err=%.2e #coeff=%.2e log(err)/(n*k)=%.2e\n", 
                       n, err2, double(size), abs(log(err2)/n/k));
	}
    }

    world.gop.fence();
    if (world.rank() == 0) print("test conv OK\n\n");
}

template <typename T, int NDIM>
void test_math(World& world) {
    bool ok = true;
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("Test basic math operations - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }

    FunctionDefaults<NDIM>::k = 9;
    FunctionDefaults<NDIM>::thresh = 1e-9;
    FunctionDefaults<NDIM>::truncate_mode = 0;
    FunctionDefaults<NDIM>::refine = true;
    FunctionDefaults<NDIM>::autorefine = false;
    FunctionDefaults<NDIM>::initial_level = 2;
    for (int i=0; i<NDIM; i++) {
        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    }

    const coordT origin(0.0);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));
    functorT functsq(new Gaussian<T,NDIM>(origin, 2.0*expnt, coeff*coeff));
    
    // First make sure out of place squaring works
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);

    // Out-of-place squaring without autoref
    f.set_autorefine(false);  world.gop.fence();

    double err = f.err(*functor);
    CHECK(err, 1e-10, "err in f before squaring");
    Function<T,NDIM> fsq = square(f);
    double new_err = f.err(*functor);
    CHECK(new_err-err,1e-14*err,"err in f after squaring");
    double errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-10, "err in fsq");

    // Test same with autorefine
    f.set_autorefine(true); world.gop.fence();
    fsq = square(f);
    errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-10, "err in fsq with autoref");

    // Repeat after agressive truncating to see if autorefine really works
    f.set_autorefine(false); world.gop.fence();
    f.truncate(1e-5);
    f.verify_tree();
    err = f.err(*functor);
    CHECK(err, 1e-5, "error in f after truncating");

    fsq = square(f);
    errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-5, "error in fsq after truncating");

    f.set_autorefine(true); world.gop.fence();
    fsq = square(f);
    errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-5, "error in fsq truncate+autoref");

    // Finally inplace squaring
    f.square();
    double new_errsq = f.err(*functsq);
    CHECK(new_errsq - errsq, 1e-14*errsq, "err in fsq trunc+auto+inplace");

    fsq.clear();

    // Test adding a constant in scaling function and wavelet bases
    double val = f(origin);
    f.reconstruct();
    f.add_scalar(3.0);
    double val2 = f(origin);

    f.compress();
    f.add_scalar(5.0);
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-(val+8),1e-12,"add scalar in place compressed");

    // Test in-place scaling by a constant in scaling function and wavelet bases
    f.reconstruct();
    f.scale(3.0);
    val2 = f(origin);
    CHECK(val2-3.0*(val+8),1e-12,"in-place scaling reconstructed");

    f.compress();
    f.scale(4.0);
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-12.0*(val+8),1e-12,"in-place scaling compressed");

    // Same but using operator notation
    f.reconstruct();
    f *= 7.0;
    val2 = f(origin);
    CHECK(val2-7.0*12.0*(val+8),1e-11,"in-place scaling (op) recon");

    f.compress();
    f *= 7.0;
    f.reconstruct();
    val2 = f(origin);
    CHECK(val2-7.0*7.0*12.0*(val+8),1e-10,"in-place scaling (op) comp");


    // Test squaring a function by multiplication
    f = Function<T,NDIM>(FunctionFactory<T,NDIM>(world).functor(functor));
    fsq = f*f;
    errsq = fsq.err(*functsq);
    CHECK(errsq, 1e-8, "err in fsq by multiplication");

    // Test composing operations using general expression(s)
    f.compress();
    err = f.err(*functor);
    f.compress();
    Function<T,NDIM> f6 = f*3.0 + 4.0*f - f;
    new_err = f.err(*functor);
    CHECK(new_err-err,1e-14,"general op unchanged input");
    new_err = (f6 - f.scale(6.0)).norm2();
    CHECK(new_err,1e-13,"general op output");


    if (world.rank() == 0) print("\nTest multiplying random functions");
    for (int i=0; i<10; i++) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &product<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        Function<T,NDIM> c = a*b;
        a.verify_tree();
        b.verify_tree();
        c.verify_tree();
        double err1 = a.err(*f1);
        double err2 = b.err(*f2);
        double err3 = c.err(*f3);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
        CHECK(err3,1e-8,"err3");
    }      

    if (world.rank() == 0) print("\nTest adding random functions out of place");
    for (int i=0; i<10; i++) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        Function<T,NDIM> c = a+b;
        a.verify_tree();
        b.verify_tree();
        c.verify_tree();
        double err1 = a.err(*f1);
        double err2 = b.err(*f2);
        double err3 = c.err(*f3);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
        CHECK(err3,1e-8,"err3");
    }      

    if (world.rank() == 0) print("\nTest adding random functions in place");
    for (int i=0; i<10; i++) {
        functorT f1(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussian<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOp<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
        Function<T,NDIM> a = FunctionFactory<T,NDIM>(world).functor(f1);
        Function<T,NDIM> b = FunctionFactory<T,NDIM>(world).functor(f2);
        a.verify_tree();
        b.verify_tree();
        a += b;
        double err1 = a.err(*f3);
        double err2 = b.err(*f2);
        if (world.rank() == 0) print("  test ",i);
        CHECK(err1,1e-8,"err1");
        CHECK(err2,1e-8,"err2");
    }      

    // Test basic math operations

    //MADNESS_ASSERT(ok);

    world.gop.fence();
}


template <typename T, int NDIM>
void test_diff(World& world) {
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("\nTest differentiation - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.0);
    //for (int i=0; i<NDIM; i++) origin[i] = i/31.4;
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::k = 10;
    FunctionDefaults<NDIM>::thresh = 1e-10;
    FunctionDefaults<NDIM>::refine = true;
    FunctionDefaults<NDIM>::initial_level = 2;
    FunctionDefaults<NDIM>::truncate_mode = 1;
    for (int i=0; i<NDIM; i++) {
        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    }
    
    START_TIMER; 
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT

    START_TIMER;
    f.compress();
    END_TIMER("compress");

    START_TIMER;
    f.truncate();
    END_TIMER("truncate");

    START_TIMER;
    f.reconstruct();
    END_TIMER("reconstruct");

    for (int axis=0; axis<NDIM; axis++) {
        if (world.rank() == 0) print("doing axis", axis);
        DerivativeGaussian<T,NDIM> df(origin,expnt,coeff,axis);

        START_TIMER;
        Function<T,NDIM> dfdx = diff(f,axis);
        END_TIMER("diff");

//         coordT p(0.0);
//         if (world.rank() == 0) {
//             for (int i=0; i<=40; i++) {
//                 p[axis] = (i-20.0)*0.1;
//                 print("     x, analytic, err",p[axis],df(p), dfdx(p)-df(p));
//             }
//         }
//         world.gop.fence();
 
        START_TIMER;
        double err = dfdx.err(df);
        END_TIMER("err");

        if (world.rank() == 0) print("    error", err);
    }
    world.gop.fence();
}



namespace madness {
    extern bool test_rnlp();
}

template <typename T, int NDIM>
void test_op(World& world) {

    test_rnlp();

    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<T,NDIM> > functorT;

    if (world.rank() == 0) {
        print("\nTest separated operators - type =", archive::get_type_name<T>(),", ndim =",NDIM,"\n");
    }
    const coordT origin(0.5);
    const double expnt = 1.0;
    const double coeff = pow(2.0/PI,0.25*NDIM);
    functorT functor(new Gaussian<T,NDIM>(origin, expnt, coeff));

    FunctionDefaults<NDIM>::k = 10;
    FunctionDefaults<NDIM>::thresh = 1e-12;
    FunctionDefaults<NDIM>::refine = true;
    FunctionDefaults<NDIM>::initial_level = 2;
    FunctionDefaults<NDIM>::truncate_mode = 1;
    for (int i=0; i<NDIM; i++) {
        FunctionDefaults<NDIM>::cell(i,0) = -10.0;
        FunctionDefaults<NDIM>::cell(i,1) =  10.0;
    }
    
    START_TIMER; 
    Function<T,NDIM> f = FunctionFactory<T,NDIM>(world).functor(functor);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT


    f.reconstruct();
    print("         f norm is", f.norm2());
    print("     f total error", f.err(*functor));


    START_TIMER;
    f.nonstandard();
    END_TIMER("nonstandard");


    // Convolution exp(-a*x^2) with exp(-b*x^2) is
    // exp(-x^2*a*b/(a+b))* (Pi/(a+b))^(NDIM/2)

    Tensor<double> coeffs(1), exponents(1);
    exponents(0L) = 1.0;
    coeffs(0L) = pow(exponents(0L)/PI, 0.5*NDIM);
    SeparatedConvolution<T,NDIM> op(world, FunctionDefaults<NDIM>::k, coeffs, exponents);
    START_TIMER;
    Function<T,NDIM> r = apply(op,f);
    END_TIMER("apply");
    r.reconstruct();
    r.verify_tree();

    double newexpnt = expnt*exponents(0L)/(expnt+exponents(0L));
    double newcoeff = pow(PI/(expnt+exponents(0L)),0.5*NDIM)*coeff*coeffs(0L);
    functorT fexact(new Gaussian<T,NDIM>(origin, newexpnt, newcoeff));


    print(" numeric at origin", r(origin));
    print("analytic at origin", (*fexact)(origin));
    print("      op*f norm is", r.norm2());
    print("  op*f total error", r.err(*fexact));
    for (int i=0; i<=100; i++) {
        coordT c(i*0.01);
        print("           ",i,r(c),(*fexact)(c));
    }
}

/// Computes the electrostatic potential due to a Gaussian charge distribution
class GaussianPotential : public FunctionFunctorInterface<double,3> {
public:
    typedef Vector<double,3> coordT;
    const coordT center;
    const double exponent;
    const double coefficient;

    GaussianPotential(const coordT& center, double expnt, double coefficient) 
        : center(center)
        , exponent(sqrt(expnt))
        , coefficient(coefficient*pow(PI/exponent,1.5)) {}
        
    double operator()(const coordT& x) const {
        double sum = 00;
        for (int i=0; i<3; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        double r = sqrt(sum);
        return coefficient*erf(exponent*r)/r;
    }
};

void test_coulomb(World& world) {
    typedef Vector<double,3> coordT;
    typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;

    if (world.rank() == 0) {
        print("\nTest Coulomb operator - type =", archive::get_type_name<double>(),", ndim = 3 (only)\n");
    }

    // Normalized Gaussian exponent a produces potential erf(sqrt(a)*r)/r
    const coordT origin(0.5);
    const double expnt = 1.0;
    const double coeff = pow(1.0/PI*expnt,0.5*3);
    functorT functor(new Gaussian<double,3>(origin, expnt, coeff));

    double thresh = 1e-12;

    FunctionDefaults<3>::k = 12;
    FunctionDefaults<3>::thresh = thresh;
    FunctionDefaults<3>::refine = true;
    FunctionDefaults<3>::initial_level = 2;
    FunctionDefaults<3>::truncate_mode = 0;
    for (int i=0; i<3; i++) {
        FunctionDefaults<3>::cell(i,0) = -10.0;
        FunctionDefaults<3>::cell(i,1) =  10.0;
    }
    
    START_TIMER; 
    Function<double,3> f = FunctionFactory<double,3>(world).functor(functor).thresh(1e-14);
    END_TIMER("project");

    //f.print_info();  <--------- This is not scalable and might crash the XT

    f.reconstruct();
    double norm = f.norm2(), err = f.err(*functor);
    if (world.rank() == 0) {
        print("         f norm is", norm);
        print("     f total error", err);
        //print(" truncating");
    }

//     START_TIMER;
//     f.truncate();
//     END_TIMER("truncate");
//     START_TIMER;
//     f.reconstruct();
//     END_TIMER("reconstruct");
//     norm = f.norm2();
//     err = f.err(*functor);
//     if (world.rank() == 0) {
//         print("         f norm is", norm);
//         print("     f total error", err);
//     }

    f.reconstruct();
    START_TIMER;
    f.nonstandard();
    END_TIMER("nonstandard");

    SeparatedConvolution<double,3> op = CoulombOperator<double,3>(world, FunctionDefaults<3>::k, 1e-12, thresh);
    START_TIMER;
    Function<double,3> r = apply_only(op,f);
    END_TIMER("apply");
    const double* nflop = op.get_nflop();
    if (world.rank() == 0) {
        double totalflops = 0.0;
        for (int i=0; i<64; i++) {
            if (nflop[i]) print(i,nflop[i]);
            totalflops += nflop[i];
        }
        print("total flops", totalflops);
    }
    START_TIMER;
    r.reconstruct();
    END_TIMER("reconstruct result");
    r.verify_tree();

    functorT fexact(new GaussianPotential(origin, expnt, coeff));

    double numeric=r(origin);
    double analytic=(*fexact)(origin);
    double rnorm = r.norm2();
    double rerr = r.err(*fexact);
    if (world.rank() == 0) {
        print(" numeric at origin", numeric);
        print("analytic at origin", analytic);
        print("      op*f norm is", rnorm);
        print("  op*f total error", rerr);
//         for (int i=0; i<=100; i++) {
//             coordT c(i*0.01);
//             print("           ",i,r(c),(*fexact)(c));
//         }
    }
}


#define TO_STRING(s) TO_STRING2(s)
#define TO_STRING2(s) #s

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    if (world.rank() == 0) {
        print("");
        print("--------------------------------------------");
        print("   MADNESS",MADNESS_VERSION, "multiresolution testsuite");
        print("--------------------------------------------");
        print("");
        print("   number of processors ...", world.size());
        print("    processor frequency ...", cpu_frequency());
        print("            host system ...", TO_STRING(HOST_SYSTEM));
        print("             byte order ...", TO_STRING(MADNESS_BYTE_ORDER));
        print("          configured by ...", MADNESS_CONFIGURATION_USER);
        print("          configured on ...", MADNESS_CONFIGURATION_HOST);
        print("          configured at ...", MADNESS_CONFIGURATION_DATE);
        print("                    CXX ...", MADNESS_CONFIGURATION_CXX);
        print("               CXXFLAGS ...", MADNESS_CONFIGURATION_CXXFLAGS);
#ifdef WORLD_WATCHDOG
        print("               watchdog ...", WATCHDOG_BARK_INTERVAL, WATCHDOG_TIMEOUT);
#endif
#ifdef OPTERON_TUNE
        print("             tuning for ...", "opteron");
#elif defined(CORE_DUO_TUNE)
        print("             tuning for ...", "core duo");
#else
        print("             tuning for ...", "core2");
#endif
#ifdef BOUNDS_CHECKING
        print(" tensor bounds checking ...", "enabled");
#endif
#ifdef TENSOR_INSTANCE_COUNT
        print("  tensor instance count ...", "enabled");
#endif
        print(" ");
    }        

    try {
        startup(world,argc,argv);
        if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
//           test_basic<double,1>(world);
//           test_conv<double,1>(world);
//           test_math<double,1>(world);
//           test_diff<double,1>(world);
//          test_op<double,1>(world);

//          test_basic<double,2>(world);
//          test_conv<double,2>(world);
//          test_math<double,2>(world);
//          test_diff<double,2>(world);
//          test_op<double,2>(world);

//         test_basic<double,3>(world);
//         test_conv<double,3>(world);
//         test_math<double,3>(world);
//         test_diff<double,3>(world);
         //test_op<double,3>(world);
         test_coulomb(world);

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

    if (world.rank() == 0) print("entering final fence");
    world.gop.fence();
    if (world.rank() == 0) print("done with final fence");
    if (world.rank() == 0) print("Final tensor instance count", BaseTensor::get_instance_count());
    MPI::Finalize();

    return 0;
}

