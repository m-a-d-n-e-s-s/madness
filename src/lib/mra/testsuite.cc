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
#include <cstdio>

#include <tensor/random.h>

const double PI = 3.1415926535897932384;

using namespace madness;

template <typename T, int NDIM>
class GaussianFunctor : public FunctionFunctorInterface<T,NDIM> {
private:
    typedef Vector<double,NDIM> coordT;
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
GaussianFunctor<T,NDIM>*
RandomGaussianFunctor(const Tensor<double> cell, double expntmax=1e5) {
    typedef Vector<double,NDIM> coordT;
    coordT origin;
    for (int i=0; i<NDIM; i++) {
        origin[i] = RandomNumber<double>()*(cell(i,1)-cell(i,0)) + cell(i,0);
    }
    double lo = log(0.1);
    double hi = log(expntmax);
    double expnt = exp(RandomNumber<double>()*(hi-lo) + lo);
    T coeff = pow(2.0*expnt/PI,0.25*NDIM);            
    return new GaussianFunctor<T,NDIM>(origin,expnt,coeff);
}

/// Returns a new functor combining two functors via operation op(left,right)
template <typename resultT, typename L, typename R, typename opT, int NDIM>
class BinaryOpFunctor : public FunctionFunctorInterface<resultT,NDIM> {
    typedef Vector<double,NDIM> coordT;
    typedef SharedPtr< FunctionFunctorInterface<L,NDIM> > functorL;
    typedef SharedPtr< FunctionFunctorInterface<R,NDIM> > functorR;

    functorL left;
    functorR right;
    opT op;

public:
    BinaryOpFunctor(functorL& left, functorR& right, opT& op) 
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

    functorT functor(new GaussianFunctor<T,NDIM>(origin, expnt, coeff));

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
    functorT functor(new GaussianFunctor<T,NDIM>(origin, expnt, coeff));

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
    functorT functor(new GaussianFunctor<T,NDIM>(origin, expnt, coeff));
    functorT functsq(new GaussianFunctor<T,NDIM>(origin, 2.0*expnt, coeff*coeff));
    
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
        functorT f1(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &product<T,T,T>;
        functorT f3(new BinaryOpFunctor<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
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
        functorT f1(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOpFunctor<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
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
        functorT f1(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        functorT f2(RandomGaussianFunctor<T,NDIM>(FunctionDefaults<NDIM>::cell,100.0));
        T (*p)(T,T) = &sum<T,T,T>;
        functorT f3(new BinaryOpFunctor<T,T,T,T(*)(T,T),NDIM>(f1,f2,p));
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


int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);

    try {
        startup(world,argc,argv);
        if (world.rank() == 0) print("Initial tensor instance count", BaseTensor::get_instance_count());
        test_basic<double,1>(world);
        test_conv<double,1>(world);
        test_math<double,1>(world);

        test_basic<double,2>(world);
        test_conv<double,2>(world);
        test_math<double,2>(world);

        test_basic<double,3>(world);
        test_conv<double,3>(world);
        test_math<double,3>(world);

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
