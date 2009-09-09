#ifndef WAVEF_H
#define WAVEF_H
//\file wavef.cc
//\brief The hydrogenic bound and continuum states
/************************************************************************
 * Here is a madness representation of the hydrogenic wave functions.
 * The bound states come from the Gnu Scientific Library. The unbound
 * states are generated with the confluent hypergeometric function which
 * uses gmp and mpfr for extended precision
 * 
 * Using: Gnu Scientific Library          http://www.gnu.org/software/gsl/
 *        GNU Multiple Precision library  http://gmplib.org/
 *        mpfr                            http://www.mpfr.org/
 * By:    Nick Vence
 ************************************************************************/
#include <mra/mra.h>
#include "interp.h"
#include <complex>
#include <iostream>
#include <stdio.h>
using std::printf;

using namespace madness;

const int NDIM  = 3;
typedef std::complex<double> complexd;
typedef SharedPtr< FunctionFunctorInterface<complexd,NDIM> > functorT;
typedef Vector<double,NDIM> vector3D;
const complexd I(0,1);
const double PI = M_PI;
const complexd one(1,0);

int      fact(int);
void     testFact(World&);
complexd gamma(complexd AA);
complexd gamma(double re, double im);
void     testGamma(World&);
complexd pochhammer(complexd AA,int n);
void     testPochhammer(World&);
void     testWF(World&);
void     debug1F1(World&);
void test1F1(World&, complexd (*func1F1)(complexd,complexd,complexd), const char* fileChar);

/******************************************
 * Scattering WaveFunction
 ******************************************/
class ScatteringWF : public FunctionFunctorInterface<complexd,NDIM> { 
public:
    typedef Vector<double,NDIM> vector3D;
    ScatteringWF( double Z, const vector3D& kVec );
    ~ScatteringWF();
    complexd operator()(const vector3D& x) const;
    complexd aFormNew3(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm1(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm2(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm3(complexd ZZ) const;
    complexd f11(double x) const;
    complexd approx1F1(double x) const;
    double   diffR(double x) const;
    double   diffI(double x) const;
    double   toX(int i) const;
    int      fromX(double x) const;
    void     taylor3Fit();
    void     splineFit();
    void     polyFit();
    double   toX(double s) const;
    double   toS(double x) const;
    double   r1F1(double x) const;
    double   i1F1(double x) const;
    int      n1;
    double   Z;
    vector3D kVec;
    double   k;
    double   domain;
    double   xi;
    double   ra;
private:
    CubicInterpolationTable<complex<double> > fit1F1;
    complexd expmPI_k;
    complexd expPI_2k;
    complexd expPI_2kXgamma1pI_k;
    complexd gamma1pI_k;
    complexd gammamI_k;
    complexd one;
    double   alpha;
    double   beta;
    double   boundary;
    double   TOL;
    double   dx;
    double   dx_2n1;
    double   dxn1_2;
    int n;
// splineFit
//     double*  aR;
//     double*  bR;
//     double*  cR;
//     double*  dR;
//     double*  aI;
//     double*  bI;
//     double*  cI;
//     double*  dI;
    double*   fR;
    double*   fI;
    double*  fpR;
    double*  fpI;
    double* fppR;
    double* fppI;
    double* fpppR;
    double* fpppI;
    double*    h;
    double*    x;
    struct MemFuncPtr {
        ScatteringWF* obj;
        MemFuncPtr(ScatteringWF* obj) : obj(obj) {}
        complex<double> operator()(double x) {return obj->f11(x);}
    };
};

/******************************************
 * Bound WaveFunction
 ******************************************/
class BoundWF : public FunctionFunctorInterface<complexd,NDIM> {
public:
    typedef Vector<double,NDIM> vector3D;
    BoundWF(double Z, int nn, int ll, int mm );
    complexd operator()(const vector3D& x) const;
private:
    double Z;
    int n;
    int l;
    int m;
};

/******************************************
 *Exp[ I*(k.r) ]
 ******************************************/
class Expikr : public FunctionFunctorInterface<complexd,NDIM>
{
public:
    typedef Vector<double,NDIM> vector3D;
    Expikr(const vector3D& kVec);
    complexd operator()(const vector3D& r) const;
private:
    vector3D kVec;
    double k;
    double costhK;
};

/******************************************
 *Exp[ -I*(kr + k.r) ]
 ******************************************/
class Expikr2 : public FunctionFunctorInterface<complexd,NDIM>
{
public:
    typedef Vector<double,NDIM> vector3D;
    Expikr2(const vector3D& kVec);
    complexd operator()(const vector3D& r) const;
private:
    vector3D kVec;
    double k;
    double costhK;
};
#endif
