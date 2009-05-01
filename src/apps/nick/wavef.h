#ifndef WAVEF_H
#define WAVEF_H
/***********************************************************************
 * Here are some useful hydrogenic wave functions represented as madness 
 * functors. The bound states come from the Gnu Scientific Library. The
 * scattering states are generated with the confluent hypergeometric 
 * function. 
 * 
 * Using: Gnu Scientific Library
 *        http://www.netlib.org/toms/707
 * By:    Nick Vence
 **********************************************************************/

#include <mra/mra.h>
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
//Fortran functions
extern "C" complexd hypergf_(complexd* AA, complexd* BB, complexd* X, 
			     double* EPS, int* LIMIT, int* KIN, double* ERR, 
			     int* NITS, double* FPMAX, double* ACC8,
			     double* ACC16);

/******************************************
 * Scattering WaveFunction
 ******************************************/
class ScatteringWF : public FunctionFunctorInterface<complexd,NDIM> { 
public:
    typedef Vector<double,NDIM> vector3D;
    ScatteringWF( double Z, const vector3D& kVec );
    complexd operator()(const vector3D& x) const;
    complexd aFormNew3(complexd AA, complexd BB, complexd ZZ) const;
    complexd hypergf(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm1(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm2(complexd AA, complexd BB, complexd ZZ) const;
    complexd aForm3(complexd ZZ) const;
    complexd f11(complexd AA, complexd BB, complexd ZZ) const;
private:
    vector3D kVec;
    double   Z;
    double   k;
    double   TOL;
    double   costhK;
    complexd expmPI_k;
    complexd expPI_2k;
    complexd gamma1pI_k;
    complexd gammamI_k;
    complexd one;
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
