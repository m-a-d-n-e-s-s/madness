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

complexd hypergf(complexd AA, complexd BB, complexd ZZ);
complexd conHyp(complexd AA, complexd BB, complexd ZZ);
complexd aForm(complexd AA, complexd BB, complexd ZZ);
complexd f11(complexd AA, complexd BB, complexd ZZ);
void     test1F1(World&, complexd (*func1F1)(complexd,complexd,complexd), char* fName);
int      fact(int);
void     testFact(World&);
complexd gamma(complexd AA);
complexd gamma(double re, double im);
void     testGamma(World&);
complexd pochhammer(complexd AA,int n);
void     testPochhammer(World&);
void f11Tester(World&);

//Fortran functions
extern "C" complexd hypergf_(complexd* AA, complexd* BB, complexd* X, 
			     double* EPS, int* LIMIT, int* KIN, double* ERR, 
			     int* NITS, double* FPMAX, double* ACC8,
			     double* ACC16);
extern "C" complexd conhyp_(complexd* AA, complexd* BB, complexd* ZZ, 
			    int* LNCHF, int* IP);

/******************************************
 * Virtual class for all wave functions 
 ******************************************/
class WaveFunction : public FunctionFunctorInterface<complexd,NDIM>
{
 public:
  typedef Vector<double,NDIM> vector3D;
  WaveFunction(double Z);
  virtual complexd operator()(const vector3D& x) const=0; 
 protected:
  double Z;
};


/******************************************
 * Scattering WaveFunction
 ******************************************/
class ScatteringWF : public WaveFunction
{ 
 public:
  typedef Vector<double,NDIM> vector3D;
    ScatteringWF( double Z, const vector3D& kVec );
  virtual complexd operator()(const vector3D& x) const;
 private:
  vector3D kVec;
  double k;
  double costhK;
};

/******************************************
 * Bound WaveFunction
 ******************************************/
class BoundWF : public WaveFunction
{
 public:
  typedef Vector<double,NDIM> vector3D;
  BoundWF(double Z, int nn, int ll, int mm );
  virtual complexd operator()(const vector3D& x) const;
 private:
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

#endif
