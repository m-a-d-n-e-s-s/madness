#ifndef WAVEF_H
#define WAVEF_H
/**********************
 * By:Nick Vence
 * The Coulomb WaveFunction
 * ********************/

#include <mra/mra.h>
#include <complex>
#include <iostream>
#include <stdio.h>
using std::printf;
#include <string>
using std::string;

using namespace madness;

const int NDIM  = 3;
typedef std::complex<double> complexd;
typedef SharedPtr< FunctionFunctorInterface<complexd,NDIM> > functorT;
typedef Vector<double,NDIM> vector3D;

extern "C" complexd hypergf_(complexd* AA, complexd* BB, complexd* X, 
			     double* EPS, int* LIMIT, int* KIN, double* ERR, 
			     int* NITS, double* FPMAX, double* ACC8,
			     double* ACC16);
extern "C" complexd conhyp_(complexd* AA, complexd* BB, complexd* ZZ, 
			    int* LNCHF, int* IP);

int fact(int);
complexd gamma(complexd AA);
complexd gamma(double re, double im);
complexd pochhammer(complexd AA,int n);
complexd hypergf(complexd AA, complexd BB, complexd ZZ);
complexd conHyp(complexd AA, complexd BB, complexd ZZ);
complexd aForm(complexd AA, complexd BB, complexd ZZ);
complexd f11(complexd AA, complexd BB, complexd ZZ);
void test1F1(complexd (*func1F1)(complexd,complexd,complexd), char* fName);
void f11Tester(World&);

const complexd I(0,1);
const double PI = M_PI;

/******************************************
 * Virtual class for all wave functions 
 ******************************************/
class WaveFunction : public FunctionFunctorInterface<complexd,NDIM>
{
 public:
  typedef Vector<double,NDIM> vector3D;
  WaveFunction(double M, double Z);
  virtual complexd operator()(const vector3D& x) const=0; 
 protected:
  double _M;
  double _Z;
};


/******************************************
 * Scattering WaveFunction
 ******************************************/
class ScatteringWF : public WaveFunction
{ 
 public:
  typedef Vector<double,NDIM> vector3D;
  ScatteringWF(double M, double Z, const vector3D& kVec );
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
  BoundWF(double M, double Z, int nn, int ll, int mm );
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
