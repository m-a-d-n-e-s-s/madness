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
typedef Vector<double,NDIM> coordT;

extern "C" complexd hypergf_(complexd* AA, complexd* BB, complexd* X, 
			     double* EPS, int* LIMIT, int* KIN, double* ERR, 
			     int* NITS, double* FPMAX, double* ACC8,
			     double* ACC16);
extern "C" complexd conhyp_(complexd* AA, complexd* BB, complexd* ZZ, 
			    int* LNCHF, int* IP);

int fact(int);
complexd gamma(complexd AA);
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
  typedef Vector<double,NDIM> coordT;
  WaveFunction(double M, double Z);
  virtual complexd operator()(const coordT& x) const=0; 
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
  typedef Vector<double,NDIM> coordT;
  ScatteringWF(double M, double Z, const coordT& kVec );
  virtual complexd operator()(const coordT& x) const;
 private:
  coordT kVec;
  double k;
  double costhK;
};

/******************************************
 * Bound WaveFunction
 ******************************************/
class BoundWF : public WaveFunction
{
 public:
  typedef Vector<double,NDIM> coordT;
  BoundWF(double M, double Z, int nn, int ll, int mm );
  virtual complexd operator()(const coordT& x) const;
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
  typedef Vector<double,NDIM> coordT;
  Expikr(const coordT& kVec);
  complexd operator()(const coordT& r) const;
 private:
  coordT kVec;
  double k;
  double costhK;
};

#endif
