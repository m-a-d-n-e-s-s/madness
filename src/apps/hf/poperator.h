#include <constants.h>
#include "mra/operator.h"

#define WST_PI madness::constants::pi
//#define WST_PI 3.14159265358979323846264338328

namespace madness
{
  //***************************************************************************
  const double acut1e_6 = 0.12; //0.6450626287524907;
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM> PeriodicCoulombOp(World& world, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(0.0, lo, 16.0*L[0], eps/(4.0 * WST_PI), &coeff, &expnt, false); //eps /(4.0*pi)
    coeff.scale(4.0*WST_PI);

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< SharedPtr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim[0]) && (expnt[i]*L[0]*L[0] >= acut1e_6); ++i)
    {
        coeff[i]=pow(coeff[i], 1.0/double(NDIM));
        ops.push_back(SharedPtr< Convolution1D<Q> >(new PeriodicGaussianConvolution1D<double>(k, 16, L[0]*coeff[i], expnt[i]*L[0]*L[0])));
    }

    return SeparatedConvolution<Q, NDIM>(world, k, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM>* PeriodicCoulombOpPtr(World& world, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(0.0, lo, 10.0*L[0], eps/(4.0 * WST_PI), &coeff, &expnt, false); //eps /(4.0*pi)
    coeff.scale(4.0*WST_PI);

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< SharedPtr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim[0]) && (expnt[i]*L[0]*L[0] >= acut1e_6); ++i)
    {
        coeff[i]=pow(coeff[i], 1.0/double(NDIM));
        ops.push_back(SharedPtr< Convolution1D<Q> >(new PeriodicGaussianConvolution1D<double>(k, 16, L[0]*coeff[i], expnt[i]*L[0]*L[0])));
    }

    return new SeparatedConvolution<Q, NDIM>(world, k, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM> PeriodicBSHOp(World& world, double mu, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(mu, lo, 10.0*L[0], eps, &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< SharedPtr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim[0]); ++i)
    {
      coeff[i]=pow(coeff[i], 1.0/double(NDIM));
      ops.push_back(SharedPtr< Convolution1D<Q> >(new PeriodicGaussianConvolution1D<double>(k, 16, L[0]*coeff[i], expnt[i]*L[0]*L[0])));
    }

    return SeparatedConvolution<Q, NDIM>(world, k, ops);
  }
  //***************************************************************************

  //***************************************************************************
  template<typename Q, int NDIM>
  SeparatedConvolution<Q, NDIM>* PeriodicBSHOpPtr(World& world, double mu, long k, double lo, double eps, Tensor<double> L)
  {
    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    Tensor<double> coeff, expnt;
    //if (mu==0) eps /= 4.0*pi;
    bsh_fit(mu, lo, 10.0*L[0], eps, &coeff, &expnt, false); //eps /(4.0*pi)

    // Scale coefficients according to the dimensionality and add to the list of operators
    std::vector< SharedPtr< Convolution1D<Q> > > ops;
    for (int i=0; (i < coeff.dim[0]); ++i)
    {
      coeff[i]=pow(coeff[i], 1.0/double(NDIM));
      ops.push_back(SharedPtr< Convolution1D<Q> >(new PeriodicGaussianConvolution1D<double>(k, 16, L[0]*coeff[i], expnt[i]*L[0]*L[0])));
    }

    return new SeparatedConvolution<Q, NDIM>(world, k, ops);
  }
  //***************************************************************************

};
