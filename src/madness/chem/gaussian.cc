/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file Basis/gaussian.cc
 * \brief Implementation of Gaussian basis function API and routines.
 */

#include "../../madness.h"
#include "gaussian.h"
#include <cmath>
#include <limits>
#include <stdexcept>

const static double pi = 4.*atan(1.);

namespace slymer {

///////////////////////////////////////////////////////////////////////////
// CartesianPrimitive implementations
///////////////////////////////////////////////////////////////////////////
PrimitiveGaussian::PrimitiveGaussian(const GaussianType &type,
    const std::array<double, 3> &center, const double ec)
    : prefactor(0), exppoly(2) {

  // make sure the decay constant is positive
  if(ec <= 0.)
    throw std::invalid_argument("Gaussian decay constant must be positive.");

  // set the exponent array; expand out (x-x0)^2 and similarly for y, z
  exppoly[{{2,0,0}}] = -ec;
  exppoly[{{1,0,0}}] = 2.*ec*center[0];
  exppoly[{{0,0,0}}] = -ec*center[0]*center[0];
  exppoly[{{0,2,0}}] = -ec;
  exppoly[{{0,1,0}}] = 2.*ec*center[1];
  exppoly[{{0,0,0}}] += -ec*center[1]*center[1];
  exppoly[{{0,0,2}}] = -ec;
  exppoly[{{0,0,1}}] = 2.*ec*center[2];
  exppoly[{{0,0,0}}] += -ec*center[2]*center[2];

  if(type == GaussianType::s) {
    prefactor[{{0,0,0}}] = pow(2.*ec/pi, 0.75);
  }
  else if(type == GaussianType::px) {
    prefactor = PolynomialCoeffs(1);
    prefactor[{{1,0,0}}] = 1.;
    prefactor[{{0,0,0}}] = -center[0];
    prefactor *= pow(2.*ec/pi, 0.5) * 2.*pow(2.*ec*ec*ec/pi, 0.25);
  }
  else if(type == GaussianType::py) {
    prefactor = PolynomialCoeffs(1);
    prefactor[{{0,1,0}}] = 1.;
    prefactor[{{0,0,0}}] = -center[1];
    prefactor *= pow(2.*ec/pi, 0.5) * 2.*pow(2.*ec*ec*ec/pi, 0.25);
  }
  else if(type == GaussianType::pz) {
    prefactor = PolynomialCoeffs(1);
    prefactor[{{0,0,1}}] = 1.;
    prefactor[{{0,0,0}}] = -center[2];
    prefactor *= pow(2.*ec/pi, 0.5) * 2.*pow(2.*ec*ec*ec/pi, 0.25);
  }

  // Cartesian and shared d-types
  else if(type == GaussianType::dxx) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{2,0,0}}] = 1.;
    prefactor[{{1,0,0}}] = -2. * center[0];
    prefactor[{{0,0,0}}] = center[0] * center[0];
    prefactor *= pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::dxy) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{1,1,0}}] = 1.;
    prefactor[{{1,0,0}}] = -center[1];
    prefactor[{{0,1,0}}] = -center[0];
    prefactor[{{0,0,0}}] = center[0] * center[1];
    prefactor *= pow(2.*ec/pi, 0.25) * 4.*pow(2.*ec*ec*ec/pi, 0.5);
  }
  else if(type == GaussianType::dxz) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{1,0,1}}] = 1.;
    prefactor[{{1,0,0}}] = -center[2];
    prefactor[{{0,0,1}}] = -center[0];
    prefactor[{{0,0,0}}] = center[0] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 4.*pow(2.*ec*ec*ec/pi, 0.5);
  }
  else if(type == GaussianType::dyy) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{0,2,0}}] = 1.;
    prefactor[{{0,1,0}}] = -2. * center[1];
    prefactor[{{0,0,0}}] = center[1] * center[1];
    prefactor *= pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::dyz) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{0,1,1}}] = 1.;
    prefactor[{{0,1,0}}] = -center[2];
    prefactor[{{0,0,1}}] = -center[1];
    prefactor[{{0,0,0}}] = center[1] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 4.*pow(2.*ec*ec*ec/pi, 0.5);
  }
  else if(type == GaussianType::dzz) {
    prefactor = PolynomialCoeffs(2);
    prefactor[{{0,0,2}}] = 1.;
    prefactor[{{0,0,1}}] = -2. * center[2];
    prefactor[{{0,0,0}}] = center[2] * center[2];
    prefactor *= pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }

  else if(type == GaussianType::fxxx) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{3,0,0}}] = 1.;
    prefactor[{{2,0,0}}] = -3. * center[0];
    prefactor[{{1,0,0}}] = 3. * center[0] * center[0];
    prefactor[{{0,0,0}}] = -center[0] * center[0] * center[0];
    prefactor *= pow(2.*ec/pi, 0.5) * 8.*ec*pow(2.*ec*ec*ec/pi, 0.25)/sqrt(15.);
  }
  else if(type == GaussianType::fxxy) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{2,1,0}}] = 1.;
    prefactor[{{2,0,0}}] = -center[1];
    prefactor[{{1,1,0}}] = -2. * center[0];
    prefactor[{{1,0,0}}] = 2. * center[0] * center[1];
    prefactor[{{0,1,0}}] = center[0] * center[0];
    prefactor[{{0,0,0}}] = -center[0] * center[0] * center[1];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fxxz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{2,0,1}}] = 1.;
    prefactor[{{2,0,0}}] = -center[2];
    prefactor[{{1,0,1}}] = -2. * center[0];
    prefactor[{{1,0,0}}] = 2. * center[0] * center[2];
    prefactor[{{0,0,1}}] = center[0] * center[0];
    prefactor[{{0,0,0}}] = -center[0] * center[0] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fxyy) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{1,2,0}}] = 1.;
    prefactor[{{0,2,0}}] = -center[0];
    prefactor[{{1,1,0}}] = -2. * center[1];
    prefactor[{{0,1,0}}] = 2. * center[0] * center[1];
    prefactor[{{1,0,0}}] = center[1] * center[1];
    prefactor[{{0,0,0}}] = -center[0] * center[1] * center[1];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fxyz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{1,1,1}}] = 1.;
    prefactor[{{1,1,0}}] = -center[2];
    prefactor[{{1,0,1}}] = -center[1];
    prefactor[{{0,1,1}}] = -center[0];
    prefactor[{{1,0,0}}] = center[1] * center[2];
    prefactor[{{0,1,0}}] = center[0] * center[2];
    prefactor[{{0,0,1}}] = center[0] * center[1];
    prefactor[{{0,0,0}}] = -center[0] * center[1] * center[2];
    prefactor *= 8.*pow(2.*ec*ec*ec/pi, 0.75);
  }
  else if(type == GaussianType::fxzz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{1,0,2}}] = 1.;
    prefactor[{{0,0,2}}] = -center[0];
    prefactor[{{1,0,1}}] = -2. * center[2];
    prefactor[{{0,0,1}}] = 2. * center[0] * center[2];
    prefactor[{{1,0,0}}] = center[2] * center[2];
    prefactor[{{0,0,0}}] = -center[0] * center[2] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fyyy) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{0,3,0}}] = 1.;
    prefactor[{{0,2,0}}] = -3. * center[1];
    prefactor[{{0,1,0}}] = 3. * center[1] * center[1];
    prefactor[{{0,0,0}}] = -center[1] * center[1] * center[1];
    prefactor *= pow(2.*ec/pi, 0.5) * 8.*ec*pow(2.*ec*ec*ec/pi, 0.25)/sqrt(15.);
  }
  else if(type == GaussianType::fyyz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{0,2,1}}] = 1.;
    prefactor[{{0,2,0}}] = -center[2];
    prefactor[{{0,1,1}}] = -2. * center[1];
    prefactor[{{0,1,0}}] = 2. * center[1] * center[2];
    prefactor[{{0,0,1}}] = center[1] * center[1];
    prefactor[{{0,0,0}}] = -center[1] * center[1] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fyzz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{0,1,2}}] = 1.;
    prefactor[{{0,0,2}}] = -center[1];
    prefactor[{{0,1,1}}] = -2. * center[2];
    prefactor[{{0,0,1}}] = 2. * center[1] * center[2];
    prefactor[{{0,1,0}}] = center[2] * center[2];
    prefactor[{{0,0,0}}] = -center[1] * center[2] * center[2];
    prefactor *= pow(2.*ec/pi, 0.25) * 2.*pow(2.*ec*ec*ec/pi, 0.25) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
  }
  else if(type == GaussianType::fzzz) {
    prefactor = PolynomialCoeffs(3);
    prefactor[{{0,0,3}}] = 1.;
    prefactor[{{0,0,2}}] = -3. * center[2];
    prefactor[{{0,0,1}}] = 3. * center[2] * center[2];
    prefactor[{{0,0,0}}] = -center[2] * center[2] * center[2];
    prefactor *= pow(2.*ec/pi, 0.5) * 8.*ec*pow(2.*ec*ec*ec/pi, 0.25)/sqrt(15.);
  }

  // spherical d orbitals
  else if(type == GaussianType::dzzmrr) {
    // setup as dzz - dxx/2 - dyy/2
    PolynomialCoeffs predxx(2), predyy(2), predzz(2);
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= -pow(2.*ec/pi, 0.5) * 2.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 2.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
    predzz[{{0,0,2}}] = 1.;
    predzz[{{0,0,1}}] = -2. * center[2];
    predzz[{{0,0,0}}] = center[2] * center[2];
    predzz *= pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.);
    prefactor = predxx + predyy + predzz;
  }
  else if(type == GaussianType::dxxmyy) {
    // setup as sqrt(3)*dxx/2 - sqrt(3)*dyy/2
    PolynomialCoeffs predxx(2), predyy(2);
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= pow(2.*ec/pi, 0.5) * 2.*ec*pow(2.*ec/pi, 0.25);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 2.*ec*pow(2.*ec/pi, 0.25);
    prefactor = predxx + predyy;
  }

  // spherical f orbitals
  else if(type == GaussianType::fxyymxxx) {
    // setup dxx and dyy, get (3dyy - dxx), normalize, and multiply by x
    PolynomialCoeffs xfact(1), predxx(2), predyy(2);
    xfact[{{1,0,0}}] = 1.;
    xfact[{{0,0,0}}] = -center[0];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)/sqrt(3.) * sqrt(0.5*ec);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= sqrt(3.) * pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25) * sqrt(0.5*ec);
    prefactor = xfact * (predxx + predyy);
  }
  else if(type == GaussianType::fxxzmyyz) {
    // setup dxx and dyy, get (dxx - dyy), normalize, and multiply by z
    PolynomialCoeffs zfact(1), predxx(2), predyy(2);
    zfact[{{0,0,1}}] = 1.;
    zfact[{{0,0,0}}] = -center[2];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec);
    prefactor = zfact * (predxx + predyy);
  }
  else if(type == GaussianType::fxzzmrrx) {
    // setup dxx, dyy, and dzz; get (4dzz-dxx-dyy); normalize; and multiply by x
    PolynomialCoeffs xfact(1), predxx(2), predyy(2), predzz(2);
    xfact[{{1,0,0}}] = 1.;
    xfact[{{0,0,0}}] = -center[0];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    predzz[{{0,0,2}}] = 1.;
    predzz[{{0,0,1}}] = -2. * center[2];
    predzz[{{0,0,0}}] = center[2] * center[2];
    predzz *= 4.*pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    prefactor = xfact * (predxx + predyy + predzz);
  }
  else if(type == GaussianType::fzzzmrrz) {
    // setup dxx, dyy, and dzz; get (2dzz-3dxx-3dyy); normalize; and multiply by z
    PolynomialCoeffs zfact(1), predxx(2), predyy(2), predzz(2);
    zfact[{{0,0,1}}] = 1.;
    zfact[{{0,0,0}}] = -center[2];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= -pow(2.*ec/pi, 0.5) * 12.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec/15.);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 12.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec/15.);
    predzz[{{0,0,2}}] = 1.;
    predzz[{{0,0,1}}] = -2. * center[2];
    predzz[{{0,0,0}}] = center[2] * center[2];
    predzz *= pow(2.*ec/pi, 0.5) * 8.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec/15.);
    prefactor = zfact * (predxx + predyy + predzz);
  }
  else if(type == GaussianType::fyzzmrry) {
    // setup dxx, dyy, and dzz; get (4dzz-dxx-dyy); normalize; and multiply by y
    PolynomialCoeffs yfact(1), predxx(2), predyy(2), predzz(2);
    yfact[{{0,1,0}}] = 1.;
    yfact[{{0,0,0}}] = -center[1];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    predzz[{{0,0,2}}] = 1.;
    predzz[{{0,0,1}}] = -2. * center[2];
    predzz[{{0,0,0}}] = center[2] * center[2];
    predzz *= pow(2.*ec/pi, 0.5) * 16.*ec*pow(2.*ec/pi, 0.25)*sqrt(0.1*ec);
    prefactor = yfact * (predxx + predyy + predzz);
  }
  else if(type == GaussianType::fxxymyyy) {
    // setup dxx and dyy, get (3dxx-dyy), normalize, and multiply by y
    PolynomialCoeffs yfact(1), predxx(2), predyy(2);
    yfact[{{0,1,0}}] = 1.;
    yfact[{{0,0,0}}] = -center[1];
    predxx[{{2,0,0}}] = 1.;
    predxx[{{1,0,0}}] = -2. * center[0];
    predxx[{{0,0,0}}] = center[0] * center[0];
    predxx *= pow(2.*ec/pi, 0.5) * 12.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec/6.);
    predyy[{{0,2,0}}] = 1.;
    predyy[{{0,1,0}}] = -2. * center[1];
    predyy[{{0,0,0}}] = center[1] * center[1];
    predyy *= -pow(2.*ec/pi, 0.5) * 4.*ec*pow(2.*ec/pi, 0.25)*sqrt(ec/6);
    prefactor = yfact * (predxx + predyy);
  }

  // Cartesian g orbitals
  else if(type == GaussianType::gxxxx) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{4,0,0}}] =  1.;
    prefactor[{{3,0,0}}] = -4. * center[0];
    prefactor[{{2,0,0}}] =  6. * center[0] * center[0];
    prefactor[{{1,0,0}}] = -4. * pow(center[0], 3); 
    prefactor[{{0,0,0}}] =  pow(center[0], 4);
    prefactor *= 16./sqrt(105.) * pow(ec, 2.75) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::gxxxy) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{3,1,0}}] =  1.;
    prefactor[{{2,1,0}}] = -3. * center[0];
    prefactor[{{1,1,0}}] =  3. * center[0] * center[0];
    prefactor[{{0,1,0}}] = -pow(center[0], 3);
    prefactor[{{3,0,0}}] = -center[1];
    prefactor[{{2,0,0}}] =  3. * center[0] * center[1];
    prefactor[{{1,0,0}}] = -3. * center[0] * center[0] * center[1];
    prefactor[{{0,0,0}}] =  pow(center[0], 3) * center[1];
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxxxz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{3,0,1}}] =  1.;
    prefactor[{{2,0,1}}] = -3. * center[0];
    prefactor[{{1,0,1}}] =  3. * center[0] * center[0];
    prefactor[{{0,0,1}}] = -pow(center[0], 3);
    prefactor[{{3,0,0}}] = -center[2];
    prefactor[{{2,0,0}}] =  3. * center[0] * center[2];
    prefactor[{{1,0,0}}] = -3. * center[0] * center[0] * center[2];
    prefactor[{{0,0,0}}] =  pow(center[0], 3) * center[2];
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxxyy) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{2,2,0}}] =  1.;
    prefactor[{{1,2,0}}] = -2. * center[0];
    prefactor[{{0,2,0}}] =  center[0] * center[0];
    prefactor[{{2,1,0}}] = -2. * center[1];
    prefactor[{{1,1,0}}] =  4. * center[0] * center[1];
    prefactor[{{0,1,0}}] = -2. * center[0] * center[0] * center[1];
    prefactor[{{2,0,0}}] =  center[1] * center[1];
    prefactor[{{1,0,0}}] = -2. * center[0] * center[1] * center[1];
    prefactor[{{0,0,0}}] =  center[0] * center[0] * center[1] * center[1];
    prefactor *= 16./3. * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxxyz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{2,1,1}}] =  1.;
    prefactor[{{1,1,1}}] = -2. * center[0];
    prefactor[{{0,1,1}}] =  center[0] * center[0];
    prefactor[{{2,0,1}}] = -center[1];
    prefactor[{{1,0,1}}] =  2. * center[0] * center[1];
    prefactor[{{0,0,1}}] = -center[0] * center[0] * center[1];
    prefactor[{{2,1,0}}] = -center[2];
    prefactor[{{1,1,0}}] =  2. * center[0] * center[2];
    prefactor[{{0,1,0}}] = -center[0] * center[0] * center[2];
    prefactor[{{2,0,0}}] =  center[1] * center[2];
    prefactor[{{1,0,0}}] = -2. * center[0] * center[1] * center[2];
    prefactor[{{0,0,0}}] =  center[0] * center[0] * center[1] * center[2];
    prefactor *= 16./sqrt(3.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxxzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{2,0,2}}] =  1.;
    prefactor[{{1,0,2}}] = -2. * center[0];
    prefactor[{{0,0,2}}] =  center[0] * center[0];
    prefactor[{{2,0,1}}] = -2. * center[2];
    prefactor[{{1,0,1}}] =  4. * center[0] * center[2];
    prefactor[{{0,0,1}}] = -2. * center[0] * center[0] * center[2];
    prefactor[{{2,0,0}}] =  center[2] * center[2];
    prefactor[{{1,0,0}}] = -2. * center[0] * center[2] * center[2];
    prefactor[{{0,0,0}}] =  center[0] * center[0] * center[2] * center[2];
    prefactor *= 16./3. * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxyyy) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{1,3,0}}] =  1.;
    prefactor[{{0,3,0}}] = -center[0];
    prefactor[{{1,2,0}}] = -3. * center[1];
    prefactor[{{0,2,0}}] =  3. * center[0] * center[1];
    prefactor[{{1,1,0}}] =  3. * center[1] * center[1];
    prefactor[{{0,1,0}}] = -3. * center[0] * center[1] * center[1];
    prefactor[{{1,0,0}}] = -pow(center[1], 3);
    prefactor[{{0,0,0}}] =  center[0] * pow(center[1], 3);
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxyyz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{1,2,1}}] =  1.;
    prefactor[{{0,2,1}}] = -center[0];
    prefactor[{{1,1,1}}] = -2. * center[1];
    prefactor[{{0,1,1}}] =  2. * center[0] * center[1];
    prefactor[{{1,0,1}}] =  center[1] * center[1];
    prefactor[{{0,0,1}}] = -center[0] * center[1] * center[1];
    prefactor[{{1,2,0}}] = -center[2];
    prefactor[{{0,2,0}}] =  center[0] * center[2];
    prefactor[{{1,1,0}}] =  2. * center[1] * center[2];
    prefactor[{{0,1,0}}] = -2. * center[0] * center[1] * center[2];
    prefactor[{{1,0,0}}] = -center[1] * center[1] * center[2];
    prefactor[{{0,0,0}}] =  center[0] * center[1] * center[1] * center[2];
    prefactor *= 16./sqrt(3.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxyzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{1,1,2}}] =  1.;
    prefactor[{{0,1,2}}] = -center[0]; 
    prefactor[{{1,0,2}}] = -center[1];
    prefactor[{{0,0,2}}] =  center[0] * center[1];
    prefactor[{{1,1,1}}] = -2. * center[2];
    prefactor[{{0,1,1}}] =  2. * center[0] * center[2];
    prefactor[{{1,0,1}}] =  2. * center[1] * center[2];
    prefactor[{{0,0,1}}] = -2. * center[0] * center[1] * center[2];
    prefactor[{{1,1,0}}] =  center[2] * center[2];
    prefactor[{{0,1,0}}] = -center[0] * center[2] * center[2];
    prefactor[{{1,0,0}}] = -center[1] * center[2] * center[2];
    prefactor[{{0,0,0}}] =  center[0] * center[1] * center[2] * center[2];
    prefactor *= 16./sqrt(3.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gxzzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{1,0,3}}] =  1.;
    prefactor[{{0,0,3}}] = -center[0];
    prefactor[{{1,0,2}}] = -3. * center[2];
    prefactor[{{0,0,2}}] =  3. * center[0] * center[2];
    prefactor[{{1,0,1}}] =  3. * center[2] * center[2];
    prefactor[{{0,0,1}}] = -3. * center[0] * center[2] * center[2];
    prefactor[{{1,0,0}}] = -pow(center[2], 3);
    prefactor[{{0,0,0}}] =  center[0] * pow(center[2], 3);
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gyyyy) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{0,4,0}}] = 1.;
    prefactor[{{0,3,0}}] = -4. * center[1];
    prefactor[{{0,2,0}}] =  6. * center[1] * center[1];
    prefactor[{{0,1,0}}] = -4. * pow(center[1], 3); 
    prefactor[{{0,0,0}}] = + pow(center[1], 4);
    prefactor *= 16./sqrt(105.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gyyyz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{0,3,1}}] =  1.;
    prefactor[{{0,2,1}}] = -3. * center[1];
    prefactor[{{0,1,1}}] =  3. * center[1] * center[1];
    prefactor[{{0,0,1}}] = -pow(center[1], 3);
    prefactor[{{0,3,0}}] = -center[2];
    prefactor[{{0,2,0}}] =  3. * center[1] * center[2];
    prefactor[{{0,1,0}}] = -3. * center[1] * center[1] * center[2];
    prefactor[{{0,0,0}}] =  pow(center[1], 3) * center[2];
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gyyzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{0,2,2}}] =  1.;
    prefactor[{{0,1,2}}] = -2. * center[1];
    prefactor[{{0,0,2}}] =  center[1] * center[1];
    prefactor[{{0,2,1}}] = -2. * center[2];
    prefactor[{{0,1,1}}] =  4. * center[1] * center[2];
    prefactor[{{0,0,1}}] = -2. * center[1] * center[1] * center[2];
    prefactor[{{0,2,0}}] =  center[2] * center[2];
    prefactor[{{0,1,0}}] = -2. * center[1] * center[2] * center[2];
    prefactor[{{0,0,0}}] =  center[1] * center[1] * center[2] * center[2];
    prefactor *= 16./3. * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gyzzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{0,1,3}}] =  1.;
    prefactor[{{0,0,3}}] = -center[1];
    prefactor[{{0,1,2}}] = -3. * center[2];
    prefactor[{{0,0,2}}] =  3. * center[1] * center[2];
    prefactor[{{0,1,1}}] =  3. * center[2] * center[2];
    prefactor[{{0,0,1}}] = -3. * center[1] * center[2] * center[2];
    prefactor[{{0,1,0}}] = -pow(center[2], 3);
    prefactor[{{0,0,0}}] =  center[1] * pow(center[2], 3);
    prefactor *= 16./sqrt(15.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::gzzzz) {
    prefactor = PolynomialCoeffs(4);
    prefactor[{{0,0,4}}] = 1.;
    prefactor[{{0,0,3}}] = -4. * center[2];
    prefactor[{{0,0,2}}] =  6. * center[2] * center[2];
    prefactor[{{0,0,1}}] = -4. * pow(center[2], 3); 
    prefactor[{{0,0,0}}] = pow(center[2], 4);
    prefactor *= 16./sqrt(105.) * pow(ec, 2.75) * pow(2./pi, 0.75);
  }

  // spherical g orbitals
  // (l,m) = (4,-4) 
  else if(type == GaussianType::gxydx2my2 ) {
    // get (xx-yy), multiply by xy, normalize
    PolynomialCoeffs xx(2), yy(2), x(1), y(1);
    x[{{1,0,0}}] = 1.;
    x[{{0,0,0}}] = -center[0];
    y[{{0,1,0}}] = 1.;
    y[{{0,0,0}}] = -center[1];
    xx[{{2,0,0}}] = 1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] = center[0] * center[0];
    yy[{{0,2,0}}] = 1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] = center[1] * center[1];
    double norm = 8./sqrt(3.) * pow(ec, 2.75) * pow(2./pi, 0.75);
    prefactor = norm * x * y * (xx + (-1. * yy));
  }
  // (l,m) = (4,-3)
  else if(type == GaussianType::gyzdx2my2) {
    // get (3xx - yy), multiply by y*z, normalize
    PolynomialCoeffs xx(2), yy(2), y(1), z(1);
    y[{{0,1,0}}] = 1.;
    y[{{0,0,0}}] = -center[1];
    z[{{0,0,1}}] = 1.;
    z[{{0,0,0}}] = -center[2];
    xx[{{2,0,0}}] = 1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] = center[0] * center[0];
    yy[{{0,2,0}}] = 1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] = center[1] * center[1];
    double norm = 8./sqrt(3.) * pow(ec, 2.75) * pow(2., 0.25) * pow(1./pi, 0.75);
    prefactor = norm * y * z * (3. * xx + (-1. * yy));
  }
  // (l,m) = (4,-2)
  else if(type == GaussianType::gxydz2mr2) {
    // get (7zz - rr), multiply by x*y, normalize
    PolynomialCoeffs zz(2), rr(2), x(1), y(1);
    x[{{1,0,0}}] = 1.;
    x[{{0,0,0}}] = -center[0];
    y[{{0,1,0}}] = 1.;
    y[{{0,0,0}}] = -center[1];
    zz[{{0,0,2}}] = 1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] = center[2] * center[2];    
    rr[{{2,0,0}}] = 1.;
    rr[{{1,0,0}}] = -2. * center[0];
    rr[{{0,2,0}}] = 1.;
    rr[{{0,1,0}}] = -2. * center[1];
    rr[{{0,0,2}}] = 1.;
    rr[{{0,0,1}}] = -2. * center[2];
    rr[{{0,0,0}}] = center[0] * center[0] + center[1] * center[1] + center[2] * center[2];
    double norm = 8./sqrt(21.) * pow(ec, 2.75) * pow(2./pi, 0.75);
    prefactor = norm * x * y * (7. * zz + (-1. * rr));
  }
  // (l,m) = (4,-1)
  else if(type == GaussianType::gyzdz2mr2) {
    // get (7zz - 3rr), multiply by y*z, normalize
    PolynomialCoeffs zz(2), rr(2), y(1), z(1);
    y[{{0,1,0}}] = 1.;
    y[{{0,0,0}}] = -center[1];
    z[{{0,0,1}}] = 1.;
    z[{{0,0,0}}] = -center[2];
    zz[{{0,0,2}}] = 1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] = center[2] * center[2];
    rr[{{2,0,0}}] = 1.;
    rr[{{1,0,0}}] = -2. * center[0];
    rr[{{0,2,0}}] = 1.;
    rr[{{0,1,0}}] = -2. * center[1];
    rr[{{0,0,2}}] = 1.;
    rr[{{0,0,1}}] = -2. * center[2];
    rr[{{0,0,0}}] = center[0] * center[0] + center[1] * center[1] + center[2] * center[2];
    double norm = 8./sqrt(21.) * pow(2., 0.25) * pow(ec, 2.75) * pow(1./pi, 0.75);
    prefactor = norm * y * z * (7. * zz + (-3. * rr));
  }
  // (l,m) = (4,0)
  else if(type == GaussianType::gzero) {
    // we want: 3x^4 + 6x^2y^2 + 3y^4 - 24x^2z^2 - 24y^2z^2 + 8z^4
    PolynomialCoeffs x4(4), y4(4), z4(4), xx(2), yy(2), zz(2);
    x4[{{4,0,0}}] =  1.;
    x4[{{3,0,0}}] = -4. * center[0];
    x4[{{2,0,0}}] =  6. * pow(center[0], 2);
    x4[{{1,0,0}}] = -4. * pow(center[0], 3);
    x4[{{0,0,0}}] =  pow(center[0], 4);
    y4[{{0,4,0}}] =  1.;
    y4[{{0,3,0}}] = -4. * center[1];
    y4[{{0,2,0}}] =  6. * pow(center[1], 2);
    y4[{{0,1,0}}] = -4. * pow(center[1], 3);
    y4[{{0,0,0}}] =  pow(center[1], 4);
    z4[{{0,0,4}}] =  1.;
    z4[{{0,0,3}}] = -4. * center[2];
    z4[{{0,0,2}}] =  6. * pow(center[2], 2);
    z4[{{0,0,1}}] = -4. * pow(center[2], 3);
    z4[{{0,0,0}}] =  pow(center[2], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  pow(center[0], 2);
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  pow(center[1], 2);
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  pow(center[2], 2);
    double norm = 2./sqrt(105.) * pow(ec, 2.75) * pow(2./pi, 0.75); 
    prefactor = norm * (3.*x4 + 3.*y4 + 8.*z4 + 6.*xx*yy + (-24.*xx*zz) + (-24.*yy*zz)); 
  }
  // (l,m) = (4,1)
  else if(type == GaussianType::gxzdz2mr2) {
    // get (7zz - 3rr), multiply by x*z, normalize
    PolynomialCoeffs zz(2), rr(2), x(1), z(1);
    x[{{1,0,0}}] = 1.;
    x[{{0,0,0}}] = -center[0];
    z[{{0,0,1}}] = 1.;
    z[{{0,0,0}}] = -center[2];
    zz[{{0,0,2}}] = 1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] = center[2] * center[2];    
    rr[{{2,0,0}}] = 1.;
    rr[{{1,0,0}}] = -2. * center[0];
    rr[{{0,2,0}}] = 1.;
    rr[{{0,1,0}}] = -2. * center[1];
    rr[{{0,0,2}}] = 1.;
    rr[{{0,0,1}}] = -2. * center[2];
    rr[{{0,0,0}}] = center[0] * center[0] + center[1] * center[1] + center[2] * center[2];
    double norm = 8./sqrt(21.) * pow(2, 0.25) * pow(ec, 2.75) * pow(1./pi, 0.75);
    // Testing shows nwchem uses the negative
    prefactor = (-1.) * norm * x * z * (7.*zz + (-3.*rr)); 
  }
  // (l,m) = (4,2)
  else if(type == GaussianType::gx2my2dz2mr2) {
    // get (7zz - rr) and (xx - yy), multiply together, normalize
    PolynomialCoeffs xx(2), yy(2), zz(2), rr(2);
    xx[{{2,0,0}}] = 1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] = center[0] * center[0];
    yy[{{0,2,0}}] = 1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] = center[1] * center[1];
    zz[{{0,0,2}}] = 1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] = center[2] * center[2]; 
    rr[{{2,0,0}}] = 1.;
    rr[{{1,0,0}}] = -2. * center[0];
    rr[{{0,2,0}}] = 1.;
    rr[{{0,1,0}}] = -2. * center[1];
    rr[{{0,0,2}}] = 1.;
    rr[{{0,0,1}}] = -2. * center[2];
    rr[{{0,0,0}}] = center[0] * center[0] + center[1] * center[1] + center[2] * center[2];
    double norm = 4./sqrt(21.) * pow(ec, 2.75) * pow(2./pi, 0.75);
    prefactor = norm * (xx + (-1*yy)) * (7.*zz + (-1.*rr));  
  }
  // (l,m) = (4,3)
  else if(type == GaussianType::gxzdx2my2) {
    // get (xx - 3yy), multiply by xz, normalize
    PolynomialCoeffs xx(2), yy(2), x(1), z(1);
    x[{{1,0,0}}] = 1.;
    x[{{0,0,0}}] = -center[0];
    z[{{0,0,1}}] = 1.;
    z[{{0,0,0}}] = -center[2];
    xx[{{2,0,0}}] = 1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] = center[0] * center[0];
    yy[{{0,2,0}}] = 1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] = center[1] * center[1];
    double norm = 8./sqrt(3.) * pow(2., 0.25) * pow(ec, 2.75) * pow(1./pi, 0.75);
    // Testing shows nwchem uses the negative
    prefactor = norm * x * z * (3.*yy + -1.*xx); 
  }
  // (l,m) = (4,4)
  else if(type == GaussianType::gx4mx2y2py4) {
    // we want: x^4 - 6x^2y^2 + y^4, then normalize 
    PolynomialCoeffs x4(4), y4(4), xx(2), yy(2);
    x4[{{4,0,0}}] =  1.;
    x4[{{3,0,0}}] = -4. * center[0];
    x4[{{2,0,0}}] =  6. * center[0] * center[0];
    x4[{{1,0,0}}] = -4. * pow(center[0], 3);
    x4[{{0,0,0}}] =  pow(center[0], 4);
    y4[{{0,4,0}}] =  1.;
    y4[{{0,3,0}}] = -4 * center[1];
    y4[{{0,2,0}}] =  6. * center[1] * center[1];
    y4[{{0,1,0}}] = -4. * pow(center[1], 3);
    y4[{{0,0,0}}] =  pow(center[1], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.; 
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    double norm = 2./sqrt(3.) * pow(ec, 2.75) * pow(2./pi, 0.75);     
    prefactor = norm * (x4 + (-6.*xx*yy) + y4); 
  }

  // Cartesian h orbitals
  else if(type == GaussianType::hxxxxx) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{5,0,0}}] =  1.;
    prefactor[{{4,0,0}}] = -5. * center[0];
    prefactor[{{3,0,0}}] =  10. * center[0] * center[0];
    prefactor[{{2,0,0}}] = -10. * pow(center[0], 3);
    prefactor[{{1,0,0}}] =  5. * pow(center[0], 4); 
    prefactor[{{0,0,0}}] = -pow(center[0], 5);
    prefactor *= 32./(3.*sqrt(105.)) * pow(ec, 3.25) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::hxxxxy) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{4,1,0}}] =  1.;
    prefactor[{{3,1,0}}] = -4. * center[0];
    prefactor[{{2,1,0}}] =  6. * center[0] * center[0];
    prefactor[{{1,1,0}}] = -4. * pow(center[0], 3);
    prefactor[{{0,1,0}}] =  pow(center[0], 4);
    prefactor[{{4,0,0}}] = -1. * center[1];
    prefactor[{{3,0,0}}] =  4. * center[0] * center[1];
    prefactor[{{2,0,0}}] = -6. * center[0] * center[0] * center[1];
    prefactor[{{1,0,0}}] =  4. * pow(center[0], 3) * center[1];
    prefactor[{{0,0,0}}] = -pow(center[0], 4) * center[1];
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxxxz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{4,0,1}}] =  1.;
    prefactor[{{3,0,1}}] = -4. * center[0];
    prefactor[{{2,0,1}}] =  6. * center[0] * center[0];
    prefactor[{{1,0,1}}] = -4. * pow(center[0], 3);
    prefactor[{{0,0,1}}] =  pow(center[0], 4);
    prefactor[{{4,0,0}}] = -center[2];
    prefactor[{{3,0,0}}] =  4. * center[0] * center[2];
    prefactor[{{2,0,0}}] = -6. * center[0] * center[0] * center[2];
    prefactor[{{1,0,0}}] =  4. * pow(center[0], 3) * center[2];
    prefactor[{{0,0,0}}] = -pow(center[0], 4) * center[2];
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75);    
  }
  else if(type == GaussianType::hxxxyy) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{3,2,0}}] =  1.;
    prefactor[{{2,2,0}}] = -3. * center[0];
    prefactor[{{1,2,0}}] =  3. * center[0] * center[0];
    prefactor[{{0,2,0}}] = -pow(center[0], 3);
    prefactor[{{3,1,0}}] = -2. * center[1];
    prefactor[{{2,1,0}}] =  6. * center[0] * center[1];
    prefactor[{{1,1,0}}] = -6. * center[0] * center[0] * center[1];
    prefactor[{{0,1,0}}] =  2. * pow(center[0], 3) * center[1];
    prefactor[{{3,0,0}}] =  center[1] * center[1];
    prefactor[{{2,0,0}}] = -3. * center[0] * center[1] * center[1];
    prefactor[{{1,0,0}}] =  3. * center[0] * center[0] * center[1] * center[1];
    prefactor[{{0,0,0}}] = -pow(center[0], 3) * center[1] * center[1];
    prefactor *= 32./(3.*sqrt(5.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxxyz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{3,1,1}}] =  1.;
    prefactor[{{2,1,1}}] = -3. * center[0];
    prefactor[{{1,1,1}}] =  3. * center[0] * center[0];
    prefactor[{{0,1,1}}] = -pow(center[0], 3);
    prefactor[{{3,0,1}}] = -center[1];
    prefactor[{{2,0,1}}] =  3. * center[0] * center[1];
    prefactor[{{1,0,1}}] = -3. * center[0] * center[0] * center[1];
    prefactor[{{0,0,1}}] =  pow(center[0], 3) * center[1];
    prefactor[{{3,1,0}}] = -center[2];
    prefactor[{{2,1,0}}] =  3. * center[0] * center[2];
    prefactor[{{1,1,0}}] = -3. * center[0] * center[0] * center[2];
    prefactor[{{0,1,0}}] =  pow(center[0], 3) * center[2];
    prefactor[{{3,0,0}}] =  center[1] * center[2];
    prefactor[{{2,0,0}}] = -3. * center[0] * center[1] * center[2];
    prefactor[{{1,0,0}}] =  3. * center[0] * center[0] * center[1] * center[2];
    prefactor[{{0,0,0}}] = -pow(center[0], 3) * center[1] * center[2];
    prefactor *= 32./sqrt(15.) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxxzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{3,0,2}}] =  1.;
    prefactor[{{2,0,2}}] = -3. * center[0];
    prefactor[{{1,0,2}}] =  3. * center[0] * center[0];
    prefactor[{{0,0,2}}] = -pow(center[0], 3);
    prefactor[{{3,0,1}}] = -2. * center[2];
    prefactor[{{2,0,1}}] =  6. * center[0] * center[2];
    prefactor[{{1,0,1}}] = -6. * center[0] * center[0] * center[2];
    prefactor[{{0,0,1}}] =  2. * pow(center[0], 3) * center[2];
    prefactor[{{3,0,0}}] =  center[2] * center[2];
    prefactor[{{2,0,0}}] = -3. * center[0] * center[2] * center[2];
    prefactor[{{1,0,0}}] =  3. * center[0] * center[0] * center[2] * center[2];
    prefactor[{{0,0,0}}] = -pow(center[0], 3) * center[2] * center[2];
    prefactor *= 32./(3.*sqrt(5.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxyyy) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{2,3,0}}] =  1.;
    prefactor[{{1,3,0}}] = -2. * center[0];
    prefactor[{{0,3,0}}] =  center[0] * center[0];
    prefactor[{{2,2,0}}] = -3. * center[1];
    prefactor[{{1,2,0}}] =  6. * center[0] * center[1];
    prefactor[{{0,2,0}}] = -3. * center[0] * center[0] * center[1];
    prefactor[{{2,1,0}}] =  3. * center[1] * center[1];
    prefactor[{{1,1,0}}] = -6. * center[0] * center[1] * center[1];
    prefactor[{{0,1,0}}] =  3. * center[0] * center[0] * center[1] * center[1];
    prefactor[{{2,0,0}}] = -pow(center[1], 3);
    prefactor[{{1,0,0}}] =  2. * center[0] * pow(center[1], 3);
    prefactor[{{0,0,0}}] = -center[0] * center[0] * pow(center[1], 3);
    prefactor *= 32./(3.*sqrt(5)) * pow(ec, 13./4) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxyyz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{2,2,1}}] =  1.;
    prefactor[{{1,2,1}}] = -2. * center[0];
    prefactor[{{0,2,1}}] =  center[0] * center[0];
    prefactor[{{2,1,1}}] = -2. * center[1];
    prefactor[{{1,1,1}}] =  4. * center[0] * center[1];
    prefactor[{{0,1,1}}] = -2. * center[0] * center[0] * center[1];
    prefactor[{{2,0,1}}] =  center[1] * center[1];
    prefactor[{{1,0,1}}] = -2. * center[0] * center[1] * center[1];
    prefactor[{{0,0,1}}] =  center[0] * center[0] * center[1] * center[1];
    prefactor[{{2,2,0}}] = -center[2];
    prefactor[{{1,2,0}}] =  2. * center[0] * center[2];
    prefactor[{{0,2,0}}] = -center[0] * center[0] * center[2];
    prefactor[{{2,1,0}}] =  2. * center[1] * center[2];
    prefactor[{{1,1,0}}] = -4. * center[0] * center[1] * center[2];
    prefactor[{{0,1,0}}] =  2. * center[0] * center[0] * center[1] * center[2];
    prefactor[{{2,0,0}}] = -center[1] * center[1] * center[2];
    prefactor[{{1,0,0}}] =  2. * center[0] * center[1] * center[1] * center[2];
    prefactor[{{0,0,0}}] = -center[0] * center[0] * center[1] * center[1] * center[2];
    prefactor *= 32./3. * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxxyzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{2,1,2}}] =  1.;
    prefactor[{{1,1,2}}] = -2. * center[0];
    prefactor[{{0,1,2}}] =  center[0] * center[0];
    prefactor[{{2,0,2}}] = -center[1];
    prefactor[{{1,0,2}}] =  2. * center[0] * center[1];
    prefactor[{{0,0,2}}] = -center[0] * center[0] * center[1];
    prefactor[{{2,1,1}}] = -2. * center[2];
    prefactor[{{1,1,1}}] =  4. * center[0] * center[2];
    prefactor[{{0,1,1}}] = -2. * center[0] * center[0] * center[2];
    prefactor[{{2,0,1}}] =  2. * center[1] * center[2];
    prefactor[{{1,0,1}}] = -4. * center[0] * center[1] * center[2];
    prefactor[{{0,0,1}}] =  2. * center[0] * center[0] * center[1] * center[2];
    prefactor[{{2,1,0}}] =  center[2] * center[2];
    prefactor[{{1,1,0}}] = -2. * center[0] * pow(center[2], 2);
    prefactor[{{0,1,0}}] =  center[0] * center[0] * center[2] * center[2];
    prefactor[{{2,0,0}}] = -center[1] * center[2] * center[2];
    prefactor[{{1,0,0}}] =  2. * center[0] * center[1] * center[2] * center[2];
    prefactor[{{0,0,0}}] = -center[0] * center[0] * center[1] * center[2] * center[2];
    prefactor *= 32./3. * pow(ec, 3.25) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::hxxzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{2,0,3}}] =  1.;
    prefactor[{{1,0,3}}] = -2. * center[0];
    prefactor[{{0,0,3}}] =  center[0] * center[0];
    prefactor[{{2,0,2}}] = -3. * center[2];
    prefactor[{{1,0,2}}] =  6. * center[0] * center[2];
    prefactor[{{0,0,2}}] = -3. * center[0] * center[0] * center[2];
    prefactor[{{2,0,1}}] =  3. * center[2] * center[2];
    prefactor[{{1,0,1}}] = -6. * center[0] * center[2] * center[2];
    prefactor[{{0,0,1}}] =  3. * center[0] * center[0] * center[2] * center[2];
    prefactor[{{2,0,0}}] = -pow(center[2], 3);
    prefactor[{{1,0,0}}] =  2. * center[0] * pow(center[2], 3);
    prefactor[{{0,0,0}}] = -center[0] * center[0] * pow(center[2], 3);
    prefactor *= 32./(3.*sqrt(5.)) * pow(ec, 3.25) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::hxyyyy) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{1,4,0}}] =  1.;
    prefactor[{{0,4,0}}] = -center[0];
    prefactor[{{1,3,0}}] = -4. * center[1];
    prefactor[{{0,3,0}}] =  4. * center[0] * center[1];
    prefactor[{{1,2,0}}] =  6. * center[1] * center[1];
    prefactor[{{0,2,0}}] = -6. * center[0] * center[1] * center[1];
    prefactor[{{1,1,0}}] = -4. * pow(center[1], 3);
    prefactor[{{0,1,0}}] =  4. * center[0] * pow(center[1], 3);
    prefactor[{{1,0,0}}] =  pow(center[1], 4);
    prefactor[{{0,0,0}}] = -center[0] * pow(center[1], 4);
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::hxyyyz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{1,3,1}}] =  1.;
    prefactor[{{0,3,1}}] = -center[0];
    prefactor[{{1,2,1}}] = -3. * center[1];
    prefactor[{{0,2,1}}] =  3. * center[0] * center[1];
    prefactor[{{1,1,1}}] =  3. * center[1] * center[1];
    prefactor[{{0,1,1}}] = -3. * center[0] * center[1] * center[1];
    prefactor[{{1,0,1}}] = -pow(center[1], 3);
    prefactor[{{0,0,1}}] =  center[0] * pow(center[1], 3);
    prefactor[{{1,3,0}}] = -center[2];
    prefactor[{{0,3,0}}] =  center[0] * center[2];
    prefactor[{{1,2,0}}] =  3. * center[1] * center[2];
    prefactor[{{0,2,0}}] = -3. * center[0] * center[1] * center[2];
    prefactor[{{1,1,0}}] = -3. * center[1] * center[1] * center[2];
    prefactor[{{0,1,0}}] =  3. * center[0] * center[1] * center[1] * center[2];
    prefactor[{{1,0,0}}] =  pow(center[1], 3) * center[2];
    prefactor[{{0,0,0}}] = -center[0] * pow(center[1], 3) * center[2];
    prefactor *= 32./sqrt(15.) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxyyzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{1,2,2}}] =  1.;
    prefactor[{{0,2,2}}] = -center[0];
    prefactor[{{1,1,2}}] = -2. * center[1];
    prefactor[{{0,1,2}}] =  2. * center[0] * center[1];
    prefactor[{{1,0,2}}] =  center[1] * center[1];
    prefactor[{{0,0,2}}] = -center[0] * pow(center[1], 2);
    prefactor[{{1,2,1}}] = -2. * center[2];
    prefactor[{{0,2,1}}] =  2. * center[0] * center[2];
    prefactor[{{1,1,1}}] =  4. * center[1] * center[2];
    prefactor[{{0,1,1}}] = -4. * center[0] * center[1] * center[2];
    prefactor[{{1,0,1}}] = -2. * center[1] * center[1] * center[2];
    prefactor[{{0,0,1}}] =  2. * center[0] * center[1] * center[1] * center[2];
    prefactor[{{1,2,0}}] =  center[2] * center[2];
    prefactor[{{0,2,0}}] = -center[0] * center[2] * center[2];
    prefactor[{{1,1,0}}] = -2. * center[1] * center[2] * center[2];
    prefactor[{{0,1,0}}] =  2. * center[0] * center[1] * center[2] * center[2];
    prefactor[{{1,0,0}}] =  center[1] * center[1] * center[2] * center[2];
    prefactor[{{0,0,0}}] = -center[0] * center[1] * center[1] * center[2] * center[2];
    prefactor *= 32./3. * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxyzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{1,1,3}}] =  1.;
    prefactor[{{0,1,3}}] = -center[0];
    prefactor[{{1,0,3}}] = -center[1];
    prefactor[{{0,0,3}}] =  center[0] * center[1];
    prefactor[{{1,1,2}}] = -3. * center[2];
    prefactor[{{0,1,2}}] =  3. * center[0] * center[2];
    prefactor[{{1,0,2}}] =  3. * center[1] * center[2];
    prefactor[{{0,0,2}}] = -3. * center[0] * center[1] * center[2];
    prefactor[{{1,1,1}}] =  3. * center[2] * center[2];
    prefactor[{{0,1,1}}] = -3. * center[0] * center[2] * center[2];
    prefactor[{{1,0,1}}] = -3. * center[1] * center[2] * center[2];
    prefactor[{{0,0,1}}] =  3. * center[0] * center[1] * center[2] * center[2];
    prefactor[{{1,1,0}}] = -pow(center[2], 3);
    prefactor[{{0,1,0}}] =  center[0] * pow(center[2], 3);
    prefactor[{{1,0,0}}] =  center[1] * pow(center[2], 3);
    prefactor[{{0,0,0}}] = -center[0] * center[1] * pow(center[2], 3);
    prefactor *= 32./sqrt(15.) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }
  else if(type == GaussianType::hxzzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{1,0,4}}] =  1.;
    prefactor[{{0,0,4}}] = -center[0];
    prefactor[{{1,0,3}}] = -4. * center[2];
    prefactor[{{0,0,3}}] =  4. * center[0] * center[2];
    prefactor[{{1,0,2}}] =  6. * center[2] * center[2];
    prefactor[{{0,0,2}}] = -6. * center[0] * center[2] * center[2];
    prefactor[{{1,0,1}}] = -4. * pow(center[2], 3);
    prefactor[{{0,0,1}}] =  4. * center[0] * pow(center[2], 3);
    prefactor[{{1,0,0}}] =  pow(center[2], 4);
    prefactor[{{0,0,0}}] = -center[0] * pow(center[2], 4);
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75);  
  }
  else if(type == GaussianType::hyyyyy) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,5,0}}] =  1.;
    prefactor[{{0,4,0}}] = -5. * center[1];
    prefactor[{{0,3,0}}] =  10. * center[1] * center[1];
    prefactor[{{0,2,0}}] = -10. * pow(center[1], 3);
    prefactor[{{0,1,0}}] =  5. * pow(center[1], 4); 
    prefactor[{{0,0,0}}] = -pow(center[1], 5);
    prefactor *= 32./(3.*sqrt(105.)) * pow(ec, 3.25) * pow(2./pi, 0.75);  
  }
  else if(type == GaussianType::hyyyyz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,4,1}}] =  1.;
    prefactor[{{0,3,1}}] = -4. * center[1];
    prefactor[{{0,2,1}}] =  6. * center[1] * center[1];
    prefactor[{{0,1,1}}] = -4. * pow(center[1], 3);
    prefactor[{{0,0,1}}] =  pow(center[1], 4);
    prefactor[{{0,4,0}}] = -center[2];
    prefactor[{{0,3,0}}] =  4. * center[1] * center[2];
    prefactor[{{0,2,0}}] = -6. * center[1] * center[1] * center[2];
    prefactor[{{0,1,0}}] =  4. * pow(center[1], 3) * center[2];
    prefactor[{{0,0,0}}] = -pow(center[1], 4) * center[2];
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75);     
  }
  else if(type == GaussianType::hyyyzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,3,2}}] =  1.;
    prefactor[{{0,2,2}}] = -3. * center[1];
    prefactor[{{0,1,2}}] =  3. * center[1] * center[1];
    prefactor[{{0,0,2}}] = -pow(center[1], 3);
    prefactor[{{0,3,1}}] = -2. * center[2];
    prefactor[{{0,2,1}}] =  6. * center[1] * center[2];
    prefactor[{{0,1,1}}] = -6. * center[1] * center[1] * center[2];
    prefactor[{{0,0,1}}] =  2. * pow(center[1], 3) * center[2];
    prefactor[{{0,3,0}}] =  center[2] * center[2];
    prefactor[{{0,2,0}}] = -3. * center[1] * center[2] * center[2];
    prefactor[{{0,1,0}}] =  3. * center[1] * center[1] * center[2] * center[2];
    prefactor[{{0,0,0}}] = -pow(center[1], 3) * center[2] * center[2];
    prefactor *= 32./(3.*sqrt(5.)) * pow(ec, 3.25) * pow(2./pi, 0.75); 
  }
  else if(type == GaussianType::hyyzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,2,3}}] =  1.;
    prefactor[{{0,1,3}}] = -2. * center[1];
    prefactor[{{0,0,3}}] =  center[1] * center[1];
    prefactor[{{0,2,2}}] = -3. * center[2];
    prefactor[{{0,1,2}}] =  6. * center[1] * center[2];
    prefactor[{{0,0,2}}] = -3. * center[1] * center[1] * center[2];
    prefactor[{{0,2,1}}] =  3. * center[2] * center[2];
    prefactor[{{0,1,1}}] = -6. * center[1] * pow(center[2], 2);
    prefactor[{{0,0,1}}] =  3. * center[1] * center[1] * center[2] * center[2];
    prefactor[{{0,2,0}}] = -pow(center[2], 3);
    prefactor[{{0,1,0}}] =  2. * center[1] * pow(center[2], 3);
    prefactor[{{0,0,0}}] = -center[1] * center[1] * pow(center[2], 3);
    prefactor *= 32./(3.*sqrt(5.)) * pow(ec, 3.25) * pow(2./pi, 0.75);  
  }
  else if(type == GaussianType::hyzzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,1,4}}] =  1.;
    prefactor[{{0,0,4}}] = -center[1];
    prefactor[{{0,1,3}}] = -4. * center[2];
    prefactor[{{0,0,3}}] =  4. * center[1] * center[2];
    prefactor[{{0,1,2}}] =  6. * center[2] * center[2];
    prefactor[{{0,0,2}}] = -6. * center[1] * center[2] * center[2];
    prefactor[{{0,1,1}}] = -4. * pow(center[2], 3);
    prefactor[{{0,0,1}}] =  4. * center[1] * pow(center[2], 3);
    prefactor[{{0,1,0}}] =  pow(center[2], 4);
    prefactor[{{0,0,0}}] = -center[1] * pow(center[2], 4);
    prefactor *= 32./sqrt(105.) * pow(ec, 3.25) * pow(2./pi, 0.75);   
  }
  else if(type == GaussianType::hzzzzz) {
    prefactor = PolynomialCoeffs(5);
    prefactor[{{0,0,5}}] =  1.;
    prefactor[{{0,0,4}}] = -5. * center[2];
    prefactor[{{0,0,3}}] =  10. * center[2] * center[2];
    prefactor[{{0,0,2}}] = -10. * pow(center[2], 3);
    prefactor[{{0,0,1}}] =  5. * pow(center[2], 4); 
    prefactor[{{0,0,0}}] = -pow(center[2], 5);
    prefactor *= 32./(3.*sqrt(105.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
  }

  // h sphericals
  // (l,m) = (5,-5) 
  else if(type == GaussianType::hm5) {
    // (5x^4-10x^2y^2+y^4), multiply by y, normalize
    PolynomialCoeffs y(1), x4(4), xx(2), yy(2),  y4(4);
    y[{{0,1,0}}]  =  1.;
    y[{{0,0,0}}]  = -center[1];
    x4[{{4,0,0}}]  =  1.;
    x4[{{3,0,0}}]  = -4. * center[0];
    x4[{{2,0,0}}]  =  6. * center[0] * center[0];
    x4[{{1,0,0}}]  = -4. * pow(center[0], 3);
    x4[{{0,0,0}}]  =  pow(center[0], 4);
    y4[{{0,4,0}}]  =  1.;
    y4[{{0,3,0}}]  = -4. * center[1];
    y4[{{0,2,0}}]  =  6. * center[1] * center[1];
    y4[{{0,1,0}}]  = -4. * pow(center[1], 3);
    y4[{{0,0,0}}]  =  pow(center[1], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.; 
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    double norm =  4./sqrt(15.) * pow(2., 0.25) * pow(ec, 3.25) * pow(1./pi, 0.75);
    prefactor = norm * y * (5.*x4 + (-10.*xx*yy) + y4);
    }
  // (l,m) = (5,-4)
  else if(type == GaussianType::hm4) {
    // get (x^2-y^2), multiply by xyz, normalize
    PolynomialCoeffs xx(2), yy(2), x(1), y(1), z(1);
    x[{{1,0,0}}] =  1.;
    x[{{0,0,0}}] = -center[0];
    y[{{0,1,0}}] =  1.;
    y[{{0,0,0}}] = -center[1];
    z[{{0,0,1}}] =  1.;
    z[{{0,0,0}}] = -center[2];
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    double norm = 16./sqrt(3.) * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * x * y * z * (xx + (-1.*yy));
  }
  // (l,m) = (5,-3)
  else if(type == GaussianType::hm3) {
    // get y, (3x^2-y^2) and (8z^2-x^2-y^2), multiply together, normalize
    PolynomialCoeffs xx(2), yy(2), zz(2), y(1);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    y[{{0,1,0}}] =  1.;
    y[{{0,0,0}}] = -center[1];
    double norm = 4./(3.*sqrt(3.)) * pow(2., 0.25) * pow(ec, 3.25) * pow(1./pi, 0.75);
    prefactor = norm * y * (3.*xx + (-1.*yy)) * (8.*zz + (-1.*xx) + (-1.*yy));
  }
  // (l,m) = (5,-2)
  else if(type == GaussianType::hm2) {
    // get (2z^2-x^2-y^2), multiply by xyz, normalize
    PolynomialCoeffs xx(2), yy(2), zz(2), x(1), y(1), z(1);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    x[{{1,0,0}}] =  1.;
    x[{{0,0,0}}] = -center[0];
    y[{{0,1,0}}] =  1.; 
    y[{{0,0,0}}] = -center[1];
    z[{{0,0,1}}] =  1.;
    z[{{0,0,0}}] = -center[2];
    double norm = 16./3. * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * x * y * z * (2.*zz + (-1.*xx) + (-1.*yy)); 
  }
  // (l,m) = (5,-1)
  else if(type == GaussianType::hm1) {
    // get x^4+2x^2y^2-12x^2z^2+y^4-12y^2z^2+8z^4
    // multiply by y, normalize
    PolynomialCoeffs x4(4), xx(2), yy(2), y4(4), zz(2), z4(4), y(1);
    x4[{{4,0,0}}] =  1.;
    x4[{{3,0,0}}] = -4. * center[0];
    x4[{{2,0,0}}] =  6. * center[0] * center[0];
    x4[{{1,0,0}}] = -4. * pow(center[0], 3);
    x4[{{0,0,0}}] =  1. * pow(center[0], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    y4[{{0,4,0}}] =  1.;
    y4[{{0,3,0}}] = -4. * center[1];
    y4[{{0,2,0}}] =  6. * center[1] * center[1];
    y4[{{0,1,0}}] = -4. * pow(center[1], 3);
    y4[{{0,0,0}}] =  1. * pow(center[1], 4);
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    z4[{{0,0,4}}] =  1.;
    z4[{{0,0,3}}] = -4. * center[2];
    z4[{{0,0,2}}] =  6. * center[2] * center[2];
    z4[{{0,0,1}}] = -4. * pow(center[2], 3);
    z4[{{0,0,0}}] =  1. * pow(center[2], 4);
    y[{{0,1,0}}] =  1.;
    y[{{0,0,0}}] = -center[1];
    double norm = 4./(3.*sqrt(7.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * y * (x4 + 2.*xx*yy + (-12.*xx*zz) + y4 + (-12.*yy*zz) + 8.*z4);
  }
  // (l,m) = (5,0)
  else if(type == GaussianType::hzero) {
    // get 15x^4+30x^2y^2-40x^2z^2+15y^4-40y^2z^2+8z^4
    // multiply by z, normalize
    PolynomialCoeffs x4(4), xx(2), zz(2), y4(4), yy(2), z4(4), z(1);
    x4[{{4,0,0}}] =  1.;
    x4[{{3,0,0}}] = -4. * center[0];
    x4[{{2,0,0}}] =  6. * center[0] * center[0];
    x4[{{1,0,0}}] = -4. * pow(center[0], 3);
    x4[{{0,0,0}}] =  1. * pow(center[0], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    y4[{{0,4,0}}] =  1.;
    y4[{{0,3,0}}] = -4. * center[1];
    y4[{{0,2,0}}] =  6. * center[1] * center[1];
    y4[{{0,1,0}}] = -4. * pow(center[1], 3);
    y4[{{0,0,0}}] =  1. * pow(center[1], 4);
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  pow(center[1], 2);
    z4[{{0,0,4}}] =  1.;
    z4[{{0,0,3}}] = -4. * center[2];
    z4[{{0,0,2}}] =  6. * center[2] * center[2];
    z4[{{0,0,1}}] = -4. * pow(center[2], 3);
    z4[{{0,0,0}}] =  1. * pow(center[2], 4);
    z[{{0,0,1}}] =  1.;
    z[{{0,0,0}}] = -center[2];
    double norm = 4./(3.*sqrt(105.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * z * (15.*x4 + 30.*xx*yy + (-40.*xx*zz) + 15.*y4 + (-40.*yy*zz) + 8.*z4); 
  }
  // (l,m) = (5,1)
  else if(type == GaussianType::hp1) {
    // get x^4+2x^2y^2-12x^2z^2+y^4-12y^2z^2+8z^4
    // multiply by x, normalize
    PolynomialCoeffs x4(4), xx(2), zz(2), y4(4), yy(2), z4(4), x(1);
    x4[{{4,0,0}}] =  1.;
    x4[{{3,0,0}}] = -4. * center[0];
    x4[{{2,0,0}}] =  6. * center[0] * center[0];
    x4[{{1,0,0}}] = -4. * pow(center[0], 3);
    x4[{{0,0,0}}] =  1. * pow(center[0], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    y4[{{0,4,0}}] =  1.;
    y4[{{0,3,0}}] = -4. * center[1];
    y4[{{0,2,0}}] =  6. * center[1] * center[1];
    y4[{{0,1,0}}] = -4. * pow(center[1], 3);
    y4[{{0,0,0}}] =  1. * pow(center[1], 4);
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    z4[{{0,0,4}}] =  1.;
    z4[{{0,0,3}}] = -4. * center[2];
    z4[{{0,0,2}}] =  6. * center[2] * center[2];
    z4[{{0,0,1}}] = -4. * pow(center[2], 3);
    z4[{{0,0,0}}] =  1. * pow(center[2], 4);
    x[{{1,0,0}}] =  1.;
    x[{{0,0,0}}] = -center[0];
    double norm = 4./(3.*sqrt(7.)) * pow(ec, 3.25) * pow(2./pi, 0.75);
    // Testing shows nwchem uses the negative of this
    prefactor = (-1.) * norm * x * (x4 + 2.*xx*yy + (-12.*xx*zz) + y4 + (-12.*yy*zz) + 8.*z4);
  }
  // (l,m) = (5,2)
  else if(type == GaussianType::hp2) {
    // get z * (x^2-y^2) * (2z^2 - x^2 - y^2), and normalize
    PolynomialCoeffs xx(2), yy(2), zz(2), z(1);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    z[{{0,0,1}}] =  1.;
    z[{{0,0,0}}] = -center[2];
    double norm = 8./3. * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * z * (xx + (-1.*yy)) * (2.*zz + (-1.*xx) + (-1.*yy)); 
  }
  // (l,m) = (5,3)
  else if(type == GaussianType::hp3) {
    // get x, (x^2-3y^2) and (8z^2-x^2-y^2), 
    // multiply and normalize
    PolynomialCoeffs xx(2), yy(2), zz(2), x(1);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.;
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    zz[{{0,0,2}}] =  1.;
    zz[{{0,0,1}}] = -2. * center[2];
    zz[{{0,0,0}}] =  center[2] * center[2];
    x[{{1,0,0}}] =  1.;
    x[{{0,0,0}}] = -center[0];
    double norm = 4./(3.*sqrt(3.)) * pow(2., 0.25) * pow(ec, 3.25) * pow(1./pi, 0.75);
    // Testing shows nwchem uses the negative of this
    prefactor = (-1.) * norm * x * (xx + (-3.*yy)) * (8.*zz + (-1.*xx) + (-1.*yy)); 
  } 
  // (l,m) = (5,4)
  else if(type == GaussianType::hp4) {
    // setup x^4, x^2y^2, y^4, get (x^4-6x^2y^2+y^4), normalize, multiply by z
    PolynomialCoeffs z(1), x4(4), xx(2), yy(2), y4(4);
    z[{{0,0,1}}]  =  1.;
    z[{{0,0,0}}]  = -center[2];
    x4[{{4,0,0}}]  =  1.;
    x4[{{3,0,0}}]  = -4. * center[0];
    x4[{{2,0,0}}]  =  6. * center[0] * center[0];
    x4[{{1,0,0}}]  = -4. * pow(center[0], 3);
    x4[{{0,0,0}}]  =  pow(center[0], 4);
    y4[{{0,4,0}}]  =  1.;
    y4[{{0,3,0}}]  = -4. * center[1];
    y4[{{0,2,0}}]  =  6. * center[1] * center[1];
    y4[{{0,1,0}}]  = -4. * pow(center[1], 3);
    y4[{{0,0,0}}]  =  pow(center[1], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.; 
    yy[{{0,1,0}}] = -2. *center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    double norm = 4./sqrt(3.) * pow(ec, 3.25) * pow(2./pi, 0.75);
    prefactor = norm * z * (x4 + (-6.*xx*yy) + y4); 
  }
  // (l,m) = (5,5)
  else if(type == GaussianType::hp5) {
    // get x, (x^4-10x^2y^2+5y^4), multiply, normalize
    PolynomialCoeffs x(1), x4(4), xx(2), yy(2),  y4(4);
    x[{{1,0,0}}]  =  1.;
    x[{{0,0,0}}]  = -center[0];
    x4[{{4,0,0}}]  =  1.;
    x4[{{3,0,0}}]  = -4. * center[0];
    x4[{{2,0,0}}]  =  6. * center[0] * center[0];
    x4[{{1,0,0}}]  = -4. * pow(center[0], 3);
    x4[{{0,0,0}}]  =  pow(center[0], 4);
    y4[{{0,4,0}}]  =  1.;
    y4[{{0,3,0}}]  = -4. * center[1];
    y4[{{0,2,0}}]  =  6. * center[1] * center[1];
    y4[{{0,1,0}}]  = -4. * pow(center[1], 3);
    y4[{{0,0,0}}]  =  pow(center[1], 4);
    xx[{{2,0,0}}] =  1.;
    xx[{{1,0,0}}] = -2. * center[0];
    xx[{{0,0,0}}] =  center[0] * center[0];
    yy[{{0,2,0}}] =  1.; 
    yy[{{0,1,0}}] = -2. * center[1];
    yy[{{0,0,0}}] =  center[1] * center[1];
    double norm = 4./sqrt(15.) * pow(2., 0.25) * pow(ec, 3.25) * pow(1./pi, 0.75);
    // Testing shows nwchem uses the negative of this
    prefactor = (-1.) * norm * x * (x4 + (-10.*xx*yy) + 5.*y4);
  }
}

PrimitiveGaussian::PrimitiveGaussian(const PolynomialCoeffs &exppoly_,
    const PolynomialCoeffs &prefactor_)
    : prefactor(prefactor_), exppoly(exppoly_) {

  const unsigned degree = exppoly_.get_degree();
  if(degree > 2)
    throw std::invalid_argument("Exponential in a Gaussian must have degree at most 2.");
  else if(degree < 2) {
    // it will be assumed in other functions that exppoly's degree is 2
    exppoly = PolynomialCoeffs(2);
    exppoly[{{0,0,0}}] = exppoly_[{{0,0,0}}];
    if(degree == 1) {
      exppoly[{{1,0,0}}] = exppoly_[{{1,0,0}}];
      exppoly[{{0,1,0}}] = exppoly_[{{0,1,0}}];
      exppoly[{{0,0,1}}] = exppoly_[{{0,0,1}}];
    }
  }
}

double PrimitiveGaussian::operator() (const std::array<double, 3> &x) const {
  return exp(exppoly(x)) * prefactor(x);
}

PrimitiveGaussian PrimitiveGaussian::operator*(const PrimitiveGaussian &rhs) const {
  PrimitiveGaussian ret;
  ret.prefactor = prefactor * rhs.prefactor;
  ret.exppoly = exppoly + rhs.exppoly;

  return ret;
}

///////////////////////////////////////////////////////////////////////////
// Implementation of the GaussianFunction class.
///////////////////////////////////////////////////////////////////////////
GaussianFunction::GaussianFunction(const GaussianType type,
    const std::array<double, 3> &center, const std::vector<double> &expcoeff,
    const std::vector<double> &coeff) {

  if(expcoeff.size() != coeff.size())
    throw std::invalid_argument("A Gaussian function must have the same number of linear expansion and exponential coefficients.");

  // go through each of the primitives and construct the Gaussian function
  for(auto iterexpcoeff = expcoeff.cbegin(),
           itercoeff = coeff.cbegin(),
           end = expcoeff.cend();
      iterexpcoeff != end; 
      iterexpcoeff++, itercoeff++) {

    terms.emplace_front(*itercoeff, PrimitiveGaussian(type, center, *iterexpcoeff));
  }

  // prune any small terms in the expansion
  removeSmallTerms();

  // store the center
  this->center = center;

  // store the exponents
  this->expcoeff = expcoeff;
}

double GaussianFunction::operator() (const std::array<double, 3> &x) const {
  double ret = 0.;
  for(const TermT &term : terms)
    ret += term.first * term.second(x);
  return ret;
}

GaussianFunction GaussianFunction::operator- () const {
  GaussianFunction ret(*this);

  for(TermT &term : ret.terms)
    term.first *= -1.;

  return ret;
}

GaussianFunction GaussianFunction::operator+(const GaussianFunction &rhs) const {
  GaussianFunction ret;
  ret.terms = terms; // deep copy
  for(const TermT &term : rhs.terms)
    ret.terms.push_front(term);
  return ret;
}

GaussianFunction &GaussianFunction::operator+=(const GaussianFunction &rhs) {
  if(this == &rhs)
    for(TermT &term : terms)
      term.first *= 2.;
  else
    for(const TermT &term : rhs.terms)
      terms.push_front(term);
  return *this;
}

GaussianFunction &GaussianFunction::operator*=(const double rhs) {
  for(TermT &term : terms)
    term.first *= rhs;
  removeSmallTerms();
  return *this;
}

GaussianFunction GaussianFunction::operator*(const GaussianFunction &rhs) const {
  GaussianFunction ret;
  for(const TermT &term : terms)
    for(const TermT &rterm : rhs.terms) {
      ret.terms.emplace_front(
        term.first * rterm.first,
        term.second * rterm.second
      );
    }
  ret.removeSmallTerms();
  return ret;
}

} // namespace slymer
