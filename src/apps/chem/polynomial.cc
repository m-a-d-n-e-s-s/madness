/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file Basis/polynomial.cc
 * \brief Implementation of polynomial API and routines.
 */

#include <cmath>
#include <vector>
#include "polynomial.h"

namespace slymer {

double PolynomialCoeffs::operator() (const std::array<double, 3> &pt) const {
  double ret = 0.;

  for(unsigned j = 0; j <= degree; ++j)
    for(unsigned k = 0; k <= degree - j; ++k)
      for(unsigned l = 0; l <= degree - j - k; ++l)
      ret += coeffs[polyIndex({{j,k,l}})]
        * ::pow(pt[0], j) * ::pow(pt[1], k) * ::pow(pt[2], l);

  return ret;
}

PolynomialCoeffs PolynomialCoeffs::operator+ (const PolynomialCoeffs &rhs) const
{
  if(degree >= rhs.degree) {
    PolynomialCoeffs ret(degree);

    unsigned size = rhs.coeffs.size();
    for(unsigned j = 0; j < size; ++j)
      ret.coeffs[j] = coeffs[j] + rhs.coeffs[j];
    size = coeffs.size();
    for(unsigned j = rhs.coeffs.size(); j < size; ++j)
      ret.coeffs[j] = coeffs[j];

    return ret;
  }
  else {
    PolynomialCoeffs ret(rhs.degree);

    unsigned size = coeffs.size();
    for(unsigned j = 0; j < size; ++j)
      ret.coeffs[j] = coeffs[j] + rhs.coeffs[j];
    size = rhs.coeffs.size();
    for(unsigned j = coeffs.size(); j < size; ++j)
      ret.coeffs[j] = rhs.coeffs[j];

    return ret;
  }
}

PolynomialCoeffs &PolynomialCoeffs::operator+= (const PolynomialCoeffs &rhs) {
  if(degree >= rhs.degree) {
    const unsigned size = rhs.coeffs.size();
    for(unsigned j = 0; j < size; ++j)
      coeffs[j] += rhs.coeffs[j];
  }
  else {
    PolynomialCoeffs newpoly(rhs.degree);

    unsigned size = coeffs.size();
    for(unsigned j = 0; j < size; ++j)
      newpoly.coeffs[j] = coeffs[j] + rhs.coeffs[j];
    size = rhs.coeffs.size();
    for(unsigned j = coeffs.size(); j < size; ++j)
      newpoly.coeffs[j] = rhs.coeffs[j];

    std::swap(coeffs, newpoly.coeffs);
    degree = rhs.degree;
  }

  return *this;
}

PolynomialCoeffs PolynomialCoeffs::operator* (const double c) const {
  PolynomialCoeffs ret(*this);

  for(double &coeff : ret.coeffs)
    coeff *= c;

  return ret;
}

PolynomialCoeffs &PolynomialCoeffs::operator*= (const double c) {
  for(double &coeff : coeffs)
    coeff *= c;

  return *this;
}

PolynomialCoeffs PolynomialCoeffs::operator* (const PolynomialCoeffs &rhs) const
{
  // make an output vector with the appropriate size for a multinomial of degree j
  const unsigned deg = degree + rhs.degree;
  PolynomialCoeffs ret(deg);

  // go through all combinations for the product
  for(unsigned x1 = 0; x1 <= degree; ++x1)
    for(unsigned y1 = 0; y1 <= degree - x1; ++y1)
      for(unsigned z1 = 0; z1 <= degree - x1 - y1; ++z1) {

        for(unsigned x2 = 0; x2 <= rhs.degree; ++x2)
          for(unsigned y2 = 0; y2 <= rhs.degree - x2; ++y2)
            for(unsigned z2 = 0; z2 <= rhs.degree - x2 - y2; ++z2) {
              ret[{{x1+x2, y1+y2, z1+z2}}] +=
                operator[]({{x1, y1, z1}}) * rhs[{{x2, y2, z2}}];
            }
      }

  return ret;
}

PolynomialCoeffs PolynomialCoeffs::pow(const unsigned j) const {
  // this function is only intended to work with polynomials of degree 1
  if(degree != 1)
    throw std::runtime_error("PolynomialCoeffs::pow can only be applied to polynomials with degree 1.");

  // array with factorials
  const std::array<unsigned, 11> fact{{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800}};

  if(j >= fact.size())
    throw std::invalid_argument("Degree is too big in PolynomialCoeffs::pow.");

  // setup a PolynomialCoeffs for the return values
  PolynomialCoeffs ret(j);

  // cycle through all exponent orders that add up to j
  for(unsigned ex1 = 0; ex1 <= j; ++ex1)
    for(unsigned ex2 = 0; ex2 <= j - ex1; ++ex2)
      for(unsigned ex3 = 0; ex3 <= j - ex1 - ex2; ++ex3) {
        unsigned ex4 = j - ex1 - ex2 - ex3; // the ex coefficients add up to j

        // calculate the coefficient on x^ex1 y^ex2 z^ex3, per the multinomial theorem
        ret[{{ex1, ex2, ex3}}] =
          (fact[j] / fact[ex1] / fact[ex2] / fact[ex3] / fact[ex4])
          * ::pow(operator[]({{1,0,0}}), ex1) * ::pow(operator[]({{0,1,0}}), ex2)
          * ::pow(operator[]({{0,0,1}}), ex3) * ::pow(operator[]({{0,0,0}}), ex4);
      }

  return ret;
}

} // namespace slymer
