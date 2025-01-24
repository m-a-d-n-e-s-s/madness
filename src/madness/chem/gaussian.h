/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file Basis/gaussian.h
 * \brief Gaussian basis function API and routines.
 *
 * Handles Gaussian basis functions as special cases (class inheritance). Both
 * primitive and contracted Gaussians from electronic structure codes are
 * implemented, but not using those labels.
 */

#ifndef __Basis_gaussian_h__
#define __Basis_gaussian_h__

#include <forward_list>
#include <utility>
#include<madness/chem/basis.h>
#include<madness/chem/polynomial.h>
#include<madness/chem/ESInterface.h>
#include <madness/mra/mra.h>
#include <madness/constants.h>

namespace slymer {

/// Implemented types of Gaussian orbitals (from quantum chemistry codes).
enum class GaussianType {
  // both Cartesian and spherical
  s, px, py, pz,

  // Cartesian (and possibly spherical)
  dxx, dxy, dxz, dyy, dyz, dzz,
  fxxx, fxxy, fxxz, fxyy, fxyz, fxzz, fyyy, fyyz, fyzz, fzzz,

  // spherical only ('m' denotes "minus")
  dzzmrr, dxxmyy,
  fxyymxxx, fxxzmyyz, fxzzmrrx, fzzzmrrz, fyzzmrry, fxxymyyy,

  // g
  // cartesian
  gxxxx, gxxxy, gxxxz, gxxyy, gxxyz, gxxzz, gxyyy, gxyyz, gxyzz,
  gxzzz, gyyyy, gyyyz, gyyzz, gyzzz, gzzzz,
  // spherical
  // ('d' denotes "dot" for multiplication)
  // (numbers in names indicate powers of preceding variable)
  gxydx2my2, gyzdx2my2, gxydz2mr2, gyzdz2mr2, gzero,
  gxzdz2mr2, gx2my2dz2mr2, gxzdx2my2, gx4mx2y2py4, 

  // h
  // cartesian
  hxxxxx, hxxxxy, hxxxxz, hxxxyy, hxxxyz, hxxxzz, hxxyyy, hxxyyz,
  hxxyzz, hxxzzz, hxyyyy, hxyyyz, hxyyzz, hxyzzz, hxzzzz, hyyyyy,
  hyyyyz, hyyyzz, hyyzzz, hyzzzz, hzzzzz,
  // spherical
  // at this point, just giving names corresponding to value
  // quantum number m
  hm5, hm4, hm3, hm2, hm1, hzero, hp1, hp2, hp3, hp4, hp5

};

// forward declaration
class GaussianFunction;

/**
 * \brief A primitive Gaussian function.
 *
 * This basis function has the form \f[
 * f(x,y,z) = p(x,y,z) \exp[(x-x_0)^T \Sigma (x-x_0)],
 * \f]
 * where \f$p\f$ is a polynomial and \f$x_0\f$ is the ``center'' of the
 * Gaussian and \f$\Sigma\f$ is the variance/covariance matrix.
 *
 * This class is designed to work with Gaussian basis functions in quantum
 * chemistry codes, but is made to be a bit more general for other
 * applications.
 *
 * The wikipedia page on multivariable normal distributions may be helpful
 * (https://en.wikipedia.org/wiki/Multivariate_normal_distribution).
 */
class PrimitiveGaussian {

protected:
  /// Polynomial prefactor (\f$p\f$ in the class description).
  PolynomialCoeffs prefactor;

  /// Coefficients in the exponent's quadratic function.
  PolynomialCoeffs exppoly;

  /**
   * \brief Return type for the function's canonical form.
   *
   * The first element is the polynomial prefactor, converted to use the principal coordinates.
   * The second element is the array of decay values (eigenvalues of the covariance matrix).
   * The third element is the \f$\nu\f$ vector (essentially the center of the Gaussian).
   * The fourth element are the canonical axes in the original coordinates (for debugging purposes).
   */
  using CF = std::tuple<PolynomialCoeffs, std::array<double, 3>, std::array<double, 3>, std::array<std::array<double, 3>, 3>>;

public:
  /// Default constructor: Make the function 1 (\f$p=1\f$, \f$f(\vec{x})=0\f$).
  PrimitiveGaussian()
    : prefactor(0), exppoly(2)
  {
    prefactor[{{0,0,0}}] = 1.;
  }

  /**
   * \brief Create a primitive Gaussian function centered on the specific point,
   *    with the specified decay constant and type.
   *
   * \throw std::invalid_argument if the decay constant is non-positive.
   *
   * \param[in] type The orbital type of the Gaussian.
   * \param[in] center The center of the Gaussian.
   * \param[in] ec The decay constant in the exponential.
   */
  PrimitiveGaussian(const GaussianType &type,
      const std::array<double, 3> &center, const double ec);

  /**
   * \brief Create a primitive Gaussian function with the specified exponential
   *    polynomial (quadratic) and polynomial prefactor.
   *
   * \throw std::invalid_argument if the exponential polynomial is not
   *    quadratic (or constant or linear as subsets).
   *
   * \param[in] exppoly_ The exponential polynomial.
   * \param[in] prefactor_ The degree of the prefactor.
   */
  PrimitiveGaussian(const PolynomialCoeffs &exppoly_,
      const PolynomialCoeffs &prefactor_);

  /**
   * \brief Evaluate the Gaussian primitive at the specified point.
   *
   * \param[in] x The point.
   * \return The Gaussian primitive evaluated at the point x.
   */
  double operator() (const std::array<double, 3> &x) const;

  /**
   * \brief Multiply two primitive Gaussians. The result is another.
   *
   * \param[in] rhs The right-hand object.
   * \result The product.
   */
  PrimitiveGaussian operator*(const PrimitiveGaussian &rhs) const;
};


/**
 * \brief A Gaussian basis function used by chemistry electronic structure
 *   codes.
 *
 * This basis function has the form \f[
 * \sum_j c_j x^{l_j} y^{m_j} z^{n_j} \exp[f_j(\vec{x})],
 * \f]
 * where \f$f_j(\vec{x})\f$ is a quadratic function. Note that this functional
 * form is a bit more general than that used by standard electronic structure
 * codes; we do not require the prinicpal axes of the Gaussian function to be
 * coincident with the x-, y-, and z-axes. This class can also hold/handle any
 * linear combination of primitive Gaussian functions (the exponents and
 * centers do not have to be the same).
 */
class GaussianFunction : public BasisFunction {

protected:
  /// Helper class for storing terms in the sum of Gaussian primitives.
  using TermT = std::pair<double, PrimitiveGaussian>;

  /// The list of terms in the sum of Gaussian primitives.
  std::forward_list<TermT> terms;

  /**
   * \brief Removes terms with very small linear expansion coefficients.
   *
   * \param[in] eps The threshold for "small". Defaults to 1.e-6.
   */
  void removeSmallTerms(const double eps = 1.e-6) {
    terms.remove_if(
      [eps](const TermT &t) {
        return std::abs(t.first) < eps;
      }
    );
  }

public:

  /// The center of the gaussian
  std::array<double, 3> center;

  /// The exponent coefficients
  std::vector<double> expcoeff;

  /// Default constructor: Make the function 0.
  GaussianFunction() : terms(0), center({{0,0,0}}), expcoeff(std::vector<double>{0.0}) {}

  /**
   * \brief Create a Gaussian orbital (uncontracted or contracted) from
   *    quantum chemistry codes.
   *
   * Both Cartesian and spherical orbitals are supported. The expcoeff and
   * coeff parameters must be arrays. If the orbital is uncontracted, the
   * arrays will be length 1.
   *
   * All uncontracted primitive Gaussian orbitals used underneath the
   * GaussianFunction are normalized.
   *
   * \throw invalid_argument If a non-positive exponential coefficient is
   *    passed or if the lengths of expcoeff and coeff are not equal.
   *
   * \param[in] type The orbital type (from enum OrbitalType).
   * \param[in] center The center of the Gaussian orbital.
   * \param[in] expcoeff Array of coefficients in the exponents of the
   *    primitive Gaussians.
   * \param[in] coeff Array of linear expansion coefficients of the primitive
   *    Gaussians.
   */
  GaussianFunction(const GaussianType type, const std::array<double, 3> &center,
    const std::vector<double> &expcoeff, const std::vector<double> &coeff);

  /**
   * \brief Create an uncontracted Gaussian orbital (shortcut to the more
   *    general constructor for Gaussian orbitals).
   *
   * \param[in] type The orbital type (from enum OrbitalType).
   * \param[in] center The center of the Gaussian orbital.
   * \param[in] expcoeff The coefficient in the exponent of the primitive
   *    Gaussian.
   */
  GaussianFunction(const GaussianType type, const std::array<double, 3> &center,
    const double expcoeff)
    : GaussianFunction(type, center, {expcoeff}, {1.}) {}

  /**
   * \brief Evaluate the Gaussian function at the specified point.
   *
   * \param[in] x The point.
   * \return The Gaussian function evaluated at the point x.
   */
  virtual double operator() (const std::array<double, 3> &x) const override;

  /**
   * \brief Evaluate the Gaussian function at the specified point.
   *
   * \param[in] x The point.
   * \return The Gaussian function evaluated at the point x.
   */
  double operator() (const madness::coord_3d& r) const {
     return operator()(std::array<double, 3>{{r[0], r[1], r[2]}});
  };

  /**
   * \brief Additive inverse of the GaussianFunction (deep copy).
   *
   * \return The negated function.
   */
  GaussianFunction operator- () const;

  /**
   * \brief Add two GaussianFunction objects, resulting in a new one.
   *
   * \param[in] rhs The right-hand GaussianFunction.
   * \return The sum of this GaussianFunction and rhs.
   */
  GaussianFunction operator+(const GaussianFunction &rhs) const;

  /**
   * \brief Add a GaussianFunction to this GaussianFunction.
   *
   * \param[in] rhs The right-hand GaussianFunction.
   * \return The sum of this GaussianFunction and rhs, stored in *this.
   */
  GaussianFunction &operator+=(const GaussianFunction &rhs);

  /**
   * \brief Subtract a GaussianFunction from this one, returning the difference.
   *
   * \param[in] rhs The GaussianFunction to be subtracted.
   * \return The difference of the two GaussianFunction objects.
   */
  GaussianFunction operator-(const GaussianFunction &rhs) const {
    return *this + (-rhs);
  }

  /**
   * \brief Subtract a GaussianFunction from this one.
   *
   * \param[in] rhs The GaussianFunction to be subtracted.
   * \return This GaussianFunction (the difference).
   */
  GaussianFunction &operator-=(const GaussianFunction &rhs) {
    return operator+=(-rhs);
  }

  /**
   * \brief Multiply the GaussianFunction by a scalar.
   *
   * \param[in] rhs The scalar multiplicative factor.
   * \return This GaussianFunction, which has been scaled.
   */
  GaussianFunction &operator*=(const double rhs);

  /**
   * \brief Multiply the GaussianFunction by another.
   *
   * \param[in] rhs The right-hand GaussianFunction object.
   * \return This GaussianFunction, which is now the product.
   */
  GaussianFunction &operator*=(const GaussianFunction &rhs) {
    GaussianFunction ret = (*this) * rhs;
    terms.swap(ret.terms);
    return *this;
  }

  /**
   * \brief Multiply the GaussianFunction by a scalar.
   *
   * \param[in] rhs The scalar multiplicative factor.
   * \return The scaled GaussianFunction.
   */
  GaussianFunction operator*(const double rhs) const {
    GaussianFunction ret(*this);
    ret *= rhs;
    return ret;
  }

  /**
   * \brief Multiply two GaussianFunctions.
   *
   * \param[in] rhs The right-hand GaussianFunction object.
   * \return The product.
   */
  GaussianFunction operator*(const GaussianFunction &rhs) const;

  /**
   * \brief Divide the GaussianFunction by a scalar.
   *
   * Does not check that the scalar is nonzero.
   *
   * \param[in] rhs The scalar factor.
   * \return This GaussianFunction, which has been scaled.
   */
  GaussianFunction &operator/=(const double rhs) {
    operator*=(1./rhs);
    return *this;
  }

  /**
   * \brief Divide the GaussianFunction by a scalar.
   *
   * \param[in] rhs The scalar factor.
   * \return The scaled GaussianFunction.
   */
  GaussianFunction operator/(const double rhs) const {
    return (*this) * (1./rhs);
  }
};


// helper functions

/**
 * \brief Multiply a GaussianFunction by a scalar.
 *
 * \param[in] lhs The scalar multiplicative factor.
 * \param[in] rhs The GaussianFunction.
 * \return The scaled GaussianFunction.
 */
inline GaussianFunction operator*(const double lhs, const GaussianFunction &rhs) {
  return rhs * lhs;
}

class Gaussian_Functor : public madness::FunctionFunctorInterface<double, 3> {
private:
    GaussianFunction func;
    std::vector<madness::coord_3d> centers;
public:
    Gaussian_Functor(GaussianFunction func, std::vector<madness::coord_3d> centers) : func(func), centers(centers) {}

    double operator()(const madness::coord_3d& r) const {
        return func(r);
    }
   
    std::vector<madness::coord_3d> special_points() const override final {
       return centers;       
    }

    madness::Level special_level() const override final {
       // From Robert: 
       // Pick initial level such that average gap between quadrature points
       // will find a significant value
       const int N = 6; // looking for where exp(-a*x^2) < 10**-N
       const int K = 6; // typically the lowest order of the polyn
       const double log10 = std::log(10.0);
       const double log2 = std::log(2.0); 
       const double L = madness::FunctionDefaults<3>::get_cell_min_width(); 
       const double expnt = *(std::max_element(func.expcoeff.begin(), func.expcoeff.end()));
       const double a = expnt*L*L;
       //const double fac = pow(2.0/(expnt*madness::constants::pi), 3.0/2.0); 
       const double fac = 100.0; 
       double n = std::log(a/(4*K*K*(N*log10+std::log(fac))))/(2*log2);
       return (n < 2 ? 2.0 : std::ceil(n)+1); // Added in the +1 for paranoia's sake
    }
};

} // namespace slymer
#endif
