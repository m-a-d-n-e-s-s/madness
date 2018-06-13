/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file Basis/polynomial.h
 * \brief API and helper functions for describing polynomial functions of
 *   three variables.
 *
 * The array is indexed by the polyIndex function. Each element of the array
 * is a coefficient on the polynomial.
 */

#ifndef __Basis_polynomial_h__
#define __Basis_polynomial_h__

#include <array>
#include <vector>

namespace slymer {

/**
 * \brief Array for storing coefficients of a polynomial of three variables
 *    with specified degree.
 *
 * The general polynomial of three variables is written as
 * \f[ \sum_{i,j,k} c_{i,j,k} x^i y^j z^k, \f]
 * where \f$0\le i+j+k\le N\f$. This data structure stores the coefficients
 * in a vector and provides access routines in terms of the exponents \f$i\f$,
 * \f$j\f$, and \f$k\f$.
 *
 * \note This class is not designed for polynomials with a large degree.
 */
class PolynomialCoeffs {

protected:
  /// The degree of the polynomial.
  unsigned degree;

  /**
   * \brief Array storing the coefficients.
   * 
   * For a given degree n, there are \f$(n+1)(n+2)/2\f$ monomials. The total
   * array thus has \f$(n+1)*(n+2)*(n+3)/6\f$ terms. The array stores all terms
   * with degree 0, then degree 1, then degree 2, etc.
   */
  std::vector<double> coeffs;

  /**
   * \brief Gets the array index for the coefficient of a particular term.
   *
   * Terms of degree 0 are stored before those of degree 1, etc.
   *
   * Within a particular degree, think of the problem in the combinatorial "stars
   * and bars" idea: there are n stars and two bars dividing them into three
   * parts. (Each part corresponds to the number of x, y, or z factors.) The first
   * array element will be bars in the 0 and 1 positions, then 0 and 2, ..., 0 and
   * n+1, then 1 and 2, etc.
   *
   * \param[in] pows The powers of x, y, and z for this term.
   * \return The array offset, as described above.
   */
  static inline unsigned polyIndex(const std::array<unsigned, 3> &pows) {
    const unsigned n = pows[0] + pows[1] + pows[2];

    // calculate the offset within the level of the degree...
    // x-y offset. All terms with pows[0] == 0 are first; there are n+1 of them.
    // Then pows[0] == 1; there are n, etc.
    unsigned ret = (n+1)*pows[0] - (pows[0]-1)*pows[0] / 2;
    // y-z offset. how many y factors are there? pows[1]. This value being 0
    // comes first, then 1, then 2, etc.
    ret += pows[1];

    // calculate the offset for smaller degrees. There are (n+1)(n+2)/2 elements
    // for degree n; add up all the smaller degrees.
    ret += n*(n+1)*(n+2)/6;

    return ret;
  }

public:
  /**
   * \brief Get the degree of the polynomial.
   *
   * \return The degree of the polynomial.
   */
  unsigned get_degree() const {
    return degree;
  }

  PolynomialCoeffs() = delete;

  /**
   * \brief Create a linear function with the specified coefficients.
   *
   * \param[in] constant The constant term.
   * \param[in] x The coefficient on x.
   * \param[in] y The coefficient on y.
   * \param[in] z The coefficient on z.
   */
  PolynomialCoeffs(const double constant, const double x, const double y,
      const double z) : degree(1), coeffs(4) {

    coeffs[polyIndex({{0,0,0}})] = constant;
    coeffs[polyIndex({{1,0,0}})] = x;
    coeffs[polyIndex({{0,1,0}})] = y;
    coeffs[polyIndex({{0,0,1}})] = z;
  }

  /**
   * \brief Create a polynomial of degree N.
   *
   * \param[in] N The degree of the polynomial.
   */
  PolynomialCoeffs(const unsigned N)
      : degree(N), coeffs((N+1)*(N+2)*(N+3)/6) {

    for(double &coeff : coeffs)
      coeff = 0.;
  }

  /**
   * \brief Access the coefficient for the specified term.
   *
   * \param[in] pows The powers of x, y, and z for this term.
   * \return Reference to the array term.
   */
  double &operator[] (const std::array<unsigned, 3> &pows) {
    return coeffs[polyIndex(pows)];
  }

  /**
   * \brief Access the coefficient for the specified term.
   *
   * \param[in] pows The powers of x, y, and z for this term.
   * \return The array term.
   */
  double operator[] (const std::array<unsigned, 3> &pows) const {
    return coeffs[polyIndex(pows)];
  }

  /**
   * \brief Evaluate the polynomial at the specified point.
   *
   * \note This implemention is a bit naive; it is not intended for
   *    polynomials with a large degree.
   *
   * \param[in] pt The point at which the polynomial is evaluated.
   * \return The value of the polynomial at the point.
   */
  double operator() (const std::array<double, 3> &pt) const;

  /**
   * \brief Add two polynomials together.
   *
   * *this is one of the polynomials.
   *
   * \param[in] rhs The other polynomial.
   * \return The sum of *this and rhs.
   */
  PolynomialCoeffs operator+ (const PolynomialCoeffs &rhs) const;

  /**
   * \brief Adds another polynomial to this one.
   *
   * \param[in] rhs The other polynomial.
   * \return The sum of *this and rhs, in *this.
   */
  PolynomialCoeffs &operator+= (const PolynomialCoeffs &rhs);

  /**
   * \brief Multiply a polynomial by a constant value.
   *
   * *this is one of the polynomials.
   *
   * \param[in] c The constant value.
   * \return The product of c and the polynomial.
   */
  PolynomialCoeffs operator* (const double c) const;

  /**
   * \brief Scale the polynomial by a constant.
   *
   * \param[in] c The scale factor.
   * \return The scaled polynomial (*this).
   */
  PolynomialCoeffs &operator*= (const double c);

  /**
   * \brief Multiply two polynomials together.
   *
   * *this is one of the polynomials.
   *
   * \param[in] rhs The other polynomial.
   * \return The product of *this and rhs.
   */
  PolynomialCoeffs operator* (const PolynomialCoeffs &rhs) const;

  /**
   * \brief Expands a linear polynomial raised to a power.
   *
   * Essentially expands the polynomial
   * \f[ (coeffs[0] + coeffs[1]*x + coeffs[2]*y + coeffs[3]*z)^j \f]
   * into a new PolynomialCoeffs object. That is, output[{{p,q,r}}] is the
   * coefficient on the \f$x^p y^q z^r\f$ term.
   *
   * This function requires multinomial coefficients, which involve
   * factorials. We assume that the arguments to the factorials will not
   * be large and hard-code small values. An error will be thrown if
   * the maximum value is exceeded.
   *
   * \throw std::runtime_error if the polynomial does not have degree 1.
   *
   * \param[in] j The exponent of the expansion.
   * \return The polynomial, ordered using the polyIndex scheme.
   */
  PolynomialCoeffs pow(const unsigned j) const;

}; // class PolynomialCoeffs

/**
 * \brief Multiply a polynomial by a constant value.
 *
 * \param[in] c The constant value.
 * \param[in] poly The polynomial.
 * \return The product of c and the polynomial.
 */
inline PolynomialCoeffs operator* (const double c, const PolynomialCoeffs &poly) {
  return poly * c;
}

} // namespace slymer

#endif
