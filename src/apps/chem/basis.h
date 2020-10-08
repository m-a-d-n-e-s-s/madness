/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file Basis/basis.h
 * \brief Basis function API and routines.
 *
 * Sets up the interface for general basis functions and related calculations.
 */

#ifndef __Basis_basis_h__
#define __Basis_basis_h__

#include <array>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

namespace slymer {

/// Abstract base class for generic basis functions.
class BasisFunction {
public:
  virtual ~BasisFunction() = default;
  
  /**
   * \brief Evaluate the basis function at the specified point.
   *
   * \param[in] x The point.
   * \return The basis function evaluated at the point x.
   */
  virtual double operator() (const std::array<double, 3> &x) const = 0;

};

/// Type for a basis set (collection of basis functions).
using BasisSet = std::vector<std::reference_wrapper<BasisFunction>>;

/**
 * \brief Convert a generic basis set to basis functions with the specific
 *    type.
 *
 * \throw std::runtime_error if a basis function in the basis set is not
 *    of the specified type.
 *
 * \tparam T The intended type of each basis function.
 * \param[in] bset The basis set.
 * \return A vector with cast references to the basis functions.
 */
template<typename T>
std::vector<std::reference_wrapper<const T>> cast_basis(const BasisSet &bset) {
  const unsigned size = bset.size();
  std::vector<std::reference_wrapper<const T>> casts;

  for(unsigned j = 0; j < size; ++j) {
    const T *ptr = dynamic_cast<const T*>(&bset[j].get());
    if(ptr == nullptr)
      throw std::runtime_error("A basis function is not of the specified type.");
    casts.emplace_back(std::cref(*ptr));
  }

  return casts;
}

} // namespace slymer

#endif
