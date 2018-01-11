/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file ESInterface/ESInterface.h
 * \brief API and helper routines for interfacing with electronic structure codes.
 *
 * Provides a common (abstract) framework for extracting information from
 * external electronic structure calculations.
 */

#ifndef __ESInterface_ESInterface_h__
#define __ESInterface_ESInterface_h__

#include <madness/mra/mra.h>
#include <array>
#include <functional>
#include <ostream>
#include <string>
#include "basis.h"

namespace slymer {

/// An atom (symbol and position).
struct Atom {
  std::string symbol; ///< The atom's symbol.
  std::array<double, 3> position; ///< The atom's location, in angstroms.
};

/// A set of atoms.
using Atoms = std::vector<Atom>;

/// Abstract base class for interfacing with electronic structure codes.

/// \todo Create a copy and move constructor, since they're deleted (by default) due to the presence of the references.
class ES_Interface {
public:
  /**
   * \brief Different properties that can be read from electronic structure codes.
   *
   * Some properties might require reading others, and this framework is designed
   * to facilitate reading multiple properties in one go through the output file(s).
   */
  enum class Properties : unsigned {
    None = 0,
    Basis = 1,
    Atoms = 2,
    Energies = 4,
    MOs = 8,
    Occupancies = 16
  };

protected:
  Properties my_properties; ///< The properties that have been read.
  BasisSet my_basis_set; ///< The basis set.
  Atoms my_atoms; ///< The atoms (symbols and positions, in angstroms).
  madness::Tensor<double> my_energies; ///< Molecular orbital energies (in eV).
  madness::Tensor<double> my_MOs; ///< Molecular orbital expansions coefficients. Column is the MO, row is the basis function.
  madness::Tensor<double> my_occupancies; ///< Molecular orbital occupancies.

public:
  std::reference_wrapper<std::ostream> err; ///< Output stream for messages.
  const Properties &properties; ///< Publically accessible list of read properties.
  const BasisSet &basis_set; ///< Publicly accessible basis set.
  const Atoms &atoms; ///< Publically accessible list of atoms.
  const madness::Tensor<double> &energies; ///< Publically accessible list of MO energies (in eV).
  const madness::Tensor<double> &MOs; ///< Publically accessible MO expansions coefficients. Column is the MO, row is the basis function.
  const madness::Tensor<double> &occupancies; ///< Publically accessible list of MO occupancies (in eV).

  ES_Interface() = delete;

  /** 
   * \brief Constructor that sets the error/warning stream and the references.
   *
   * \param[in,out] err_ Output stream for messages. This can be updated later.
   */
  ES_Interface(std::ostream &err_)
    : my_properties{Properties::None}, my_energies(1), my_MOs(1, 1),
      my_occupancies(1), err(err_), properties(my_properties), 
      basis_set(my_basis_set), atoms(my_atoms), energies(my_energies), 
      MOs(my_MOs), occupancies(my_occupancies)
  {}

  virtual ~ES_Interface() = default;

protected:
  /// Reset the interface.
  void reset() {
    my_properties = Properties::None;
    my_basis_set.clear();
    my_atoms.clear();
    my_energies.reshape(1);
    my_MOs.reshape(1,1);
    my_occupancies.reshape(1);
  }

public:
  /** 
   * \brief Read the specified properties and store them in the member variables.
   *
   * \param[in] props The properties to be read, using a bit flag combination.
   */
  virtual void read(const Properties props) = 0;

};

/// \cond nodoc
inline ES_Interface::Properties operator| (ES_Interface::Properties lhs, ES_Interface::Properties rhs) {
  return static_cast<ES_Interface::Properties>(static_cast<unsigned>(lhs) | static_cast<unsigned>(rhs));
}

inline ES_Interface::Properties operator& (ES_Interface::Properties lhs, ES_Interface::Properties rhs) {
  return static_cast<ES_Interface::Properties>(static_cast<unsigned>(lhs) & static_cast<unsigned>(rhs));
}
/// \endcond

} // namespace slymer

#endif
