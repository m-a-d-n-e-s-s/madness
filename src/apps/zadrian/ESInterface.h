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
#include <bitset>
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

namespace Properties {
  /**
   * \brief Different properties that can be read from electronic structure codes.
   *
   * C-style bitflags via (\c std::bitset) are used for specifying the properties.
   *
   * Some properties might require reading others, and this framework is designed
   * to facilitate reading multiple properties in one go through the output file(s).
   */
  using Properties = std::bitset<5>;

  constexpr Properties None = 0; ///< No properties.
  constexpr Properties Basis = 1 << 0; ///< The basis set.
  constexpr Properties Atoms = 1 << 1; ///< The atoms & positions.
  constexpr Properties Energies = 1 << 2; ///< The MO energies.
  constexpr Properties MOs = 1 << 3; ///< The MO vector coefficients.
  constexpr Properties Occupancies = 1 << 4; ///< MO occupancies.
} // namespace Properties

/// Abstract base class for interfacing with electronic structure codes.
class ES_Interface {
protected:
  Properties::Properties my_properties; ///< The properties that have been read.
  BasisSet my_basis_set; ///< The basis set.
  Atoms my_atoms; ///< The atoms (symbols and positions, in angstroms).
  unsigned int my_lineardeps; ///< Number of linear dependencies in the basis
  madness::Tensor<double> my_energies; ///< Alpha molecular orbital energies 
  madness::Tensor<double> my_MOs; ///< Alpha molecular orbital expansions coefficients. Column is the MO, row is the basis function.
  madness::Tensor<double> my_occupancies; ///< Alpha molecular orbital occupancies.
  madness::Tensor<double> my_beta_energies; ///< Beta molecular orbital energies 
  madness::Tensor<double> my_beta_MOs; ///< Beta molecular orbital expansions coefficients. Column is the MO, row is the basis function.
  madness::Tensor<double> my_beta_occupancies; ///< Beta molecular orbital occupancies.
 
public:
  std::reference_wrapper<std::ostream> err; ///< Output stream for messages.
  const Properties::Properties &properties; ///< Publically accessible list of read properties.
  const BasisSet &basis_set; ///< Publicly accessible basis set.
  const Atoms &atoms; ///< Publically accessible list of atoms.
  const unsigned int &lineardeps; ///< Publically accessible number of linear dependencies
  const madness::Tensor<double> &energies; ///< Publically accessible list of alpha MO energies.
  const madness::Tensor<double> &MOs; ///< Publically accessible alpha MO expansions coefficients. Column is the MO, row is the basis function.
  const madness::Tensor<double> &occupancies; ///< Publically accessible list of alpha MO occupancies (in eV).
  const madness::Tensor<double> &beta_energies; ///< Publically accessible list of beta MO energies (in eV).
  const madness::Tensor<double> &beta_MOs; ///< Publically accessible beta MO expansions coefficients. Column is the MO, row is the basis function.
  const madness::Tensor<double> &beta_occupancies; ///< Publically accessible list of beta MO occupancies (in eV).


  /// No default constructor.
  ES_Interface() = delete;

  /**
   * \brief Move constructor.
   *
   * \param[in] es The existing interface to move.
   */
  ES_Interface(ES_Interface &&es)
    : my_properties{std::move(es.my_properties)}, my_lineardeps{std::move(es.lineardeps)}, my_energies{std::move(es.my_energies)},
      my_MOs{std::move(es.my_MOs)}, my_occupancies{std::move(es.my_occupancies)},
      my_beta_energies{std::move(es.my_beta_energies)}, my_beta_MOs{std::move(es.my_beta_MOs)},
      my_beta_occupancies{std::move(es.my_beta_occupancies)}, err(es.err), properties(my_properties), 
      basis_set(my_basis_set), atoms(my_atoms), lineardeps(my_lineardeps), energies(my_energies), MOs(my_MOs), occupancies(my_occupancies),
      beta_energies(my_beta_energies), beta_MOs(my_beta_MOs), beta_occupancies(my_beta_occupancies)
  {}

  /**
   * \brief Copy constructor.
   *
   * \param[in] es The existing interface to copy.
   */
  ES_Interface(const ES_Interface &es)
    : my_properties{es.my_properties}, my_lineardeps{es.lineardeps}, my_energies{es.my_energies},
      my_MOs{es.my_MOs}, my_occupancies{es.my_occupancies}, my_beta_energies{es.my_beta_energies},
      my_beta_MOs{es.my_beta_MOs}, my_beta_occupancies{es.my_occupancies},
      err(es.err), properties(my_properties), basis_set(my_basis_set), atoms(my_atoms), lineardeps(my_lineardeps),
      energies(my_energies), MOs(my_MOs), occupancies(my_occupancies), beta_energies(my_beta_energies),
      beta_MOs(my_beta_MOs), beta_occupancies(my_beta_occupancies)
  {}

  /** 
   * \brief Constructor that sets the error/warning stream and the references.
   *
   * \param[in,out] err_ Output stream for messages. This can be updated later.
   */
  ES_Interface(std::ostream &err_)
    : my_properties{Properties::None}, my_lineardeps(0), my_energies(1), my_MOs(1, 1),
      my_occupancies(1), my_beta_energies(1), my_beta_MOs(1,1), my_beta_occupancies(1),
      err(err_), properties(my_properties), basis_set(my_basis_set), atoms(my_atoms), 
      lineardeps(my_lineardeps), energies(my_energies), MOs(my_MOs), occupancies(my_occupancies), 
      beta_energies(my_beta_energies), beta_MOs(my_beta_MOs), beta_occupancies(my_beta_occupancies)
  {}

  virtual ~ES_Interface() = default;

protected:
  /// Reset the interface.
  void reset() {
    my_properties = Properties::None;
    my_basis_set.clear();
    my_atoms.clear();
    my_lineardeps = 0;
    my_energies.reshape(1);
    my_MOs.reshape(1, 1);
    my_occupancies.reshape(1);
    my_beta_energies.reshape(1);
    my_beta_MOs.reshape(1, 1);
    my_beta_occupancies.reshape(1);
  }

public:
  /** 
   * \brief Read the specified properties and store them in the member variables.
   *
   * \param[in] props The properties to be read, using a bit flag combination.
   */
  virtual void read(Properties::Properties props) = 0;

};

} // namespace slymer

#endif
