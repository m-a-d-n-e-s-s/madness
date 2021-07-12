/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file ESInterface/NWChem.h
 * \brief API and helper routines for interfacing with NWChem.
 */

#ifndef __ESInterface_NWChem_h__
#define __ESInterface_NWChem_h__

#include <iostream>
#include <memory>
#include "ESInterface.h"
#include "gaussian.h"

namespace slymer {

/// Class for interfacing with NWChem (tested on version 6.6).
class NWChem_Interface : public ES_Interface {
protected:
  /// The base file name of the NWChem output.
  std::string my_fname;

  /// Storage for the actual basis functions.
  std::vector<std::unique_ptr<GaussianFunction>> gaussians;

public:
  NWChem_Interface() = delete;

  /// Publically-accessible version of the file name.
  const std::string &fname;

  /**
   * \brief Wrap the output of a NWChem computation.
   *
   * The parameter is the base file name (potentially including path) for the
   * NWChem output files. For example, the NWChem log file will be fname.out,
   * the molecular orbitals file will be fname.movecs, etc.
   *
   * \param[in] fname_ Base file name for the NWChem computation.
   * \param[in,out] err_ Output stream for error or warning messages.
   */
  NWChem_Interface(const std::string &fname_, std::ostream &err_)
    : ES_Interface(err_), my_fname(fname_), gaussians(0), fname(my_fname) {}

  /**
   * \brief Changes the base file name.
   *
   * \param[in] fname_ The new base file name for the NWChem computation.
   */
  void reset(const std::string &fname_) {
    ES_Interface::reset();
    my_fname = fname_;
    gaussians.clear();
  }

  /** 
   * \brief Read the specified properties and store them in the member variables.
   *
   * \throw std::invalid_argument if the NWChem logfile (fname.out) cannot be
   *    opened.
   * \throw std::runtime_error if there is an error reading the NWChem output
   *    file.
   *
   * \param[in] props The properties to be read, using a bit flag combination.
   */
  virtual void read(Properties::Properties props) override;

protected:
  /**
   * \brief Extract and store the atom types and positions.
   *
   * \param[in,out] in The stream containing the NWChem output file.
   */
  void read_atoms(std::istream &in);

  /**
   * \brief Extract and store the basis set.
   *
   * \todo Reading Cartesian-type d and f shells has not been tested.
   * \note Only shells of type s, p, d, f, g and h are implemented.
   *
   * \param[in,out] in The stream containing the NWChem output file.
   */
  void read_basis_set(std::istream &in);

  /**
   * \brief Read the NWChem movecs file containing occupation numbers,
   *    MO energies, and MO coefficients.
   *
   * \param[in] props Properties to store from the read (from the list above).
   * \param[in,out] in The stream for the NWChem output movecs file.
   */
  void read_movecs(const Properties::Properties props, std::istream &in);
};

} // namespace slymer

#endif
