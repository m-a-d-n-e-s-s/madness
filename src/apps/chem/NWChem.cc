/* This file is a part of Slymer, which is distributed under the Creative
   Commons Attribution-NonCommercial 4.0 International Public License.

   (c) 2017 Stony Brook University. */

/**
 * \file ESInterface/NWChem.cc
 * \brief Implementation of the interface for NWChem.
 */

#include <madness/mra/mra.h>
#include "NWChem.h"
#include <cstdint>
#include <fstream>
#include <list>
#include <map>
#include <iostream>
#include <regex>
#include <stdexcept>
#include "gaussian.h"

namespace slymer {

// helper functions

void NWChem_Interface::read(Properties::Properties props) {
  // figure out which properties we need to read:
  // 1) Work out dependencies (for instance, basis needs atoms).
  // 2) Omit ones that are already read. That is, we need properties that are
  //    specified in props (post step 1), but not already read (in my_properties).
  //    ... bitwise and props with (not my_properties) -- has to be requested and
  //    not read.

  // basis requires atoms
  if((props & Properties::Basis).any())
    props |= Properties::Atoms;

  // if we're doing the MO expansion coefficients, we might as well do the
  // energies
  if((props & Properties::MOs).any())
    props |= Properties::Energies;

  // skip things we've already read
  props &= ~my_properties;

  // do we need to open and read fname.out?
  if((props & Properties::Atoms).any() ||
     (props & Properties::Basis).any()) {

    // open the NWChem log file: fname.out
    std::ifstream in(fname + ".out");
    if(!in)
      throw std::invalid_argument("Cannot open " + fname + ".out for reading.");

    if((props & Properties::Atoms).any())
      read_atoms(in);
    if((props & Properties::Basis).any())
      read_basis_set(in);

    in.close();
  }

  // do we need to open and read fname.movecs?
  if((props & Properties::Energies).any() ||
     (props & Properties::MOs).any() ||
     (props & Properties::Occupancies).any()) {

    // open the NWChem movecs file: fname.movecs
    // this file is binary
    std::ifstream in(fname + ".movecs", std::ios::binary);
    if(!in)
      throw std::invalid_argument("Cannot open " + fname + ".movecs for reading.");
    in.unsetf(std::ios::skipws);

    read_movecs(props, in);

    in.close();
  }
}

void NWChem_Interface::read_atoms(std::istream &in) {
  my_atoms.clear();

  // state variables
  std::string line; // current line of the NWChem output
  bool reading = false; // are we in the block of text where we can read basis functions?
  double unitcf = 1.; // Conversion factor to go from units in the NWChem file to angstroms

  // regex search strings
  std::regex startline{R"(Geometry)"},
             unitline{R"(Output coordinates in (a.u.|angstroms|nanometer|picometer|.+) \(scale)"},
             atomline{R"([[:digit:]]+[[:space:]]+([[:alpha:]]+)[[:space:]]+[^[:space:]]+[[:space:]]+([^[:space:]]+)[[:space:]]+([^[:space:]]+)[[:space:]]+([^[:space:]]+))"},
             doneline{R"(Atomic Mass)"};
  std::smatch matches;

  // note that NWChem outputs the orbitals by atom type, listing ang. momentum, exponents, and coefficients
  while(std::getline(in, line)) {
    // is this the start of the basis function output?
    if(!reading && std::regex_search(line, startline)) {
      reading = true;
    }

    // in what units are the atomic positions?
    else if(reading && std::regex_search(line, matches, unitline)) {
      bool doprint = true;
      if(matches[1] == "angstroms")
        unitcf = 1./0.529177;
      else if(matches[1] == "a.u.")
        unitcf = 1.0;
      else if(matches[1] == "nanometer")
        unitcf = 10./0.529177;
      else if(matches[1] == "picometer")
        unitcf = 0.01/0.529177;
      else {
        err.get() << "Unknown units: " << matches[1] << ". Assuming angstroms." << std::endl;
        unitcf = 1.;
        doprint = false;
      }
      if(doprint)
        err.get() << "NWChem uses " << matches[1] << ". Conversion to a.u. is "
          << unitcf << '.' << std::endl;
    }

    // does this line have an atomic position?
    else if(reading && std::regex_search(line, matches, atomline)) {
      err.get() << matches[1] << " at (" << matches[2] << ", " << matches[3] << ", " << matches[4] << ")." << std::endl;
      Atom addme;
      addme.symbol = matches[1];
      addme.position = {{std::stod(matches[2]) * unitcf, std::stod(matches[3]) * unitcf, std::stod(matches[4]) * unitcf}};
      my_atoms.emplace_back(std::move(addme));
    }

    // are we done?
    else if(reading && std::regex_search(line, doneline)) {
      break;
    }
  }

  my_properties = my_properties | Properties::Atoms;
}



// helper functions and data types for reading the basis set

/**
 * \brief Store information on a shell of basis functions as it's being read.
 *
 * Used when reading the basis functions.
 */
struct BasisShell {
  char type; ///< Type of basis function.
  std::vector<double> coeffs, ///< Expansion coefficient for each primitive.
                      exps; ///< Exponent for each primitive.
};

static std::ostream &operator<< (std::ostream &out, const BasisShell &bs) {
  out << "Type " << bs.type << ", exps = [ ";
  for(auto ex : bs.exps)
    out << ex << ' ';
  out << "], coeffs = [ ";
  for(auto coeff : bs.coeffs)
    out << coeff << ' ';
  out << "]";
  return out;
}

void NWChem_Interface::read_basis_set(std::istream &in) {
  my_basis_set.clear();
  // First we have to read the basis functions (and types) that are on each atom

  // state variables
  std::string line; // current line of the NWChem output
  bool reading = false; // are we in the block of text where we can read basis functions?
  bool spherical = true; // are we using spherical orbitals (true) or Cartesian (false)?
  std::string curatom; // Symbol of the current atom being read.
  unsigned curnum = 1; // The current basis function number (for the symbol) being read.
  BasisShell curfunc; // The current basis function being read.
  std::map<std::string, std::list<BasisShell>> basisperatom; // List of basis functions for the given symbol (map key).
  std::list<BasisShell> basiscuratom; // List of basis functions for the current atom.

  // regex search strings
  std::regex startline{R"(Basis.*(cartesian|spherical))"},
             newatomline{R"(([[:alpha:]]+)[[:space:]]+\(([[:alpha:]]+)\))"},
             orbitalline{R"(([[:digit:]]+)[[:space:]]+([[:alpha:]]+)[[:space:]]+([^[:space:]]+)[[:space:]]+([^[:space:]]+))"},
             doneline{R"(Summary)"};
  std::smatch matches;

  // note that NWChem outputs the orbitals by atom type, listing ang. momentum, exponents, and coefficients
  while(std::getline(in, line)) {
    // is this the start of the basis function output?
    if(!reading && std::regex_search(line, matches, startline)) {
      reading = true; // we're now in the basis function block

      // determine the type of basis function being used
      if(matches[1] == "cartesian")
        spherical = false;
      else if(matches[1] == "spherical")
        spherical = true;

      err.get() << "Using " << (spherical ? "spherical" : "Cartesian") << " orbitals." << std::endl; 
    }

    // Are we done reading orbitals?
    else if(reading && std::regex_search(line, doneline)) {
      if(curatom.length() > 0) {
        err.get() << "   Function " << curnum << ": " << curfunc << std::endl;
        basiscuratom.emplace_back(std::move(curfunc));
        basisperatom.emplace(std::move(curatom), std::move(basiscuratom));
      }
      break; // nothing more to read
    }

    // Is this the start of a new atom?
    else if(reading && std::regex_search(line, matches, newatomline)) {
      // if this isn't the first symbol, store the previous symbol's basis and reset
      if(curatom.length() > 0) {
        err.get() << "   Function " << curnum << ": " << curfunc << std::endl;

        basiscuratom.emplace_back(std::move(curfunc));
        basisperatom.emplace(std::move(curatom), std::move(basiscuratom));
      }

      curatom = matches[1];
      curnum = 1;
      err.get() << "Reading basis for " << curatom << " (" << matches[2] << ")." << std::endl;
    }

    else if(reading && std::regex_search(line, matches, orbitalline)) {
      // this line is an uncontracted orbital shell
      const unsigned func = std::stoul(matches[1]);
      if(func == curnum + 1) {
        // store the old basis function and prepare to start reading a new one
        err.get() << "   Function " << curnum << ": " << curfunc << std::endl;
        basiscuratom.emplace_back(std::move(curfunc));
        ++curnum;
      }

      if(func == curnum) {
        // add a new primitive (uncontracted Gaussian) to this overall shell
        const std::string type = matches[2];
        if(curfunc.exps.size() == 0)
          // first one. set the type
          curfunc.type = type[0];
        else if(curfunc.type != type[0])
          throw std::runtime_error("Inconsistent orbital types while reading NWChem basis set.");

        curfunc.exps.push_back(std::stod(matches[3]));
        curfunc.coeffs.push_back(std::stod(matches[4]));
      }
      else
        throw std::runtime_error("Inconsistent basis function counting while reading NWChem basis set.");
    }
  }

  // now we need to go through the atoms in the calculation and create the basis
  // functions centered on each.
  err.get() << "Processing basis functions for each atom." << std::endl;

  for(const Atom &atom : atoms) {
    unsigned nfunc = 0;

    err.get() << "Atom: " << atom.symbol << " at (" << atom.position[0] <<
      ", " << atom.position[1] << ", " << atom.position[2] << "): ";

    // add the basis functions for this atom type, centered at this location
    for(const BasisShell &bs : basisperatom[atom.symbol]) {
      if(bs.type == 'S') {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::s, atom.position, bs.exps,
            bs.coeffs));
        ++nfunc;
      }
      else if(bs.type == 'P') {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::px, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::py, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::pz, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 3;
      }
      else if(bs.type == 'D' && !spherical) {
        throw std::runtime_error("Cartesian D shell requested: These are not tested or verified.");
        
        // THESE ROUTINES HAVE NOT BEEN TESTED
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dxx, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dxy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dxz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dzz, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 6;
      }
      else if(bs.type == 'D' && spherical) {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dxy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dzzmrr, atom.position, bs.exps,
            bs.coeffs));

        std::unique_ptr<GaussianFunction> dxz(
          new GaussianFunction(GaussianType::dxz, atom.position, bs.exps,
            bs.coeffs));
        dxz->operator*=(-1.);
        gaussians.emplace_back(std::move(dxz));

        gaussians.emplace_back(
          new GaussianFunction(GaussianType::dxxmyy, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 5;
      }
      else if(bs.type == 'F' && !spherical) {
        throw std::runtime_error("Cartesian F shell requested: These are not tested or verified.");
        
        // THESE ROUTINES HAVE NOT BEEN TESTED
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxxx, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxxy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxxz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fzzz, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 10;
      }
      else if(bs.type == 'F' && spherical) {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxxymyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fyzzmrry, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fzzzmrrz, atom.position, bs.exps,
            bs.coeffs));

        std::unique_ptr<GaussianFunction> fxzzmrrx(
          new GaussianFunction(GaussianType::fxzzmrrx, atom.position, bs.exps,
            bs.coeffs));
        fxzzmrrx->operator*=(-1.);
        gaussians.emplace_back(std::move(fxzzmrrx));

        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxxzmyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::fxyymxxx, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 7;
      }
      else if(bs.type == 'G' && !spherical) {
        throw std::runtime_error("Cartesian G shell requested: These are not tested or verified.");

        // THESE FUNCTIONS HAVE NOT BEEN TESTED OR VERIFIED!!!
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxxx, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxxy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxxz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxxzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gzzzz, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 15;
      }
      else if(bs.type == 'G' && spherical) {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxydx2my2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyzdx2my2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxydz2mr2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gyzdz2mr2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gzero, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxzdz2mr2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gx2my2dz2mr2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gxzdx2my2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::gx4mx2y2py4, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 9;
      }
      else if(bs.type == 'H' && !spherical) {
        throw std::runtime_error("Cartesian H shell requested: These are not tested or verified.");

        // THESE FUNCTIONS HAVE NOT BEEN TESTED OR VERIFIED!!! 
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxxx, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxxy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxxz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxxzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxxzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxyyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxyyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxyyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxyzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hxzzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hyyyyy, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hyyyyz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hyyyzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hyyzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hyzzzz, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hzzzzz, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 21;
      }
      else if(bs.type == 'H' && spherical) {
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hm5, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hm4, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hm3, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hm2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hm1, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hzero, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hp1, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hp2, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hp3, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hp4, atom.position, bs.exps,
            bs.coeffs));
        gaussians.emplace_back(
          new GaussianFunction(GaussianType::hp5, atom.position, bs.exps,
            bs.coeffs));
        nfunc += 11;
      }
      else
        throw std::runtime_error(std::string("Unknown or unimplemented shell type: ") + bs.type);
    }

    err.get() << nfunc << " basis functions." << std::endl;
  }

  unsigned size = gaussians.size();
  for(unsigned j = 0; j < size; ++j) {
    // add a {reference to the basis function} to the basis set vector
    my_basis_set.emplace_back(std::ref(*gaussians[j].get()));
  }

  my_properties = my_properties | Properties::Basis;
}



// helper routine for read_movecs

/**
 * \brief Read bytes from a binary file into the specified type.
 *
 * Attempts to account for the endian-ness of the binary data. This should
 * hopefully make the routine more robust to running NWChem jobs on one
 * machine and processing the files on another.
 *
 * \tparam T Type of data to be read (only really need the size).
 * \param [in,out] in The input stream.
 * \param [in] swap Change the endian-ness of the read data?
 * \return The read data, in the requested type.
 */
template<typename T>
static T read_endian(std::istream &in, const bool swap) {
  T ret;
  char * const ptr = reinterpret_cast<char *>(&ret);
  const unsigned bytes = sizeof(T);

  if(swap) {
    // need to swap the endian-ness of the data
    // read byte-by-byte, in reverse order
    for(unsigned j = 0; j < bytes; ++j)
      in.read(ptr + (bytes - j - 1), 1);
  }
  else
    // use the data as-is
    in.read(ptr, bytes);

  return ret;
}

void NWChem_Interface::read_movecs(const Properties::Properties props,
  std::istream &in)
{
  int32_t num[2];
  bool swap_endian = false;
  unsigned nsets, nbasis;
  const std::runtime_error errmess("Error reading NWChem movecs file.");

  // what properties are we storing? (all will be read from the movecs file)
  const bool do_occupancies = (props & Properties::Occupancies).any();
  const bool do_energies = (props & Properties::Energies).any();
  const bool do_MOs = (props & Properties::MOs).any();

  // temporary storage places (in case errors pop up while/after reading)
  madness::Tensor<double> temp_occupancies, temp_energies;
  madness::Tensor<double> temp_MOs;

  // this function is based on the mov2asc program
  // https://github.com/jeffhammond/nwchem/blob/master/contrib/mov2asc/mov2asc.F

  // every write statement in fortran puts a 4-byte integer before and after
  // the data. this integer contains the size of the data. the num array here
  // is used to read these values.

  // read the header information
  // first up is a FORTRAN integer (4 bytes) telling the number of characters
  // this should be 142 = 3*32 + 20 + 26 or 110 = 2*32 + 20 + 26
  num[0] = read_endian<int32_t>(in, false);
  if(num[0] != 142 && num[0] != 110) {
    // try the other endian-ness
    err.get() << "Unable to read header. Trying opposite endian-ness." << std::endl;
    char * const ptr = reinterpret_cast<char *>(num);
    std::swap(ptr[0], ptr[3]);
    std::swap(ptr[1], ptr[2]); // all four of these chars make up num[0]
    swap_endian = true;
    if(num[0] != 142 && num[0] != 110) {
      err.get() << "Still unable to read header. Abort." << std::endl;
      throw errmess;
    }
    else
      err.get() << "Using opposite endian-ness for the remainder of this read." << std::endl;
  }
  // basissum, geomsum, bqsum, scftype20, and date. bqsum may be skipped
  in.seekg(num[0], std::ios_base::cur);
  // bookend for the length of the above strings combined
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[1] != num[0]) {
    err.get() << "Error reading header (basissum, ..., date)." << std::endl;
    throw errmess;
  }

  // for whatever reason the scftype20 is repeated, with bookends (now = 20)
  num[0] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != 20) {
    err.get() << "Error reading header (scftype20)." << std::endl;
    throw errmess;
  }
  // scftype20
  in.seekg(20, std::ios_base::cur);
  num[0] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != 20) {
    err.get() << "Error reading header (scftype20)." << std::endl;
    throw errmess;
  }
  
  // size of the title (lentit)
  num[0] = read_endian<int32_t>(in, swap_endian);
  // lentit. don't check size explicitly, it's probably either 4 or 8
  // depending on the fortran compiler/environment
  in.seekg(num[0], std::ios_base::cur);
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading header (lentit)." << std::endl;
    throw errmess;
  }

  // read the title itself
  num[0] = read_endian<int32_t>(in, swap_endian);
  // the title. varying length
  in.seekg(num[0], std::ios_base::cur);
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading title of NWChem job in header." << std::endl;
    throw errmess;
  }

  // length of the basis set name (lenbas)
  num[0] = read_endian<int32_t>(in, swap_endian);
  // lentit. don't check size explicitly, it's probably either 4 or 8
  // depending on the fortran compiler/environment
  in.seekg(num[0], std::ios_base::cur);
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading header (lenbas)." << std::endl;
    throw errmess;
  }

  // read the basis name itself
  num[0] = read_endian<int32_t>(in, swap_endian);
  // the title. varying length
  in.seekg(num[0], std::ios_base::cur);
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading basis set name in header." << std::endl;
    throw errmess;
  }

  // the number of MO sets (nsets)
  num[0] = read_endian<int32_t>(in, swap_endian);
  if(num[0] == 4) {
    int32_t temp;
    temp = read_endian<int32_t>(in, swap_endian);
    nsets = static_cast<unsigned>(temp);
  }
  else if(num[0] == 8) {
    int64_t temp;
    temp = read_endian<int64_t>(in, swap_endian);
    nsets = static_cast<unsigned>(temp);
  }
  else {
    err.get() << "Unknown or unimplemented integer type for the number of sets in header." << std::endl;
    throw errmess;
  }
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading number of sets in header." << std::endl;
    throw errmess;
  }
  if(nsets != 1) {
    err.get() << "The read_MOs basis has only been tested for movecs files that have a single set.\n" <<
      "The results from this function need to be verified before use (" << nsets << " sets).\n"
      "Only the last set will be kept." << std::endl;
  }

  // size of the basis set (NBF in the mov2asc script)
  num[0] = read_endian<int32_t>(in, swap_endian);
  if(num[0] == 4) {
    int32_t temp;
    temp = read_endian<int32_t>(in, swap_endian);
    nbasis = static_cast<unsigned>(temp);
  }
  else if(num[0] == 8) {
    int64_t temp;
    temp = read_endian<int64_t>(in, swap_endian);
    nbasis = static_cast<unsigned>(temp);
  }
  else {
    err.get() << "Unknown or unimplemented integer type for basis set size in header." << std::endl;
    throw errmess;
  }
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != num[1]) {
    err.get() << "Error reading basis set size in header." << std::endl;
    throw errmess;
  }
  // make sure this agrees with any information we've already read
  if((properties & Properties::Basis) != Properties::None && nbasis != basis_set.size()) {
    err.get() << "The number of basis functions in read_MOs does not match the number of " <<
      "basis functions read\nfrom the NWChem log file. Use the following results with caution."
      << std::endl;
  }

  // read the number of vectors in each set (nmo)
  std::vector<unsigned> nmo(nsets);
  for(unsigned j = 0; j < nsets; ++j) {
    num[0] = read_endian<int32_t>(in, swap_endian);
    if(num[0] == 4) {
      int32_t temp;
      temp = read_endian<int32_t>(in, swap_endian);
      nmo[j] = static_cast<unsigned>(temp);
    }
    else if(num[0] == 8) {
      int64_t temp;
      temp = read_endian<int64_t>(in, swap_endian);
      nmo[j] = static_cast<unsigned>(temp);
    }
    else {
      err.get() << "Unknown or unimplemented integer type for nmo[" << j << "] in header." << std::endl;
      throw errmess;
    }
    num[1] = read_endian<int32_t>(in, swap_endian);
    if(num[0] != num[1]) {
      err.get() << "Error reading nmo[" << j << "] in header." << std::endl;
      throw errmess;
    }
  }

  // allocate space to store the occupation numbers, the eigenvalues, and the
  // eigenvectors (MO vectors), as desired by the request
  if(do_occupancies) {
    madness::Tensor<double> one(nmo[nsets-1]);
    temp_occupancies = copy(one);
  }
  if(do_energies) {
    madness::Tensor<double> two(nmo[nsets - 1]);
    temp_energies = copy(two);
  }
  if(do_MOs) {
    madness::Tensor<double> three(nmo[nsets - 1], nmo[nsets - 1]);
    temp_MOs = copy(three);
  }

  // go through the sets
  for(unsigned set = 0; set < nsets; ++set) {
    // first read the occupancies
    // number of bits bookend (8 for double * nmo[set]);
    num[0] = read_endian<int32_t>(in, swap_endian);
    if(num[0] != static_cast<int>(8*nmo[set])) {
      err.get() << "Unexpected number of occupancies for set " << set << '.'
        << std::endl;
      throw errmess;
    }
    if(do_occupancies && set == nsets - 1)
      for(unsigned j = 0; j < nmo[set]; ++j)
        temp_occupancies[j] = read_endian<double>(in, swap_endian);
    else
      in.seekg(num[0], std::ios_base::cur); // just buzz past them.
    num[1] = read_endian<int32_t>(in, swap_endian);
    if(num[0] != num[1]) {
      err.get() << "Error reading occupancies for set " << set << '.' << std::endl;
      throw errmess;
    }


    // next up are the eigenvalues (energies)
    // number of bits bookend (8 for double * nmo[set]);
    num[0] = read_endian<int32_t>(in, swap_endian);
    if(num[0] != static_cast<int>(8*nmo[set])) {
      err.get() << "Unexpected number of energies for set " << set << '.'
        << std::endl;
      throw errmess;
    }
    if(do_energies && set == nsets - 1)
      for(unsigned j = 0; j < nmo[set]; ++j)
        // NWChem reports energies in Hartrees, we want eV.
        temp_energies[j] = read_endian<double>(in, swap_endian) * 27.21138602;
    else
      in.seekg(num[0], std::ios_base::cur); // just buzz past them.
    num[1] = read_endian<int32_t>(in, swap_endian);
    if(num[0] != num[1]) {
      err.get() << "Error reading energies for set " << set << '.' << std::endl;
      throw errmess;
    }


    // finally, read the MO vectors, which were written vector-by-vector
    for(unsigned mo = 0; mo < nmo[set]; ++mo) {
      // bookend size of the vector
      num[0] = read_endian<int32_t>(in, swap_endian);
      if(num[0] != static_cast<int>(8*nbasis)) {
        err.get() << "Unexpected number of coefficients in MO " << mo << " of set "
          << set << '.' << std::endl;
        throw errmess;
      }
      if(do_MOs && set == nsets - 1)
        for(unsigned coeff = 0; coeff < nbasis; ++coeff)
          temp_MOs(coeff, mo) = read_endian<double>(in, swap_endian);
      else
        in.seekg(num[0], std::ios_base::cur); // just buzz past the MO
      // bookend size of the vector
      num[1] = read_endian<int32_t>(in, swap_endian);
      if(num[0] != num[1]) {
        err.get() << "Error reading coefficients of MO " << mo << " of set "
          << set << '.' << std::endl;
        throw errmess;
      }
    }
  }

  // at the very end of the movecs file is the total energy and nuclear repulsion
  // energy of the system
  num[0] = read_endian<int32_t>(in, swap_endian);
  if(num[0] != 16) {
    err.get() << "Unexpected size reading the footer." << std::endl;
    throw errmess;
  }
  read_endian<double>(in, swap_endian); // total energy
  read_endian<double>(in, swap_endian); // nuclear repulsion
  num[1] = read_endian<int32_t>(in, swap_endian);
  if(num[1] != 16) {
    err.get() << "Error reading the footer." << std::endl;
    throw errmess;
  }

  // move/swap the placeholders into the class's storage space
  if(do_occupancies) {
    my_occupancies = std::move(temp_occupancies);
    my_properties = my_properties | Properties::Occupancies;
  }
  if(do_energies) {
    my_energies = std::move(temp_energies);
    my_properties = my_properties | Properties::Energies;
  }
  if(do_MOs) {
    my_MOs = std::move(temp_MOs);
    my_properties = my_properties | Properties::MOs;
  }
}

} // namespace slymer
