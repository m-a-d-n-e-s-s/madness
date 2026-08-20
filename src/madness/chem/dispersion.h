/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680
 */

#ifndef SRC_MADNESS_CHEM_DISPERSION_H_
#define SRC_MADNESS_CHEM_DISPERSION_H_

#include <string>

#include <madness/tensor/tensor.h>
#include <madness/world/world.h>
#include <madness/chem/molecule.h>

namespace madness {

/// interface class to simple-dftd3, Grimme's D3 empirical dispersion correction

/// The D3 correction is a closed-form function of the nuclear coordinates
/// alone: it touches neither the density nor the Fock matrix, so it enters a
/// calculation as one additive number in the total energy and one additive
/// tensor in the nuclear gradient.
///
/// References
///  - S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
///    J. Chem. Phys. 132, 154104 (2010), doi:10.1063/1.3382344   [the D3 model]
///  - S. Grimme, S. Ehrlich, L. Goerigk,
///    J. Comput. Chem. 32, 1456 (2011), doi:10.1002/jcc.21759    [BJ damping]
///  - implementation: https://github.com/dftd3/simple-dftd3
///
/// Requires MADNESS to be configured against simple-dftd3
/// (`-DENABLE_DFTD3=ON`, the default, plus a discoverable install); a default
/// constructed object is inactive and returns zeros, but asking for an actual
/// correction in a build without the library throws.
class DispersionCorrection {
public:

    /// default ctor -- inactive, energy() is 0.0 and gradient() is all zeros
    DispersionCorrection() = default;

    /// construct from the `dft` input group

    /// @param[in]  spec        "none", "d3bj" or "d3zero" (the `dispersion` keyword)
    /// @param[in]  xc_line     the `xc` input line; its first token is the fallback
    ///                         source for the damping parameter set
    /// @param[in]  functional  explicit functional name (the `dispersion_functional`
    ///                         keyword); empty means derive it from @c xc_line
    /// @param[in]  atm         include the three-body Axilrod-Teller-Muto term
    DispersionCorrection(const std::string& spec, const std::string& xc_line,
                         const std::string& functional, const bool atm);

    /// true if a correction is actually being applied
    bool active() const { return damping != none; }

    /// short human-readable tag, e.g. "D3(BJ)/pbe0" -- "none" if inactive
    std::string description() const;

    /// print the method and its primary references; a no-op if inactive

    /// Prints at most once per object, so callers on a path that repeats per
    /// geometry step (Nemo::value) need no guard of their own.
    void print_citation(World& world) const;

    /// the dispersion energy in Hartree; 0.0 if inactive
    double energy(World& world, const Molecule& mol) const;

    /// d E_disp / d R in Hartree/bohr, length 3*natom, laid out [3*atom + axis]

    /// Same layout as SCF::derivatives and NemoBase::compute_gradient, so the
    /// result can be added to those tensors directly. All zeros if inactive.
    Tensor<double> gradient(World& world, const Molecule& mol) const;

    /// throw if `spec` asks for a correction, naming the engine that cannot apply it

    /// For engines whose energy expression has no dispersion term. Silently
    /// dropping the correction would report an uncorrected energy as the answer
    /// to a deck that asked for a corrected one. A no-op for "none"/"", so
    /// call sites need no guard of their own.
    ///
    /// Takes the raw `dispersion` keyword because not every such engine owns an
    /// SCF to ask -- Znemo carries a bare CalculationParameters.
    static void reject(const std::string& spec, const char* engine);

    /// same, for an already-constructed correction
    void reject(const char* engine) const { reject(description(), engine); }

    /// true if this build can compute a correction at all
    static bool available();

    /// linked simple-dftd3 version as "major.minor.patch", or "" if unavailable
    static std::string library_version();

private:

    enum Damping { none, rational, zero };

    Damping damping = none;

    /// the method name handed to simple-dftd3's parameter tables
    std::string method;

    /// include the three-body ATM term
    bool atm = false;

    /// evaluate on rank 0 and broadcast; fills the cache below
    void compute(World& world, const Molecule& mol) const;

    /// memoized result -- SCF::solve asks once per iteration for a quantity
    /// that only changes when the optimizer moves the nuclei
    mutable Tensor<double> cached_coords;
    mutable Tensor<double> cached_gradient;
    mutable double cached_energy = 0.0;
    mutable bool cache_valid = false;
    mutable bool citation_printed = false;
};

} // namespace madness

#endif /* SRC_MADNESS_CHEM_DISPERSION_H_ */
