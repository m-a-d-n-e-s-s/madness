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

  $Id$
*/
#ifndef MADNESS_CHEM_POTENTIALMANAGER_H__INCLUDED
#define MADNESS_CHEM_POTENTIALMANAGER_H__INCLUDED

/**
 * @file potentialmanager.h
 * @brief Declaration of molecule-related classes and functions.
 */

#include <madness/chem/corepotential.h>
#include <madness/constants.h>
#include <madness/chem/atomutil.h>
#include <madness/chem/molecule.h>
#include <madness/mra/kernelrange.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <cmath>
#include <array>
#include <unordered_map>
#include <madness/tensor/tensor.h>
#include <madness/misc/misc.h>
#include <madness/mra/interpolation_1d.h>
#include <madness/mra/mra.h>

namespace madness {
/**
 * @class MolecularPotentialFunctor
 * @brief Functor for evaluating the nuclear attraction potential of a molecule at a given point.
 *
 * This class implements the FunctionFunctorInterface for 3D coordinates and provides
 * an interface to evaluate the nuclear attraction potential (smoothed Coulomb potential)
 * of a molecule at a specified point in space. It also provides access to special points,
 * the coordinates of all nuclei in the molecule.
 *
 */
class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
public:
    using FunctionFunctorInterface<double,3>::operator();

    MolecularPotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coord_3d& r) const {
        return molecule.nuclear_attraction_potential(r[0], r[1], r[2]);
    }

    std::vector<coord_3d> special_points() const {return molecule.get_all_coords_vec();}
};

/**
 * @class MolecularCorePotentialFunctor
 * @brief Functor for evaluating the molecular core potential at a given point in space.
 *
 * This class implements the FunctionFunctorInterface to provide a callable object
 * that computes the molecular core potential for a given 3D coordinate using the
 * associated Molecule instance.
 *
 */
class MolecularCorePotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
public:
    using FunctionFunctorInterface<double,3>::operator();

    MolecularCorePotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coord_3d& r) const {
        return molecule.molecular_core_potential(r[0], r[1], r[2]);
    }

    std::vector<coord_3d> special_points() const {return molecule.get_all_coords_vec();}
};

/// The Gaussian-smoothed nuclear charges \f$ s_i = Z_i (a/\pi)^{3/2} e^{-a r_i^2} \f$ of the
/// smoothed-density Coulomb energy: centres and effective charges
/// (Molecule::effective_nuclear_charge) of the atoms that carry one.
struct SmoothedNuclearCharges {
    double a;
    std::vector<coord_3d> centres;
    std::vector<double> charges;
    std::vector<size_t> atoms;    ///< the index of each centre in the molecule

    SmoothedNuclearCharges(const Molecule& molecule, double a) : a(a) {
        MADNESS_CHECK_THROW(a > 0.0, "the exponent of the smoothed nuclear charge density must be positive");
        for (size_t i = 0; i < molecule.natom(); ++i) {
            const double z = molecule.effective_nuclear_charge(i);
            if (z == 0.0) continue;
            const Atom& atom = molecule.get_atom(i);
            centres.push_back(coord_3d{atom.x, atom.y, atom.z});
            charges.push_back(z);
            atoms.push_back(i);
        }
    }
};

/// Centres sorted into cubic cells of a given edge. Every centre within one edge length of a
/// point lies in the point's own cell or in one of the 26 adjacent cells, so summing short-ranged
/// atomic terms costs O(1) per point, and projecting their sum costs O(natom).
class NuclearCellList {
private:
    using cell_index = std::array<long,3>;
    struct cell_hash {
        std::size_t operator()(const cell_index& c) const {
            const std::hash<long> h;
            return h(c[0]) ^ (h(c[1]) << 21) ^ (h(c[2]) << 42);
        }
    };

    double edge;
    std::unordered_map<cell_index, std::vector<size_t>, cell_hash> cells;

    cell_index cell_of(const coord_3d& x) const {
        return {long(std::floor(x[0] / edge)), long(std::floor(x[1] / edge)), long(std::floor(x[2] / edge))};
    }

public:
    NuclearCellList(const std::vector<coord_3d>& centres, double edge) : edge(edge) {
        for (size_t i = 0; i < centres.size(); ++i) cells[cell_of(centres[i])].push_back(i);
    }

    /// call \p f(i) for every centre i in the cell of \p x and in the 26 adjacent cells
    template <typename funcT>
    void for_each_near(const coord_3d& x, funcT&& f) const {
        const cell_index c = cell_of(x);
        for (long i = -1; i <= 1; ++i) {
            for (long j = -1; j <= 1; ++j) {
                for (long k = -1; k <= 1; ++k) {
                    const auto it = cells.find(cell_index{c[0] + i, c[1] + j, c[2] + k});
                    if (it == cells.end()) continue;
                    for (const size_t icentre : it->second) f(icentre);
                }
            }
        }
    }
};

/// Potential \f$ u(r;s) = \sum_i Z_i\,\mathrm{erf}(\sqrt{a}\,r_i)/r_i \f$ of the smoothed nuclear
/// charge density. Like the nuclear potential, every point sums over all atoms.
class SmoothedMolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    SmoothedNuclearCharges s;
    double sqrta;
    double peak;    ///< the r -> 0 limit of erf(sqrt(a) r)/r, 2 sqrt(a/pi)
public:
    using FunctionFunctorInterface<double,3>::operator();

    SmoothedMolecularPotentialFunctor(const Molecule& molecule, double a)
        : s(molecule, a), sqrta(std::sqrt(a)), peak(2.0 * std::sqrt(a / madness::constants::pi)) {}

    double operator()(const coord_3d& x) const override {
        double sum = 0.0;
        for (size_t i = 0; i < s.centres.size(); ++i) {
            const double dx = x[0] - s.centres[i][0];
            const double dy = x[1] - s.centres[i][1];
            const double dz = x[2] - s.centres[i][2];
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            sum += s.charges[i] * (r > 0.0 ? std::erf(sqrta * r) / r : peak);
        }
        return sum;
    }

    std::vector<coord_3d> special_points() const override {return s.centres;}
};

/// The smoothed nuclear charge density \f$ s(r) = \sum_i Z_i (a/\pi)^{3/2} e^{-a r_i^2} \f$.
/// Each Gaussian is cut off where it has decayed to 1e-18 of its peak. A cell list with that
/// edge length restricts a point to the atoms of its own and the 26 adjacent cells, so that
/// projecting s costs O(natom), and Gaussians that overlap in short bonds are summed.
class SmoothedNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    SmoothedNuclearCharges s;
    double norm;    ///< (a/pi)^{3/2}
    double rcut;    ///< exp(-a rcut^2) = 1e-18
    NuclearCellList near;

public:
    using FunctionFunctorInterface<double,3>::operator();

    SmoothedNuclearChargeDensityFunctor(const Molecule& molecule, double a)
        : s(molecule, a), norm(std::pow(a / madness::constants::pi, 1.5)),
          rcut(std::sqrt(18.0 * std::log(10.0) / a)), near(s.centres, rcut) {}

    double operator()(const coord_3d& x) const override {
        const double rcutsq = rcut * rcut;
        double sum = 0.0;
        near.for_each_near(x, [&](const size_t iatom) {
            const double dx = x[0] - s.centres[iatom][0];
            const double dy = x[1] - s.centres[iatom][1];
            const double dz = x[2] - s.centres[iatom][2];
            const double rsq = dx * dx + dy * dy + dz * dz;
            if (rsq < rcutsq) sum += s.charges[iatom] * norm * std::exp(-s.a * rsq);
        });
        return sum;
    }

    std::vector<coord_3d> special_points() const override {return s.centres;}
};

/// The local remainder c = Vnuc + u(r;s) of the nuclear potential, so that Vnuc = -G*s + c:
/// \f[ c(r) = \sum_i Z_i \left[ \mathrm{erf}(\sqrt{a}\,r_i)/r_i - u(r_i/c_i)/c_i \right] + F\cdot r, \f]
/// with u the smoothed 1/r and c_i the smoothing radius of Molecule::nuclear_attraction_potential,
/// and F the external field. Beyond 7 c_i, where u(x) = 1/x exactly, a term is
/// -Z_i erfc(sqrt(a) r_i)/r_i. It is cut off where sqrt(a) r_i > sqrt(18 ln 10), below 1e-19 Z_i,
/// and a cell list with that edge makes projecting c O(natom).
///
/// All-electron atoms only. s carries the effective nuclear charges; an atom whose charge in
/// Vnuc differs from it (pseudo-atom, core potential, ghost point charge) would leave a
/// long-ranged term in c, so the constructor refuses it.
class NuclearPotentialRemainderFunctor : public FunctionFunctorInterface<double,3> {
private:
    SmoothedNuclearCharges s;
    std::vector<double> rsmooth;    ///< 1/c_i of each centre, Molecule::get_rcut()
    double sqrta;
    double peak;    ///< the r -> 0 limit of erf(sqrt(a) r)/r, 2 sqrt(a/pi)
    double rcut;    ///< cutoff radius of every term, and the edge of the cell list
    NuclearCellList near;
    coord_3d field;

    static std::vector<double> smoothing_radii(const Molecule& molecule, const SmoothedNuclearCharges& s) {
        MADNESS_CHECK_THROW(molecule.all_electron_nuclei(),
                            "the nuclear potential by convolution supports all-electron atoms only "
                            "(no pseudopotentials, core potentials or ghost point charges): use "
                            "nuclear_potential projection");
        const std::vector<double> rcut_all = molecule.get_rcut();
        std::vector<double> r;
        r.reserve(s.atoms.size());
        for (const size_t i : s.atoms) r.push_back(rcut_all[i]);
        return r;
    }

    static double cutoff(double a, const std::vector<double>& rsmooth) {
        double r = std::sqrt(18.0 * std::log(10.0) / a);
        for (const double rs : rsmooth) r = std::max(r, 7.0 / rs);
        return r;
    }

public:
    using FunctionFunctorInterface<double,3>::operator();

    NuclearPotentialRemainderFunctor(const Molecule& molecule, double a)
        : s(molecule, a), rsmooth(smoothing_radii(molecule, s)), sqrta(std::sqrt(a)),
          peak(2.0 * std::sqrt(a / madness::constants::pi)), rcut(cutoff(a, rsmooth)),
          near(s.centres, rcut),
          field{molecule.get_field()[0], molecule.get_field()[1], molecule.get_field()[2]} {}

    double operator()(const coord_3d& x) const override {
        const double rcutsq = rcut * rcut;
        double sum = field[0] * x[0] + field[1] * x[1] + field[2] * x[2];
        near.for_each_near(x, [&](const size_t i) {
            const double dx = x[0] - s.centres[i][0];
            const double dy = x[1] - s.centres[i][1];
            const double dz = x[2] - s.centres[i][2];
            const double rsq = dx * dx + dy * dy + dz * dz;
            if (rsq >= rcutsq) return;
            const double r = std::sqrt(rsq);
            if (r * rsmooth[i] > 7.0) {
                sum -= s.charges[i] * std::erfc(sqrta * r) / r;
            } else {
                const double u_s = (r > 0.0) ? std::erf(sqrta * r) / r : peak;
                sum += s.charges[i] * (u_s - smoothed_potential(r * rsmooth[i]) * rsmooth[i]);
            }
        });
        return sum;
    }

    std::vector<coord_3d> special_points() const override {return s.centres;}
};

/**
 * @class CoreOrbitalFunctor
 * @brief Functor for evaluating a core orbital of a specific atom in a molecule.
 *
 * This class implements the FunctionFunctorInterface for evaluating the value of a core orbital
 * at a given 3D coordinate. It holds references to a molecule, the atom index, the core orbital index,
 * and the magnetic quantum number m.
 *
 */
class CoreOrbitalFunctor : public FunctionFunctorInterface<double,3> {
    const Molecule molecule;
    const int atom;
    const unsigned int core;
    const int m;
public:
    using FunctionFunctorInterface<double,3>::operator();

    CoreOrbitalFunctor(Molecule& molecule, int atom, unsigned int core, int m)
        : molecule(molecule), atom(atom), core(core), m(m) {};
    double operator()(const coord_3d& r) const {
        return molecule.core_eval(atom, core, m, r[0], r[1], r[2]);
    };
};

/**
 * @class CoreOrbitalDerivativeFunctor
 * @brief Functor for evaluating the derivative of a core orbital for a given atom in a molecule.
 *
 * It encapsulates the logic to compute the derivative of a specified core orbital with respect to a given axis
 * for a particular atom in a molecule. It holds references to a molecule, the atom index, the axis index,
 * the core orbital index, and the magnetic quantum number m.
 *
 */
class CoreOrbitalDerivativeFunctor : public FunctionFunctorInterface<double,3> {
    const Molecule molecule;
    const int atom, axis;
    const unsigned int core;
    const int m;
public:
    using FunctionFunctorInterface<double,3>::operator();

    CoreOrbitalDerivativeFunctor(Molecule& molecule, int atom, int axis, unsigned int core, int m)
        : molecule(molecule), atom(atom), axis(axis), core(core), m(m) {};
    double operator()(const coord_3d& r) const {
        return molecule.core_derivative(atom, axis, core, m, r[0], r[1], r[2]);
    };
};

/**
 * @class NuclearDensityFunctor
 * @brief Default functor for evaluating nuclear density at a given point in space.
 *
 * This class implements a functor that computes the nuclear density for a given molecule,
 * supporting both open and periodic boundary conditions. It can be used to evaluate the
 * nuclear density at any point in 3D space, and provides special points and refinement
 * level information for adaptive algorithms.
 *
 */
class NuclearDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
  const Molecule& atoms;
  BoundaryConditions<3> bc_;
  Tensor<double> cell;
  std::vector<coord_3d> special_points_;
  int maxR;
  int special_level_ = 15;
  double rscale = 1.0;
public:
  using FunctionFunctorInterface<double,3>::operator();

  /**
  * @brief Constructs a NuclearDensityFunctor for evaluating nuclear densities.
  *
  * This constructor can handle both open and periodic boundary conditions.
  *
  * @param atoms Reference to the molecule containing the atoms.
  * @param bc Boundary conditions for the simulation (default: open boundaries).
  * @param cell Simulation cell tensor (unit cell, if periodic; default: identity).
  * @param special_level The initial refinement level for special points (default: 15).
  * @param rscale Scaling factor for the nuclear radius. Setting rscale > 1 increases the effective size of a nucleus by this factor (i.e., rcut is divided by rscale).
  */
  NuclearDensityFunctor(const Molecule& atoms,
                        const BoundaryConditions<3>& bc = FunctionDefaults<3>::get_bc(),
                        const Tensor<double>& cell = FunctionDefaults<3>::get_cell(),
                        int special_level = 15,
                        double rscale = 1.0);

  double operator()(const coord_3d& x) const final;

  std::vector<coord_3d> special_points() const final;

  Level special_level() const final;

  NuclearDensityFunctor& set_rscale(double rscale);

};

/**
 * @class GaussianNuclearDensityPotentialFunctor
 * @brief Functor for evaluating the Coulomb potential of all nuclei of a molecule; nuclei are represented by primitive spherical (l=0) Gaussians.
 *
 * It computes the Gaussian potential generated by a given molecule at a specified position.
 *
 * @note DOI 10.1006/adnd.1997.0751
 */
class GaussianNuclearDensityPotentialFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const Molecule& molecule;
  int special_level_ = 15;

 public:
  using FunctionFunctorInterface<double,3>::operator();

  GaussianNuclearDensityPotentialFunctor(const madness::Molecule& molecule, int special_level = 15)
      : molecule(molecule), special_level_(special_level) {}

  double operator()(const madness::coord_3d& R) const final;

  madness::Level special_level() const final { return special_level_; };

  std::vector<coord_3d> special_points() const { return molecule.get_all_coords_vec(); }
};

/**
 * @class FermiNuclearDensityFunctor
 * @brief Functor representing the Fermi nuclear density distribution for a given atom.
 *
 * This class implements a functor that evaluates the two-parameter charge distribution for
 * the nuclear density at a given 3D coordinate. The density is significant only in a small
 * region around the atomic center.
 *
 * @note DOI 10.1006/adnd.1997.0751
 */
class FermiNuclearDensityFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const Atom& atom;
  int special_level_ = 18;

 public:
  using FunctionFunctorInterface<double,3>::operator();

  FermiNuclearDensityFunctor(const Atom& atom, int special_level = 18)
      : atom(atom), special_level_(special_level) {}

  double operator()(const madness::coord_3d& R) const final;

  madness::Level special_level() const final { return special_level_; };

  std::vector<coord_3d> special_points() const final { return {atom.get_coords()}; }
};

/**
 * @class WignerSeitzPotentialFunctor
 * @brief Functor for evaluating the Wigner-Seitz potential in a simulation cell.
 *
 * This class implements a functor that evaluates the electrostatic potential
 * in a simulation cell due to a set of point charges, with optional periodic
 * boundary conditions and configurable lattice summation range.
 *
 * The potential is computed by summing contributions from point charges
 * in the simulation cell and their periodic images, as determined by the
 * specified boundary conditions and kernel range.
 *
 * @note The lattice summation range can be overridden by the user, or
 *       determined automatically based on the boundary conditions and kernel range.
 */
class WignerSeitzPotentialFunctor : public FunctionFunctorInterface<double,3> {
public:
  using FunctionFunctorInterface<double,3>::operator();

  /**
   * @brief Constructs a WignerSeitzPotentialFunctor evaluating the potential
   *        in a simulation cell due to point charges, optionally with periodic
   *        boundary conditions and a specified lattice summation range.
   *
   * @tparam Int Integer type for the lattice summation range.
   * @param atoms List of point charges in the simulation cell.
   * @param c The simulation cell dimensions.
   * @param b The boundary conditions.
   * @param r The kernel range along each Cartesian direction.
   * @param lattice_sum_range Overrides the default lattice summation range along each axis.
   *        By default, this is determined by the number of cells in each direction with
   *        nonzero contributions to the simulation cell.
   */
  template <typename Int>
  WignerSeitzPotentialFunctor(const Molecule &atoms, Tensor<double> c,
                              BoundaryConditions<3> b, std::array<KernelRange, 3> r,
                              std::array<Int, 3> lattice_sum_range)
      : atoms(atoms), cell(std::move(c)), bc(std::move(b)), range(std::move(r)),
        cell_width{cell(0, 1) - cell(0, 0), cell(1, 1) - cell(1, 0),
                   cell(2, 1) - cell(2, 0)},
        rcell_width{1. / cell_width[0], 1. / cell_width[1],
                    1. / cell_width[2]}, lattice_sum_range{lattice_sum_range[0], lattice_sum_range[1], lattice_sum_range[2]} {
    for (int d = 0; d != 3; ++d)
      MADNESS_ASSERT(lattice_sum_range[d] >= 0);
  }

  /**
   * @brief Constructs a WignerSeitzPotentialFunctor with default lattice sum range.
   *
   * This constructor initializes the WignerSeitzPotentialFunctor using the provided molecule,
   * coefficients tensor, boundary conditions, and kernel ranges. It automatically computes
   * the default lattice sum range using the given boundary conditions and kernel ranges.
   *
   * @param atoms The molecule containing the atomic positions and properties.
   * @param c The tensor of coefficients for the potential calculation.
   * @param b The boundary conditions for the simulation cell.
   * @param r The kernel ranges for each spatial dimension.
   *
   * @note This constructor delegates to the main constructor, passing the default lattice sum range.
   * In this case, the move is a cast, and calling make_default_lattice_sum_range like this is OK.
   */
  WignerSeitzPotentialFunctor(const Molecule &atoms, Tensor<double> c,
                              BoundaryConditions<3> b, std::array<KernelRange, 3> r) :
  WignerSeitzPotentialFunctor(atoms, std::move(c), std::move(b), std::move(r), make_default_lattice_sum_range(b,r)) {}

  double operator()(const coord_3d &x) const final;

  std::vector<coord_3d> special_points() const final;

  static std::array<std::int64_t, 3> make_default_lattice_sum_range(const BoundaryConditions<3>& bc, const std::array<KernelRange, 3>& range) {
    std::array<std::int64_t, 3> result;
    for (int d = 0; d != 3; ++d) {
      result[d] = bc.is_periodic()[d] ? (range[d].iextent_x2() + 1) / 2 : 0;
    }
    return result;
  }

private:
  const Molecule &atoms;
  const Tensor<double> cell;
  const BoundaryConditions<3> bc;
  const std::array<KernelRange, 3> range;
  const std::array<double, 3> cell_width;
  const std::array<double, 3> rcell_width;
  const std::array<std::int64_t, 3> lattice_sum_range;
};

/**
 * @class SAPFunctor
 * @brief Functor for evaluating a smoothed atomic potential, supporting open and periodic boundary conditions.
 *
 * This class implements the FunctionFunctorInterface for a 3D double-valued function,
 * representing a smoothed interpolated atomic potential centered on a given atom. It supports
 * both open and periodic boundary conditions, and allows customization of the smoothing
 * parameter, simulation cell, and initial refinement level.
 *
 */
class SAPFunctor : public FunctionFunctorInterface<double,3> {
 private:
  const Atom& atom;
  double smoothing_param;
  BoundaryConditions<3> bc_;
  Tensor<double> cell;
  Level special_level_;
 public:
  using FunctionFunctorInterface<double,3>::operator();

  /**
   * @brief Constructs a SAPFunctor for evaluating a smoothed 1/r potential.
   *
   * This constructor initializes the SAPFunctor with a given atom, smoothing parameter,
   * boundary conditions, simulation cell, and an initial refinement level. It supports
   * both open and periodic boundary conditions.
   *
   * @param atom The atom for which the potential is evaluated.
   * @param smoothing_param Controls the smoothness of the 1/r potential.
   * @param bc Boundary conditions for the simulation (default: open or as specified by FunctionDefaults).
   * @param cell The simulation cell tensor (default: as specified by FunctionDefaults).
   * @param special_level The initial refinement level (default: 15).
   */
  SAPFunctor(const Atom& atom,
             double smoothing_param,
             const BoundaryConditions<3>& bc = FunctionDefaults<3>::get_bc(),
             const Tensor<double>& cell = FunctionDefaults<3>::get_cell(),
             int special_level = 15);

  double operator()(const coord_3d& x) const final;

  Level special_level() const final;

  std::vector<coord_3d> special_points() const final;
};

/**
 * @class PotentialManager
 * @brief Manages molecular potentials and core projections for quantum chemistry calculations.
 *
 * This class encapsulates the management of nuclear and core potentials for a given molecule,
 * including the construction of nuclear potentials, application of nonlocal core projectors,
 * and calculation of core projector derivatives. It provides interfaces to access the molecule,
 * core type, and nuclear potential, as well as to perform core projections and apply nonlocal potentials.
 *
 */
class PotentialManager {
private:
Molecule mol;
real_function_3d vnuc;
real_function_3d unuc;
real_function_3d sdens;
real_function_3d cnuc;          ///< c = Vnuc + u(r;s), the local remainder of the nuclear potential
double ss_coulomb = 0.0;        ///< <s; G~ s> from make_nuclear_potential_by_convolution
std::string core_type_;

public:
    PotentialManager(const Molecule& molecule, const std::string& core_type)
     : mol(molecule), core_type_(core_type) {}

    const Molecule& molecule() const {
      return this->mol;
    }

    const std::string& core_type() const {
      return this->core_type_;
    }

    const real_function_3d& vnuclear() {
        return vnuc;
    }

    const real_function_3d& unuclear() {
        return unuc;
    }

    const real_function_3d& smooth_density() {
        return sdens;
    }

    /// c = Vnuc + u(r;s), the local remainder of the nuclear potential (make_nuclear_potential_remainder)
    const real_function_3d& nuclear_potential_remainder() {
        return cnuc;
    }

    /// <s; G~ s>, the convolved part of the consistent nuclear repulsion
    /// E_nn~ = 1/2 <s; G~ s> - sqrt(a/2pi) sum_i Z_i^2 + E_ov (make_nuclear_potential_by_convolution)
    double smoothed_density_coulomb() const {
        return ss_coulomb;
    }

    /**
     * @brief Projects the input wavefunctions onto the atomic core orbitals.
     *
     * This function computes the projection of the given set of wavefunctions (`psi`)
     * onto the core orbitals of each atom in the molecule. The projection is performed
     * for each atom and each of its core orbitals, accumulating the result in the
     * returned vector of functions. Optionally, the projection can include the core
     * boundary condition factor (`Bc`).
     *
     * @param world The MADNESS World object for parallel execution and data management.
     * @param psi The input vector of real 3D functions (wavefunctions) to be projected.
     * @param include_Bc If true, includes the core boundary condition factor in the projection (default: true).
     * @return A vector of real 3D functions representing the projection of `psi` onto the core orbitals.
     *
     */
    vector_real_function_3d core_projection(World & world, const vector_real_function_3d& psi, const bool include_Bc = true)
    {
        int npsi = psi.size();
        if (npsi == 0) return psi;
        int natom = mol.natom();
        vector_real_function_3d proj = zero_functions_compressed<double,3>(world, npsi);
        real_tensor overlap_sum(static_cast<long>(npsi));

        for (int i=0; i<natom; ++i) {
            Atom at = mol.get_atom(i);
            unsigned int atn = at.atomic_number;
            unsigned int nshell = mol.n_core_orb(atn);
            if (nshell == 0) continue;
            for (unsigned int c=0; c<nshell; ++c) {
                unsigned int l = mol.get_core_l(atn, c);
                int max_m = (l+1)*(l+2)/2;
                nshell -= max_m - 1;
                for (int m=0; m<max_m; ++m) {
                    real_function_3d core = real_factory_3d(world).functor(real_functor_3d(new CoreOrbitalFunctor(mol, i, c, m)));
                    real_tensor overlap = inner(world, core, psi);
                    overlap_sum += overlap;
                    for (int j=0; j<npsi; ++j) {
                        if (include_Bc) overlap[j] *= mol.get_core_bc(atn, c);
                        proj[j] += core.scale(overlap[j]);
                    }
                }
            }
            world.gop.fence();
        }
        if (world.rank() == 0) print("sum_k <core_k|psi_i>:", overlap_sum);
        return proj;
    }

    /**
     * @brief Computes the derivative of the core projector operator with respect to a given axis for a specified atom.
     *
     * This function projects the core orbitals and their derivatives onto the molecular orbitals,
     * then evaluates the sum:
     * \f[
     * \sum_i \mathrm{occ}_i \langle \psi_i | \left( \sum_c B_c \frac{d}{dx} | \mathrm{core} \rangle \langle \mathrm{core} | \right) | \psi_i \rangle
     * \f]
     * where \f$ \psi_i \f$ are molecular orbitals, \f$ \mathrm{occ}_i \f$ are their occupations,
     * and the sum over \f$ c \f$ runs over the core orbitals of the specified atom.
     *
     * @param world The MADNESS World object for parallel computation.
     * @param mo The vector of molecular orbitals as real-valued 3D functions.
     * @param occ The occupation numbers for each molecular orbital.
     * @param atom The index of the atom for which the core projector derivative is computed.
     * @param axis The spatial axis (0=x, 1=y, 2=z) along which the derivative is taken.
     * @return The computed derivative value as a double.
     */
    double core_projector_derivative(World & world, const vector_real_function_3d& mo, const real_tensor& occ, int atom, int axis)
    {
        vector_real_function_3d cores, dcores;
        std::vector<double> bc;
        unsigned int atn = mol.get_atom(atom).atomic_number;
        unsigned int ncore = mol.n_core_orb(atn);

        // projecting core & d/dx core
        for (unsigned int c=0; c<ncore; ++c) {
            unsigned int l = mol.get_core_l(atn, c);
            int max_m = (l+1)*(l+2)/2;
            for (int m=0; m<max_m; ++m) {
                real_functor_3d func = real_functor_3d(new CoreOrbitalFunctor(mol, atom, c, m));
                cores.push_back(real_function_3d(real_factory_3d(world).functor(func).truncate_on_project()));
                func = real_functor_3d(new CoreOrbitalDerivativeFunctor(mol, atom, axis, c, m));
                dcores.push_back(real_function_3d(real_factory_3d(world).functor(func).truncate_on_project()));
                bc.push_back(mol.get_core_bc(atn, c));
            }
        }

        // calc \sum_i occ_i <psi_i|(\sum_c Bc d/dx |core><core|)|psi_i>
        double r = 0.0;
        for (unsigned int c=0; c<cores.size(); ++c) {
            double rcore= 0.0;
            real_tensor rcores = inner(world, cores[c], mo);
            real_tensor rdcores = inner(world, dcores[c], mo);
            for (unsigned int i=0; i<mo.size(); ++i) {
                rcore += rdcores[i] * rcores[i] * occ[i];
            }
            r += 2.0 * bc[c] * rcore;
        }

        return r;
    }

    void apply_nonlocal_potential(World& world, const vector_real_function_3d& amo, vector_real_function_3d Vpsi) {
        if (core_type_.substr(0,3) == "mcp") {
         //   START_TIMER(world);
            gaxpy(world, 1.0, Vpsi, 1.0, core_projection(world, amo));
         //   END_TIMER(world, "MCP Core Projector");
        }
    }

    void make_nuclear_potential(World& world) {
        double safety = 0.1;
        double vtol = FunctionDefaults<3>::get_thresh() * safety;
        vnuc = real_factory_3d(world).functor(real_functor_3d(new MolecularPotentialFunctor(mol))).thresh(vtol).truncate_on_project();
        vnuc.set_thresh(FunctionDefaults<3>::get_thresh());
        vnuc.reconstruct();
        //     "" is  legacy core_type value for all-electron (also be used by CorePotentialManager)
        // "none" is current core_type value for all-electron
        if (core_type_ != "" && core_type_ != "none") {
            real_function_3d c_pot = real_factory_3d(world).functor(real_functor_3d(new MolecularCorePotentialFunctor(mol))).thresh(vtol).initial_level(4);
            c_pot.set_thresh(FunctionDefaults<3>::get_thresh());
            c_pot.reconstruct();
            vnuc += c_pot;
            vnuc.truncate();
        }
    }

    /// project u(r;s), the potential of the smoothed nuclear charge density with exponent \p a
    void make_smoothed_nuclear_potential(World& world, double a) {
        double safety = 0.1;
        double utol = FunctionDefaults<3>::get_thresh() * safety;
        unuc = real_factory_3d(world).functor(real_functor_3d(new SmoothedMolecularPotentialFunctor(mol, a))).thresh(utol).truncate_on_project();
        unuc.set_thresh(FunctionDefaults<3>::get_thresh());
        unuc.reconstruct();
    }

    /// project the smoothed nuclear charge density s with Gaussian exponent \p a
    void make_smoothed_nuclear_charge_density(World& world, double a) {
        double safety = 0.1;
        double utol = FunctionDefaults<3>::get_thresh() * safety;
        sdens = real_factory_3d(world).functor(real_functor_3d(new SmoothedNuclearChargeDensityFunctor(mol, a))).thresh(utol).truncate_on_project();
        sdens.set_thresh(FunctionDefaults<3>::get_thresh());
        sdens.reconstruct();
        //if (world.rank() == 0) print("smooth density trace",sdens.trace());
    }

    /// project c = Vnuc + u(r;s), the local remainder of the nuclear potential, for the smoothed
    /// density with exponent \p a: O(natom), all-electron atoms only (NuclearPotentialRemainderFunctor)
    void make_nuclear_potential_remainder(World& world, double a) {
        double safety = 0.1;
        double vtol = FunctionDefaults<3>::get_thresh() * safety;
        cnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearPotentialRemainderFunctor(mol, a))).thresh(vtol).truncate_on_project();
        cnuc.set_thresh(FunctionDefaults<3>::get_thresh());
        cnuc.reconstruct();
    }

    /// The nuclear potential by convolution of the smoothed nuclear charges, Vnuc~ = c - G~ s, with
    /// \p G the Coulomb operator the caller's electronic Coulomb potential uses. Replaces the
    /// O(natom^2) projection of make_nuclear_potential by an O(natom) one and one convolution.
    /// Needs s and c (make_smoothed_nuclear_charge_density, make_nuclear_potential_remainder).
    ///
    /// Vnuc~ carries the operator error of G~ s. Only the consistent nuclear repulsion
    /// E_nn~ = 1/2 <s; G~ s> - sqrt(a/2pi) sum_i Z_i^2 + E_ov (smoothed_density_coulomb()) cancels
    /// it, together with the error of the electron repulsion; the analytic E_nn does not.
    void make_nuclear_potential_by_convolution(const real_convolution_3d& G) {
        MADNESS_CHECK_THROW(sdens.is_initialized() && cnuc.is_initialized(),
                            "the nuclear potential by convolution needs the smoothed density and the remainder c");
        real_function_3d Gs = apply(G, sdens);
        ss_coulomb = inner(sdens, Gs);
        vnuc = cnuc - Gs;
        vnuc.truncate(FunctionDefaults<3>::get_thresh() * 0.1);
        vnuc.reconstruct();
    }
};
}

#endif
