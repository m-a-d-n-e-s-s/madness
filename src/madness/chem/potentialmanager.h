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
#include <madness/tensor/tensor.h>
#include <madness/misc/misc.h>
#include <madness/misc/interpolation_1d.h>
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
    MolecularCorePotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coord_3d& r) const {
        return molecule.molecular_core_potential(r[0], r[1], r[2]);
    }

    std::vector<coord_3d> special_points() const {return molecule.get_all_coords_vec();}
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
 * @class GaussianNuclearDensityFunctor
 * @brief Functor for evaluating the Coulomb potential of all nuclei of a molecule; nuclei are represented by primitive spherical (l=0) Gaussians.
 *
 * It computes the Gaussian potential generated by a given molecule at a specified position.
 *
 * @note DOI 10.1006/adnd.1997.0751
 */
class GaussianNuclearDensityFunctor : public FunctionFunctorInterface<double, 3> {
 private:
  const Molecule& molecule;
  int special_level_ = 15;

 public:
  GaussianNuclearDensityFunctor(const madness::Molecule& molecule, int special_level = 15)
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
};
}

#endif
