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

#ifndef MADNESS_CHEM_MOLECULE_H__INCLUDED
#define MADNESS_CHEM_MOLECULE_H__INCLUDED

/// \file moldft/molecule.h
/// \brief Declaration of molecule related classes and functions

#include <chem/corepotential.h>
#include <chem/atomutil.h>
#include <madness/world/vector.h>
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

namespace madness {

class Atom {
public:
    double x, y, z, q;          ///< Coordinates and charge in atomic units
    unsigned int atomic_number; ///< Atomic number
    double mass;                ///< Mass
    bool pseudo_atom;           ///< Indicates if this atom uses a pseudopotential

    explicit Atom(double x, double y, double z, double q, unsigned int atomic_number, bool pseudo_atom)
            : x(x), y(y), z(z), q(q), atomic_number(atomic_number), pseudo_atom(pseudo_atom) {
        mass= get_atomic_data(atomic_number).mass;

        if (mass==-1.0) MADNESS_EXCEPTION("faulty element in Atom",1);

        // unstable elements are indicated by negative masses, the mass
        // is taken from the longest-living element
        if (mass<0.0) mass*=-1.0;

    }

    explicit Atom(double x, double y, double z, double q, unsigned int atomic_number)
            : x(x), y(y), z(z), q(q), atomic_number(atomic_number) {
        mass= get_atomic_data(atomic_number).mass;

        if (mass==-1.0) MADNESS_EXCEPTION("faulty element in Atom",1);

        // unstable elements are indicated by negative masses, the mass
        // is taken from the longest-living element
        if (mass<0.0) mass*=-1.0;

        pseudo_atom = false;

    }

    Atom(const Atom& a) : x(a.x), y(a.y), z(a.z), q(a.q),
            atomic_number(a.atomic_number), mass(a.mass), pseudo_atom(a.pseudo_atom) {}

    /// Default construct makes a zero charge ghost atom at origin
    Atom() : x(0), y(0), z(0), q(0), atomic_number(0), mass(0.0), pseudo_atom(false) {}

   int get_atomic_number() const {return atomic_number;}

    madness::Vector<double,3> get_coords() const {
        return madness::Vector<double,3>{x, y, z};
    }

    /// return the mass in atomic units (electron mass = 1 a.u.)
    double get_mass_in_au() const {return constants::atomic_mass_in_au * mass;}

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & x & y & z & q & atomic_number & mass & pseudo_atom;
    }
};

std::ostream& operator<<(std::ostream& s, const Atom& atom);

class Molecule {
private:
    // If you add more fields don't forget to serialize them
    std::vector<Atom> atoms;
    std::vector<double> rcut;  // Reciprocal of the smoothing radius
    double eprec;              // Error in energy/atom due to smoothing
    enum {atomic, angstrom} units;
    CorePotentialManager core_pot;
    madness::Tensor<double> field;

    void swapaxes(int ix, int iy);

    template <typename opT>
    bool test_for_op(double xaxis, double yaxis, double zaxis, opT op) const;

    bool test_for_c2(double xaxis, double yaxis, double zaxis) const;

    bool test_for_sigma(double xaxis, double yaxis, double zaxis) const;

    bool test_for_inverse() const;

public:

    /// The molecular point group
    /// is automatically assigned in the identify_pointgroup function
    std::string pointgroup_;

    /// Makes a molecule with zero atoms
    Molecule() : atoms(), rcut(), eprec(1e-4), core_pot(), field(3L) {};

    Molecule(const std::string& filename);

    void read_file(const std::string& filename);

    void read_core_file(const std::string& filename);

    std::string guess_file() const { return core_pot.guess_file(); };

    unsigned int n_core_orb_all() const ;

    unsigned int n_core_orb(unsigned int atn) const {
        if (core_pot.is_defined(atn))
            return core_pot.n_core_orb_base(atn);
        else
            return 0;
    };

    unsigned int get_core_l(unsigned int atn, unsigned int c) const {
        return core_pot.get_core_l(atn, c);
    }

    double get_core_bc(unsigned int atn, unsigned int c) const {
        return core_pot.get_core_bc(atn, c);
    }

    double core_eval(int atom, unsigned int core, int m, double x, double y, double z) const;

    double core_derivative(int atom, int axis, unsigned int core, int m, double x, double y, double z) const;

    bool is_potential_defined(unsigned int atn) const { return core_pot.is_defined(atn); };

    bool is_potential_defined_atom(int i) const { return core_pot.is_defined(atoms[i].atomic_number); };

    void add_atom(double x, double y, double z,  double q, int atn);

    void add_atom(double x, double y, double z,  double q, int atn, bool psat);

    int natom() const {
        return atoms.size();
    };

    void set_atom_charge(unsigned int i, double zeff);

    unsigned int get_atom_charge(unsigned int i) const;

    unsigned int get_atom_number(unsigned int i);

    void set_pseudo_atom(unsigned int i, bool psat);

    bool get_pseudo_atom(unsigned int i);

    void set_atom_coords(unsigned int i, double x, double y, double z);

    madness::Tensor<double> get_all_coords() const;

    std::vector< madness::Vector<double,3> > get_all_coords_vec() const;

    std::vector<double> atomic_radii;

    void set_all_coords(const madness::Tensor<double>& newcoords);

    void set_eprec(double value);

    void set_rcut(double value);

    std::vector<double> get_rcut() const {return rcut;}

    void set_core_eprec(double value) {
        core_pot.set_eprec(value);
    }

    void set_core_rcut(double value) {
        core_pot.set_rcut(value);
    }

    double get_eprec() const {
        return eprec;
    }

    double bounding_cube() const;

    const Atom& get_atom(unsigned int i) const;

    const std::vector<Atom> & get_atoms()const{return atoms;}

    void print() const;

    double inter_atomic_distance(unsigned int i,unsigned int j) const;

    double nuclear_repulsion_energy() const;

    double nuclear_repulsion_derivative(int i, int j) const;

    /// compute the nuclear-nuclear contribution to the second derivatives

    /// @param[in]  iatom   the i-th atom (row of the hessian)
    /// @param[in]  jatom   the j-th atom (column of the hessian)
    /// @param[in]  iaxis   the xyz axis of the i-th atom
    /// @param[in]  jaxis   the xyz axis of the j-th atom
    /// return the (3*iatom + iaxis, 3*jatom + jaxis) matix element of the hessian
    double nuclear_repulsion_second_derivative(int iatom, int jatom,
            int iaxis, int jaxis) const;

    /// return the hessian matrix of the second derivatives d^2/dxdy V

    /// no factor 0.5 included
    Tensor<double> nuclear_repulsion_hessian() const;

    /// compute the dipole moment of the nuclei

    ///  @param[in] axis the axis (x, y, z)
    double nuclear_dipole(int axis) const;

    /// compute the derivative of the nuclear dipole wrt a nuclear displacement

    /// @param[in]  atom    the atom which will be displaced
    /// @param[in]  axis    the axis where the atom will be displaced
    /// @return     a vector which all 3 components of the dipole derivative
    Tensor<double> nuclear_dipole_derivative(const int atom, const int axis) const;

    double nuclear_charge_density(double x, double y, double z) const;

    double mol_nuclear_charge_density(double x, double y, double z) const;

    double smallest_length_scale() const;

    void identify_point_group();

    /// Moves the center of nuclear charge to the origin
    void center();

    /// rotates the molecule and the external field

    /// @param[in]  D   the rotation matrix
    void rotate(const Tensor<double>& D);

    /// translate the molecule
    void translate(const Tensor<double>& translation);

    Tensor<double> center_of_mass() const;

    /// compute the mass-weighting matrix for the hessian

    /// use as
    /// mass_weighted_hessian=inner(massweights,inner(hessian,massweights));
    Tensor<double> massweights() const {

        Tensor<double> M(3*natom(),3*natom());
        for (int i=0; i<natom(); i++) {
            const double sqrtmass=1.0/sqrt(get_atom(i).get_mass_in_au());
            M(3*i  ,3*i  )=sqrtmass;
            M(3*i+1,3*i+1)=sqrtmass;
            M(3*i+2,3*i+2)=sqrtmass;
        }
        return M;
    }



    Tensor<double> moment_of_inertia() const;

    void orient();

    double total_nuclear_charge() const;

    /// nuclear attraction potential for the whole molecule
    double nuclear_attraction_potential(double x, double y, double z) const;

    /// nuclear attraction potential for a specific atom in the molecule
    double atomic_attraction_potential(int iatom, double x, double y, double z) const;

    double molecular_core_potential(double x, double y, double z) const;

    double core_potential_derivative(int atom, int axis, double x, double y, double z) const;

    double nuclear_attraction_potential_derivative(int atom, int axis, double x, double y, double z) const;

    double nuclear_attraction_potential_second_derivative(int atom, int iaxis,
            int jaxis, double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & atoms & rcut & eprec & core_pot;
    }
};
}

#endif
