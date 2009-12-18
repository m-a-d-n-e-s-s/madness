#ifndef MADNESS_MOLECULE_H
#define MADNESS_MOLECULE_H

/// \file moldft/molecule.h
/// \brief Declaration of molecule related classes and functions

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <cmath>
#include <misc/misc.h>

struct AtomicData {
    // !!! The order of declaration here must match the order in the initializer !!!

    // Nuclear info from L. Visscher and K.G. Dyall, Dirac-Fock
    // atomic electronic structure calculations using different
    // nuclear charge distributions, Atom. Data Nucl. Data Tabl., 67,
    // (1997), 207.
    //
    // http://dirac.chem.sdu.dk/doc/FiniteNuclei/FiniteNuclei.shtml
    const char* const symbol;
    const char* const symbol_lowercase;
    const unsigned int atomic_number;
    const int isotope_number;
    const double nuclear_radius;     ///< Radius of the nucleus for the finite nucleus models (in atomic units).
    const double nuclear_half_charge_radius; ///< Half charge radius in the Fermi Model (in atomic units).
    const double nuclear_gaussian_exponent; ///< Exponential parameter in the Gaussian Model (in atomic units).

    /// Covalent radii stolen without shame from NWChem
    const double covalent_radius;
};

const AtomicData& get_atomic_data(unsigned int atn);

unsigned int symbol_to_atomic_number(const std::string& symbol);


class Atom {
public:
    double x, y, z, q;          ///< Coordinates and charge in atomic units
    unsigned int atomic_number; ///< Atomic number

    explicit Atom(double x, double y, double z, double q, unsigned int atomic_number)
            : x(x), y(y), z(z), q(q), atomic_number(atomic_number) {}

    Atom(const Atom& a)
            : x(a.x), y(a.y), z(a.z), q(a.q), atomic_number(a.atomic_number) {}

    /// Default construct makes a zero charge ghost atom at origin
    Atom()
            : x(0), y(0), z(0), q(0), atomic_number(0) {}

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & x & y & z & q & atomic_number;
    }
};

std::ostream& operator<<(std::ostream& s, const Atom& atom);

class Molecule {
private:
    // If you add more fields don't forget to serialize them
    std::vector<Atom> atoms;
    std::vector<double> rcut;  // Reciprocal of the smoothing radius
    double eprec;              // Error in energy/atom due to smoothing

    void swapaxes(int ix, int iy);

    template <typename opT>
    bool test_for_op(double xaxis, double yaxis, double zaxis, opT op) const;

    bool test_for_c2(double xaxis, double yaxis, double zaxis) const;

    bool test_for_sigma(double xaxis, double yaxis, double zaxis) const;

    bool test_for_inverse() const;

public:
    /// Makes a molecule with zero atoms
    Molecule() : atoms(), rcut(), eprec(1e-4) {};

    Molecule(const std::string& filename);

    void read_file(const std::string& filename);

    void add_atom(double x, double y, double z,  double q, int atn);

    int natom() const {
        return atoms.size();
    };

    void set_atom_coords(unsigned int i, double x, double y, double z);

    void set_eprec(double value) {
        eprec = value;
    }

    double get_eprec() const {
        return eprec;
    }

    double bounding_cube() const;

    const Atom& get_atom(unsigned int i) const;

    void print() const;

    double inter_atomic_distance(unsigned int i,unsigned int j) const;

    double nuclear_repulsion_energy() const;

    double nuclear_repulsion_derivative(int i, int j) const;

    double smallest_length_scale() const;

    void identify_point_group();

    void center();

    void orient();

    double total_nuclear_charge() const;

    double nuclear_attraction_potential(double x, double y, double z) const;

    double nuclear_attraction_potential_derivative(int atom, int axis, double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & atoms & rcut & eprec;
    }
};


#endif
