#ifndef MENTITY_H_

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

    Atom(double x, double y, double z, double q, unsigned int atomic_number)
        : x(x), y(y), z(z), q(q), atomic_number(atomic_number)
    {}

    Atom(const Atom& a)
        : x(a.x), y(a.y), z(a.z), q(a.q), atomic_number(a.atomic_number)
    {}

    /// Default construct makes a zero charge ghost atom at origin
    Atom()
        : x(0), y(0), z(0), q(0), atomic_number(0)
    {}

    template <typename Archive>
    void serialize(Archive& ar) {ar & x & y & z & q & atomic_number;}
};

std::ostream& operator<<(std::ostream& s, const Atom& atom);

class MolecularEntity {
private:
    // If you add more fields don't forget to serialize them
    std::vector<Atom> atoms;
    std::vector<double> rcut;  // Reciprocal of the smoothing radius

public:
    /// Makes a MolecularEntity with zero atoms
    MolecularEntity() : atoms() {};

    MolecularEntity(const std::string& filename, bool fractional);

    void read_file(const std::string& filename, bool fractional);

    void add_atom(double x, double y, double z, int atn, double q);

    int natom() const {return atoms.size();};

    void set_atom_coords(unsigned int i, double x, double y, double z);

    double bounding_cube() const;

    const Atom& get_atom(unsigned int i) const;

    void print() const;

    double inter_atomic_distance(unsigned int i,unsigned int j) const;

    double nuclear_repulsion_energy() const;

    double smallest_length_scale() const;

    void center();

    double total_nuclear_charge() const;

    double nuclear_attraction_potential(double x, double y, double z) const;

    double nuclear_charge_density(double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {ar & atoms & rcut;}
};

#define MENTITY_H_


#endif /* MENTITY_H_ */
