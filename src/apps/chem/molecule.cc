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


  $Id: test.cc 257 2007-06-25 19:09:38Z HartmanBaker $
*/


/// \file moldft/molecule.cc
/// \brief Simple management of molecular information and potential

#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/constants.h>
#include <chem/molecule.h>
#include <chem/atomutil.h>
#include <madness/misc/misc.h>
#include <iomanip>
#include <set>

namespace madness {
static inline double distance(double x1, double y1, double z1, double x2, double y2, double z2) {
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return sqrt(xx*xx + yy*yy + zz*zz);
}

static inline double distance_sq(double x1, double y1, double z1, double x2, double y2, double z2) {
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return xx*xx + yy*yy + zz*zz;
}


/// read molecule from the input file and return part of the header for
/// a Gaussian cube file.
/// @param[in]  filename input file name (usually "input")
std::vector<std::string> cubefile_header(std::string filename, const bool& no_orient=false) {
    Molecule molecule;
    molecule.read_file(filename);
    if(no_orient==false) molecule.orient();
    return molecule.cubefile_header();
}


std::ostream& operator<<(std::ostream& s, const Atom& atom) {
    s << "Atom([" << atom.x << ", " << atom.y << ", " << atom.z << "], " << atom.q << "," << atom.atomic_number << ")";
    return s;
}

/// Read coordinates from a file

/// Scans the file for the first geometry block in the format
/// \code
///    geometry
///       tag x y z
///       ...
///    end
/// \endcode
/// The charge \c q is inferred from the tag which is
/// assumed to be the standard symbol for an element.
/// Same as the simplest NWChem format.  For ghost
/// atoms (\c bq ) the  charge is read as a fifth field
/// on the line.
///
/// This code is just for the examples ... don't trust it!
Molecule::Molecule(const std::string& filename) :
		atoms(), rcut(), eprec(1e-4), core_pot(), field(3L) {
    read_file(filename);
}

std::vector<std::string> Molecule::cubefile_header() const {
	std::vector<std::string> molecular_info;
	for (unsigned int i = 0; i < natom(); ++i) {
		std::stringstream ss;
		const int charge = get_atom(i).get_atomic_number();
		ss << charge << " " << charge << " ";
		ss << std::fixed;
		ss.precision(8);
		const Vector<double, 3> coord = get_atom(i).get_coords();
		ss << coord[0] << " " << coord[1] << " " << coord[2] << " \n";
		molecular_info.push_back(ss.str());
	}
	return molecular_info;
}

void Molecule::read_file(const std::string& filename) {
    std::ifstream f(filename.c_str());
    if(f.fail()) {
        std::string errmsg = std::string("Failed to open file: ") + filename;
        MADNESS_EXCEPTION(errmsg.c_str(), 0);
    }
    this->read(f);
}

void Molecule::read(std::istream& f) {
    atoms.clear();
    rcut.clear();
    eprec = 1e-4;
    units = atomic;
    madness::position_stream(f, "geometry");
    double scale = 1.0; // Default is atomic units

    std::string s;
    while (std::getline(f,s)) {
        std::istringstream ss(s);
        std::string tag;
        ss >> tag;
        if (tag == "end") {
            goto finished;
        }
        else if (tag == "units") {
            if (natom()) throw "Molecule: read_file: presently units must be the first line of the geometry block";
            ss >> tag;
            if (tag == "a.u." || tag == "au" || tag == "atomic") {
                std::cout << "\nAtomic units being used to read input coordinates\n\n";
                scale = 1.0;
            }
            else if (tag == "angstrom" || tag == "angs") {
                units = angstrom;
                scale = 1e-10 / madness::constants::atomic_unit_of_length;
                printf("\nAngstrom being used to read input coordinates (1 Bohr = %.8f Angstrom)\n\n", scale);
            }
            else {
                throw "Molecule: read_file: unknown units requested";
            }
        }
        else if (tag == "eprec") {
            ss >> eprec;
            // adapt nuclear smoothing factor rcut
            this->set_eprec(eprec);
        }
        else if (tag == "field") {
            ss >> field[0] >> field[1] >> field[2];
        }
        else {
            double xx, yy, zz;
            ss >> xx >> yy >> zz;
            xx *= scale;
            yy *= scale;
            zz *= scale;
            int atn = symbol_to_atomic_number(tag);
            double qq = atn;
            if (atn == 0) ss >> qq; // Charge of ghost atom
            //check if pseudo-atom or not
            bool psat = check_if_pseudo_atom(tag);
            add_atom(xx,yy,zz,qq,atn,psat);
        }
    }
    throw "No end to the geometry in the input file";
finished:
    ;
}

//version without pseudo-atoms
void Molecule::add_atom(double x, double y, double z, double q, int atomic_number) {
    atoms.push_back(Atom(x,y,z,q,atomic_number));
    double c = smoothing_parameter(q, eprec); // eprec is error per atom
    //printf("smoothing param %.6f\n", c);
    double radius = get_atomic_data(atomic_number).covalent_radius;//Jacob added
    atomic_radii.push_back(radius*1e-10/madness::constants::atomic_unit_of_length);// Jacob added
    rcut.push_back(1.0/c);
}

//version specifying pseudo-atoms
void Molecule::add_atom(double x, double y, double z, double q, int atomic_number, bool pseudo_atom) {
    atoms.push_back(Atom(x,y,z,q,atomic_number,pseudo_atom));
    double c = smoothing_parameter(q, eprec); // eprec is error per atom
    //printf("smoothing param %.6f\n", c);
    double radius = get_atomic_data(atomic_number).covalent_radius;//Jacob added
    atomic_radii.push_back(radius*1e-10/madness::constants::atomic_unit_of_length);// Jacob added
    rcut.push_back(1.0/c);
}

void Molecule::set_atom_charge(unsigned int i, double zeff) {
  if (i>=atoms.size()) throw "trying to set charge of invalid atom";
  atoms[i].q = zeff;
}

unsigned int Molecule::get_atom_charge(unsigned int i) const {
  if (i>=atoms.size()) throw "trying to get charge of invalid atom";
  return atoms[i].q;
}

unsigned int Molecule::get_atom_number(unsigned int i) {
  if (i>=atoms.size()) throw "trying to get number of invalid atom";
  return atoms[i].atomic_number;
}

void Molecule::set_pseudo_atom(unsigned int i, bool psat) {
  if (i>=atoms.size()) throw "trying to set charge of invalid atom";
  atoms[i].pseudo_atom = psat;
}

bool Molecule::get_pseudo_atom(unsigned int i) const {
  if (i>=atoms.size()) throw "trying to get pseudo atom for invalid atom";
  return atoms[i].pseudo_atom;
}

void Molecule::set_atom_coords(unsigned int i, double x, double y, double z) {
    if (i>=atoms.size()) throw "trying to set coords of invalid atom";
    atoms[i].x = x;
    atoms[i].y = y;
    atoms[i].z = z;
}

madness::Tensor<double> Molecule::get_all_coords() const {
    madness::Tensor<double> c(natom(),3);
    for (size_t i=0; i<natom(); ++i) {
        const Atom atom = get_atom(i);
        c(i,0) = atom.x;
        c(i,1) = atom.y;
        c(i,2) = atom.z;
    }
    return c;
}

std::vector< madness::Vector<double,3> > Molecule::get_all_coords_vec() const {
  std::vector< madness::Vector<double,3> > c(natom());
  for (size_t i=0; i<natom(); ++i) {
    const Atom atom = get_atom(i);
    c[i][0] = atom.x;
    c[i][1] = atom.y;
    c[i][2] = atom.z;
  }
  return c;
}

void Molecule::set_all_coords(const madness::Tensor<double>& c) {
    MADNESS_ASSERT(c.ndim()==2u && size_t(c.dims()[0])==natom() && c.dims()[1]==3);
    for (size_t i=0; i<natom(); ++i) {
        atoms[i].x = c(i,0);
        atoms[i].y = c(i,1);
        atoms[i].z = c(i,2);
    }
}

/// updates rcuts with given eprec
void Molecule::set_eprec(double value) {
    eprec = value;
    for (size_t i=0; i<atoms.size(); ++i) {
        rcut[i] = 1.0 / smoothing_parameter(atoms[i].q, eprec);
    }
    core_pot.set_eprec(value);
}

void Molecule::set_rcut(double value) {
    for (size_t i=0; i<atoms.size(); ++i) {
        rcut[i] = (value<=0.0) ? 1.0 : value;
    }
}

const Atom& Molecule::get_atom(unsigned int i) const {
    if (i>=atoms.size()) throw "trying to get coords of invalid atom";
    return atoms[i];
}

void Molecule::print() const {
    std::cout.flush();
    printf(" geometry\n");
    printf("   eprec %.1e\n", eprec);
    printf("   units atomic\n");
    //printf("   Finite Field %11.8f %20.8f %20.8f\n", field[0], field[1], field[2]);
    for (size_t i=0; i<natom(); ++i) {
        printf("   %-2s  %20.8f %20.8f %20.8f", get_atomic_data(atoms[i].atomic_number).symbol,
               atoms[i].x, atoms[i].y, atoms[i].z);
        if (atoms[i].atomic_number == 0) printf("     %20.8f", atoms[i].q);
        printf("\n");
    }
    printf(" end\n");
}

double Molecule::inter_atomic_distance(unsigned int i,unsigned int j) const {
    if (i>=atoms.size()) throw "trying to compute distance with invalid atom";
    if (j>=atoms.size()) throw "trying to compute distance with invalid atom";
    return distance(atoms[i].x, atoms[i].y, atoms[i].z,
                    atoms[j].x, atoms[j].y, atoms[j].z);
}

double Molecule::nuclear_repulsion_energy() const {
    double sum = 0.0;
    unsigned int z1, z2;
    for (size_t i=0; i<atoms.size(); ++i) {
        if (atoms[i].pseudo_atom){
            z1 = atoms[i].q;}
        else{
            z1 = atoms[i].atomic_number;}
        if (core_pot.is_defined(z1)) z1 -= core_pot.n_core_orb(z1) * 2;
        for (size_t j=i+1; j<atoms.size(); ++j) {
            if (atoms[j].pseudo_atom){
                z2 = atoms[j].q;}
            else{
                z2 = atoms[j].atomic_number;}
            if (core_pot.is_defined(z2)) z2 -= core_pot.n_core_orb(z2) * 2;
            sum += z1 * z2 / inter_atomic_distance(i,j);
        }
    }
    return sum;
}

double Molecule::nuclear_dipole(int axis) const {
    double sum = 0.0;
    for (size_t atom = 0; atom < atoms.size(); ++atom) {
        unsigned int z;
        if (atoms[atom].pseudo_atom){
            z = atoms[atom].q;}
        else{
            z = atoms[atom].atomic_number;}
        if (core_pot.is_defined(z)) z -= core_pot.n_core_orb(z) * 2;
        double r;
        switch (axis) {
            case 0: r = atoms[atom].x; break;
            case 1: r = atoms[atom].y; break;
            case 2: r = atoms[atom].z; break;
            default: MADNESS_EXCEPTION("invalid axis", 0);
        }
        sum += r*z;
    }
    return sum;
}

Tensor<double> Molecule::nuclear_dipole_derivative(const int atom, const int axis) const{
    Tensor<double> mu_x(3);
    mu_x(axis)=atoms[atom].q;
    return mu_x;
}


double Molecule::nuclear_repulsion_derivative(size_t i, int axis) const {
    double sum = 0.0;
    unsigned int z1 = atoms[i].atomic_number;
    if (core_pot.is_defined(z1)) z1 -= core_pot.n_core_orb(z1) * 2;
    for (size_t j=0; j<atoms.size(); ++j) {
        if (j != i) {
            size_t z2 = atoms[j].atomic_number;
            if (core_pot.is_defined(z2)) z2 -= core_pot.n_core_orb(z2) * 2;
            double r = inter_atomic_distance(i,j);
            double xx;
            if (axis == 0) xx = atoms[i].x - atoms[j].x;
            else if (axis == 1) xx = atoms[i].y - atoms[j].y;
            else xx = atoms[i].z - atoms[j].z;
            sum -= xx * z1 * z2/ (r * r * r);
        }
    }
    return sum;
}

/// compute the nuclear-nuclear contribution to the molecular hessian
Tensor<double> Molecule::nuclear_repulsion_hessian() const {

    Tensor<double> hessian(3*natom(),3*natom());
    for (size_t iatom=0; iatom<natom(); ++iatom) {
        for (int iaxis=0; iaxis<3; ++iaxis) {
            for (size_t jatom=0; jatom<natom(); ++jatom) {
                for (int jaxis=0; jaxis<3; ++jaxis) {
                    hessian(3*iatom+iaxis, 3*jatom+jaxis)=
                            nuclear_repulsion_second_derivative(iatom,jatom,iaxis,jaxis);
                }
            }
        }
    }
    return hessian;
}


/// compute the nuclear-nuclear contribution to the second derivatives

/// @param[in]  iatom   the i-th atom (row of the hessian)
/// @param[in]  jatom   the j-th atom (column of the hessian)
/// @param[in]  iaxis   the xyz axis of the i-th atom
/// @param[in]  jaxis   the xyz axis of the j-th atom
/// return the (3*iatom + iaxis, 3*jatom + jaxis) matix element
double Molecule::nuclear_repulsion_second_derivative(int iatom, int jatom,
        int iaxis, int jaxis) const {

    double sum = 0.0;
    unsigned int ZA = atoms[iatom].atomic_number;
    unsigned int ZB = atoms[jatom].atomic_number;

    Tensor<double> RA(3), RB(3);
    RA(0l)=atoms[iatom].x; RA(1)=atoms[iatom].y; RA(2)=atoms[iatom].z;
    RB(0l)=atoms[jatom].x; RB(1)=atoms[jatom].y; RB(2)=atoms[jatom].z;

    if (core_pot.is_defined(ZA)) MADNESS_EXCEPTION("no core potentials in the hessian",1);
    if (core_pot.is_defined(ZB)) MADNESS_EXCEPTION("no core potentials in the hessian",1);

    // first term is (for A\neq B, i.e. iatom/=jatom):
    // \frac{\partial^2}{\partial X_A\partial Y_B}\frac{Z_AZ_B}{R_{AB}}
    if (iatom != jatom) {
        const double RAB = inter_atomic_distance(iatom,jatom);
        if (iaxis==jaxis) {
            sum+=(RAB*RAB - 3*std::pow(RA(iaxis)-RB(iaxis),2.0))/std::pow(RAB,5.0);
        } else {
            sum-=3*(RA(iaxis)-RB(iaxis))*(RA(jaxis)-RB(jaxis)) /std::pow(RAB,5.0);
        }
        sum*=ZA*ZB;
    }

    // second term is (for A==B, i.e. iatom==jatom):
    // \sum_{C\neq A}\frac{\partial^2}{\partial X_A\partial X_B}Z_AZ_B\frac{Z_CZ_B}{R_{AC}}
    if (iatom==jatom) {
        for (unsigned int katom=0; katom<atoms.size(); ++katom) {
            double RAC = inter_atomic_distance(iatom,katom);
            Tensor<double> RC(3);
            RC(0l)=atoms[katom].x; RC(1)=atoms[katom].y; RC(2)=atoms[katom].z;
            const unsigned int ZC=atoms[katom].atomic_number;

            if (katom != (unsigned int)(iatom)) {
                if (iaxis==jaxis) {
                    sum-=ZA*ZC*(RAC*RAC - 3.0*std::pow(RA(iaxis)-RC(iaxis),2.0))/std::pow(RAC,5.0);
                } else {
                    sum+=ZA*ZC*3.0*(RA(iaxis)-RC(iaxis))*(RA(jaxis)-RC(jaxis)) /std::pow(RAC,5.0);
                }

            }
        }
    }
    return sum;
}


double Molecule::smallest_length_scale() const {
    double rcmax = 0.0;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        rcmax = std::max(rcmax,rcut[i]);
    }
    return 1.0/rcmax;
}


/// Moves the center of nuclear charge to the origin
void Molecule::center() {
    double xx=0.0, yy=0.0, zz=0.0, qq=0.0;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        xx += atoms[i].x*atoms[i].q;
        yy += atoms[i].y*atoms[i].q;
        zz += atoms[i].z*atoms[i].q;
        qq += atoms[i].q;
    }
    xx /= qq;
    yy /= qq;
    zz /= qq;
    Tensor<double> translation(3);
    translation(0l)=-xx;
    translation(1l)=-yy;
    translation(2l)=-zz;
    translate(translation);
}

/// translate the molecule
 void Molecule::translate(const Tensor<double>& translation) {
    for (unsigned int i=0; i<atoms.size(); ++i) {
        atoms[i].x += translation(0l);
        atoms[i].y += translation(1l);
        atoms[i].z += translation(2l);
    }
}

/// Apply to (x,y,z) a C2 rotation about an axis thru the origin and (xaxis,yaxis,zaxis)
static void apply_c2(double xaxis, double yaxis, double zaxis, double& x, double& y, double& z) {
    double raxissq = xaxis*xaxis + yaxis*yaxis + zaxis*zaxis;
    double dx = x*xaxis*xaxis/raxissq;
    double dy = y*yaxis*yaxis/raxissq;
    double dz = z*zaxis*zaxis/raxissq;

    x = 2.0*dx - x;
    y = 2.0*dy - y;
    z = 2.0*dz - z;
}

/// Apply to (x,y,z) a reflection through a plane containing the origin with normal (xaxis,yaxis,zaxis)
static void apply_sigma(double xaxis, double yaxis, double zaxis, double& x, double& y, double& z) {
    double raxissq = xaxis*xaxis + yaxis*yaxis + zaxis*zaxis;
    double dx = x*xaxis*xaxis/raxissq;
    double dy = y*yaxis*yaxis/raxissq;
    double dz = z*zaxis*zaxis/raxissq;

    x = x - 2.0*dx;
    y = y - 2.0*dy;
    z = z - 2.0*dz;
}

static void apply_inverse(double xjunk, double yjunk, double zjunk, double& x, double& y, double& z) {
    x = -x;
    y = -y;
    z = -z;
}

template <typename opT>
bool Molecule::test_for_op(double xaxis, double yaxis, double zaxis, opT op) const {
    const double symtol = 1e-2;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        double x=atoms[i].x, y=atoms[i].y, z=atoms[i].z;
        op(xaxis, yaxis, zaxis, x, y, z);
        bool found = false;
        for (unsigned int j=0; j<atoms.size(); ++j) {
            double r = distance(x, y, z, atoms[j].x, atoms[j].y, atoms[j].z);
            if (r < symtol) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

bool Molecule::test_for_c2(double xaxis, double yaxis, double zaxis) const {
    return test_for_op(xaxis, yaxis, zaxis, apply_c2);
}

bool Molecule::test_for_sigma(double xaxis, double yaxis, double zaxis) const {
    return test_for_op(xaxis, yaxis, zaxis, apply_sigma);
}

bool Molecule::test_for_inverse() const {
    return test_for_op(0.0, 0.0, 0.0, apply_inverse);
}

void Molecule::swapaxes(int ix, int iy) {
    for (unsigned int i=0; i<atoms.size(); ++i) {
        double r[3] = {atoms[i].x, atoms[i].y, atoms[i].z};
        std::swap(r[ix],r[iy]);
        atoms[i].x=r[0]; atoms[i].y=r[1]; atoms[i].z=r[2];
    }

    // field rotation
    double r[3] = {field[0], field[1], field[2]};
    std::swap(r[ix],r[iy]);
    field[0]=r[0]; field[1]=r[1]; field[2]=r[2];
}

void Molecule::identify_point_group() {
    // C2 axes must be along the Cartesian axes and
    // mirror planes must be orthogonal to them

    bool x_is_c2 = test_for_c2(1.0,0.0,0.0);
    bool y_is_c2 = test_for_c2(0.0,1.0,0.0);
    bool z_is_c2 = test_for_c2(0.0,0.0,1.0);
    bool xy_is_sigma = test_for_sigma(0.0,0.0,1.0);
    bool xz_is_sigma = test_for_sigma(0.0,1.0,0.0);
    bool yz_is_sigma = test_for_sigma(1.0,0.0,0.0);
    bool inverse = test_for_inverse();

    /*
      .   (i,c,s)
      Ci  (1,0,0) --- i
      C2h (1,1,1) --- c2z, sxy, i
      D2h (1,3,3) --- c2z, c2y, c2x, sxy, sxz, syz, i

      C1  (0,0,0) ---
      Cs  (0,0,1) --- sxy
      C2  (0,1,0) --- c2z
      C2v (0,1,2) --- c2z, sxz, syz
      D2  (0,3,0) --- c2z, c2y, c2x

     */

    int nc2 = int(x_is_c2) + int(y_is_c2) + int(z_is_c2);
    int nsi = int(xy_is_sigma) + int(xz_is_sigma) + int(yz_is_sigma);

    std::string pointgroup;
    if (inverse && nc2==0 && nsi==0) {
        pointgroup = "Ci";
    }
    else if (inverse && nc2==1 && nsi==1) {
        pointgroup = "C2h";
        if (x_is_c2) swapaxes(0,2);
        if (y_is_c2) swapaxes(1,2);
    }
    else if (inverse && nc2==3 && nsi==3) {
        pointgroup = "D2h";
    }
    else if (!inverse && nc2==0 && nsi==0) {
        pointgroup = "C1";
    }
    else if (!inverse && nc2==0 && nsi==1) {
        pointgroup = "Cs";
        if (xz_is_sigma) swapaxes(1,2);
        if (yz_is_sigma) swapaxes(0,2);
    }
    else if (!inverse && nc2==1 && nsi==0) {
        pointgroup = "C2";
        if (x_is_c2) swapaxes(0,2);
        if (y_is_c2) swapaxes(1,2);
    }
    else if (!inverse && nc2==1 && nsi==2) {
        pointgroup = "C2v";
        if (x_is_c2) swapaxes(0,2);
        if (y_is_c2) swapaxes(1,2);
    }
    else if (!inverse && nc2==3 && nsi==0) {
        pointgroup = "D2";
    }
    else {
        madness::print("Not-quite-symmetric geometry (clean up to fix), will assume C1");
        pointgroup = "C1";
    }

//    madness::print("\n The point group is", pointgroup);
    // assign to member variable
    pointgroup_ = pointgroup;
}

/// compute the center of mass
Tensor<double> Molecule::center_of_mass() const {
    Tensor<double> com(3);
    double xx=0.0, yy=0.0, zz=0.0, qq=0.0;
    for (unsigned int i=0; i<natom(); ++i) {
        xx += get_atom(i).x*get_atom(i).mass;
        yy += get_atom(i).y*get_atom(i).mass;
        zz += get_atom(i).z*get_atom(i).mass;
        qq += get_atom(i).mass;
    }
    com(0l)=xx/qq;
    com(1l)=yy/qq;
    com(2l)=zz/qq;
    return com;
}

// Align molecule with axes of inertia
Tensor<double> Molecule::moment_of_inertia() const {
    madness::Tensor<double> I(3L,3L);
    for (unsigned int i=0; i<atoms.size(); ++i) {
        double q = atoms[i].mass, x[3] = {atoms[i].x, atoms[i].y, atoms[i].z};
        for (int j=0; j<3; ++j) {
            I(j,j)+=q*(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
            for (int k=0; k<3; ++k)
                I(j,k) -= q*x[j]*x[k];
        }
    }
    return I;
}

/// Centers and orients the molecule in a standard manner
void Molecule::orient(bool verbose) {

    center();

    // Align molecule with axes of charge inertia
    madness::Tensor<double> I(3L,3L);
    for (unsigned int i=0; i<atoms.size(); ++i) {
        double q = atoms[i].atomic_number, x[3] = {atoms[i].x, atoms[i].y, atoms[i].z};
        for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k)
                I(j,k) += q*x[j]*x[k];
    }
    madness::Tensor<double> U, e;
    madness::syev(I, U, e);
    // madness::print("Moment of inertia eigenvalues and tensor\n");
    // madness::print(I);
    // madness::print(U);
    // madness::print(e);

    // rotate the molecule and the external field
    rotate(U);


    if (verbose) {
		// Try to resolve degenerate rotations
		double symtol = 1e-2;
		if (fabs(e[0]-e[1])<symtol && fabs(e[1]-e[2])<symtol) {
			madness::print("Cubic point group");
		}
		else if (fabs(e[0]-e[1])<symtol) {
			madness::print("XY degenerate");
		}
		else if (fabs(e[0]-e[2])<symtol) {
			madness::print("XZ degenerate");
		}
		else if (fabs(e[1]-e[2])<symtol) {
			madness::print("YZ degenerate");
		}
		else {
			madness::print("Abelian pointgroup");
		}
    }

    // Now hopefully have mirror planes and C2 axes correctly oriented
    // Figure out what elements are actually present and enforce
    // conventional ordering
    identify_point_group();
}

/// rotates the molecule and the external field

/// @param[in]  D   the rotation matrix
void Molecule::rotate(const Tensor<double>& D) {
    madness::Tensor<double> r(3L), rU;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        r[0]=atoms[i].x; r[1]=atoms[i].y, r[2]= atoms[i].z;
        rU = inner(r,D);
        atoms[i].x=rU[0]; atoms[i].y=rU[1]; atoms[i].z=rU[2];
    }
    // field rotation
    rU = inner(field,D);
    field[0] = rU[0];
    field[1] = rU[1];
    field[2] = rU[2];
}

/// Returns the half width of the bounding cube

/// The molecule will be contained in the cube [-L,+L].
double Molecule::bounding_cube() const {
    double L = 0.0;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        L = std::max(L, fabs(atoms[i].x));
        L = std::max(L, fabs(atoms[i].y));
        L = std::max(L, fabs(atoms[i].z));
    }
    return L;
}

double Molecule::total_nuclear_charge() const {
    double sum = 0.0;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        sum += atoms[i].q;
    }
    return sum;
}
//Nuclear charge density of the molecule
double Molecule::mol_nuclear_charge_density(double x, double y, double z) const {
    // Only one atom will contribute due to the short range of the nuclear            
    // charge density                                                                                                                                        
    for (unsigned int i=0; i<atoms.size(); i++) {
        double r = distance(x, y, z, atoms[i].x, atoms[i].y, atoms[i].z)*rcut[i];
        if (r < 6.0) {
            return atoms[i].atomic_number*smoothed_density(r)*rcut[i]*rcut[i]*rcut[i];
        }
    }
    return  0.0;
}
double Molecule::nuclear_attraction_potential(double x, double y, double z) const {
    // This is very inefficient since it scales as O(ngrid*natom)
    // ... we can easily make an O(natom) version using
    // the integral operator and sparse projection of an effective
    // density ... its potential can be evaluated at the same
    // time as the electronic Coulomb potential so it will be
    // essentially free.

    double sum = 0.0;
    for (unsigned int i=0; i<atoms.size(); ++i) {
        //make sure this isn't a pseudo-atom
        if (atoms[i].pseudo_atom) continue;

        double r = distance(atoms[i].x, atoms[i].y, atoms[i].z, x, y, z);
        //sum -= atoms[i].q/(r+1e-8);
        sum -= atoms[i].q * smoothed_potential(r*rcut[i])*rcut[i];
    }

    // field contribution
    sum += field[0] * x + field[1] * y + field[2] * z;

    return sum;
}

double Molecule::atomic_attraction_potential(int iatom, double x, double y,
        double z) const {

    //make sure this isn't a pseudo-atom
    if (atoms[iatom].pseudo_atom) return 0.0;

    double r = distance(atoms[iatom].x, atoms[iatom].y, atoms[iatom].z, x, y, z);
    double sum =- atoms[iatom].q * smoothed_potential(r*rcut[iatom])*rcut[iatom];

    return sum;
}


double Molecule::nuclear_attraction_potential_derivative(int atom, int axis, double x, double y, double z) const {
    double r = distance(atoms[atom].x, atoms[atom].y, atoms[atom].z, x, y, z);
    double rc = rcut[atom];
    double coord;
    if (axis == 0) coord = x-atoms[atom].x;
    else if (axis == 1) coord = y-atoms[atom].y;
    else coord = z-atoms[atom].z;

    //double q = atoms[atom].q;
    //if (core_pot.is_defined(atoms[atom].atomic_number)) {
    //    q -= core_pot.n_core_orb(atoms[atom].atomic_number) * 2;
    //    rc = 1.0/smoothing_parameter(q, eprec);
    //}

    double dv = atoms[atom].q * (coord / r) * dsmoothed_potential(r * rc) * (rc * rc);
    //double dv = q * (coord / r) * dsmoothed_potential(r * rc) * (rc * rc);
    double df = field[axis];
    return dv + df;
}

/// the second derivative of the (smoothed) nuclear potential Z/r

/// \f[
/// V(R,r_{el}) -V(r) =\frac{Z}{|R-r_{el}|} \approx Z u(r) \f]
/// with
/// \f[
/// \frac{\partial^2 V}{\partial X_i\partial X_j}
///   =  Z \left(\frac{\partial^2 u(r)}{\partial r^2} \frac{\partial r}{\partial X_i}
///      \frac{\partial r}{\partial X_j} + \frac{\partial u}{\partial r}
///      \frac{\partial^2 r}{\partial X_i \partial X_j}\right)
/// \f]
double Molecule::nuclear_attraction_potential_second_derivative(int atom,
        int iaxis, int jaxis, double x, double y, double z) const {

    const Vector<double,3> rr={x-atoms[atom].x,y-atoms[atom].y,z-atoms[atom].z};
    double r = rr.normf();
    double rc = rcut[atom];

    double u=smoothed_potential(r*rc)*rc;
    double d2u=d2smoothed_potential(r * rc) * (rc * rc * rc);
//    double rreg=r+1.e-5;
//    double rinv=1./(rreg);
//    double r3inv = 1.0/(rreg * rreg* rreg);
    double rinv=u;
    double r3inv = 0.5*d2u;

    double di=rr[iaxis]*rinv;
    double dj=rr[jaxis]*rinv;
    double term1=3.0*r3inv*di*dj;

//    double term2=-di*dj*r3inv;
//    if (iaxis==jaxis) term2+=rinv;

    double result=-atoms[atom].q * (term1);
    if (iaxis==jaxis) return 0.0;
    return result;

}


double Molecule::nuclear_charge_density(double x, double y, double z) const {
  // Only one atom will contribute due to the short range of the nuclear charge density

  for (unsigned int i=0; i<atoms.size(); i++) {
      double rsq = distance_sq(atoms[i].x, atoms[i].y, atoms[i].z, x, y, z)*rcut[i]*rcut[i];
      if (rsq < 36.0) {
          double r = sqrt(rsq);
          return atoms[i].q * smoothed_density(r)*rcut[i]*rcut[i]*rcut[i];
      }
  }
  return 0.0;
}


unsigned int Molecule::n_core_orb_all() const {
    int natom = atoms.size();
    unsigned int sum = 0;

    for (int i=0; i<natom; ++i) {
        unsigned int atn = atoms[i].atomic_number;
        if (core_pot.is_defined(atn)) sum += core_pot.n_core_orb(atn);
    }

    return sum;
}

double Molecule::core_eval(int atom, unsigned int core, int m, double x, double y, double z) const {
    unsigned int atn = atoms[atom].atomic_number;
    double xx = x - atoms[atom].x;
    double yy = y - atoms[atom].y;
    double zz = z - atoms[atom].z;
    double rsq = xx*xx + yy*yy + zz*zz;
    return core_pot.core_eval(atn, core, m, rsq, xx, yy, zz);
}

double Molecule::core_derivative(int atom, int axis, unsigned int core, int m, double x, double y, double z) const {
    unsigned int atn = atoms[atom].atomic_number;
    double xx = x - atoms[atom].x;
    double yy = y - atoms[atom].y;
    double zz = z - atoms[atom].z;
    double rsq = xx*xx + yy*yy + zz*zz;
    double xi;
    if (axis == 0) xi = xx;
    else if (axis == 1) xi = yy;
    else xi = zz;
    return core_pot.core_derivative(atn, core, m, axis, xi, rsq, xx, yy, zz);
}

double Molecule::molecular_core_potential(double x, double y, double z) const {
    int natom = atoms.size();
    double sum = 0.0;

    for (int i=0; i<natom; ++i) {
        unsigned int atn = atoms[i].atomic_number;
        if (core_pot.is_defined(atn)) {
            double r = distance(atoms[i].x, atoms[i].y, atoms[i].z, x, y, z);
            sum += core_pot.potential(atn, r);
        }
    }

    return sum;
}

double Molecule::core_potential_derivative(int atom, int axis, double x, double y, double z) const {
    int natom = atoms.size();
    if (natom <= atom) return 0.0;

    unsigned int atn = atoms[atom].atomic_number;
    //if (!core_pot.is_defined(atn)) return 0.0;

    double xi;
    if (axis == 0) xi = x-atoms[atom].x;
    else if (axis == 1) xi = y-atoms[atom].y;
    else xi = z-atoms[atom].z;
    double r = distance(atoms[atom].x, atoms[atom].y, atoms[atom].z, x, y, z);
    return core_pot.potential_derivative(atn, xi, r);
}

void Molecule::read_core_file(const std::string& filename) {
    std::set<unsigned int> atomset;
    int natom = atoms.size();
    for (int i=0; i<natom; ++i) {
        if (atomset.count(atoms[i].atomic_number) == 0)
            atomset.insert(atoms[i].atomic_number);
    }

    core_pot.read_file(filename, atomset, eprec);

    //return;

    // rcut update
    for (int i=0; i<natom; ++i) {
        unsigned int atn = atoms[i].atomic_number;
        if (core_pot.is_defined(atn)) {
            double q = atoms[i].q - core_pot.n_core_orb(atn) * 2;
            if (q == 0.0) {
                rcut[i] = 1.0;
                continue;
            }
            double r = rcut[i];
            rcut[i] = 1.0/smoothing_parameter(q, eprec);
            //rcut[i] = 1.0/smoothing_parameter(q, 1.0);
            madness::print("rcut update", i, r, "to", rcut[i]);
        }
    }

    return;
}

}
