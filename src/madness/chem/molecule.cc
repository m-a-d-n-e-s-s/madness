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
#include <madness/world/world.h>
#include<madness/chem/molecule.h>
#include<madness/chem/gth_pseudopotential.h>
#include<madness/chem/atomutil.h>
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

Molecule::Molecule(std::vector<Atom> atoms, double eprec, CorePotentialManager core_pot, madness::Tensor<double> field) : atoms(std::move(atoms)), core_pot(std::move(core_pot)), field(std::move(field)) {
  atomic_radii.reserve(this->atoms.size());
  for(auto&& atom: this->atoms) {
    double radius =
        get_atomic_data(atom.get_atomic_number()).covalent_radius;
    atomic_radii.emplace_back(
        radius * 1e-10 /
        madness::constants::atomic_unit_of_length);
  }
  this->update_rcut_with_eprec(eprec);
}

  Molecule::Molecule(World& world, const commandlineparser& parser) :atoms(), rcut(), core_pot(), field(3L), parameters(world,parser)
{
    try {
        if (world.rank()==0) get_structure();
        world.gop.broadcast_serializable(*this, 0);
        MADNESS_CHECK(parameters.field().size()==3);
        for (int i=0; i<3; ++i) field(i)=parameters.field()[i];
    } catch (...) {
        std::cout << "\n\nsomething went wrong in the geometry input" << std::endl;
        std::cout << "geometry parameters " << std::endl << std::endl;
        parameters.print("geometry","end");
        throw madness::MadnessException("faulty geometry input",0,1,__LINE__,__FUNCTION__,__FILE__);
    }
}


void Molecule::print_parameters() {
    GeometryParameters param;
    madness::print("default parameters for the geometry input:\n");
    param.print("geometry","end");
    madness::print("");
    madness::print("");


    madness::print("If the molecular geometry is provided in the input file you need to specify");
    madness::print("the coordinates inside the geometry block\n");
    madness::print("Example:\n");
    madness::print("geometry");
    madness::print("  units  atomic ");
    madness::print("  O                     0                   0          0.21300717 ");
    madness::print("  H                     0           1.4265081         -0.85202867 ");
    madness::print("  H                     0          -1.4265081         -0.85202867 ");
    madness::print("end\n");
}

void Molecule::get_structure() {

    // read input parameters from the input file
    std::string sourcetype=parameters.source_type();
    std::string sourcename=parameters.source_name();

    if (sourcetype == "inputfile") {
        try {
            std::ifstream ifile(sourcename);
            read(ifile);
        } catch (const std::exception& err) {
            std::cout << "caught runtime error: " << err.what() << std::endl;
            std::cout << "failed to load geometry from input file" << std::endl;
            MADNESS_EXCEPTION("failed to load geometry from input file",1);
        }


    } else if (sourcetype == "library") {
        read_structure_from_library(sourcename);
    } else if (sourcetype == "xyz") {
        try {
            read_xyz(sourcename);
        } catch (std::exception& err) {
            std::cout << "could not read from xyz-file " << sourcename <<  std::endl;
            std::cout << "source " << sourcetype << "  " << sourcename << std::endl;
            MADNESS_EXCEPTION("failed to get geometry",1);
        }

    } else {
        std::cout << "could not determine molecule" << std::endl;
        std::cout << " source_type " << sourcetype << std::endl;
        std::cout << " source_name " << sourcename << std::endl;
        MADNESS_EXCEPTION("failed to get geometry",1);
    }


    // set derived parameters for the molecule

    //if psp_calc is true, set all atoms to PS atoms
    //if not, check whether some atoms are PS atoms or if this a pure AE calculation
    if (parameters.get<bool>("psp_calc")) {
        for (size_t iatom = 0; iatom < natom(); iatom++) {
            set_pseudo_atom(iatom, true);
        }
    }

    //modify atomic charge for complete PSP calc or individual PS atoms
    for (size_t iatom = 0; iatom < natom(); iatom++) {
        if (get_pseudo_atom(iatom)) {
            unsigned int an = get_atomic_number(iatom);
            double zeff = madness::get_charge_from_file("gth.xml", an);
            set_atom_charge(iatom, zeff);
        }
    }

    if (parameters.core_type() != "none") {
        read_core_file(parameters.core_type());
    }

    if (not parameters.no_orient()) orient();
    if (natom()==0) {
        MADNESS_EXCEPTION("no molecule was given",1);
    }


};

std::string Molecule::get_structure_library_path() {
    std::string chemdata_dir(MRA_CHEMDATA_DIR);
    if (getenv("MRA_CHEMDATA_DIR")) chemdata_dir=std::string(getenv("MRA_CHEMDATA_DIR"));
    return chemdata_dir+"/structure_library";
}

std::istream& Molecule::position_stream_in_library(std::ifstream& f, const std::string& name) {
    // get the location of the structure library
    std::string library=get_structure_library_path();

    f.open(library);
    std::string errmsg;

    if(f.fail()) {
        errmsg=std::string("Failed to open structure library: ") + library;
    } else {
        try {
            std::string full_line="structure="+name;
            madness::position_stream_to_word(f, full_line,'#',true,true);
        } catch (...) {
            errmsg = "could not find structure " + name + " in the library\n\n";
        }
    }
    MADNESS_CHECK_THROW(errmsg.empty(),errmsg.c_str());
    return f;
}

void Molecule::read_structure_from_library(const std::string& name) {

    try {
        std::ifstream f;
        position_stream_in_library(f,name);
        this->read(f);
    } catch (...) {
        std::string errmsg = "could not find structure " + name + " in the library\n\n";
        MADNESS_EXCEPTION(errmsg.c_str(), 0);
    }
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
    madness::position_stream(f, "geometry",false);		// do not rewind
    double scale = 1.0; // Default is atomic units
    if (parameters.units()=="angstrom") scale = 1e-10 / madness::constants::atomic_unit_of_length;

    std::string s, tag;
    while (std::getline(f,s)) {
        std::istringstream ss(s);
        ss >> tag;
        // check if tag is a keyword to be ignored
        bool ignore=false;
        for (auto& p : parameters.get_all_parameters()) {
            if (tag==p.first) ignore=true;
        }
        if (ignore) continue;
        if (tag == "end") {
            goto finished;
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
    update_rcut_with_eprec(parameters.eprec());
}


void Molecule::read_xyz(const std::string filename) {
    std::ifstream f(filename.c_str());
    if(f.fail()) {
        std::string errmsg = std::string("Failed to open file: ") + filename;
        MADNESS_EXCEPTION(errmsg.c_str(), 0);
    }

    atoms.clear();
    rcut.clear();
    MADNESS_CHECK(parameters.units()=="angstrom"); // xyz is always in angs
    const double scale = 1e-10 / madness::constants::atomic_unit_of_length;

    long natom_expected=1;
    long current_line=0;

    std::string s, tag;
    while (std::getline(f,s)) {
        std::istringstream ss(s);
        current_line++;
        if (current_line==1) {
            ss >> natom_expected;
            MADNESS_CHECK(natom_expected>0);
            continue;
        }
        if (current_line==2) continue;      // ignore comment line
        double xx, yy, zz;
        if (not (ss >> tag >>  xx >> yy >> zz)) {
            MADNESS_EXCEPTION(std::string("error reading the xyz input file"+filename).c_str(),1);
        };
        xx *= scale;
        yy *= scale;
        zz *= scale;
        int atn = symbol_to_atomic_number(tag);
        double qq = atn;
        if (atn == 0) ss >> qq; // Charge of ghost atom
        //check if pseudo-atom or not
        bool psat = check_if_pseudo_atom(tag);
        add_atom(xx,yy,zz,qq,atn,psat);
        if (current_line==natom_expected+2) break;
    }
    MADNESS_CHECK(size_t(natom_expected)==natom());
    update_rcut_with_eprec(parameters.eprec());
}

//version without pseudo-atoms
void Molecule::add_atom(double x, double y, double z, double q, int atomic_number) {
    atoms.push_back(Atom(x,y,z,q,atomic_number));
    double c = smoothing_parameter(q, get_eprec()); // eprec is error per atom
    //printf("smoothing param %.6f\n", c);
    double radius = get_atomic_data(atomic_number).covalent_radius;//Jacob added
    atomic_radii.push_back(radius*1e-10/madness::constants::atomic_unit_of_length);// Jacob added
    rcut.push_back(1.0/c);
}

//version specifying pseudo-atoms
void Molecule::add_atom(double x, double y, double z, double q, int atomic_number, bool pseudo_atom) {
    atoms.push_back(Atom(x,y,z,q,atomic_number,pseudo_atom));
    double c = smoothing_parameter(q, get_eprec()); // eprec is error per atom
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

unsigned int Molecule::get_atomic_number(unsigned int i) const {
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
void Molecule::update_rcut_with_eprec(double value) {
  if (value != get_eprec()) {
    parameters.set_user_defined_value("eprec", value);
  }
  rcut.clear();
  rcut.reserve(atoms.size());
  for (auto &&atom : atoms) {
    rcut.emplace_back(1.0 / smoothing_parameter(atom.q, value));
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

// Returns molecule in qc-schema format
// symbols (nat,) atom symbols in title case. array[string]
// geometry (3*nat,) vector of xyz coordinates [a0] of the atoms.  array[number]
// There are optional parameters yet to be implemented
// https://molssi-qc-schema.readthedocs.io/en/latest/auto_topology.html
nlohmann::json Molecule::to_json() const {
    nlohmann::json mol_schema;
    mol_schema["symbols"] = {};
    mol_schema["geometry"] = {};

//    get_atomic_data(atoms[0].atomic_number).symbol;
    for (size_t i = 0; i < natom(); ++i) {
        mol_schema["symbols"].push_back(get_atomic_data(atoms[i].atomic_number).symbol);
        mol_schema["geometry"].push_back({atoms[i].x, atoms[i].y, atoms[i].z});
    }
    return mol_schema;
}



void Molecule::print() const {
    std::string p =parameters.print_to_string();
    std::cout.flush();
    std::stringstream sstream;
    sstream << "geometry" << std::endl;
    sstream << p << std::endl;
//    sstream << "   eprec  " << std::scientific << std::setw(1) << get_eprec()  << std::endl << std::fixed;
//    sstream << "   units atomic" << std::endl;
    for (size_t i=0; i<natom(); ++i) {
        sstream << std::setw(5) << get_atomic_data(atoms[i].atomic_number).symbol << "  ";
        sstream << std::setw(20) << std::setprecision(8) << atoms[i].x
                << std::setw(20) << atoms[i].y
                << std::setw(20) << atoms[i].z;
        if (atoms[i].atomic_number == 0) sstream << "     " << atoms[i].q;
        sstream << std::endl;
    }
    sstream << "end" << std::endl;
    std::cout << sstream.str();
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


template <typename opT>
bool Molecule::test_for_op(opT op, const double symtol) const {
    // all atoms must have a symmetry-equivalent partner
    for (unsigned int i=0; i<atoms.size(); ++i) {
        if (find_symmetry_equivalent_atom(i,op,symtol)==-1) return false;
    }
    return true;
}


template <typename opT>
int Molecule::find_symmetry_equivalent_atom(int iatom, opT op, const double symtol) const  {
    double x=atoms[iatom].x, y=atoms[iatom].y, z=atoms[iatom].z;
    op(x, y, z);
    for (unsigned int j=0; j<atoms.size(); ++j) {
        double r = distance(x, y, z, atoms[j].x, atoms[j].y, atoms[j].z);
        if (r < symtol) return j;
    }
    return -1;
}

template <typename opT>
void Molecule::symmetrize_for_op(opT op, const double symtol) {

    MADNESS_CHECK(test_for_op(op, symtol));

    for (unsigned int iatom = 0; iatom < atoms.size(); ++iatom) {
        int jatom = find_symmetry_equivalent_atom(iatom, op, symtol);

        // symmetrize
        double x=atoms[iatom].x,y=atoms[iatom].y,z=atoms[iatom].z;
        op(x,y,z);
        double r = distance(atoms[jatom].x, atoms[jatom].y, atoms[jatom].z, x,y,z);
        MADNESS_CHECK(r<symtol);

        x=0.5*(x+atoms[jatom].x);
        y=0.5*(y+atoms[jatom].y);
        z=0.5*(z+atoms[jatom].z);
        atoms[iatom].x=x;
        atoms[iatom].y=y;
        atoms[iatom].z=z;
        op(x,y,z);
        atoms[jatom].x=x;
        atoms[jatom].y=y;
        atoms[jatom].z=z;

        // check result
        int jatom2=find_symmetry_equivalent_atom(iatom, op, 1.e-12);
        MADNESS_CHECK(jatom==jatom2);
    }
}


bool Molecule::test_for_c2(double xaxis, double yaxis, double zaxis, const double symtol) const {
    return test_for_op(apply_c2(xaxis,yaxis,zaxis), symtol);
}

bool Molecule::test_for_sigma(double xaxis, double yaxis, double zaxis, const double symtol) const {
    return test_for_op(apply_sigma(xaxis, yaxis, zaxis), symtol);
}

bool Molecule::test_for_inverse(const double symtol) const {
    return test_for_op(apply_inverse(0.0, 0.0, 0.0), symtol);
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

std::string Molecule::symmetrize_and_identify_point_group(const double symtol) {
    // C2 axes must be along the Cartesian axes and
    // mirror planes must be orthogonal to them

    bool x_is_c2 = test_for_c2(1.0,0.0,0.0,fabs(symtol));
    bool y_is_c2 = test_for_c2(0.0,1.0,0.0,fabs(symtol));
    bool z_is_c2 = test_for_c2(0.0,0.0,1.0,fabs(symtol));
    bool xy_is_sigma = test_for_sigma(0.0,0.0,1.0,fabs(symtol));
    bool xz_is_sigma = test_for_sigma(0.0,1.0,0.0,fabs(symtol));
    bool yz_is_sigma = test_for_sigma(1.0,0.0,0.0,fabs(symtol));
    bool inverse = test_for_inverse(fabs(symtol));

    // this is stupid, but necessary for backwards compatibility
    // FIXME SYMTOL
    if (symtol>0.0) {
        if (x_is_c2) symmetrize_for_op(apply_c2(1.0,0.0,0.0), symtol);
        if (y_is_c2) symmetrize_for_op(apply_c2(0.0,1.0,0.0), symtol);
        if (z_is_c2) symmetrize_for_op(apply_c2(0.0,0.0,1.0), symtol);
        if (xy_is_sigma) symmetrize_for_op(apply_sigma(0.0,0.0,1.0),symtol);
        if (xz_is_sigma) symmetrize_for_op(apply_sigma(0.0,1.0,0.0),symtol);
        if (yz_is_sigma) symmetrize_for_op(apply_sigma(1.0,0.0,0.0),symtol);
        if (inverse) symmetrize_for_op(apply_inverse(0.0,0.0,0.0),symtol);
    }

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
    return pointgroup;
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
		double symtol = fabs(parameters.symtol());
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
    pointgroup_= symmetrize_and_identify_point_group(parameters.symtol());
    if (parameters.symtol()>0.0) {
        std::string pointgroup_tight= symmetrize_and_identify_point_group(1.e-12);
        MADNESS_CHECK(pointgroup_tight==pointgroup_);
    }
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

    core_pot.read_file(filename, atomset, get_eprec());

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
            rcut[i] = 1.0/smoothing_parameter(q, get_eprec());
            //rcut[i] = 1.0/smoothing_parameter(q, 1.0);
            madness::print("rcut update", i, r, "to", rcut[i]);
        }
    }

    return;
}

}
