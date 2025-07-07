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

#include<madness/chem/molecularbasis.h>
#include "NWChem.h"

namespace madness {

std::ostream& operator<<(std::ostream& s, const ContractedGaussianShell& c) {
    static const char* tag[] = {"s","p","d","f","g"};
    const std::size_t maxbufsize=32768;
    char buf[maxbufsize];
    std:size_t bufsize= maxbufsize;
    char* p = buf;
    const std::vector<double>& coeff = c.get_coeff();
    const std::vector<double>& expnt = c.get_expnt();

    p += snprintf(p,bufsize,"%s [",tag[c.angular_momentum()]);
    bufsize-=5;
    for (int i=0; i<c.nprim(); ++i) {
        p += snprintf(p,bufsize, "%.6f(%.6f)",coeff[i],expnt[i]);
        bufsize-=14;
        if (i != (c.nprim()-1)) {
            p += snprintf(p,bufsize, ", ");
            bufsize-=2;
        }
    }
    p += snprintf(p,bufsize, "]");
    s << buf;
    return s;
}

std::ostream& operator<<(std::ostream& s, const AtomicBasis& c) {
    const std::vector<ContractedGaussianShell>& shells = c.get_shells();
    for (int i=0; i<c.nshell(); ++i) {
        s << "     " << shells[i] << std::endl;
    }
    if (c.has_guess_info()) {
        s << "     " << "Guess density matrix" << std::endl;
        s << c.get_dmat();
    }
    if (c.has_guesspsp_info()) {
        s << "     " << "Guess density matrix (psp)" << std::endl;
        s << c.get_dmatpsp();
    }

    return s;
}

std::ostream& operator<<(std::ostream& s, const AtomicBasisFunction& a) {
    a.print_me(s);
    return s;
}

void AtomicBasisFunction::print_me(std::ostream& s) const {
    s << "atomic basis function: center " << xx << " " << yy << " " << zz << " : ibf " << ibf << " nbf " << nbf << " : shell " << shell << std::endl;
}

/// Print basis info for atoms in the molecule (once for each unique atom type)
void AtomicBasisSet::print(const Molecule& molecule) const {
    molecule.print();
    std::cout << "\n " << name << " atomic basis set" << std::endl;
    for (size_t i=0; i<molecule.natom(); ++i) {
        const Atom& atom = molecule.get_atom(i);
        const unsigned int atn = atom.atomic_number;
        for (size_t j=0; j<i; ++j) {
            if (molecule.get_atom(j).atomic_number == atn)
                goto doneitalready;
        }
        std::cout << std::endl;
        std::cout << "   " <<  get_atomic_data(atn).symbol << std::endl;
        std::cout << ag[atn];
doneitalready:
        ;
    }
}

/// Print basis info for all supported atoms
void AtomicBasisSet::print_all() const {
    std::cout << "\n " << name << " atomic basis set" << std::endl;
    for (unsigned int i=0; i<ag.size(); ++i) {
        if (ag[i].nbf() > 0) {
            std::cout << "   " <<  get_atomic_data(i).symbol << std::endl;
            std::cout << ag[i];
        }
    }
}

void AtomicBasisSet::read_file(std::string filename) {
    static const bool debug = false;
    const char* data_dir = MRA_CHEMDATA_DIR;

    std::string full_filename(data_dir);
    full_filename+="/"+filename;

    // override default location for the basis set
    if (getenv("MRA_CHEMDATA_DIR")) {
    	char* chemdata_dir=getenv("MRA_CHEMDATA_DIR");
        full_filename=std::string(chemdata_dir)+"/"+filename;
    }


    TiXmlDocument doc(full_filename);

    // try to read the AO basis from current directory, otherwise from
    // the environment variable MRA_DATA_DIR
    if (!doc.LoadFile()) {

    	std::cout << "AtomicBasisSet: Failed loading from file " << filename
    			<< " : ErrorDesc  " << doc.ErrorDesc()
    			<< " : Row " << doc.ErrorRow()
    			<< " : Col " << doc.ErrorCol() << std::endl;
    	MADNESS_EXCEPTION("AtomicBasisSet: Failed loading basis set",0);
    }
    for (TiXmlElement* node=doc.FirstChildElement(); node; node=node->NextSiblingElement()) {
        if (strcmp(node->Value(),"name") == 0) {
            name = node->GetText();
            if (debug) std::cout << "Loading basis set " << name << std::endl;
        }
        else if (strcmp(node->Value(), "basis") == 0) {
            const char* symbol = node->Attribute("symbol");
            if (debug) std::cout << "  found basis set for " << symbol << std::endl;
            int atn = symbol_to_atomic_number(symbol);
            std::vector<ContractedGaussianShell> g;
            for (TiXmlElement* shell=node->FirstChildElement(); shell; shell=shell->NextSiblingElement()) {
                const char* type = shell->Attribute("type");
                int nprim=-1;
                shell->Attribute("nprim",&nprim);
                if (debug) std::cout << "      found shell " << type << " " << nprim << std::endl;
                std::vector<double> expnt = load_tixml_vector<double>(shell, nprim, "exponents");
                if (strcmp(type,"L") == 0) {
                    std::vector<double> scoeff = load_tixml_vector<double>(shell, nprim, "scoefficients");
                    std::vector<double> pcoeff = load_tixml_vector<double>(shell, nprim, "pcoefficients");
                    g.push_back(ContractedGaussianShell(0,scoeff,expnt));
                    g.push_back(ContractedGaussianShell(1,pcoeff,expnt));
                }
                else {
                    static const char* tag[] = {"S","P","D","F","G"};
                    int i;
                    for (i=0; i<5; ++i) {
                        if (strcmp(type,tag[i]) == 0) goto foundit;
                    }
                    MADNESS_EXCEPTION("Loading atomic basis set: bad shell type?",0);
foundit:
                    std::vector<double> coeff = load_tixml_vector<double>(shell, nprim, "coefficients");
                    g.push_back(ContractedGaussianShell(i, coeff, expnt));
                }
            }
            ag[atn] = AtomicBasis(g);
        }
        else if (strcmp(node->Value(), "atomicguess") == 0) {
            const char* symbol = node->Attribute("symbol");
            if (debug) std::cout << "  atomic guess info for " << symbol << std::endl;
            int atn = symbol_to_atomic_number(symbol);
            MADNESS_ASSERT(is_supported(atn));
            int nbf = ag[atn].nbf();
            Tensor<double> dmat = load_tixml_matrix<double>(node, nbf, nbf, "guessdensitymatrix");
            Tensor<double> dmatpsp = dmat;
            Tensor<double> aocc = load_tixml_matrix<double>(node, nbf, 1, "alphaocc");
            Tensor<double> bocc = load_tixml_matrix<double>(node, nbf, 1, "betaocc");
	    Tensor<double> aeps = load_tixml_matrix<double>(node, nbf, 1, "alphaeps");
	    Tensor<double> beps = load_tixml_matrix<double>(node, nbf, 1, "betaeps");
            Tensor<double> aoccpsp = aocc;
	    Tensor<double> boccpsp = bocc;
            Tensor<double> avec = load_tixml_matrix<double>(node, nbf, nbf, "alphavectors");
            Tensor<double> bvec = load_tixml_matrix<double>(node, nbf, nbf, "betavectors");
            ag[atn].set_guess_info(dmat, dmatpsp, avec, bvec, aocc, bocc, aeps, beps, aoccpsp, boccpsp);
        }
        else {
            MADNESS_EXCEPTION("Loading atomic basis set: unexpected XML element", 0);
        }
    }

}

void AtomicBasisSet::read_nw_file(std::string filename) {

    // Construct the slymer object that contains interface
    std::ostream bad(nullptr);
    slymer::NWChem_Interface nwchem(filename, bad);

    // Read in the molecule info
    nwchem.read(slymer::Properties::Atoms);

    // Let madness know a basis exists on each atom...
    for(const slymer::Atom &atom : nwchem.atoms) {
        int atn = symbol_to_atomic_number(atom.symbol);
 
        // We need to add to ag[atn] so madness doesn't baulk
        // These functions will not be used in anyway. Need 
        // to add at least the number of orbitals to work.
        // Adding in 2 basis functions per electron, just to 
        // be safe.
        if (ag[atn].nbf() == 0) {
            std::vector<ContractedGaussianShell> g;
            for(int i = 0; i < atn; i++) {
               g.push_back(ContractedGaussianShell(0,{1},{1}));
               g.push_back(ContractedGaussianShell(1,{2},{2}));
            }
            ag[atn] = AtomicBasis(g);
        }
    }
}



void AtomicBasisSet::modify_dmat_psp(int atn, double zeff){
    static const bool debug = false;

    // number of core states to be eliminated
    int zcore = atn - round(zeff);

    Tensor<double> aocc = ag[atn].get_aoccpsp();
    Tensor<double> bocc = ag[atn].get_boccpsp();

    double occ_sum = aocc.sum() + bocc.sum();

    if (debug) std::cout << "before: atn, zeff, occ_sum  " << atn << " " << zeff << "  " << occ_sum << std::endl;
    if (debug) std::cout << "aocc:" << std::endl;
    if (debug) std::cout << aocc << std::endl;
    if (debug) std::cout << "bocc:" << std::endl;
    if (debug) std::cout << bocc << std::endl;

    // return immediately if we already modified this atom type or if there are no core states
    // allow for noise in occupancy sum
    double tol=1e-4;
    if (zcore==0 || (occ_sum < zeff+tol && occ_sum > zeff-tol)) return;

    // otherwise check that total occupancy matches atomic number within tolerance
    if (occ_sum > atn+tol || occ_sum < atn-tol){
        MADNESS_EXCEPTION("Problem with occupancy of initial guess", 0);
    }

    // set occupancies for relevant core states to zero
    // assuming that there is an even number of core states with occupancy 1
    for (int i=0;i<zcore/2;++i){
        aocc[i] = 0.0;
        bocc[i] = 0.0;
    }

    if (debug) std::cout << "after: atn, zeff, occ_sum  " << atn << " " << zeff << "  " << occ_sum << std::endl;
    if (debug) std::cout << "aocc:" << std::endl;
    if (debug) std::cout << aocc << std::endl;
    if (debug) std::cout << "bocc:" << std::endl;
    if (debug) std::cout << bocc << std::endl;

    ag[atn].set_aoccpsp(aocc);
    ag[atn].set_boccpsp(bocc);

    // recalculate the density matrix with new occupancies
    Tensor<double> avec = ag[atn].get_avec();
    Tensor<double> bvec = ag[atn].get_bvec();

    Tensor<double> aovec = transpose(avec);
    Tensor<double> bovec = transpose(bvec);

    // multiply vectors by occupancies
    for (int i=0;i<aocc.size();++i){
        for (int j=0;j<aocc.size();++j){
            aovec(i,j) *= aocc[i];
            bovec(i,j) *= bocc[i];
        }
    }

    Tensor<double> dmata = inner(avec, aovec);
    Tensor<double> dmatb = inner(bvec, bovec);
    Tensor<double> dmat = dmata + dmatb;

    if (debug) std::cout << "dmat:" << std::endl;
    if (debug) std::cout << dmat << std::endl;

    ag[atn].set_dmatpsp(dmat);

}


}
