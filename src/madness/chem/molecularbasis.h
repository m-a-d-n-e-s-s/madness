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

#ifndef MADNESS_CHEM_MOLECULAR_BASIS_H__INCLUDED
#define MADNESS_CHEM_MOLECULAR_BASIS_H__INCLUDED

#include <madness/madness_config.h>
#include <madness/constants.h>
#include<madness/chem/molecule.h>
#include<madness/chem/atomutil.h>
#include <madness/external/tinyxml/tinyxml.h>
#include <madness/tensor/tensor.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>

namespace madness {
/// Represents a single shell of contracted, Cartesian, Gaussian primitives
class ContractedGaussianShell {
public:
    static const int maxtype=6; ///< Maximum angular momentum supported
    static const int maxbf=(maxtype+1)*(maxtype+2)/2; ///< Maximum number of basis functions in a shell
private:
    int type;  ///< Angular momentum = 0, 1, 2, ...
    std::vector<double> coeff;
    std::vector<double> expnt;
    double rsqmax;
    int numbf;  ///< Number of basis functions in shell (type+1)*(type+2)/2

    void normalize() {
        // nwcchem cartesian normalization conventions
        // translation of nmcoeff.F into python and thence to c++
        int np = coeff.size();
        if (np == 1) coeff[0] = 1.0e0;

        double pi32=pow(madness::constants::pi,1.5);
        int l_lim = 2*type - 1;
        double f = 1.0e00;
        for (int n=l_lim; n>1; n-=2) f *= n;

        for (int n=0; n<np; ++n)
            coeff[n] *= pow(2.e0*expnt[n]/madness::constants::pi,0.75e0)*pow(4.e0*expnt[n],0.5E0*type)/sqrt(f);

        double sum = 0.0;
        for (int n1=0; n1<np; ++n1) {
            for (int n2=0; n2<np; ++n2) {
                double S =pi32/pow(expnt[n1]+expnt[n2],1.5e0+type)/pow(2e0,type);
                sum = sum + coeff[n1]*coeff[n2]*S;
            }
        }
        sum *= f;

        f = 1e0/sqrt(sum);
        for (int n=0; n<np; ++n) coeff[n] *= f;
    }

public:
    ContractedGaussianShell()
            : type(-1), coeff(), expnt(), rsqmax(0.0), numbf(0) {};

    ContractedGaussianShell(int type,
                            const std::vector<double>& coeff,
                            const std::vector<double>& expnt,
                            bool donorm=true)
            : type(type), coeff(coeff), expnt(expnt), numbf((type+1)*(type+2)/2) {
        if (donorm) normalize();
        double minexpnt = expnt[0];
        for (unsigned int i=1; i<expnt.size(); ++i)
            minexpnt = std::min(minexpnt,expnt[i]);
        rsqmax = 27.6/minexpnt;  // 27.6 = log(1e12)
    }


    /// Returns square of the distance beyond which function is less than 1e-8.
    double rangesq() const {
        return rsqmax;
    }


    /// Evaluates the radial part of the contracted function
    double eval_radial(double rsq) const {
        if (rsq > rsqmax) return 0.0;
        double sum = 0.0;
        for (unsigned int i=0; i<coeff.size(); ++i) {
            double ersq = expnt[i]*rsq;
            if (ersq < 27.6) sum += coeff[i]*exp(-ersq); // 27.6 = log(1e12)
        }
        return sum;
    }


    /// Evaluates the entire shell returning the incremented result pointer
    double* eval(double rsq, double x, double y, double z, double* bf) const {
        double R = eval_radial(rsq);
        if (fabs(R) < 1e-12) {
            for (int i=0; i<numbf; ++i) bf[i] = 0.0;

        }
        else {
            switch (type) {
            case 0:
                bf[0] =  R;
                break;
            case 1:
                bf[0] =  R*x;
                bf[1] =  R*y;
                bf[2] =  R*z;
                break;
            case 2:
              { // braces need by some compilers to limit scope of fac
                static const double fac = 1.0; //sqrt(3.0);
                bf[0] = R*x*x;
                bf[1] = R*x*y*fac;
                bf[2] = R*x*z*fac;
                bf[3] = R*y*y;
                bf[4] = R*y*z*fac;
                bf[5] = R*z*z;
              }
                break;
            case 3:
                bf[0] = R*x*x*x;
                bf[1] = R*x*x*y;
                bf[2] = R*x*x*z;
                bf[3] = R*x*y*y;
                bf[4] = R*x*y*z;
                bf[5] = R*x*z*z;
                bf[6] = R*y*y*y;
                bf[7] = R*y*y*z;
                bf[8] = R*y*z*z;
                bf[9] = R*z*z*z;
                break;

            default:
                throw "UNKNOWN ANGULAR MOMENTUM";
            }
        }
        return bf+numbf;
    }


    /// Returns the shell angular momentum
    int angular_momentum() const {
        return type;
    }

    /// Returns the number of basis functions in the shell
    int nbf() const {
        return numbf;
    }

    /// Returns the number of primitives in the contraction
    int nprim() const {
        return coeff.size();
    }

    /// Returns a const reference to the coefficients
    const std::vector<double>& get_coeff() const {
        return coeff;
    }

    /// Returns a const reference to the exponents
    const std::vector<double>& get_expnt() const {
        return expnt;
    }

    /// Returns a string description of the basis function type
    const char* get_desc(int ibf) const {
        static const char* tags[4][10] = {
            {"s"   ,""    ,""    ,""    ,""    ,""    ,""    ,""    ,""    ,""    } ,
            {"px"  ,"py"  ,"pz"  ,""    ,""    ,""    ,""    ,""    ,""    ,""    } ,
            {"dxx" ,"dxy" ,"dxz" ,"dyy" ,"dyz" ,"dzz" ,""    ,""    ,""    ,""    } ,
            {"fxxx","fxxy","fxxz","fxyy","fxyz","fxzz","fxzz","fyyy","fyzz","fzzz"}
        };
        MADNESS_ASSERT(ibf<numbf && ibf >= 0);
        return tags[type][ibf];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & type & coeff & expnt & rsqmax & numbf;
    }
};

/// Represents multiple shells of contracted gaussians on a single center
class AtomicBasis {
    std::vector<ContractedGaussianShell> g;
    double rmaxsq;
    int numbf;
    Tensor<double> dmat, dmatpsp, avec, bvec, aocc, bocc, aeps, beps, aoccpsp, boccpsp;

public:
    AtomicBasis() : g(), rmaxsq(0.0), numbf(0) {};

    AtomicBasis(const std::vector<ContractedGaussianShell>& g)
            : g(g) {
        rmaxsq = 0.0;
        numbf = 0;
        for (unsigned int i=0; i<g.size(); ++i) {
            rmaxsq = std::max(rmaxsq, g[i].rangesq());
            numbf += g[i].nbf();
        }
    }

    void set_guess_info(const Tensor<double>& dmat, const Tensor<double>& dmatpsp,
                        const Tensor<double>& avec, const Tensor<double>& bvec,
                        const Tensor<double>& aocc, const Tensor<double>& bocc,
                        const Tensor<double>& aeps, const Tensor<double>& beps,
                        const Tensor<double>& aoccpsp, const Tensor<double>& boccpsp) {
        this->dmat = copy(dmat);
        this->dmatpsp = copy(dmatpsp);
        this->avec = copy(avec);
        this->bvec = copy(bvec);
        this->aocc = copy(aocc);
        this->beps = copy(beps);
        this->aeps = copy(aeps);
        this->bocc = copy(bocc);
        this->aoccpsp = copy(aoccpsp);
        this->boccpsp = copy(boccpsp);
    }

    /// Returns the number of basis functions on the center
    int nbf() const {
        return numbf;
    }

    /// Returns the number of shells on the center
    int nshell() const {
        return g.size();
    }

    /// Returns a const reference to the shells
    const std::vector<ContractedGaussianShell>& get_shells() const {
        return g;
    };

    /// Evaluates the basis functions at point x, y, z relative to atomic center

    /// The array bf[] must be large enough to hold nbf() values.
    ///
    /// Returned is the incremented pointer.
    double* eval(double x, double y, double z, double* bf) const {
        double rsq = x*x + y*y + z*z;
        if (rsq > rmaxsq) {
            for (int i=0; i<numbf; ++i) bf[i] = 0.0;
            return bf+numbf;
        }

        double* bfstart = bf;
        for (unsigned int i=0; i<g.size(); ++i) {
            bf = g[i].eval(rsq, x, y, z, bf);
        }
        // paranoia is good
        MADNESS_CHECK(bf-bfstart == numbf);
        return bf;
    }

    /// Evaluates the guess atomic density at point x, y, z relative to atomic center
    double eval_guess_density(double x, double y, double z, bool pspat) const {
        MADNESS_ASSERT(has_guess_info());
        double rsq = x*x + y*y + z*z;
        if (rsq > rmaxsq) return 0.0;

        double bf[ContractedGaussianShell::maxbf];
        eval(x, y, z, bf);
        const double* p;
        // check if pseudo-atom
        if (pspat){
            p = dmatpsp.ptr();}
        else{
            p = dmat.ptr();}
        double sum = 0.0;
        for (int i=0; i<numbf; ++i, p+=numbf) {
            double sumj = 0.0;
            for (int j=0; j<numbf; ++j)
                sumj += p[j]*bf[j];
            sum += bf[i]*sumj;
        }
        return sum;
    }

    /// Return shell that contains basis function ibf and also return index of function in the shell
    const ContractedGaussianShell& get_shell_from_basis_function(int ibf, int& ibf_in_shell) const {
        int n=0;
        for (unsigned int i=0; i<g.size(); ++i) {
            int nbf_in_shell = g[i].nbf();
            if (ibf>=n && ibf<(n+nbf_in_shell)) {
                ibf_in_shell = ibf-n;
                return g[i];
            }
            else {
                n += g[i].nbf();
            }
        }
        MADNESS_EXCEPTION("AtomicBasis: get_shell_from_basis_function", ibf*100000 + nbf());
    }

    bool has_guess_info() const {
        return dmat.size()>0;
    }

    const Tensor<double>& get_dmat() const {
        return dmat;
    };

    void set_dmat(Tensor<double>& mat) {
       dmat = mat;
    };

    bool has_guesspsp_info() const {
        return dmatpsp.size()>0;
    }

    const Tensor<double>& get_dmatpsp() const {
        return dmatpsp;
    };

    void set_dmatpsp(Tensor<double>& mat) {
       dmatpsp = mat;
    };

    const Tensor<double>& get_avec() const {
        return avec;
    };

    const Tensor<double>& get_bvec() const {
        return bvec;
    };

    const Tensor<double>& get_aeps() const {
        return aeps;
    };

    const Tensor<double>& get_beps() const {
        return beps;
    };

    const Tensor<double>& get_aocc() const {
        return aocc;
    };

    const Tensor<double>& get_bocc() const {
        return bocc;
    };

    void set_aocc(Tensor<double>& occ) {
       aocc = occ;
    };

    void set_bocc(Tensor<double>& occ)  {
       bocc = occ;
    };

    const Tensor<double>& get_aoccpsp() const {
        return aoccpsp;
    };

    const Tensor<double>& get_boccpsp() const {
        return boccpsp;
    };

    void set_aoccpsp(Tensor<double>& occ) {
       aoccpsp = occ;
    };

    void set_boccpsp(Tensor<double>& occ)  {
       boccpsp = occ;
    };

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & g & rmaxsq & numbf & dmat & dmatpsp & avec & bvec & aocc & bocc & aeps & beps & aoccpsp & boccpsp;
    }

};

/// Used to represent one basis function from a shell on a specific center
class AtomicBasisFunction {
private:
    const double xx, yy, zz; // Coordinates of the center
    const ContractedGaussianShell& shell; // Reference to the underlying atomic shell
    const int ibf; // Index of basis function in the shell (0, 1, ...)
    const int nbf; // Number of functions in the shell

public:
    AtomicBasisFunction(double x, double y, double z,
                        const ContractedGaussianShell& shell, int ibf)
            : xx(x), yy(y), zz(z), shell(shell), ibf(ibf), nbf(shell.maxbf) {}


    AtomicBasisFunction(const AtomicBasisFunction& aofunc)
            : xx(aofunc.xx)
            , yy(aofunc.yy)
            , zz(aofunc.zz)
            , shell(aofunc.shell)
            , ibf(aofunc.ibf)
            , nbf(aofunc.nbf) {}

    double operator()(double x, double y, double z) const {
        double bf[ContractedGaussianShell::maxbf];
        x-=xx;
        y-=yy;
        z-=zz;
        double rsq = x*x + y*y + z*z;
        shell.eval(rsq, x, y, z, bf);
        return bf[ibf];
    }

    void print_me(std::ostream& s) const;

    const ContractedGaussianShell& get_shell() const {
        return shell;
    }

    int get_index() const {
        return ibf;
    }

    const char* get_desc() const {
        return shell.get_desc(ibf);
    }

    void get_coords(double& x, double& y, double& z) const {
    	x=xx; y=yy; z=zz;
        return;
    }

    madness::Vector<double,3> get_coords_vec() const {
        return madness::Vector<double,3>{xx, yy, zz};
    }

    double rangesq() const {
        return shell.rangesq();
    }
};

/// Contracted Gaussian basis
class AtomicBasisSet {
    std::string name;
    std::vector<AtomicBasis> ag;  ///< Basis associated by atomic number = 1, 2, ...; 0=Bq.

    template <typename T>
    std::vector<T> load_tixml_vector(TiXmlElement* node, int n, const char* name) {
        TiXmlElement* child = node->FirstChildElement(name);
        MADNESS_ASSERT(child);
        std::istringstream s(child->GetText());
        std::vector<T> r(n);
        for (int i=0; i<n; ++i) {
            s >> r[i];
            MADNESS_ASSERT(s);
        }
        return r;
    }

    template <typename T>
    Tensor<T> load_tixml_matrix(TiXmlElement* node, int n, int m, const char* name) {
        TiXmlElement* child = node->FirstChildElement(name);
        MADNESS_ASSERT(child);
        std::istringstream s(child->GetText());
        Tensor<T> r(n,m);
        for (int i=0; i<n; ++i) {
            for (int j=0; j<m; ++j) {
                s >> r(i,j);
                MADNESS_ASSERT(s);
            }
        }
        return r;
    }

public:
    AtomicBasisSet() : name("unknown"), ag(110) {}


    AtomicBasisSet(std::string filename) : name(""), ag(110) {
        read_file(filename);
    }

    std::string get_name() const {return name;}

    /// read the atomic basis set from file

    /// use the default location MRA_CHEMDATA_DIR as defined in the Makefile.am
    /// unless it is overridden by the environment variable MRA_CHEMDATA_DIR
    /// @param[in]	filename	the name of the basis set (sto-3g, 6-31g, etc)
    void read_file(std::string filename);

    /// read the atomic basis set from file

    /// @param[in]	filename	the base name of nwchem files (.out and .movecs)
    void read_nw_file(std::string filename);



    /// Makes map from atoms to first basis function on atom and number of basis functions on atom
    void atoms_to_bfn(const Molecule& molecule, std::vector<int>& at_to_bf, std::vector<int>& at_nbf) const {
        at_to_bf = std::vector<int>(molecule.natom());
        at_nbf   = std::vector<int>(molecule.natom());

        int n = 0;
        for (size_t i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            at_to_bf[i] = n;
            at_nbf[i] = ag[atn].nbf();
            n += at_nbf[i];
        }
    }

    /// Makes map from shells to first basis function on she and number of basis functions on sh
    void shells_to_bfn(const Molecule& molecule, std::vector<int>& sh_to_bf, std::vector<int>& sh_nbf) const {
        sh_to_bf = std::vector<int>();
        sh_nbf   = std::vector<int>();

        int nbf = 0;
        for (size_t i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            const auto& shells = ag[atn].get_shells();
            for (const auto& sh : shells) {
                int n = sh.nbf();
                sh_nbf.push_back(n);
                sh_to_bf.push_back(nbf);
                nbf += n;
            }
        }
    }

    /// Returns the atomic alpha eigenvectors for atom iat
    const Tensor<double>& get_avec(const Molecule& molecule, size_t iat) const {
      MADNESS_ASSERT(iat>=0 && iat<molecule.natom());
      const Atom& atom = molecule.get_atom(iat);
      const int atn = atom.atomic_number;
      MADNESS_ASSERT(is_supported(atn));
      return ag[atn].get_avec();
    }

    /// Returns the atomic alpha eigenvalues for atom iat
    const Tensor<double>& get_aeps(const Molecule& molecule, size_t iat) const {
      MADNESS_ASSERT(iat>=0 && iat<molecule.natom());
      const Atom& atom = molecule.get_atom(iat);
      const int atn = atom.atomic_number;
      MADNESS_ASSERT(is_supported(atn));
      return ag[atn].get_aeps();
    }

    /// Returns the number of the atom the ibf'th basis function is on
    int basisfn_to_atom(const Molecule& molecule, size_t ibf) const {
        MADNESS_ASSERT(ibf >= 0);
        size_t n = 0;
        for (size_t i=0; i<molecule.natom(); ++i) {
            // Is the desired function on this atom?
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            const int nbf_on_atom = ag[atn].nbf();
            if (ibf >= n  && (n+nbf_on_atom) > ibf) {
                return i;
            }
            else {
                n += nbf_on_atom;
            }
        }
        MADNESS_EXCEPTION("AtomicBasisSet: get_atomic_basis_function: confused?", ibf);
    }

    /// Returns the ibf'th atomic basis function
    AtomicBasisFunction get_atomic_basis_function(const Molecule& molecule, size_t ibf) const {
        MADNESS_ASSERT(ibf >= 0);
        size_t n = 0;
        for (size_t i=0; i<molecule.natom(); ++i) {
            // Is the desired function on this atom?
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            const int nbf_on_atom = ag[atn].nbf();
            if (ibf >= n  && (n+nbf_on_atom) > ibf) {
                int index;
                const ContractedGaussianShell& shell =
                    ag[atn].get_shell_from_basis_function(ibf-n, index);
                return AtomicBasisFunction(atom.x, atom.y, atom.z, shell, index);
            }
            else {
                n += nbf_on_atom;
            }
        }
        MADNESS_EXCEPTION("AtomicBasisSet: get_atomic_basis_function: confused?", ibf);
    }


    /// Given a molecule count the number of basis functions
    int nbf(const Molecule& molecule) const {
        int n = 0;
        for (size_t i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            n += ag[atn].nbf();
        }
        return n;
    }

    /// Evaluates the basis functions
    void eval(const Molecule& molecule, double x, double y, double z, double *bf) const {
        for (size_t i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            bf = ag[atn].eval(x-atom.x, y-atom.y, z-atom.z, bf);
        }
    }


    /// Evaluates the guess density
    double eval_guess_density(const Molecule& molecule, double x, double y, double z) const {
        double sum = 0.0;
        bool pspat;
        for (size_t i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            if (atom.pseudo_atom){
                pspat=true;}
            else{
                pspat=false;}
            const int atn = atom.atomic_number;
            sum += ag[atn].eval_guess_density(x-atom.x, y-atom.y, z-atom.z, pspat);
        }
        return sum;
    }

    bool is_supported(int atomic_number) const {
        return ag[atomic_number].nbf() > 0;
    }

    /// Print basis info for atoms in the molecule (once for each unique atom type)
    void print(const Molecule& molecule) const;

    /// Eliminates core orbitals from the density matrix for pseudopotential calculations
    void modify_dmat_psp(int atn, double zeff);

    template <typename T>
    class AnalysisSorter {
        const Tensor<T> v;
    public:
        AnalysisSorter(const Tensor<T>& v) : v(v) {}
        bool operator()(long i, long j) const {
            return std::abs(v[i]) > std::abs(v[j]);
        }
    };

    /// Given a vector of AO coefficients prints an analysis

    /// For each significant coeff it prints
    /// - atomic symbol
    /// - atom number
    /// - basis function type (e.g., dxy)
    /// - basis function number
    /// - MO coeff
    template <typename T>
    void print_anal(const Molecule& molecule, const Tensor<T>& v) const {
        const double thresh = 0.2*v.normf();
        if (thresh == 0.0) {
            printf("    zero vector\n");
            return;
        }
        long nbf = int(v.dim(0));
        std::vector<long> list(nbf);
        long ngot=0;
        for (long i=0; i<nbf; ++i) {
            if (std::abs(v(i)) > thresh) {
                list[ngot++] = i;
            }
        }
        std::sort(list.begin(),list.begin()+ngot,AnalysisSorter<T>(v));

        const char* format;
        if (molecule.natom() < 10) {
            format = "  %2s(%1d)%4s(%2ld)%6.3f  ";
        }
        else if (molecule.natom() < 100) {
            format = "  %2s(%2d)%4s(%3ld)%6.3f  ";
        }
        else if (molecule.natom() < 1000) {
            format = "  %2s(%3d)%4s(%4ld)%6.3f  ";
        }
        else {
            format = "  %2s(%4d)%4s(%5ld)%6.3f  ";
        }
        printf("         ");
        for (long ii=0; ii<ngot; ++ii) {
            long ibf = list[ii];

            const int iat = basisfn_to_atom(molecule, ibf);
            const Atom& atom = molecule.get_atom(iat);
            const AtomicBasisFunction ao = get_atomic_basis_function(molecule, ibf);
            const char* desc = ao.get_desc();
            const char* element = get_atomic_data(atom.atomic_number).symbol;

            // This will need wrapping in a template for a complex MO vector
            printf(format, element, iat, desc, ibf, v[ibf]);
        }
        printf("\n");
    }

    /// Print basis info for all supported atoms
    void print_all() const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & name & ag;
    }
};
}


#endif
