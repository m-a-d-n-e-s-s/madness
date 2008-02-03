#ifndef MOLECULAR_BASIS_H
#define MOLECULAR_BASIS_H

#include <examples/molecule.h>
#include <examples/constants.h>
#include <tinyxml/tinyxml.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>

/// Represents a single shell of contracted, Cartesian, Gaussian primitives
class ContractedGaussianShell {
    int type;  //< Angular momentum = 0, 1, 2, ...
    std::vector<double> coeff;
    std::vector<double> expnt;
    double rsqmax;

    void normalize() {
        // nwcchem cartesian normalization conventions
        // translation of nmcoeff.F into python and thence to c++
        int np = coeff.size();
        if (np == 1) coeff[0] = 1.0e0;

        double pi32=pow(madness::constants::pi,1.5);
        int l_lim = 2*type - 1;
        double f = 1.0e00;
        for (int n=l_lim; n>1; n-=2) f *= n;
        
        for (int n=0; n<np; n++)
            coeff[n] *= pow(2.e0*expnt[n]/madness::constants::pi,0.75e0)*pow(4.e0*expnt[n],0.5E0*type)/sqrt(f);
                
        double sum = 0.0;
        for (int n1=0; n1<np; n1++) {
            for (int n2=0; n2<np; n2++) {
                double S =pi32/pow(expnt[n1]+expnt[n2],1.5e0+type)/pow(2e0,type);
                sum = sum + coeff[n1]*coeff[n2]*S;
            }
        }
        sum *= f;

        f = 1e0/sqrt(sum);
        for (int n=0; n<np; n++) coeff[n] *= f;
    }
    
public:
    ContractedGaussianShell() 
        : type(-1), coeff(), expnt(), rsqmax(0.0)
    {};

    ContractedGaussianShell(int type, 
                            const std::vector<double>& coeff, 
                            const std::vector<double>& expnt, 
                            bool donorm=true)
        : type(type), coeff(coeff), expnt(expnt) 
    {
        if (donorm) normalize();
        double minexpnt = expnt[0];
        for (unsigned int i=1; i<expnt.size(); i++) 
            minexpnt = std::min(minexpnt,expnt[i]);
        rsqmax = 18.4/minexpnt;  // 18.4 = 8*ln(10)
    }

    /// Returns square of the distance beyond which function is less than 1e-8.
    double rangesq() const {return rsqmax;}

    /// Evaluates the radial part of the contracted function
    double operator()(double rsq) const {
        double sum = 0.0;
        for (unsigned int i=0; i<coeff.size(); i++) {
            double ersq = expnt[i]*rsq;
            if (ersq < 18.4) sum += coeff[i]*exp(-ersq);
        }
        return sum;
    }

    /// Returns the shell angular momentum
    int angular_momentum() const {return type;}

    /// Returns the number of basis functions in the shell
    int nbf() const {return (type+1)*(type+2)/2;}

    /// Returns the number of primitives in the contraction
    int nprim() const {return coeff.size();}

    /// Returns a const reference to the coefficients
    const std::vector<double>& get_coeff() const {return coeff;}

    /// Returns a const reference to the exponents
    const std::vector<double>& get_expnt() const {return expnt;}

    template <typename Archive>
    void serialize(Archive& ar) {ar & type & coeff & expnt & rsqmax;}
};


std::ostream& operator<<(std::ostream& s, const ContractedGaussianShell& c) {
    static const char* tag[] = {"s","p","d","f","g"};
    char buf[32768];
    char* p = buf;
    const std::vector<double>& coeff = c.get_coeff();
    const std::vector<double>& expnt = c.get_expnt();

    p += sprintf(p,"%s [",tag[c.angular_momentum()]);
    for (int i=0; i<c.nprim(); i++){
        p += sprintf(p, "%.6f(%.6f)",coeff[i],expnt[i]);
        if (i != (c.nprim()-1)) p += sprintf(p, ", "); 
    }
    p += sprintf(p, "]");
    s << buf;
    return s;
}


/// Represents multiple shells of contracted gaussians on a single center
class AtomicBasis {
    std::vector<ContractedGaussianShell> g;
    double rmaxsq;
    int numbf;

public:
    AtomicBasis() : g(), rmaxsq(0.0), numbf(0) {};

    AtomicBasis(const std::vector<ContractedGaussianShell>& g) 
        : g(g) 
    {
        rmaxsq = 0.0; 
        numbf = 0;
        for (unsigned int i=0; i<g.size(); i++) {
            rmaxsq = std::max(rmaxsq, g[i].rangesq());
            numbf += g[i].nbf();
        }
    }
    
    /// Returns the number of basis functions on the center
    int nbf() const {return numbf;}

    /// Returns the number of shells on the center
    int nshell() const {return g.size();}

    /// Returns a const reference to the shells
    const std::vector<ContractedGaussianShell>& get_shells() const {return g;};

    /// Evaluates the basis functions at point x, y, z

    /// The array bf[] must be large enough to hold nbf() values.
    ///
    /// Returned is the incremented pointer.
    double* operator()(double x, double y, double z, double* bf) const {
        double rsq = x*x + y*y + z*z;
        if (rsq > rmaxsq) {
            for (int i=0; i<numbf; i++) bf[i] = 0.0;
            return bf+numbf;
        }

        double* bfstart;
        for (unsigned int i=0; i<g.size(); i++) {
            if (rsq < g[i].rangesq()) {
                double R = g[i](rsq);
                switch (g[i].angular_momentum()) {
                case 0:
                    bf[0] =  R;
                    break;
                case 1:
                    bf[0] =  R*x;
                    bf[1] =  R*y;
                    bf[2] =  R*z;
                    break;
                case 2:
                    bf[0] = R*x*x;
                    bf[1] = R*x*y;
                    bf[2] = R*x*z ;
                    bf[3] = R*y*y ;
                    bf[4] = R*y*z ;
                    bf[5] = R*z*z;
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
                }
            }
            bf += g[i].nbf();
        }
        // paranoia is good
        MADNESS_ASSERT(bf-bfstart == numbf);
        return bf;
    }

    template <typename Archive>
    void serialize(Archive& ar) {ar & g & rmaxsq & numbf;}

};

std::ostream& operator<<(std::ostream& s, const AtomicBasis& c) {
    const std::vector<ContractedGaussianShell>& shells = c.get_shells();
    for (int i=0; i<c.nshell(); i++) {
        s << "   " << shells[i] << std::endl;
    }
    return s;
}


/// Contracted Gaussian basis 
class AtomicBasisSet {
    std::string name;
    std::vector<AtomicBasis> ag;  // Basis associated by atomic number = 1, 2, ...; 0=Bq.

    template <typename T>
    std::vector<T> load_tixml_vector(TiXmlElement* node, int n, const char* name) {
        TiXmlElement* child = node->FirstChildElement(name);
        MADNESS_ASSERT(child);
        std::istringstream s(child->GetText());
        std::vector<T> r(n);
        for (int i=0; i<n; i++) {
            MADNESS_ASSERT(s >> r[i]);
        }
        return r;
    }

public:
    AtomicBasisSet() : name("unknown"), ag(110) {}


    AtomicBasisSet(std::string filename) : name(""), ag(110) {
        load_from_file(filename);
    }

    void load_from_file(std::string filename) {
        static const bool debug = false;
        TiXmlDocument doc(filename);
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
                        for (i=0; i<5; i++) {
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
            else {
                MADNESS_EXCEPTION("Loading atomic basis set: unexpected XML element", 0);
            }
        }
    }

    /// Given a molecule count the number of basis functions
    int nbf(const Molecule& molecule) const {
        int n = 0;
        for (int i=0; i<molecule.natom(); i++) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_EXCEPTION("AtomicBasisSet: unsupported atom?",atn);
            n += ag[atn].nbf();
        }
        return n;
    }

    void operator()(const Molecule& molecule, double x, double y, double z, double *bf) const {
        for (int i=0; i<molecule.natom(); i++) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            bf = ag[atn](x, y, z, bf);
        }
    }

    bool is_supported(int atomic_number) {
        return ag[atomic_number].nbf() > 0;
    }

    /// Print basis info for atoms in the molecule (once for each unique atom type)
    void print(const Molecule& molecule) const {
        for (int i=0; i<molecule.natom(); i++) {
            const Atom& atom = molecule.get_atom(i);
            const unsigned int atn = atom.atomic_number;
            for (int j=0; j<i; j++) {
                if (molecule.get_atom(j).atomic_number == atn) 
                    goto doneitalready;
            }
            std::cout << std::endl;
            std::cout << "Basis functions on " << get_atomic_data(atn).symbol << std::endl;
            std::cout << ag[atn];
        doneitalready:
            ;
        }
    }

    /// Print basis info for all supported atoms
    void print_all() const {
        for (unsigned int i=0; i<ag.size(); i++) {
            if (ag[i].nbf() > 0) {
                std::cout << std::endl;
                std::cout << "Basis functions on " << get_atomic_data(i).symbol << std::endl;
                std::cout << ag[i];
            }
        }
    }

    template <typename Archive>
    void serialize(Archive& ar) {ar & name & ag;}

        
    


};



#endif
