//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <tinyxml/tinyxml.h>

//using namespace madness;

#include <chem/molecule.h>
#include <chem/molecularbasis.h>
#include <chem/xcfunctional.h>

namespace madness {

typedef Tensor<double> tensorT;

double get_charge_from_file(const std::string filename, unsigned int atype);

/*template <typename Q, int NDIM>
struct function_real2complex_op
{
  typedef std::complex<Q> resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = resultT(*_p0,0.0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};

Function<std::complex<double>,3> function_real2complex(const Function<double,3>& r)
{
  return unary_op_coeffs(r, function_real2complex_op<double,3>());
}*/

class VLocalFunctor : public FunctionFunctorInterface<double,3> {
private:
    double Zeff, zi, C1, C2, C3, C4;
    coord_3d center;
    std::vector<coord_3d> specialpts;

public:
    VLocalFunctor(double Zeff, double zi, 
        double C1, double C2, double C3, double C4, const coord_3d& center)
          : Zeff(Zeff), zi(zi), C1(C1), C2(C2), 
            C3(C3), C4(C4), center(center) {
        specialpts.push_back(center);
    }

    double operator()(const coord_3d& r) const {
        const double x = r[0]-center[0]; const double y = r[1]-center[1]; const double z = r[2]-center[2];
        double rr = std::sqrt(x*x + y*y + z*z);
        double rs = rr/zi;
        double rs2 = rs*rs; double rs4 = rs2*rs2; double rs6 = rs2*rs4;
        return -(Zeff/rr)*std::erf(rs/std::sqrt(2.0)) 
          + std::exp(-0.5*rs2)*
          (C1 + C2*rs2 + C3*rs4 + C4*rs6);
    }
              
    std::vector<coord_3d> special_points() const {return specialpts;}

    Level special_level() {
      return 15;
    }
};

class ProjFunctor : public FunctionFunctorInterface<double,3> {
private:
    double alpha; // radius
    int l, i; // i = 1,2,3 and l = 0,1,2,3 (angular momentum)
    coord_3d center;
    std::vector<coord_3d> specialpts;
    // gamma of half-integers starting and 1 (which means 1/2)
    /*const double gamma_data[17] = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
        945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};*/

    double sqrtPI;
    int itmp, itmp2;
    double t1;
    
public:

    static const double gamma_data[17];
    
    ProjFunctor(double alpha, int l, int i, const coord_3d& center) 
     : alpha(alpha), l(l), i(i), center(center) {
        specialpts.push_back(coord_3d(0.0));
        sqrtPI = std::sqrt(constants::pi);
        itmp = 2*l + (4*i-1);
        itmp2 = l + 2*(i-1);
        t1 = 1./std::pow(alpha, 0.5*(double)itmp)/std::sqrt(gamma_data[itmp-1]*sqrtPI);
    }

    double operator()(const coord_3d& r) const {
        double rsq = (r[0]-center[0])*(r[0]-center[0])+(r[1]-center[1])*(r[1]-center[1])+
            (r[2]-center[2])*(r[2]-center[2]);
        double rr = std::sqrt(rsq);
        double rval = t1;
        if (itmp2 == 0)
            rval *= std::sqrt(2);
        else if (itmp2 == 1)
            rval *= rr*std::sqrt(2);
        else if (itmp2 == 2)
            rval *= rsq*std::sqrt(2);
        else if (itmp2 == 3)
            rval *= rr*rsq*std::sqrt(2);
        else if (itmp2 == 4)
            rval *= rsq*rsq*std::sqrt(2);
        else if (itmp2 == 5)
            rval *= rr*rsq*rsq*std::sqrt(2);
        else if (itmp2 == 6)
            rval *= rsq*rsq*rsq*std::sqrt(2);
        else if (itmp2 == 7)
            rval *= rr*rsq*rsq*rsq*std::sqrt(2);
        rval *= std::exp(-0.5*(rsq/alpha/alpha));
        return rval;
    }
    
    std::vector<coord_3d> special_points() const {return specialpts;}

    Level special_level() {
      return 15;
    }
};

class ProjRLMFunctor : public FunctionFunctorInterface<double,3> {
private:
    double alpha; // radius
    int l, m, i; // i = 1,2,3 and m = 0,1,2,3 and l = 0,1,2,3 (angular momentum)
    coord_3d center;
    std::vector<coord_3d> specialpts;
    // gamma of half-integers starting and 1 (which means 1/2)
    /*const double gamma_data[17] = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
        945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};*/
    double sqrtPI;
    int itmp, itmp2;
    double t1;
    
public:

    static const double gamma_data[17];
    
    ProjRLMFunctor(double alpha, int l, int m, int i, const coord_3d& center) 
     : alpha(alpha), l(l), m(m), i(i), center(center) {
        specialpts.push_back(coord_3d(0.0));
        sqrtPI = std::sqrt(constants::pi);
        itmp = 2*l + (4*i-1);
        itmp2 = 2*(i-1);
        t1 = 1./std::pow(alpha, 0.5*(double)itmp)/std::sqrt(gamma_data[itmp-1]*sqrtPI);
    }

    double operator()(const coord_3d& r) const {
        double x = r[0]-center[0]; double y = r[1]-center[1]; double z = r[2]-center[2];
        double rsq = x*x + y*y + z*z;
        double rr = std::sqrt(rsq);
        double rval = t1;
        const double PI = constants::pi;
        // Radial part
        if (itmp2 == 0) {
            rval *= std::sqrt(2);
        }
        else if (itmp2 == 1) {
            rval *= rr*std::sqrt(2);
        }
        else if (itmp2 == 2) {
            rval *= rsq*std::sqrt(2);
        }
        else if (itmp2 == 3) {
            rval *= rr*rsq*std::sqrt(2);
        }
        else if (itmp2 == 4) {
            rval *= rsq*rsq*std::sqrt(2);
        }
        else if (itmp2 == 5) {
            rval *= rr*rsq*rsq*std::sqrt(2);
        }
        else if (itmp2 == 6) {
            rval *= rsq*rsq*rsq*std::sqrt(2);
        }
        else if (itmp2 == 7) {
            rval *= rr*rsq*rsq*rsq*std::sqrt(2);
        }
        // Angular part
        if (l == 0) {
            rval *= (1./2.)*std::sqrt(1./PI);
        } else if (l == 1) {
            if (m == 0) {
                rval *= std::sqrt(3./4./PI)*x;
            }
            else if (m == 1) {
                rval *= std::sqrt(3./4./PI)*y;
            }
            else if (m == 2) {
                rval *= std::sqrt(3./4./PI)*z;
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 1", 0);
            }
        } else if (l == 2) {
            if (m == 0) {
                rval *= (1./4.)*std::sqrt(5./PI)*(-x*x - y*y + 2*z*z);
            }
            else if (m == 1) {
                rval *= (1./2.)*std::sqrt(15./PI)*(y*z);
            }
            else if (m == 2) {
                rval *= (1./2.)*std::sqrt(15./PI)*(x*z);
            }
            else if (m == 3) {
                rval *= (1./2.)*std::sqrt(15./PI)*(x*y);
            }
            else if (m == 4) {
                rval *= (1./4.)*std::sqrt(15./PI)*(x*x - y*y);
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 2", 0);
            }
        }
        rval *= std::exp(-0.5*(rsq/alpha/alpha));
        return rval;
    }
    
    std::vector<coord_3d> special_points() const {return specialpts;}

    Level special_level() {
      return 15;
    }
};

class ProjStore {
private:
    int maxL;
    vector_real_function_3d projs;

public:
    ProjStore(World& world, const real_tensor& radii, const coord_3d& center) 
     : maxL(radii.dim(0)) {
       double vtol = 1e-1 * FunctionDefaults<3>::get_thresh();
       for (int l = 0; l < maxL; l++) {
           for (int i = 1; i <= 3; i++) {
               //if (world.rank() == 0) print("ProjStore constructor:  ", l, i);
               //projs.push_back(real_factory_3d(world).functor(real_functor_3d(new ProjFunctor(radii(l), l, i, center))).
               //    truncate_mode(0).truncate_on_project());
               real_function_3d f1 = 
                 real_factory_3d(world).functor(real_functor_3d(new ProjFunctor(radii(l), l, i, center))).truncate_mode(0);
               f1.truncate(vtol);
               projs.push_back(f1);
           }
       }
    }
  
    real_function_3d nlproj( int l, int i) {
        // i = 1, 2, 3
        return projs[3*l + i - 1];
    }
};

class ProjRLMStore {
private:
    int maxL;
    real_tensor radii;
    coord_3d center;

public:
    ProjRLMStore(const real_tensor& radii, const coord_3d& center) 
     : maxL(radii.dim(0)), radii(radii), center(center) {}
  
    real_function_3d nlmproj(World& world, int l, int m, int i) {
        real_function_3d f1 = (m < 2*l+1) ?
           real_factory_3d(world).functor(real_functor_3d(new ProjRLMFunctor(radii(l), l, m, i, center))).truncate_mode(0)
             : real_factory_3d(world);
        return f1;
    }
};

class YlmFunctor : public FunctionFunctorInterface<double_complex, 3> {
private:
    int l, m;
    coord_3d center;
    std::vector<coord_3d> specialpts;
    
public:
    YlmFunctor(int l, int m, const coord_3d& center) 
     : l(l), m(m), center(center) {
        specialpts.push_back(center);
    }


    double_complex operator()(const coord_3d& r) const {
        const double x = r[0]-center[0]; const double y = r[1]-center[1]; const double z = r[2]-center[2];
        double rr = std::sqrt(x*x + y*y + z*z);
        double PI = constants::pi;
        if (l == 0) {
            return (1./2.)*std::sqrt(1./PI); 
        } else if (l == 1) {
            if (m == -1) {
                return (1./2.)*std::sqrt(3./2./PI)*double_complex(x, -y)/rr;
            }
            if (m == 0) {
                return (1./2.)*std::sqrt(3./PI)*z/rr;
            }
            if (m == -1) {
                return -(1./2.)*std::sqrt(3./2./PI)*double_complex(x, y)/rr;
            }
        }
    }
    std::vector<coord_3d> special_points() const {return specialpts;}

    Level special_level() {
      return 15;
    }
};

class RlmFunctor : public FunctionFunctorInterface<double, 3> {
private:
    int l, m;
    coord_3d center;
    std::vector<coord_3d> specialpts;

public:
    RlmFunctor(int l, int m, const coord_3d& center) 
     : l(l), m(m), center(center) {
        specialpts.push_back(center);
    }

    double operator()(const coord_3d& r) const {
        const double x = r[0]-center[0]; const double y = r[1]-center[1]; const double z = r[2]-center[2];
        double rr = std::sqrt(x*x + y*y + z*z);
        double PI = constants::pi;
        if (l == 0) {
            return (1./2.)*std::sqrt(1./PI);
        } else if (l == 1) {
            if (m == 0) {
                return std::sqrt(3./4./PI)*x/rr;
            }
            else if (m == 1) {
                return std::sqrt(3./4./PI)*y/rr;
            }
            else if (m == 2) {
                return std::sqrt(3./4./PI)*z/rr;
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 1", 0);
            }
        } else if (l == 2) {
            if (m == 0) {
                return (1./4.)*std::sqrt(5./PI)*(-x*x - y*y + 2*z*z)/rr/rr;
            }
            else if (m == 1) {
                return (1./2.)*std::sqrt(15./PI)*(y*z)/rr/rr;
            }
            else if (m == 2) {
                return (1./2.)*std::sqrt(15./PI)*(x*z)/rr/rr;
            }
            else if (m == 3) {
                return (1./2.)*std::sqrt(15./PI)*(x*y)/rr/rr;
            }
            else if (m == 4) {
                return (1./4.)*std::sqrt(15./PI)*(x*x - y*y)/rr/rr;
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 2", 0);
            }
        }
    }
    std::vector<coord_3d> special_points() const {return specialpts;}

    Level special_level() {
      return 15;
    }
};

class YlmStore {
private:
    int maxL;
    vector_complex_function_3d ylms;

public:
    YlmStore(World& world, int maxL, const coord_3d& center) : maxL(maxL) {
        for (int il = 0; il <= maxL; il++) {
            for (int im = -il; im <= il; im++) {
                ylms.push_back(complex_factory_3d(world).functor(complex_functor_3d(new YlmFunctor(il, im, center))).
                    truncate_mode(0).truncate_on_project());
            }
        }
    }
    // m = -l, -l+1, ..., -1, 0, 1, ... l-1, l 
    complex_function_3d ylm(int l, int m)
    {
        if (l > maxL) MADNESS_EXCEPTION("YlmStore: l out of bounds", 0);
        if (l == 0) {
            return ylms[0];
        }
        else if (l == 1) {
            return ylms[2+m];
        } 
        else if (l == 2) {
            return ylms[6+m];
        } 
        else if (l == 3) {
            return ylms[12+m];
        } 
    }
};

class RlmStore {
private:
    int maxL;
    vector_real_function_3d rlms;

public:
    RlmStore(World& world, int maxL, const coord_3d center) : maxL(maxL) {
       double vtol = 1e-1 * FunctionDefaults<3>::get_thresh();
        for (int il = 0; il <= maxL; il++) {
            for (int im = 0; im < 2*il+1; im++) {
//                rlms.push_back(real_factory_3d(world).functor(real_functor_3d(new RlmFunctor(il, im, center))).
//                    truncate_mode(0).truncate_on_project());
                real_function_3d f1 = 
                    real_factory_3d(world).functor(real_functor_3d(new RlmFunctor(il, im, center))).truncate_mode(0);
                f1.truncate(vtol);
                rlms.push_back(f1);
            }
        }
    }
 
    // m = 0, 1, 2, ..., 2*l+1 (different that Ylm)
    real_function_3d rlm(int l, int m) {
        if (l > maxL) MADNESS_EXCEPTION("RlmStore: l out of bounds", 0);
        if (l == 0) {
            return rlms[0];
        }
        else if (l == 1) {
            return rlms[1+m];
        } 
        else if (l == 2) {
            return rlms[4+m];
        } 
        else if (l == 3) {
            return rlms[9+m];
        } 
    }
};

template <typename Q>
class GTHPseudopotential {
private:
public:
    Molecule molecule;
    std::array<real_tensor,118> localp;
    std::array<real_tensor,118> radii;
    std::array<real_tensor,118> hlij;
    std::array<real_tensor,118> klij;
    real_function_3d vlocalp;
    std::vector<unsigned int> atoms_with_projectors;


public:
    GTHPseudopotential(World& world, Molecule molecule) : molecule(molecule) {}

    void make_pseudo_potential(World& world) {
        // Load info from file
        load_pseudo_from_file(world, "gth.xml");
        atoms_with_projectors.clear();
        
        // fill list with atoms-with-projectors (i.e. not H or He)
        for (int iatom = 0; iatom < molecule.natom(); iatom++) {
            Atom atom = molecule.get_atom(iatom);
            unsigned int atype = atom.atomic_number;
            if (radii[atype-1].dim(0) > 0)
              atoms_with_projectors.push_back(iatom);
        }

        vlocalp = real_factory_3d(world);
        vlocalp.compress();
        for (int iatom = 0; iatom < molecule.natom(); iatom++) {
            // Get atom and it's associated GTH tensors
            Atom atom = molecule.get_atom(iatom);
            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            // do local part
            real_tensor atom_localp = localp[atype-1];
            if (world.rank() == 0) {
               print("Atomic PP parameters: ", atom_localp[0], atom_localp[1], atom_localp[2], atom_localp[3], atom_localp[4], atom_localp[5]);
               print("center at: ", center);
            }
            real_function_3d temp = real_factory_3d(world).functor(
                real_functor_3d(new 
                VLocalFunctor(atom_localp[0], atom_localp[1], atom_localp[2], atom_localp[3], atom_localp[4], atom_localp[5], center))).
                truncate_mode(0).truncate_on_project();
            temp.compress();
            //vlocalp += temp;
            vlocalp.gaxpy(1.0, temp, 1.0, true);
        }
    }

    real_function_3d vlocalpot() {
        return vlocalp;
    }

    void reproject(int k, double thresh) {
      vlocalp = madness::project(vlocalp, k, thresh, true);
    }

    void load_pseudo_from_file(World& world, const std::string filename) {
        bool debug = true;
       
        TiXmlDocument doc(filename);
        if (!doc.LoadFile()) {
            MADNESS_EXCEPTION("Failed to load GTH pseudopotential file", 0);
        }

        for (int iatom = 0; iatom < molecule.natom(); iatom++) {
            Atom atom = molecule.get_atom(iatom);
            unsigned int atype = atom.atomic_number;
            if (world.rank() == 0) {printf("atom atomic_number = %d\n", atype);}
    
            bool success = false;
            for (TiXmlElement* node=doc.FirstChildElement(); node && !success; node=node->NextSiblingElement()) {
                if (strcmp(node->Value(),"name") == 0) {
                    std::string name = node->GetText();
                    if (debug && world.rank() == 0) std::cout << "Loading pseudopotential file " << name << std::endl;
                }
                else if (strcmp(node->Value(), "atom") == 0) {
                    const char* symbol = node->Attribute("symbol");
                    unsigned int atn = symbol_to_atomic_number(symbol);
                    if (atype == atn) {
                        success = true;
                        if (debug && world.rank() == 0) std::cout << "  found atomic pseudopotential " << symbol << std::endl;
                        int lmax = -1;
                        node->Attribute("lmax", &lmax);
                        if (debug && world.rank() == 0) std::cout << "  maximum L is " << lmax << std::endl;
                        real_tensor t_radii((long)lmax+1); 
                        real_tensor t_hlij((long)lmax+1, (long)3, (long)3);
                        real_tensor t_klij((long)lmax+1, (long)3, (long)3);
                        // local part
                        TiXmlElement* xmlVLocal = node->FirstChildElement();
                        real_tensor t_localp((long)6);
                        double zeff = 0.0; xmlVLocal->Attribute("Zeff", &zeff); t_localp[0] = zeff;
                        double lradius = 0.0; xmlVLocal->Attribute("radius", &lradius); t_localp(1) = lradius;
                        double C1 = 0.0; xmlVLocal->Attribute("C1", &C1); t_localp[2] = C1;
                        double C2 = 0.0; xmlVLocal->Attribute("C2", &C2); t_localp[3] = C2;
                        double C3 = 0.0; xmlVLocal->Attribute("C3", &C3); t_localp[4] = C3;
                        double C4 = 0.0; xmlVLocal->Attribute("C4", &C4); t_localp[5] = C4;
                        // loop through nonlocal part
                        for (TiXmlElement* xmlLnlproj = xmlVLocal->NextSiblingElement(); 
                            xmlLnlproj; xmlLnlproj=xmlLnlproj->NextSiblingElement()) {
                            int lvalue = -1; xmlLnlproj->Attribute("l", &lvalue);
                            double radius = 0.0; xmlLnlproj->Attribute("radius", &radius); t_radii[lvalue] = radius; 
                            double h00 = 0.0; xmlLnlproj->Attribute("h00", &h00); t_hlij(lvalue, 0, 0) = h00;
                            double h11 = 0.0; xmlLnlproj->Attribute("h11", &h11); t_hlij(lvalue, 1, 1) = h11;
                            double h22 = 0.0; xmlLnlproj->Attribute("h22", &h22); t_hlij(lvalue, 2, 2) = h22;
                            double k00 = 0.0; xmlLnlproj->Attribute("k00", &k00); t_klij(lvalue, 0, 0) = k00;
                            double k11 = 0.0; xmlLnlproj->Attribute("k11", &k11); t_klij(lvalue, 1, 1) = k11;
                            double k22 = 0.0; xmlLnlproj->Attribute("k22", &k22); t_klij(lvalue, 2, 2) = k22;
                        }
                        // off-diagonal elements
                        if (lmax >= 0) {
                            t_hlij(0, 0, 1) = -1./2.*std::sqrt(3./5.)*t_hlij(0, 1, 1); 
                            t_hlij(0, 1, 0) = t_hlij(0, 0, 1);
                            t_hlij(0, 0, 2) = 1./2.*std::sqrt(5./21.)*t_hlij(0, 2, 2);
                            t_hlij(0, 2, 0) = t_hlij(0, 0, 2);
                            t_hlij(0, 1, 2) = -1./2.*std::sqrt(100./63.)*t_hlij(0, 2, 2);
                            t_hlij(0, 2, 1) = t_hlij(0, 1, 2);
                        } if (lmax >= 1) {
                            t_hlij(1, 0, 1) = 1./2.*std::sqrt(5./7.)*t_hlij(1, 1, 1);
                            t_hlij(1, 1, 0) = t_hlij(1, 0, 1);
                            t_hlij(1, 0, 2) = 1./6.*std::sqrt(35./11.)*t_hlij(1, 2, 2);
                            t_hlij(1, 2, 0) = t_hlij(1, 0, 2);
                            t_hlij(1, 1, 2) = -1./6.*14./std::sqrt(11.)*t_hlij(1, 2, 2);
                            t_hlij(1, 2, 1) = t_hlij(1, 1, 2);
                        } if (lmax >= 2) {
                            t_hlij(2, 0, 1) = -1./2.*std::sqrt(7./9.)*t_hlij(2, 1, 1);
                            t_hlij(2, 1, 0) = t_hlij(2, 0, 1);
                            t_hlij(2, 0, 2) = 1./2.*std::sqrt(63./143.)*t_hlij(2, 2, 2);
                            t_hlij(2, 2, 0) = t_hlij(2, 0, 2);
                            t_hlij(2, 1, 2) = -1./2.*18./std::sqrt(143.)*t_hlij(2, 2, 2);
                            t_hlij(2, 2, 1) = t_hlij(2, 1, 2);
                        }

                        // Copy to main array
                        localp[atype-1] = t_localp;
                        radii[atype-1] = t_radii;
                        hlij[atype-1] = t_hlij;
                        klij[atype-1] = t_klij;
                    }
                }
            }
        }
    }

    std::vector<Function<Q,3> > apply_potential_slow(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi) {
        bool debug = (world.rank() == 0) && false;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential + vlocalp), psi, vtol);
        if (world.rank() == 0) {print("multiplied local");}

        // Non-local part of potential
        unsigned int natoms = atoms_with_projectors.size();
        for (int iatom = 0; iatom < natoms; iatom++) {
            // Get atom and it's associated GTH tensors
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            real_tensor& atom_hlij = hlij[atype-1];
            if (debug && world.rank() == 0) {
              print("pseudo atom:  ", atype);
              print("center:");
              print(center);
              print("radii:");
              print(atom_radii);
              print("hlij:");
              print(atom_hlij);
            }
            // Create function stores for projectors
            int maxL = atom_radii.dim(0)-1;
            RlmStore rlmstore(world, maxL, center);
            ProjStore projstore(world, atom_radii, center);

            int norbs = psi.size();

            ProjRLMStore prlmstore(atom_radii, center);
            
            // NEW (VECTORIZED) ... hopefully
            vector_real_function_3d localproj;
            int lidx = 0;
            Tensor<int> Pilm_lookup((unsigned long) 3, (unsigned long) maxL+1, (unsigned long) 2*maxL+1);
            for (unsigned int j = 1; j <= 3; j++) {
                for (int l = 0; l <= maxL; l++) {
                    for (int m = 0; m < 2*maxL+1; m++) {
                       Pilm_lookup(j-1, l, m) = lidx++;
                       if (m < 2*l+1) localproj.push_back(prlmstore.nlmproj(world,l,m,j)); 
                       else localproj.push_back(real_factory_3d(world));
                    }
                }
            }
            if (world.rank() == 0) {printf("size of localproj: %d\n", localproj.size());}

            Tensor<Q> Pilm = matrix_inner(world, localproj, psi);
            Pilm = Pilm.reshape(3, maxL+1, 2*maxL+1, norbs);

            Tensor<Q> Qilm((unsigned long) 3, (unsigned long) maxL+1, (unsigned long) 2*maxL+1, (unsigned int) norbs);
            for (unsigned int iorb=0; iorb<psi.size(); iorb++) {
                for (unsigned int i = 1; i <= 3; i++) {
                    for (int l = 0; l <= maxL; l++) {
                        for (int m = 0; m < 2*l+1; m++) {
                            Q s = 0.0;
                            for (unsigned int j = 1; j <= 3; j++) {
                                s += atom_hlij(l,i-1,j-1)*Pilm(j-1,l,m,iorb);
                            }
                            Qilm(i-1,l,m,iorb) = s;
                        }
                    }
                }
            }
            Qilm = Qilm.reshape(3*(maxL+1)*(2*maxL+1),norbs);

            double vtol2 = 1e-4*thresh;
            double trantol = vtol2 / std::min(30.0, double(localproj.size()));
            if (world.rank() == 0) {print("localproj size: ", localproj.size(), "Qilm dims: ", Qilm.dim(0), Qilm.dim(1));}
            vector_real_function_3d dpsi = transform(world, localproj, Qilm, trantol, true);
            gaxpy(world, 1.0, vpsi, 1.0, dpsi);    
        }
        return vpsi;
    }

    std::vector<Function<Q,3> > apply_potential(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi, const tensorT & occ, Q & enl) {
        bool debug = (world.rank() == 0) && false;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential), psi, vtol);
        if (world.rank() == 0) {print("multiplied local");}

        unsigned int norbs = psi.size();
        unsigned int natoms = atoms_with_projectors.size();

        // NEW (VECTORIZED) ... hopefully
        vector_real_function_3d localproj;
        int lidx = 0;

        // Hard code maxLL 
        unsigned int maxLL = 0;
        for (unsigned int iatom = 0; iatom < natoms; iatom++) {
            // Get atom and its associated GTH tensors
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            if (atom_radii.dim(0) > 0)
              maxLL = std::max(maxLL,(unsigned int)atom_radii.dim(0)-1);
        }
        if (world.rank() == 0) printf("maxLL = %d\n", maxLL);

        Tensor<int> Pilm_lookup((unsigned int) natoms, (unsigned long) 3, (unsigned long) maxLL+1, (unsigned long) 2*maxLL+1);
        for (unsigned int iatom = 0; iatom < natoms; iatom++) {
            // Get atom and its associated GTH tensors
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            real_tensor& atom_hlij = hlij[atype-1];

            // Create function stores for projectors
            ProjRLMStore prlmstore(atom_radii, center);
            for (unsigned int j = 1; j <= 3; j++) {
                for (unsigned int l = 0; l <= maxLL; l++) {
                    for (unsigned int m = 0; m < 2*maxLL+1; m++) {
                       Pilm_lookup(iatom, j-1, l, m) = lidx++;
                       if (m < 2*l+1) localproj.push_back(prlmstore.nlmproj(world,l,m,j)); 
                       else localproj.push_back(real_factory_3d(world));
                    }
                }
            }
        }
        if (world.rank() == 0) {printf("size of localproj: %d\n", localproj.size());}

        Tensor<Q> Pilm = matrix_inner(world, localproj, psi);
        Pilm = Pilm.reshape(natoms, 3, maxLL+1, 2*maxLL+1, norbs);

        Tensor<Q> Qilm((unsigned int) natoms, (unsigned long) 3, (unsigned long) maxLL+1, (unsigned long) 2*maxLL+1, (unsigned int) norbs);
        for (unsigned int iorb=0; iorb<psi.size(); iorb++) {
            for (unsigned int iatom = 0; iatom < natoms; iatom++) {
                // Get atom and its associated GTH tensors
                Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
                unsigned int atype = atom.atomic_number;
                real_tensor& atom_radii = radii[atype-1];
                real_tensor& atom_hlij = hlij[atype-1];
                int maxL = atom_radii.dim(0)-1;
                for (unsigned int i = 1; i <= 3; i++) {
                    for (int l = 0; l <= maxL; l++) {
                        for (int m = 0; m < 2*l+1; m++) {
                            Q s = 0.0;
                            for (unsigned int j = 1; j <= 3; j++) {
                                s += atom_hlij(l,i-1,j-1)*Pilm(iatom, j-1,l,m,iorb);
                            }
                            Qilm(iatom, i-1,l,m,iorb) = s;
                        }
                    }
                }
            }
        }
        Qilm = Qilm.reshape(natoms*3*(maxLL+1)*(2*maxLL+1),norbs);

        double vtol2 = 1e-4*thresh;
        double trantol = vtol2 / std::min(30.0, double(localproj.size()));
        if (world.rank() == 0) {print("localproj size: ", localproj.size(), "Qilm dims: ", Qilm.dim(0), Qilm.dim(1));}
        vector_real_function_3d dpsi = transform(world, localproj, Qilm, trantol, true);

        // calculate non-local energy
        tensorT nlmat = matrix_inner(world, dpsi, psi, true);
        int nocc = occ.size();
        enl = 0.0;
        for(int i = 0;i < nocc;++i){
            enl += occ[i] * nlmat(i, i);
        }

        //test
        /*tensorT lmat = matrix_inner(world, vpsi, psi, true);
        Q el = 0.0;
        for(int i = 0;i < nocc;++i){
            el += occ[i] * lmat(i, i);
        }

        if(world.rank() == 0){
            printf("\n              enl, el, epot %16.8f  %16.8f  %16.8f\n", enl, el, enl+el);
        }*/

        gaxpy(world, 1.0, vpsi, 1.0, dpsi);  

        return vpsi;
    }


};




}
