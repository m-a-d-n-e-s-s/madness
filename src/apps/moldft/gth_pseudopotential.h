#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include <tinyxml/tinyxml.h>

using namespace madness;

#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>
#include <moldft/xcfunctional.h>

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

double proj_s_1_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double sqrtPI = std::sqrt(constants::pi);
    double sqroot2 = std::sqrt(2.0);
    double r00 = 0.5/sqrtPI;
    double rc = 0.179488;
    double n = sqroot2*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,1.5)*std::sqrt(sqrtPI/2.0);
    return r00*(n/d);
}

double proj_s_2_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double sqrtPI = std::sqrt(constants::pi);
    double sqroot2 = std::sqrt(2.0);
    double r00 = 0.5/sqrtPI;
    double rc = 0.179488;
    double n = sqroot2*rr*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,3.5)*std::sqrt(15*sqrtPI/8.0);
    return r00*(n/d);
}

double proj_p_1_x_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1x = std::sqrt(3.0/4.0/constants::pi)*r[0]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,2.5)*std::sqrt(3.0*sqrtPI/4.0);
    return r1x*(n/d);
}

double proj_p_2_x_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1x = std::sqrt(3.0/4.0/constants::pi)*r[0]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*rr*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,4.5)*std::sqrt(105.0*sqrtPI/16.0);
    return r1x*(n/d);
}

double proj_p_1_y_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1y = std::sqrt(3.0/4.0/constants::pi)*r[1]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,2.5)*std::sqrt(3.0*sqrtPI/4.0);
    return r1y*(n/d);
}

double proj_p_2_y_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1y = std::sqrt(3.0/4.0/constants::pi)*r[1]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*rr*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,4.5)*std::sqrt(105.0*sqrtPI/16.0);
    return r1y*(n/d);
}

double proj_p_1_z_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1z = std::sqrt(3.0/4.0/constants::pi)*r[2]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,2.5)*std::sqrt(3.0*sqrtPI/4.0);
    return r1z*(n/d);
}

double proj_p_2_z_ne(const coord_3d& r) {
    double rr = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    double r1z = std::sqrt(3.0/4.0/constants::pi)*r[2]/rr;
    double rc = 0.214913;
    double sqroot2 = std::sqrt(2.0);
    double sqrtPI = std::sqrt(constants::pi);
    double n = sqroot2*rr*rr*rr*std::exp(-0.5*(rr/rc)*(rr/rc));
    double d = std::pow(rc,4.5)*std::sqrt(105.0*sqrtPI/16.0);
    return r1z*(n/d);
}

class ProjFunctor : public FunctionFunctorInterface<double,3> {
private:
    double alpha; // radius
    int l, i; // i = 1,2,3 and l = 0,1,2,3 (angular momentum)
    coord_3d center;
    std::vector<coord_3d> specialpts;
    double gamma_data[17];
    double sqrtPI;
    int itmp, itmp2;
    double t1;
    
public:
    
    ProjFunctor(double alpha, int l, int i, const coord_3d& center) 
     : alpha(alpha), l(l), i(i), center(center) {
        specialpts.push_back(coord_3d(0.0));
        // gamma of half-integers starting and 1 (which means 1/2)
        // "whole" integers are zero
        gamma_data = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
            945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};
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
    double gamma_data[17];
    double sqrtPI;
    int itmp, itmp2;
    double t1;
    
public:
    
    ProjRLMFunctor(double alpha, int l, int m, int i, const coord_3d& center) 
     : alpha(alpha), l(l), m(m), i(i), center(center) {
        specialpts.push_back(coord_3d(0.0));
        gamma_data = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
            945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};
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

public:
    GTHPseudopotential(World& world, Molecule molecule) : molecule(molecule) {
        load_pseudo_from_file("gth.xml");

        vlocalp = real_factory_3d(world);
        vlocalp.compress();
        for (int iatom = 0; iatom < molecule.natom(); iatom++) {
            // Get atom and it's associated GTH tensors
            Atom atom = molecule.get_atom(iatom);
            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            // do local part
            real_tensor atom_localp = localp[atype-1];
            print("Atomic PP parameters: ", atom_localp[0], atom_localp[1], atom_localp[2], atom_localp[3], atom_localp[4], atom_localp[5]);
            print("center at: ", center);
            real_function_3d temp = real_factory_3d(world).functor(
                real_functor_3d(new 
                VLocalFunctor(atom_localp[0], atom_localp[1], atom_localp[2], atom_localp[3], atom_localp[4], atom_localp[5], center))).
                truncate_mode(0).truncate_on_project();
            temp.compress();
            //vlocalp += temp;
            vlocalp.gaxpy(1.0, temp, 1.0, true);
        }
    }

    void reproject(int k, double thresh) {
      vlocalp = madness::project(vlocalp, k, thresh, true);
    }

    void load_pseudo_from_file(const std::string filename) {
        bool debug = true;
       
        TiXmlDocument doc(filename);
        if (!doc.LoadFile()) {
            MADNESS_EXCEPTION("Failed to load GTH pseudopotential file", 0);
        }

        for (int iatom = 0; iatom < molecule.natom(); iatom++) {
            Atom atom = molecule.get_atom(iatom);
            unsigned int atype = atom.atomic_number;
            printf("atom atomic_number = %d\n", atype);
    
            bool success = false;
            for (TiXmlElement* node=doc.FirstChildElement(); node && !success; node=node->NextSiblingElement()) {
                if (strcmp(node->Value(),"name") == 0) {
                    std::string name = node->GetText();
                    if (debug) std::cout << "Loading pseudopotential file " << name << std::endl;
                }
                else if (strcmp(node->Value(), "atom") == 0) {
                    const char* symbol = node->Attribute("symbol");
                    int atn = symbol_to_atomic_number(symbol);
                    if (atype == atn) {
                        success = true;
                        if (debug) std::cout << "  found atomic pseudopotential " << symbol << std::endl;
                        int lmax = -1;
                        node->Attribute("lmax", &lmax);
                        if (debug) std::cout << "  maximum L is " << lmax << std::endl;
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

    std::vector<Function<Q,3> > apply_potential(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi) {
        bool debug = (world.rank() == 0) && false;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential + vlocalp), psi, vtol);
        print("multiplied local");
        for (unsigned int iorb=0; iorb<psi.size(); iorb++) {
            // Non-local part of potential
            for (int iatom = 0; iatom < molecule.natom(); iatom++) {
                // Get atom and it's associated GTH tensors
                Atom atom = molecule.get_atom(iatom);
                coord_3d center = atom.get_coords();
                unsigned int atype = atom.atomic_number;
                real_tensor& atom_radii = radii[atype-1];
                real_tensor& atom_hlij = hlij[atype-1];
                if (debug) {
                  print("orbital:    ", iorb);
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

                for (unsigned int i = 1; i <= 3; i++) {
                    if (debug) printf("\n");
                    for (unsigned int j = 1; j <= 3; j++) {
                        //unsigned int j = i;
                        for (int l = 0; l <= maxL; l++) {
                            for (int m = 0; m < 2*l+1; m++) {
                                //if (atom_hlij(l, i-1, j-1) > vtol) {
                                if (true) {
                                    if (debug) printf("atom:  %d     l:  %d    m:  %d    i:  %d    j:  %d    %15.8e  %15.8e  %15.8e\n", iatom, l, m, i, j, atom_hlij(l,i-1,j-1), inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]),atom_hlij(l,i-1,j-1)*inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]));
                                    //double_complex t1 = inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    Q t1 = inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    if (std::abs(t1) > thresh) vpsi[iorb] += atom_hlij(l,i-1,j-1)*inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb])*rlmstore.rlm(l,m)*projstore.nlproj(l,i);
                                }
                            }
                        }
                    }
                }
            }
        }
        return vpsi;
    }

    std::vector<Function<Q,3> > apply_potential2(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi) {
        bool debug = (world.rank() == 0) && true;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential + vlocalp), psi, vtol);
        print("multiplied local");
        for (unsigned int iorb=0; iorb<psi.size(); iorb++) {
            // Non-local part of potential
            for (int iatom = 0; iatom < molecule.natom(); iatom++) {
                // Get atom and it's associated GTH tensors
                Atom atom = molecule.get_atom(iatom);
                coord_3d center = atom.get_coords();
                unsigned int atype = atom.atomic_number;
                real_tensor& atom_radii = radii[atype-1];
                real_tensor& atom_hlij = hlij[atype-1];
                if (debug) {
                  print("orbital:    ", iorb);
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

                for (unsigned int i = 1; i <= 3; i++) {
                    if (debug) printf("\n");
                    //for (unsigned int j = 1; j <= 3; j++) {
                    {
                        unsigned int j = i;
                        for (int l = 0; l <= maxL; l++) {
                            for (int m = 0; m < 2*l+1; m++) {
                                //if (atom_hlij(l, i-1, j-1) > vtol) {
                                if (true) {
                                    if (debug) printf("atom:  %d     l:  %d    m:  %d    i:  %d    j:  %d    %15.8e    %15.8e\n", iatom, l, m, i, j, atom_hlij(l,i-1,j-1), atom_hlij(l,i-1,j-1)*std::abs(inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb])));
                                    //double_complex t1 = inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    Q t1 = inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    if (std::abs(t1) > thresh) vpsi[iorb] += atom_hlij(l,i-1,j-1)*inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb])*rlmstore.rlm(l,m)*projstore.nlproj(l,i);
                                }
                            }
                        }
                    }
                }
            }
        }
        return vpsi;
    }

    std::vector<Function<Q,3> > apply_potential_ne(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi) {
        bool debug = (world.rank() == 0) && false;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential + vlocalp), psi, vtol);
        print("multiplied local ne");

        real_function_3d prjs1 = real_factory_3d(world).f(proj_s_1_ne);
        real_function_3d prjs2 = real_factory_3d(world).f(proj_s_2_ne);
        real_function_3d prjpx1 = real_factory_3d(world).f(proj_p_1_x_ne);
        real_function_3d prjpx2 = real_factory_3d(world).f(proj_p_2_x_ne);
        real_function_3d prjpy1 = real_factory_3d(world).f(proj_p_1_y_ne);
        real_function_3d prjpy2 = real_factory_3d(world).f(proj_p_2_y_ne);
        real_function_3d prjpz1 = real_factory_3d(world).f(proj_p_1_z_ne);
        real_function_3d prjpz2 = real_factory_3d(world).f(proj_p_2_z_ne);

        Atom atom = molecule.get_atom(0);
        coord_3d center = atom.get_coords();
        unsigned int atype = atom.atomic_number;
        real_tensor& atom_radii = radii[atype-1];
        real_tensor& atom_hlij = hlij[atype-1];

        for (unsigned int iorb=0; iorb < psi.size(); iorb++) {
            vpsi[iorb] += atom_hlij(0,0,0)*inner(prjs1,psi[iorb])*prjs1;
            vpsi[iorb] += atom_hlij(0,1,1)*inner(prjs2,psi[iorb])*prjs2;
            vpsi[iorb] += atom_hlij(1,0,0)*inner(prjpx1,psi[iorb])*prjpx1;
            vpsi[iorb] += atom_hlij(1,1,1)*inner(prjpx2,psi[iorb])*prjpx2;
            vpsi[iorb] += atom_hlij(1,0,0)*inner(prjpy1,psi[iorb])*prjpy1;
            vpsi[iorb] += atom_hlij(1,1,1)*inner(prjpy2,psi[iorb])*prjpy2;
            vpsi[iorb] += atom_hlij(1,0,0)*inner(prjpz1,psi[iorb])*prjpz1;
            vpsi[iorb] += atom_hlij(1,1,1)*inner(prjpz2,psi[iorb])*prjpz2;
        }

        for (unsigned int iorb=0; iorb < psi.size(); iorb++) {
            //double t1 = atom_hlij(0,0,0)*inner(prjs1,psi[iorb]);
            //double t2 = atom_hlij(0,1,1)*inner(prjs2,psi[iorb]);
            //double t3 = atom_hlij(1,0,0)*inner(prjpx1,psi[iorb]);
            //double t4 = atom_hlij(1,1,1)*inner(prjpx2,psi[iorb]);
            //double t5 = atom_hlij(1,0,0)*inner(prjpy1,psi[iorb]);
            //double t6 = atom_hlij(1,1,1)*inner(prjpy2,psi[iorb]);
            //double t7 = atom_hlij(1,0,0)*inner(prjpz1,psi[iorb]);
            //double t8 = atom_hlij(1,1,1)*inner(prjpz2,psi[iorb]);

            double t1 = inner(prjs1,psi[iorb]);
            double t2 = inner(prjs2,psi[iorb]);
            double t3 = inner(prjpx1,psi[iorb]);
            double t4 = inner(prjpx2,psi[iorb]);
            double t5 = inner(prjpy1,psi[iorb]);
            double t6 = inner(prjpy2,psi[iorb]);
            double t7 = inner(prjpz1,psi[iorb]);
            double t8 = inner(prjpz2,psi[iorb]);

            printf("atom_hlij(0,0,0)*inner(prjs1,psi[iorb]) = %15.8f   %15.8f  %15.8f\n",  atom_hlij(0,0,0), t1, atom_hlij(0,0,0)*t1);
            printf("atom_hlij(0,1,1)*inner(prjs2,psi[iorb]) = %15.8f   %15.8f  %15.8f\n",  atom_hlij(0,1,1), t2, atom_hlij(0,1,1)*t2);
            printf("atom_hlij(1,0,0)*inner(prjpx1,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,0,0), t3, atom_hlij(1,0,0)*t3);
            printf("atom_hlij(1,1,1)*inner(prjpx2,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,1,1), t4, atom_hlij(1,1,1)*t4);
            printf("atom_hlij(1,0,0)*inner(prjpy1,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,0,0), t5, atom_hlij(1,0,0)*t5);
            printf("atom_hlij(1,1,1)*inner(prjpy2,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,1,1), t6, atom_hlij(1,1,1)*t6);
            printf("atom_hlij(1,0,0)*inner(prjpz1,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,0,0), t7, atom_hlij(1,0,0)*t7);
            printf("atom_hlij(1,1,1)*inner(prjpz2,psi[iorb]) = %15.8f   %15.8f  %15.8f\n", atom_hlij(1,1,1), t8, atom_hlij(1,1,1)*t8);
            printf("\n\n");
        }
        return vpsi;
    }

    vector_real_function_3d apply_potential_wsttiger(World& world, const real_function_3d& potential, const vector_complex_function_3d& psi) {
        bool debug = (world.rank() == 0) && false;
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        vector_real_function_3d rpsi(psi.size());
        for (unsigned int i = 0; i < psi.size(); i++) rpsi[i] = real(psi[i]);
        vector_real_function_3d vpsi = real(mul_sparse(world,(potential + vlocalp), rpsi, vtol));
        print("multiplied local");
        for (unsigned int iorb=0; iorb<psi.size(); iorb++) {
            // Non-local part of potential
            for (int iatom = 0; iatom < molecule.natom(); iatom++) {
                // Get atom and it's associated GTH tensors
                Atom atom = molecule.get_atom(iatom);
                coord_3d center = atom.get_coords();
                unsigned int atype = atom.atomic_number;
                real_tensor& atom_radii = radii[atype-1];
                real_tensor& atom_hlij = hlij[atype-1];
                if (debug) {
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
                YlmStore ylmstore(world, maxL, center);
                ProjStore projstore(world, atom_radii, center);

                for (unsigned int i = 1; i <= 3; i++) {
                    if (debug) printf("\n");
                    // WSTHORNTON
                    //for (unsigned int j = 1; j <= 3; j++) {
                    {
                      unsigned int j = i;
                        for (int l = 0; l <= maxL; l++) {
                            for (int m = 0; m < 2*l+1; m++) {
                                if (atom_hlij(l, i-1, j-1) > vtol) {
                                //if (true) {
                                    //if (debug) printf("atom:  %d     l:  %d    i:  %d    j:  %d    %15.8e    %15.8e\n", iatom, l, i, j, atom_hlij(l,i-1,j-1), std::abs(inner(conj(ylmstore.ylm(l,m))*projstore.nlproj(l,j), psi[iorb])));
                                    //double_complex t1 = inner(conj(rlmstore.rlm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    double_complex t1 = inner(conj(ylmstore.ylm(l,m))*projstore.nlproj(l,j), psi[iorb]);
                                    if (std::abs(t1) > thresh) vpsi[iorb] += real(atom_hlij(l,i-1,j-1)*inner(conj(ylmstore.ylm(l,m))*projstore.nlproj(l,j), psi[iorb])*ylmstore.ylm(l,m)*projstore.nlproj(l,i));
                                }
                            }
                        }
                    }
                }
            }
        }
        return vpsi;
    }
};
