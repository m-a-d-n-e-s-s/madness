/// \file gth_pseudopotential.cc
/// \brief GTH pseudopotential functionality
/// \defgroup moldft The molecular density functional and Hartree-Fock code

#ifndef MADNESS_CHEM_GTH_PSEUDOPOTENTIAL_H__INCLUDED
#define MADNESS_CHEM_GTH_PSEUDOPOTENTIAL_H__INCLUDED

#include <madness/mra/mra.h>
#include <madness/external/tinyxml/tinyxml.h>

#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/xcfunctional.h>

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
        return -(Zeff/rr)*erf(rs/std::sqrt(2.0))
          + std::exp(-0.5*rs2)*
          (C1 + C2*rs2 + C3*rs4 + C4*rs6);
    }
              
    std::vector<coord_3d> special_points() const override final {return specialpts;}

    Level special_level() const override final {
      return 6;
    }
};

class ProjRLMFunctor : public FunctionFunctorInterface<double,3> {
private:
    double alpha; // radius
    int l, m; // i = 1,2,3 and m = 0,1,2,3 and l = 0,1,2,3 (angular momentum) (i just used in constructor)
    coord_3d center;
    std::vector<coord_3d> specialpts;
    // gamma of half-integers starting and 1 (which means 1/2)
    /*const double gamma_data[17] = {1.0, 0.0, 1.0/2.0, 0.0, 3.0/4.0, 0.0, 15.0/8.0, 0.0, 105.0/16.0, 0.0, 
        945.0/32.0, 0.0, 10395.0/64.0, 0.0, 135135.0/128.0, 0.0, 2027025.0/256.0};*/
    double sqrtPI;
    int itmp, itmp2;
    double t1;
    
public:

    virtual bool supports_vectorized() const {return false;}

    static const double gamma_data[17];
    
    ProjRLMFunctor(double alpha, int l, int m, int i, const coord_3d& center) 
     : alpha(alpha), l(l), m(m), center(center) {
        specialpts.push_back(coord_3d(0.0));
        sqrtPI = std::sqrt(constants::pi);
        itmp = 2*l + (4*i-1);
        itmp2 = 2*(i-1);
        t1 = 1./std::pow(alpha, 0.5*(double)itmp)/std::sqrt(gamma_data[itmp-1]*sqrtPI);
    }

    double operator()(const coord_3d& r) const {
        double x = r[0]-center[0]; double y = r[1]-center[1]; double z = r[2]-center[2];
        double rsq = x*x + y*y + z*z;

        if (rsq > 40.0) return 0.0;

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
   
    virtual bool screened(const coord_3d& c1, const coord_3d& c2) const {
        double ftol = 1e-12;

        double x1 = c1[0]; double y1 = c1[1]; double z1 = c1[2];
        double x2 = c2[0]; double y2 = c2[1]; double z2 = c2[2];

        // if center is inside box, then return false
        // otherwise, look for the closest point and check
        bool inside = (center[0] >= x1) && (center[0] <= x2) &&
                      (center[1] >= y1) && (center[1] <= y2) &&
                      (center[2] >= z1) && (center[2] <= z2);
        if (inside) {
          return false;
        }
        else {
            //printf("GTH_pseudopotential: (point not inside)\n");
            //print("  c1: ", c1, "     c2: ", c2);
            double minr = 1e10;
            int ii = -1; int jj = -1; int kk = -1;
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    for (int k = 0; k <= 1; k++) {
                        double x = (i == 0) ? c1[0] : c2[0];
                        double y = (j == 0) ? c1[1] : c2[1];
                        double z = (k == 0) ? c1[2] : c2[2];
                        coord_3d rr = coord_3d{x, y, z} - center;
                        double rsq = rr[0]*rr[0]+rr[1]*rr[1]+rr[2]*rr[2];
             //           print("current minr:  ", minr, "     point: ", {x,y,z}, "     center: ", center, "     p: ", {x,y,z}, "     rsq: ", rsq);
                        if (minr > rsq) {
                            minr = rsq;
                            ii = i; jj = j; kk = k;
                        }
                    }
                }
            }
            //print("ii: ", ii, "jj: ", jj, "kk: ", kk);
            //printf("\n");
            if ((ii < 0) || (jj < 0) || (kk < 0)) MADNESS_EXCEPTION("GTH_Pseudopotential: failed to find suitable minimum point\n", 0);
            double x = (ii == 0) ? c1[0] : c2[0]; 
            double y = (jj == 0) ? c1[1] : c2[1]; 
            double z = (kk == 0) ? c1[2] : c2[2]; 
            double fval = this->operator()({x, y, z});
            if (fabs(fval) < ftol) return true;
            else return false;
        }
    }
 
    virtual void operator()(const Vector<double*,3>& xvals, double* MADNESS_RESTRICT fvals, int npts) const {
        
        double* x = new double[npts];
        double* y = new double[npts];
        double* z = new double[npts];
        double* rsq = new double[npts];
        double* rr = new double[npts];

        double* x1 = xvals[0]; double* x2 = xvals[1]; double* x3 = xvals[2];
        for (int i = 0; i < npts; i++) {
            x[i] = x1[i]-center[0];
            y[i] = x2[i]-center[1];
            z[i] = x3[i]-center[2];
            rsq[i] = x[i]*x[i] + y[i]*y[i] + z[i]*z[i];
            rr[i] = std::sqrt(rsq[i]);
            fvals[i] = t1;
        }

        const double PI = constants::pi;

        // Radial part
        if (itmp2 == 0) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= std::sqrt(2);
            }
        }
        else if (itmp2 == 1) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rr[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 2) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rsq[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 3) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rr[i]*rsq[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 4) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rsq[i]*rsq[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 5) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rr[i]*rsq[i]*rsq[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 6) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rsq[i]*rsq[i]*rsq[i]*std::sqrt(2);
            }
        }
        else if (itmp2 == 7) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= rr[i]*rsq[i]*rsq[i]*rsq[i];
            }
        }
        // Angular part
        if (l == 0) {
            for (int i = 0; i < npts; i++) {
                fvals[i] *= (1./2.)*std::sqrt(1./PI)*std::exp(-0.5*(rsq[i]/alpha/alpha));
            }
        } else if (l == 1) {
            if (m == 0) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= std::sqrt(3./4./PI)*x[i]*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 1) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= std::sqrt(3./4./PI)*y[i]*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 2) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= std::sqrt(3./4./PI)*z[i]*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 1", 0);
            }
        } else if (l == 2) {
            if (m == 0) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= (1./4.)*std::sqrt(5./PI)*(-x[i]*x[i] - y[i]*y[i] + 2*z[i]*z[i])*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 1) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= (1./2.)*std::sqrt(15./PI)*(y[i]*z[i])*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 2) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= (1./2.)*std::sqrt(15./PI)*(x[i]*z[i])*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 3) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= (1./2.)*std::sqrt(15./PI)*(x[i]*y[i])*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else if (m == 4) {
                for (int i = 0; i < npts; i++) {
                    fvals[i] *= (1./4.)*std::sqrt(15./PI)*(x[i]*x[i] - y[i]*y[i])*std::exp(-0.5*(rsq[i]/alpha/alpha));
                }
            }
            else {
                MADNESS_EXCEPTION("m out of range for l = 2", 0);
            }
        }
        //for (int i = 0; i < npts; i++) {
        //    fvals[i] *= std::exp(-0.5*(rsq[i]/alpha/alpha));
        //}

        delete [] x;
        delete [] y;
        delete [] z;
        delete [] rsq;
        delete [] rr;
    }

    std::vector<coord_3d> special_points() const override final {return specialpts;}

    Level special_level() const override final {
      return 6;
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
      // real_function_3d f1 = (m < 2*l+1) ?
      // 	real_factory_3d(world).functor(real_functor_3d(new ProjRLMFunctor(radii(l), l, m, i, center))).
      // 	truncate_on_project().nofence().truncate_mode(0) : real_factory_3d(world);

      real_function_3d f1;
      if (m < 2*l+1) {
	auto functor = real_functor_3d(new ProjRLMFunctor(radii(l), l, m, i, center));
	f1 = real_factory_3d(world).functor(functor).truncate_on_project().nofence().truncate_mode(0);
      }
      else {
	f1 = real_factory_3d(world);
      }
	
      return f1;
    }

    ProjRLMFunctor nlmproj_functor(World& world, int l, int m, int i) {
       return ProjRLMFunctor(radii(l), l, m, i, center);
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
        for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
            Atom atom = molecule.get_atom(iatom);

            //make sure this is actually a pseudo-atom
            if (!atom.pseudo_atom) continue;

            unsigned int atype = atom.atomic_number;
            if (radii[atype-1].dim(0) > 0)
              atoms_with_projectors.push_back(iatom);
        }

        vlocalp = real_factory_3d(world);
        vlocalp.compress();
        for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
            // Get atom and its associated GTH tensors
            Atom atom = molecule.get_atom(iatom);

            //make sure this is actually a pseudo-atom
            if (!atom.pseudo_atom) continue;

            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            // do local part
            real_tensor atom_localp = localp[atype-1];
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

        for (size_t iatom = 0; iatom < molecule.natom(); iatom++) {
            Atom atom = molecule.get_atom(iatom);
            unsigned int atype = atom.atomic_number;
            if (debug && world.rank() == 0) {printf("atom atomic_number = %d\n", atype);}
    
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
                            t_hlij(1, 0, 1) = -1./2.*std::sqrt(5./7.)*t_hlij(1, 1, 1);
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


    std::vector<Function<Q,3> > apply_potential(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi, const tensorT & occ, Q & enl) {
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential), psi, vtol);

        unsigned int norbs = psi.size();
        unsigned int natoms = atoms_with_projectors.size();

        // NEW (VECTORIZED) ... hopefully
        vector_real_function_3d localproj;
        int lidx = 0;
        unsigned int maxLL = 0;
        // loop through all of the atom types in the molecule and get the maximum L value
        // we need this because we are going to create a fixed sized tensor to store 
        // mapping between a linear index and (ias=which atom, i=which projector, l=angular momentum, 
        // m=angular momentum projection)
        for (unsigned int iatom = 0; iatom < natoms; iatom++) {
            // Get atom and its associated GTH tensors
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            if (atom_radii.dim(0) > 0)
              maxLL = std::max(maxLL,(unsigned int)atom_radii.dim(0)-1);
        }

        // Pilm_lookup is a mapping between a linear index and (ias,i,l,m)
        Tensor<int> Pilm_lookup((unsigned int) natoms, (unsigned long) 3, (unsigned long) maxLL+1, (unsigned long) 2*maxLL+1);
        for (unsigned int iatom = 0; iatom < natoms; iatom++) {
            // Get atom and its associated GTH tensors
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            coord_3d center = atom.get_coords();
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            //real_tensor& atom_hlij = hlij[atype-1];

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
            // Somehow this scares me ... thread-safety of the container localproj (???)
            world.gop.fence();
        }
        truncate(world, localproj, FunctionDefaults<3>::get_thresh());
        compress(world, localproj);
        //truncate(world, psi, FunctionDefaults<3>::get_thresh());
        compress(world, psi);
        //truncate(world, vpsi, FunctionDefaults<3>::get_thresh());
        compress(world, vpsi);

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
        vector_real_function_3d dpsi = transform(world, localproj, Qilm, trantol, true);

        // calculate non-local energy
        tensorT nlmat = matrix_inner(world, dpsi, psi, true);
        int nocc = occ.size();
        enl = 0.0;
        for(int i = 0;i < nocc;++i){
            enl += occ[i] * nlmat(i, i);
        }

        //debug printing
        /*tensorT lmat = matrix_inner(world, vpsi, psi, true);
        Q el = 0.0;
        for(int i = 0;i < nocc;++i){
            el += occ[i] * lmat(i, i);
            std::cout << "nloc/loc " << i << "  " << occ[i] << "  " << nlmat(i,i) << "  "<< lmat(i,i) << std::endl;
        }

        if(world.rank() == 0){
            printf("\n              enl, el, epot %16.8f  %16.8f  %16.8f\n", enl, el, enl+el);
        }*/

        gaxpy(world, 1.0, vpsi, 1.0, dpsi);

        return vpsi;
    }

    std::vector<Function<Q,3> > apply_potential_simple(World& world, const real_function_3d& potential, const std::vector<Function<Q,3> >& psi, const tensorT & occ, Q & enl) {
        double thresh = FunctionDefaults<3>::get_thresh();
        double vtol = 1e-2*thresh;
        std::vector<Function<Q,3> > vpsi = mul_sparse(world,(potential), psi, vtol);

        //unsigned int norbs = psi.size();
        unsigned int natoms = atoms_with_projectors.size();

        for (unsigned int iatom = 0; iatom < natoms; iatom++) {
            Atom atom = molecule.get_atom(atoms_with_projectors[iatom]);
            unsigned int atype = atom.atomic_number;
            real_tensor& atom_radii = radii[atype-1];
            real_tensor& atom_hlij = hlij[atype-1];
            coord_3d center = atom.get_coords();
            ProjRLMStore prlmstore(atom_radii, center);
            unsigned int maxLL = atom_radii.dim(0)-1;
            for (unsigned int l = 0; l <= maxLL; l++) {
                for (unsigned int m = 0; m < 2*l+1; m++) {
                    for (unsigned int i = 1; i <= 3; i++) {
                        real_function_3d fproji = prlmstore.nlmproj(world,l,m,i);
                        for (unsigned int j = 1; j <= 3; j++) {
                            real_function_3d fprojj = prlmstore.nlmproj(world,l,m,j);
                            for (unsigned int iorb = 0; iorb < psi.size(); iorb++) {
                               double val = atom_hlij(l,i-1,j-1)*(fprojj.inner(psi[iorb]));
                               if (std::abs(val) > vtol*1e-2) {
                                   vpsi[iorb] += val*fproji;
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




}

#endif // MADNESS_CHEM_GTH_PSEUDOPOTENTIAL_H__INCLUDED
