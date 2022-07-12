//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <madness/mra/mra.h>
#include "mentity.h"

using namespace madness;

static double L = 6.5;

typedef Vector<double,3> coordT;
typedef Function<double,3> rfunctionT;
typedef FunctionFactory<double,3> rfactoryT;
typedef std::vector<rfunctionT> rvecfuncT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > rfunctorT;

//*************************************************************************
double compute_volume()
{
  return L*L*L;
}
//*************************************************************************

//***************************************************************************
struct vectorLengthFunctor : public std::binary_function<Vector<double,3>, Vector<double,3>, bool>
{
    bool operator()( Vector<double,3> lhs, Vector<double,3> rhs)
    {
      double llen = sqrt(lhs[0]*lhs[0] + lhs[1]*lhs[1] + lhs[2]*lhs[2]);
      double rlen = sqrt(rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);
      return (llen < rlen);
    }
};
//***************************************************************************

//*************************************************************************
class GaussianFunctor : public FunctionFunctorInterface<double,3> {
private:
  double coeff;
  double expnt;
  std::vector<coordT> specialpts;

public:
  GaussianFunctor(double coeff, double expnt)
        : coeff(coeff), expnt(expnt)
  {
    specialpts.push_back({0.0,0.0,0.0});
  }

    virtual std::vector<coordT> special_points() const
    {
      return specialpts;
    }

    virtual Level special_level()
    {
      return 10;
    }

    double operator()(const coord_3d& r) const {
        double x = r[0]; double y = r[1]; double z = r[2];
        //return coeff*std::exp(-expnt*(x*x + y*y + z*z));
        return coeff*std::exp(-expnt*x*x)*std::exp(-expnt*y*y)*std::exp(-expnt*z*z);
    }
};
//*************************************************************************
void gen_ce(double mu, double xlo, double eps, Tensor<double>& c, Tensor<double>& e)
{
  double T = 25.0;

  double slo = 0.0;
  if (mu > 1e-12)
  {
    slo = -0.5*std::log(4*T/mu/mu);
    double xhi = 1.0;
    slo = std::log10(eps/xhi)-1.0;
  }
  double shi = 0.5*std::log(T/xlo/xlo);
  slo = -28.63102;
  double h = 1.0/(0.2-0.47*std::log10(eps));
  unsigned int npt = std::floor(std::abs((shi-slo)/h));
  h = (shi-slo)/(npt+1);
  c = Tensor<double>((int)npt);
  e = Tensor<double>((int)npt);

  for (unsigned int i = 0; i < npt; i++)
  {
    double s = slo + h*(i+1);
    c[i] = (2*h/std::sqrt(constants::pi)) * (std::exp(0.25*mu*mu*std::exp(-2.0*s)))*std::exp(s);
    e[i] = std::exp(2.0*s);
  }

  print("npt:  ", npt, "  slo:  ", slo, "  shi:  ", shi);
  for (unsigned int i = 0; i < c.dim(0); i++)
  {
    print("i:  ", i, "  c:  ", c[i], "  e:  ", e[i]);
  }
}
//*************************************************************************

//*************************************************************************

//*************************************************************************
std::vector< Vector<double,3> > generate_R_vectors(World& world, double maxRlen = 50.0)
{
  const double t1 = L;

  std::vector< Vector<double,3> > rvecs;
  std::vector< Vector<double,3> > rvecs2;
  double maxRlen20 = 1.2*maxRlen;

  int rlo = -500;
  int rhi = 501;
  for (int ir1 = rlo; ir1 < rhi; ir1++)
  {
    for (int ir2 = rlo; ir2 < rhi; ir2++)
    {
      for (int ir3 = rlo; ir3 < rhi; ir3++)
      {
        double rlen = t1*std::sqrt(ir1*ir1 + ir2*ir2 + ir3*ir3);
        if (rlen <= maxRlen)
        {
          Vector<double,3> rvec {t1*ir1, t1*ir2, t1*ir3};
          rvecs.push_back(rvec);
        }
        else if (rlen < maxRlen20)
        {
          Vector<double,3> rvec {t1*ir1, t1*ir2, t1*ir3};
          rvecs2.push_back(rvec);
        }
      }
    }
  }
  std::sort(rvecs.begin(), rvecs.end(), vectorLengthFunctor());
  std::sort(rvecs2.begin(), rvecs2.end(), vectorLengthFunctor());

//      if (_world.rank() == 0) printf("Size of vectors:  %d\n", rvecs.size());
//      if (_world.rank() == 0)
//        printf("\nR-vectors:\n");

//  for (unsigned int ir = 0; ir < rvecs.size(); ir++)
//  {
//    Vector<double,3> rvec = rvecs[ir];
//    double rlen = std::sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);
//
//        if (world.rank() == 0)
//          printf("%10.5f %10.5f %10.5f   %10.5f\n",rvec[0],rvec[1],rvec[2],rlen);
//  }

  // Get size info
  unsigned int rsize = rvecs.size();
  Vector<double,3> rvec = rvecs[rsize-1];
  double maxRlen2 = std::sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);
  unsigned int rsize2 = rvecs2.size();
  if (world.rank() == 0) print("rsize:  ", rsize, "rsize2:  ", rsize2);
  Vector<double,3> rvec2 = rvecs2[0];
  double maxRlen20_2 = std::sqrt(rvec2[0]*rvec2[0] + rvec2[1]*rvec2[1] + rvec2[2]*rvec2[2]);

  if (world.rank() == 0)
    print("R-max length requested:  ", maxRlen, "    R-max length:  ", maxRlen2,
        "    R-max length (120%):  ", maxRlen20_2);

  return rvecs;
}
//*************************************************************************

//*************************************************************************
// generate G-vectors for a SCC in terms of the real space lattice vectors
std::vector< Vector<double,3> > generate_G_vectors(World& world, double maxGlen = 15.0)
{
  const double TWO_PI = 2*constants::pi;
  const double t1 = TWO_PI/L;

  std::vector< Vector<double,3> > gvecs;

  int glo = -500;
  int ghi = 501;
  for (int ig1 = glo; ig1 < ghi; ig1++)
  {
    for (int ig2 = glo; ig2 < ghi; ig2++)
    {
      for (int ig3 = glo; ig3 < ghi; ig3++)
      {
        double glen = t1*std::sqrt(ig1*ig1 + ig2*ig2 + ig3*ig3);
        if (glen <= maxGlen)
        {
          Vector<double,3> gvec {t1*ig1, t1*ig2, t1*ig3};
          gvecs.push_back(gvec);
        }
      }
    }
  }
  std::sort(gvecs.begin(), gvecs.end(), vectorLengthFunctor());
//      if (_world.rank() == 0)
//        printf("\nG-vectors:\n");
  for (unsigned int ig = 0; ig < gvecs.size(); ig++)
  {
    Vector<double,3> gvec = gvecs[ig];
    double glen = std::sqrt(gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]);
//        if (_world.rank() == 0)
//          printf("%10.5f %10.5f %10.5f   %10.5f\n",gvec[0],gvec[1],gvec[2],glen);
  }

  unsigned int gsize = gvecs.size();
  Vector<double,3> gvec = gvecs[gsize-1];
  double maxGlen2 = std::sqrt(gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]);
  if (world.rank() == 0) print("G-max length requested:  ", maxGlen, "    G-max length:  ", maxGlen2);
  return gvecs;

}
//*************************************************************************

//*************************************************************************
class OnesFunctor : public FunctionFunctorInterface<double,3> {
private:
  double val;

public:
  OnesFunctor(double val=1.0) : val(val) {}

  virtual double operator()(const coordT& r) const
  {
    return val;
  }
};
//*************************************************************************


//*************************************************************************
class EwaldNuclearPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
  std::vector< Vector<double,3> > rvecs;
  std::vector< Vector<double,3> > gvecs;
  std::vector<double_complex> gsfactor;
  double alpha;
  MolecularEntity* mentity;

public:
  EwaldNuclearPotentialFunctor(World& world, MolecularEntity* mentity, double alpha)
   : alpha(alpha), mentity(mentity)
  {
    if (world.rank() == 0) print("EwaldNuclearPotentialFunctor [started]");

    double gthresh = 1e-6;
    double gmax = 550.0;
    double gmin = 0.01;
    bool gstop = false;
    double gval2 = 0.5*(gmax+gmin);
    for (unsigned int i = 0; i < 500 && !gstop; i++)
    {
      gval2 = 0.5*(gmax+gmin);
      double t1 = std::exp(-gval2/4/alpha/alpha)/gval2;
      if (((t1 > 0.9e-10) && (t1 < 1.1e-10)) || (gmax-gmin/2) < gthresh) gstop = true;
      if (!gstop)
      {
        if (t1 > 1e-10)
        {
          gmin = gval2;
        }
        else
        {
          gmax = gval2;
        }
      }
      if (world.rank() == 0) print("i: ", i, "  gval: ", std::sqrt(gval2), "  t1: ", t1);
    }

    double rthresh = 1e-6;
    double rmax = 550.0;
    double rmin = 0.01;
    bool rstop = false;
    double rval = 0.5*(rmax+rmin);
    for (unsigned int i = 0; i < 500 && !rstop; i++)
    {
      rval = 0.5*(rmax+rmin);
      double t1 = erfc(alpha*rval)/rval;
      if (((t1 > 0.9e-10) && (t1 < 1.1e-10)) || (rmax-rmin/2) < rthresh) rstop = true;
      if (!rstop)
      {
        if (t1 > 1e-10)
        {
          rmin = rval;
        }
        else
        {
          rmax = rval;
        }
      }
      if (world.rank() == 0) print("i: ", i, "  rval: ", rval, "  t1: ", t1);
    }

    const double TWOPI = 2.0*constants::pi;
    rvecs = generate_R_vectors(world,std::max(rval,10*L));
    gvecs = generate_G_vectors(world,std::max(std::sqrt(gval2),10*TWOPI/L));
    if (world.rank() == 0) print("gvecs size: ", gvecs.size());
    if (world.rank() == 0) print("rvecs size: ", rvecs.size());
    unsigned int natoms = mentity->natom();
    unsigned int NG = 7;
    gsfactor.push_back(0.0);
    for (unsigned int ig = 1; ig < gvecs.size(); ig++)
//    for (unsigned int ig = 1; ig < NG; ig++)
    {
      // Get G-vector from list
      Vector<double,3> gvec = gvecs[ig];
      double G2 = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
      double_complex rhon = double_complex(0.0,0.0);
//      if (world.rank() == 0)
//        print("ig:  ", ig, gvec[0], gvec[1], gvec[2]);
//      if (world.rank() == 0)
//        print("================================================");
      for (unsigned int ia = 0; ia < natoms; ia++)
      {
        Atom iatom = mentity->get_atom(ia);
        Vector<double,3> tvec {1*(-iatom.x),
                               1*(-iatom.y),
                               1*(-iatom.z)};
        double_complex t1 = std::exp(double_complex(0.0,gvec[0]*tvec[0] +
                                     gvec[1]*tvec[1] + gvec[2]*tvec[2]));
        rhon += iatom.q*t1;
//        if (world.rank() == 0)
//          print("ia:  ", ia, tvec[0], tvec[1], tvec[2]);
//        if (world.rank() == 0)
//          print("t1:  ", t1);
      }
//      if (world.rank() == 0) print("rhon:  ", rhon, "\n");
      double_complex t2 = rhon*std::exp(-G2/4.0/alpha/alpha)/G2;
      double_complex t3 = (std::abs(t2) > 1e-8) ? t2 : 0.0;
      gsfactor.push_back(t3);
//      if (world.rank() == 0)
//        print("ig:  ", ig, "  G2:  ", G2, "  gsfactor:  ", t3);
    }

    if (world.rank() == 0) print("EwaldNuclearPotentialFunctor [end]");
  }

  double_complex do_G_sum_v1(const coordT& r) const
  {
    // number of atoms in unit cell
    unsigned int natoms = mentity->natom();
    // other parameters
    const double TWOPI = 2*constants::pi;
    double v = compute_volume();

    // RECIPROCAL SPACE SUM
    double charge = mentity->total_nuclear_charge();
    double_complex s1 = -charge/alpha/alpha/4.0;
    // skip G=0
    for (unsigned int ig = 1; ig < gvecs.size(); ig++)
    {
      Vector<double,3> gvec = gvecs[ig];
      double_complex t1 = std::exp(double_complex(0.0,gvec[0]*r[0]+gvec[1]*r[1]+gvec[2]*r[2]));
      s1 += gsfactor[ig]*t1;
    }
    s1 *= 2.0*TWOPI/v;
    return s1;
  }

  double_complex do_G_sum_v2(const coordT& r) const
  {
    // number of atoms in unit cell
    unsigned int natoms = mentity->natom();
    // other parameters
    const double TWOPI = 2*constants::pi;
    double v = compute_volume();

    // RECIPROCAL SPACE SUM
    unsigned int NG = 7;
    double charge = mentity->total_nuclear_charge();
    double_complex s1 = -charge/alpha/alpha/4.0;
    for (unsigned int ia = 0; ia < natoms; ia++)
    {
      // skip G=0
      for (unsigned int ig = 1; ig < gvecs.size(); ig++)
//      for (unsigned int ig = 1; ig < NG; ig++)
      {
        Vector<double,3> gvec = gvecs[ig];
        double G2 = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
        Atom iatom = mentity->get_atom(ia);
        Vector<double,3> rr {r[0]-iatom.x,r[1]-iatom.y,r[2]-iatom.z};
        double_complex t1 = std::exp(double_complex(0.0,gvec[0]*rr[0]+gvec[1]*rr[1]+gvec[2]*rr[2]));
        s1 += iatom.q*t1*std::exp(-G2/4.0/alpha/alpha)/G2;
//        print("ig: ", ig, "  g: ", gvec[0], gvec[1], gvec[2]);
//        print("t1: ", t1, "  factor: ", std::exp(-G2/4.0/alpha/alpha)/G2);
//        print("");
      }
    }
    MADNESS_ASSERT(std::imag(s1) < 1e-10);
    s1 *= 2.0*TWOPI/v;
    return s1;
  }

  double_complex do_R_sum_v1(const coordT& r) const
  {
    // number of atoms in unit cell
    unsigned int natoms = mentity->natom();
    // other parameters
    const double TWOPI = 2*constants::pi;
    double v = compute_volume();

    // REAL SPACE SUM
    double_complex s2 = 0.0;
    for (unsigned int ia = 0; ia < natoms; ia++)
    {
      Atom iatom = mentity->get_atom(ia);
      for (unsigned int ir = 0; ir < rvecs.size(); ir++)
      {
        Vector<double,3> rvec = rvecs[ir];
        Vector<double,3> tvec {r[0]-iatom.x-rvec[0],
                               r[1]-iatom.y-rvec[1],
                               r[2]-iatom.z-rvec[2]};
        double t2 = std::sqrt(tvec[0]*tvec[0] +
            tvec[1]*tvec[1] + tvec[2]*tvec[2]);
        s2 += iatom.q*(1-erf(alpha*t2))/t2;
//        print("ia: ", ia, "  x: ", iatom.x, "  y: ", iatom.y, "  z: ", iatom.z);
      }
    }
    return s2;
  }

//  double_complex do_const_sum(const coordT& r) const
//  {
//    // number of atoms in unit cell
//    unsigned int natoms = mentity->natom();
//    // other parameters
//    const double TWOPI = 2*constants::pi;
//    double v = compute_volume();
//
//    double_complex s3 = 0.0;
//    double sqrtpi = std::sqrt(constants::pi);
//    for (unsigned int ia = 0; ia < natoms; ia++)
//    {
//      Atom iatom = mentity->get_atom(ia);
//      double_complex rr = std::abs(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
//      s3 += 2.0*iatom.q*alpha*rr/sqrtpi;
//    }
//    return s3;
//  }

  double operator()(const coordT& r) const
  {
    double_complex s1 = do_G_sum_v1(r);
    double_complex s2 = do_R_sum_v1(r);
//    print("s1: ", std::real(s1), "  s2: ", std::real(s2));
    return -std::real(s1+s2);
  }

//  double operator()(const coordT& r) const
//  {
//    // number of atoms in unit cell
//    unsigned int natoms = mentity->natom();
//    // other parameters
//    const double TWOPI = 2*constants::pi;
//    double v = compute_volume();
//
//    // RECIPROCAL SPACE SUM
//    double charge = mentity->total_nuclear_charge();
//    double_complex s1 = -charge/alpha/alpha/4.0;
//    // skip G=0
//    for (unsigned int ig = 1; ig < gvecs.size(); ig++)
//    {
//      Vector<double,3> gvec = gvecs[ig];
//      double_complex t1 = std::exp(double_complex(0.0,gvec[0]*r[0]+gvec[1]*r[1]+gvec[2]*r[2]));
//      s1 += gsfactor[ig]*t1;
//    }
//    s1 *= 2.0*TWOPI/v;
//
//    // REAL SPACE SUM
//    double_complex s2 = 0.0;
//    for (unsigned int ia = 0; ia < mentity->natom(); ia++)
//    {
//      Atom iatom = mentity->get_atom(ia);
//      Vector<double,3> tvec {iatom.x, iatom.y, iatom.z};
//      for (unsigned int ir = 0; ir < rvecs.size(); ir++)
//      {
//        Vector<double,3> rvec = rvecs[ir];
//        Vector<double,3> tvec2 {r[0]-tvec[0]-rvec[0],
//                                r[1]-tvec[1]-rvec[1],
//                                r[2]-tvec[1]-rvec[1]};
//        double t2 = std::sqrt(tvec2[0]*tvec2[0] +
//            tvec2[1]*tvec2[1] + tvec2[2]*tvec2[2]);
//        s2 += erfc(alpha*t2)/std::abs(t2);
//      }
//    }
//
//    double_complex s3 = 0.0;
//    double sqrtpi = std::sqrt(constants::pi);
//    for (unsigned int ia = 0; ia < natoms; ia++)
//    {
//      Atom iatom = mentity->get_atom(ia);
//      s3 += 2.0*iatom.q*alpha/sqrtpi;
//    }
//
//    double rvalue = std::real(s1 - s3);
////    MADNESS_ASSERT(std::imag(s1) < 1e-10);
////    double rvalue = std::real(s1);
//    return rvalue;
//  }

  virtual ~EwaldNuclearPotentialFunctor() {}
};
//*************************************************************************

//*************************************************************************
double nuclear_potential(const Tensor<double>& c, const Tensor<double>& e, double r)
{
  double val = 0.0;
  for (unsigned int mu = 0; mu < c.dim(0); mu++)
  {
    val += c[mu]*std::exp(-e[mu]*r*r);
  }
  return val;
}
//*************************************************************************

//*************************************************************************
class MolecularNuclearPotentialFunctor2 : public FunctionFunctorInterface<double,3> {
private:
  MolecularEntity mentity;
  Tensor<double> coeff;
  Tensor<double> expnt;
  BoundaryConditions<3> bc;
  std::vector< Vector<double,3> > rvecs;


public:
  MolecularNuclearPotentialFunctor2(World& world, const MolecularEntity& mentity,
      const double& R, double eps, double lo=1e-12) : mentity(mentity)
  {
    bc=FunctionDefaults<3>::get_bc();
    if (bc(0,0) == BC_PERIODIC)
    {
      rvecs = generate_R_vectors(world, 8.0*R);
    }
    else
    {
      rvecs.push_back({0.0,0.0,0.0});
    }
    const Tensor<double>& cell_width = FunctionDefaults<3>::get_cell_width();
    double hi = cell_width.normf(); // Diagonal width of cell
    if (bc(0,0) == BC_PERIODIC) hi *= 100;
    else hi *= 5; // Extend range for periodic summation
    const double pi = constants::pi;

    print("hi:  ", hi);

//    gen_ce(0.0,1e-5,eps,coeff,expnt);

    // bsh_fit generates representation for 1/4Pir but we want 1/r
    // so have to scale eps by 1/4Pi
    bsh_fit(0.0, lo, hi, eps/(4.0*pi), &coeff, &expnt, false);

    if (bc(0,0) == BC_PERIODIC) {
        truncate_periodic_expansion(coeff, expnt, cell_width.max(), true);
    }
    coeff.scale(4.0*pi);

    for (unsigned int i = 0; i < coeff.dim(0); i++)
    {
      coeff[i] = coeff[i]*R;
      expnt[i] = expnt[i]*R*R;
    }

//    for (unsigned int i = 0; i < coeff.dim(0); i++)
//    {
//      print("i:  ", i,"  coeff:  ", coeff[i], "  expnt:  ", expnt[i]);
//    }
//
//    for (unsigned int i = 0; i < rvecs.size(); i++)
//    {
//      Vector<double,3> rvec = rvecs[i];
//      double R2 = rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];
//      print("i:  ", i, "  R:  ", rvecs[i], "  R2:  ", R2);
//    }
  }

//  double operator()(const coordT& r) const
//  {
//    double rval = 0.0;
//    for (unsigned int ia = 0; ia < mentity.natom(); ia++)
//    {
//      Atom iatom = mentity.get_atom(ia);
//      for (unsigned int mu = 0; mu < coeff.dim(0); mu++)
//      {
//        for (unsigned int ir = 0; ir < rvecs.size(); ir++)
//        {
//          Vector<double,3> rvec = rvecs[ir];
//          Vector<double,3> tvec {r[0]-iatom.x-rvec[0],
//                                 r[1]-iatom.y-rvec[1],
//                                 r[2]-iatom.z-rvec[2]};
//          double T2 = tvec[0]*tvec[0] + tvec[1]*tvec[1] + tvec[2]*tvec[2];
//          rval += std::exp(-expnt[mu]*T2);
//        }
//        rval *= coeff[mu];
//      }
//      rval *= iatom.q;
//    }
//    return rval;
//  }

  double operator()(const coordT& r) const
  {
    double rval = 0.0;
    for (unsigned int ia = 0; ia < mentity.natom(); ia++)
    {
      Atom iatom = mentity.get_atom(ia);
      for (unsigned int ir = 0; ir < rvecs.size(); ir++)
      {
        Vector<double,3> rvec = rvecs[ir];
        Vector<double,3> tvec {r[0]-iatom.x-rvec[0],
                               r[1]-iatom.y-rvec[1],
                               r[2]-iatom.z-rvec[2]};
        double t1 = std::sqrt(tvec[0]*tvec[0] + tvec[1]*tvec[1] + tvec[2]*tvec[2]);
        rval += -iatom.q*nuclear_potential(coeff,expnt,t1);
      }
    }
    return rval;
  }

  Tensor<double> get_coeff() {return coeff;}
  Tensor<double> get_expnt() {return expnt;}
};
//*************************************************************************

//*************************************************************************
class MolecularNuclearChargeDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const MolecularEntity& _mentity;
    const double R;
    const bool periodic;
    const std::vector<coordT> _specialpts;
    double charge;

public:
    MolecularNuclearChargeDensityFunctor(const MolecularEntity& mentity, const double& R,
        const bool& periodic, const std::vector<coordT>& specialpts)
      : _mentity(mentity), R(R), periodic(periodic), _specialpts(specialpts) {
      charge = mentity.total_nuclear_charge();
      print("charge is: ", charge);
    }

    virtual std::vector<coordT> special_points() const
    {
      return _specialpts;
    }

    virtual Level special_level()
    {
      return 10;
    }

    double operator()(const coordT& x) const
    {
        double big = 0.5*R + 6.0*_mentity.smallest_length_scale();
//        double big = 2*R + 6.0*_mentity.smallest_length_scale();
        // Only one contribution at any point due to the short
        // range of the nuclear charge density
        //printf("big: %10.8f\n\n", big);
        double value = 0.0;
//        double value = charge;
        if (periodic)
        {
            for (int xr = -1; xr <= 1; xr += 1)
            {
                double xx = x[0] + xr*R;
                //printf("x[0]: %10.8f     xx: %10.8f\n", x[0], xx);
                if (xx < big && xx > -big)
                {
                    for (int yr = -1; yr <= 1; yr += 1)
                    {
                        double yy = x[1] + yr*R;
                        //printf("y[0]: %10.8f     yy: %10.8f\n", x[1], yy);
                        if (yy < big && yy > -big)
                        {
                            for (int zr = -1; zr <= 1; zr += 1)
                            {
                                double zz = x[2] + zr*R;
                                //printf("z[0]: %10.8f     zz: %10.8f\n", x[2], zz);
                                if (zz < big && zz > -big)
                                {
                                    double t1 = _mentity.nuclear_charge_density(xx, yy, zz);
                                    value += t1;
                                    //printf("t1: %10.8f     value: %10.8f\n", t1, value);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            value = _mentity.nuclear_charge_density(x[0], x[1], x[2]);
        }
        return value;
    }
};
//*************************************************************************

//*************************************************************************
rfunctionT make_nuclear_charge_density(World& world,
    const MolecularEntity& mentity, double thresh = 1e-6)
{
  std::vector<coordT> specialpts;
  for (int i = 0; i < mentity.natom(); i++)
  {
    coordT pt(0.0);
    Atom atom = mentity.get_atom(i);
    pt[0] = atom.x; pt[1] = atom.y; pt[2] = atom.z;
    specialpts.push_back(pt);
//    if (world.rank() == 0) print("Special point: ", pt);
  }

  rfunctionT rhon = rfactoryT(world).functor(
      rfunctorT(new MolecularNuclearChargeDensityFunctor(mentity, L, true, specialpts))).
      thresh(thresh).initial_level(6).truncate_on_project();

  return rhon;

}
//*************************************************************************

//*************************************************************************
rvecfuncT make_nuclear_charge_density_individual(World& world,
            const MolecularEntity& mentity, double thresh)
{
  rvecfuncT ndensity;
  for (int i = 0; i < mentity.natom(); i++)
  {
    // do special points for single density
    coordT pt(0.0);
    Atom atom = mentity.get_atom(i);
    pt[0] = atom.x; pt[1] = atom.y; pt[2] = atom.z;
    std::vector<coordT> specialpts;
    specialpts.push_back(pt);
    // create single density
    MolecularEntity m = mentity.get_entity(i);
    rfunctionT rho_i = rfactoryT(world).functor(
        rfunctorT(new MolecularNuclearChargeDensityFunctor(m, L, true, specialpts))).
        thresh(thresh).initial_level(6).truncate_on_project();
    ndensity.push_back(rho_i);
  }
  return ndensity;
}
//*************************************************************************

//*************************************************************************

//*************************************************************************

//***************************************************************************
void test_gaussian_num_coeffs(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(10);
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
      FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

      // Create normalized gaussian function with wide ranging exponents
      for (int in = 0; in < 20; in++)
      {
        int n = in - 8;
        double expnt = std::pow(2.0, (int) n);
        double coeff = std::pow(expnt/constants::pi, 1.5);
        rfunctionT fexp = rfactoryT(world).functor(
            rfunctorT(new GaussianFunctor(coeff,expnt))).truncate_on_project();
        // how many nodes needed to represent function
        int maxnodes = fexp.max_nodes();
        double tr = fexp.trace();
        // full-width half maximum
        double fwhm = 2*std::log(2)/expnt;
        // full-width 1/10th maximum
        double fwtm = 2*std::log(10)/expnt;
        if (world.rank() == 0) print("n:  ", n, "  coeff:  ", coeff, "  expnt:  ",
            expnt, "  trace:  ", tr, "  max nodes:  ", maxnodes, "  fwhm:  ", fwhm);
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();


}
//***************************************************************************


//*************************************************************************
void compute_madelung_energy_PWSCF(World& world, MolecularEntity mentity,
    double alpha = 8.5, double rmax = 100.0, double gmax = 100.0)
{
  // generate real and reciprocal lattice vectors
  std::vector< Vector<double,3> > rvecs = generate_R_vectors(world, rmax);
  std::vector< Vector<double,3> > gvecs = generate_G_vectors(world, gmax);
  if (world.rank() == 0)
    printf("rvecs size: %d     gvecs size: %d\n", rvecs.size(), gvecs.size());
  // number of atoms in unit cell
  unsigned int natoms = mentity.natom();
  // other parameters
  const double TWOPI = 2*constants::pi;
  double v = compute_volume();

  if (world.rank() == 0) mentity.print();

  // RECIPROCAL SPACE SUM
  int NG = 7;
  double charge = mentity.total_nuclear_charge();
  double_complex s1 = -charge*charge/alpha/alpha/4.0;
    // skip G=0
    for (unsigned int ig = 1; ig < gvecs.size(); ig++)
//    for (unsigned int ig = 1; ig < NG; ig++)
    {
      Vector<double,3> gvec = gvecs[ig];
      double G2 = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
      double_complex rhon = double_complex(0.0,0.0);
      for (unsigned int ia = 0; ia < natoms; ia++)
      {
        Atom iatom = mentity.get_atom(ia);
        Vector<double,3> tvec {1*(-iatom.x),
                               1*(-iatom.y),
                               1*(-iatom.z)};
        double_complex t1 = std::exp(double_complex(0.0,gvec[0]*tvec[0] +
                                     gvec[1]*tvec[1] + gvec[2]*tvec[2]));
        rhon += iatom.q*t1;
      }
//      print("rhon:   ", std::abs(rhon));
      s1 += std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2;
//      print("extra amount:     ", std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2);
    }
  s1 *= 2.0*TWOPI/v;

//  // REAL SPACE SUM
  double_complex s2 = 0.0;
  // loop over R-lattice vectors
  for (unsigned int ir = 0; ir < rvecs.size(); ir++)
  {
    Vector<double,3> rvec = rvecs[ir];
    // loop through the atoms in the unit cell
    for (unsigned int ia = 0; ia < natoms; ia++)
    {
      for (unsigned int ja = 0; ja < natoms; ja++)
      {
        // do not include term if R=(0.0,0.0,0.0) and i==j
        if ((ir > 0) || ((ir==0) && (ia != ja)))
        {
          Atom iatom = mentity.get_atom(ia);
          Atom jatom = mentity.get_atom(ja);
          Vector<double,3> dvec {L*(iatom.x-jatom.x),
                                 L*(iatom.y-jatom.y),
                                 L*(iatom.z-jatom.z)};
          Vector<double,3> tvec {dvec[0]+rvec[0],
                                 dvec[1]+rvec[1],
                                 dvec[2]+rvec[2]};
          double tnorm = std::sqrt(tvec[0]*tvec[0] + tvec[1]*tvec[1] +
                                   tvec[2]*tvec[2]);
//              s2 += iatom.q*jatom.q*(1.0 - erf(alpha*tnorm))/tnorm;
          s2 += iatom.q*jatom.q*(erfc(alpha*tnorm))/tnorm;
        }
      }
    }
  }
  s2 *= 0.5;
//
  double_complex s3 = 0.0;
  if (world.rank() == 0) printf("\n");
  double sqrtpi = std::sqrt(constants::pi);
  for (unsigned int ia = 0; ia < natoms; ia++)
  {
    Atom iatom = mentity.get_atom(ia);
//    printf("atom %d    charge: %8.4f\n",ia,iatom.q);
    s3 += 2.0*iatom.q*iatom.q*alpha/sqrtpi;
//    s3 += iatom.q*iatom.q*alpha*std::sqrt(8.0/TWOPI);
//    print("update energy: ", 2.0*iatom.q*iatom.q*alpha/sqrtpi);
  }
//  print("value: ", 2*alpha/sqrtpi);


  double energy = std::real(s1 + s2 - s3);

  if (world.rank()==0)
  {
    printf("alpha: %8.4f G-sum: %8.4f  R-sum: %8.4f  total energy:  %15.7f\n",
        alpha, std::real(s1-s3), std::real(s2), energy);
  }
}
//*************************************************************************

//*************************************************************************
void compute_madelung_energy(World& world, MolecularEntity mentity,
    double alpha = 1.5, double rmax = 200.0, double gmax = 200.0)
{
  // generate real and reciprocal lattice vectors
  std::vector< Vector<double,3> > rvecs = generate_R_vectors(world, rmax);
  std::vector< Vector<double,3> > gvecs = generate_G_vectors(world, gmax);
  if (world.rank() == 0)
    printf("rvecs size: %d     gvecs size: %d\n", rvecs.size(), gvecs.size());
  // number of atoms in unit cell
  unsigned int natoms = mentity.natom();
  // other parameters
  const double TWOPI = 2*constants::pi;
  double v = compute_volume();
  double charge = -mentity.total_nuclear_charge();

  // RECIPROCAL SPACE SUM
  double_complex s1 = 0.0;
  // skip G=0
  for (unsigned int ig = 1; ig < gvecs.size(); ig++)
  //for (unsigned int ig = 1; ig < 7; ig++)
  {
    Vector<double,3> gvec = gvecs[ig];
    double G2 = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
    double_complex rhon = double_complex(0.0,0.0);
    for (unsigned int ia = 0; ia < natoms; ia++)
    {
      Atom iatom = mentity.get_atom(ia);
      Vector<double,3> tvec {1*(iatom.x),
                             1*(iatom.y),
                             1*(iatom.z)};
      double_complex t1 = std::exp(double_complex(0.0,gvec[0]*tvec[0] +
                                   gvec[1]*tvec[1] + gvec[2]*tvec[2]));
      rhon += iatom.q*t1;
    }
    s1 += std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2;
//    print("abs(rhon):  ", std::abs(rhon), "  G1:  ", std::sqrt(G2), "     g-update: ",
//        std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2);

  }
  s1 *= TWOPI/v;

//  print("reciprocal space sum:  ", s1);

  // REAL SPACE SUM
  double_complex s2 = 0.0;
  // loop over R-lattice vectors
  for (unsigned int ir = 0; ir < rvecs.size(); ir++)
//    for (unsigned int ir = 0; ir < 7; ir++)
  {
    Vector<double,3> rvec = rvecs[ir];
    // loop through the atoms in the unit cell
    for (unsigned int ia = 0; ia < natoms; ia++)
    {
      for (unsigned int ja = 0; ja <= ia ; ja++)
      {
        // do not include term if R=(0.0,0.0,0.0) and i==j
        if ((ir > 0) || ((ir==0) && (ia != ja)))
        {
          Atom iatom = mentity.get_atom(ia);
          Atom jatom = mentity.get_atom(ja);
          Vector<double,3> dvec {1*(iatom.x-jatom.x),
                                 1*(iatom.y-jatom.y),
                                 1*(iatom.z-jatom.z)};
          Vector<double,3> tvec {dvec[0]+rvec[0],
                                 dvec[1]+rvec[1],
                                 dvec[2]+rvec[2]};
          double tnorm = std::sqrt(tvec[0]*tvec[0] + tvec[1]*tvec[1] +
                                   tvec[2]*tvec[2]);
//              s2 += iatom.q*jatom.q*(1.0 - erf(alpha*tnorm))/tnorm;
          s2 += 0.5*iatom.q*jatom.q*(erfc(alpha*tnorm))/tnorm;
//          print("tnorm:  ", tnorm, "     r-update: ", iatom.q*jatom.q*(erfc(alpha*tnorm))/tnorm);
        }
      }
    }
  }

  double_complex s3 = 0.0;
  double sqrtpi = std::sqrt(constants::pi);
  for (unsigned int ia = 0; ia < natoms; ia++)
  {
    Atom iatom = mentity.get_atom(ia);
    s3 += iatom.q*iatom.q*alpha/sqrtpi;
  }

  // MISC. TERM
  double s4 = -charge*charge/alpha/alpha/4.0;


  double energy = std::abs(s1) + std::abs(s2) - std::abs(s3);

  if (world.rank()==0)
  {
    printf("G-sum:  %8.4f    %8.4f\n",std::real(s1), std::imag(s1));
    printf("R-sum:  %8.4f    %8.4f\n",std::real(s2), std::imag(s2));
    printf("C-sum:  %8.4f    %8.4f\n",std::real(s3), std::imag(s3));
    printf("C2-sum: %8.4f\n",s4);
    printf("sum1:   %8.4f\n\n", std::abs(s1-s3+s2));
    printf("sum2:   %8.4f\n\n", std::abs(s1-s3+s2+s4));
    printf("\n\nalpha: %8.4f G-sum: %8.4f  C-sum: %8.4f  R-sum: %8.4f  total energy:  %15.7f\n",
        alpha, std::abs(s1), std::abs(s3), std::abs(s2), energy);
//    printf("%8.4f   %15.7f\n",
//        alpha, energy);
  }
}
//*************************************************************************

//*************************************************************************
void test_gence(int argc, char** argv)
{
  Tensor<double> c,e;
  gen_ce(0.0, 1e-8, 1e-12, c, e);
  unsigned int npts = 250;
  double LL = 100.0;
  for (unsigned int i = 1; i < npts; i++)
  {
    double r = i*LL/npts;
    double val = 0.0;
    for (unsigned int mu = 0; mu < c.dim(0); mu++)
    {
      val += c[mu]*std::exp(-e[mu]*r*r);
    }
    double val2 = 1/r;
    double rerr = std::abs(val-val2)/val;
    printf("%10.5f  %15.10e  %15.10e  %15.10e\n", r, val2, val, rerr);
  }
}
//*************************************************************************

//*************************************************************************
void test_gence2(int argc, char** argv)
{
  Tensor<double> c,e;
//  gen_ce(0.0, 1e-8, 1e-12, c, e);
  // bsh_fit generates representation for 1/4Pir but we want 1/r
  // so have to scale eps by 1/4Pi
  double lo = 1e-8;
  double eps = 1e-12;
  bsh_fit(0.0, lo, 100.0*L, eps/(4.0*constants::pi), &c, &e, false);

//  if (bc(0,0) == BC_PERIODIC) {
//      truncate_periodic_expansion(c, e, cell_width.max(), true);
//  }
  c.scale(4.0*constants::pi);

  for (unsigned int i = 0; i < c.dim(0); i++)
  {
    c[i] = c[i]*L;
    e[i] = e[i]*L*L;
  }


  unsigned int npts = 250;
  double LL = 100.0;
  for (unsigned int i = 1; i < npts; i++)
  {
    double r = i*LL/npts;
    double val = nuclear_potential(c, e, r);
    double val2 = 1/r;
    double rerr = std::abs(val-val2)/val;
    printf("%10.5f  %15.10e  %15.10e  %15.10e\n", r, val2, val, rerr);
  }
}
//*************************************************************************

//*************************************************************************
void test_nuclear_potential2(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(8);
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_FREE));
      FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

      MolecularEntity mentity;
      mentity.add_atom(0.0,0.0,0.0,1.0,1.0);
      mentity.center();

      MolecularNuclearPotentialFunctor2 npoty(world, mentity, L, 1e-16);
      if (world.rank() == 0) mentity.print();

      Tensor<double> c = npoty.get_coeff();
      Tensor<double> e = npoty.get_expnt();

      unsigned int npts = 153;
      for (unsigned int ir = 1; ir < npts; ir++)
      {
        double rpt = ir*L/npts-L/2;
        coordT pt(rpt);
        double npotpt = npoty(pt);
        rpt *= std::sqrt(3);
        double npotpt2 = -1/std::abs(rpt);
        double nerror = std::abs((npotpt-npotpt2)/npotpt2);
        if (world.rank() == 0) printf("%10.5f     %15.10e     %15.10e     %15.10e\n", rpt, npotpt, npotpt2, nerror);
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();
}
//*************************************************************************

//*************************************************************************
void test_nuclear_potential(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(8);
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
      FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

      MolecularEntity mentity;
      mentity.add_atom(0.0,0.0,0.0,7,9);
      mentity.add_atom(L/2,L/2,L/2,1,3);
      mentity.center();

//      if (world.rank() == 0) print("Making Ewald nuclear potential ...");
//      rfunctionT npot2 = rfactoryT(world).functor(
//          rfunctorT(new EwaldNuclearPotentialFunctor(world, &mentity, 0.2))).
//          thresh(1e-6).initial_level(6).truncate_on_project();
      if (world.rank() == 0) print("Making charge density ...");
      rfunctionT ndensity = make_nuclear_charge_density(world,mentity,1e-6);
      if (world.rank() == 0) print("Setting avg value of charge density to zero ...");
      double ndensav = ndensity.trace();
      if (world.rank() == 0) print("avg value of charge density: ", ndensav);
      rfunctionT avgval = rfactoryT(world).functor(rfunctorT(new OnesFunctor(ndensav)));
      //ndensity -= avgval;

      print("mentity total charge:  ", mentity.total_nuclear_charge());

      SeparatedConvolution<double,3> cop =
          CoulombOperator(world,1e-2,FunctionDefaults<3>::get_thresh() * 0.1);
      if (world.rank() == 0) print("Making normal nuclear potential ...");
      rfunctionT npot1 = apply(cop,ndensity);
      if (world.rank() == 0) print("Creating difference potential ...");
//      rfunctionT npotd = npot1-npot2;
//      double nptr = npotd.trace();
//      if (world.rank() == 0) print("The trace of the difference of potentials: ", nptr);

//      if (world.rank() == 0) print("Making *new* nuclear potential ...");
//      rfunctionT vnuc2 = rfactoryT(world).functor(
//          rfunctorT(new MolecularNuclearPotentialFunctor2(world, mentity, L, 1e-8))).
//          truncate_on_project();
//      if (world.rank() == 0) print("Done making normal nuclear potential ...");


      //      EwaldNuclearPotentialFunctor ewaldf(world, &mentity, 0.7);
      //      EwaldNuclearPotentialFunctor ewaldf1(world, &mentity, 1.1);
      //      EwaldNuclearPotentialFunctor ewaldf2(world, &mentity, 2.0);
      EwaldNuclearPotentialFunctor ewaldf(world, &mentity, 0.05);
      EwaldNuclearPotentialFunctor ewaldf1(world, &mentity, 0.2);
      EwaldNuclearPotentialFunctor ewaldf2(world, &mentity, 0.4);
      MolecularNuclearPotentialFunctor2 mnpotf(world, mentity, L , 1e-7, 1e-2);
      unsigned int npts = 250;
      for (unsigned int ir = 0; ir < npts; ir++)
      {
        double rpt = ir*L/npts-L/2;
        coordT pt(rpt);
        double npot1pt = npot1(pt);
        double npot2pt = ewaldf(pt);
        double npot3pt = ewaldf1(pt);
        double npot4pt = ewaldf2(pt);
        double npot5pt = mnpotf(pt);
        double nerror2 = std::abs(npot1pt-npot2pt);
        double nerror3 = std::abs(npot1pt-npot3pt);
        double nerror4 = std::abs(npot1pt-npot4pt);
        double nerror5 = std::abs(npot1pt-npot5pt);
//        if (world.rank() == 0) print(rpt, "  ", npot1pt, "  ", npot2pt,
//            "  ", npot3pt, "  ", npot4pt, "  ", nerror2, nerror3, nerror4);
//        if (world.rank() == 0) print(rpt, "  ", npot1pt, "  ", npot3pt, "  ",
//            npot5pt, "  ", nerror3, "  ", nerror5);
        if (world.rank() == 0)
          printf("%10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e\n",
            rpt, npot1pt, npot2pt, npot3pt, npot4pt, npot5pt,
            nerror2, nerror3, nerror4, nerror5);
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();


}
//*************************************************************************

//*************************************************************************
void test_nuclear_potential3(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(8);
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
      FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

      MolecularEntity mentity;
      mentity.add_atom(0.0,0.0,0.0,7,9);
      mentity.add_atom(L/2,L/2,L/2,1,3);
      mentity.center();

      if (world.rank() == 0) print("Making charge density ...");
      rfunctionT ndensity2 = make_nuclear_charge_density(world,mentity,1e-2);
      rfunctionT ndensity4 = make_nuclear_charge_density(world,mentity,1e-4);
      rfunctionT ndensity6 = make_nuclear_charge_density(world,mentity,1e-6);
      rfunctionT ndensity8 = make_nuclear_charge_density(world,mentity,1e-8);

      if (world.rank() == 0) print("Making normal nuclear potential (1e-2) ...");
      SeparatedConvolution<double,3> cop2 =
          CoulombOperator(world,1e-2,1e-2);
      rfunctionT npot2 = apply(cop2,ndensity2);
      if (world.rank() == 0) print("Making normal nuclear potential (1e-4) ...");
      SeparatedConvolution<double,3> cop4 =
          CoulombOperator(world,1e-2,1e-4);
      rfunctionT npot4 = apply(cop4,ndensity4);
      if (world.rank() == 0) print("Making normal nuclear potential (1e-6) ...");
      SeparatedConvolution<double,3> cop6 =
          CoulombOperator(world,1e-2,1e-6);
      rfunctionT npot6 = apply(cop6,ndensity6);
      if (world.rank() == 0) print("Making normal nuclear potential (1e-8) ...");
      SeparatedConvolution<double,3> cop8 =
          CoulombOperator(world,1e-2,1e-8);
      rfunctionT npot8 = apply(cop8,ndensity8);

      unsigned int npts = 250;
      for (unsigned int ir = 0; ir < npts; ir++)
      {
        double rpt = ir*L/npts-L/2;
        coordT pt(rpt);
        double npot2pt = npot2(pt);
        double npot4pt = npot4(pt);
        double npot6pt = npot6(pt);
        double npot8pt = npot8(pt);
        double nerror2 = std::abs(npot8pt-npot2pt);
        double nerror4 = std::abs(npot8pt-npot4pt);
        double nerror6 = std::abs(npot8pt-npot6pt);
        double nerror8 = std::abs(npot8pt-npot8pt);
        if (world.rank() == 0)
          printf("%10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e  %10.7e\n",
            rpt, npot2pt, npot4pt, npot6pt, npot8pt,
            nerror2, nerror4, nerror6, nerror8);
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();


}
//*************************************************************************

void test_nuclear_potential_big_unit_cell(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(8);
//      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_FREE));
      FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

      MolecularEntity mentity;
      mentity.add_atom(0.0,0.0,0.0,7,9);
      mentity.add_atom(L/2,L/2,L/2,1,3);
      mentity.center();

      if (world.rank() == 0) print("Making charge density ...");
      rfunctionT ndensity = make_nuclear_charge_density(world,mentity,1e-6);
      SeparatedConvolution<double,3> cop =
          CoulombOperator(world,1e-2,FunctionDefaults<3>::get_thresh() * 0.1);
      if (world.rank() == 0) print("Making normal nuclear potential ...");
      rfunctionT npot1 = apply(cop,ndensity);
      if (world.rank() == 0) print("Creating difference potential ...");
      unsigned int npts = 250;

      Atom at0 = mentity.get_atom(0);
      Atom at1 = mentity.get_atom(1);
      mentity.print();
      for (unsigned int ir = 0; ir < npts; ir++)
      {
        double rpt = ir*L/npts-L/2;
        coordT pt(rpt);
        double npot1pt = npot1(pt);
        double npot2pt = -at0.q/std::sqrt((rpt-at0.x)*(rpt-at0.x) +
            (rpt-at0.y)*(rpt-at0.y) + (rpt-at0.z)*(rpt-at0.z)) -
            at1.q/std::sqrt((rpt-at1.x)*(rpt-at1.x) +
                (rpt-at1.y)*(rpt-at1.y) + (rpt-at1.z)*(rpt-at1.z));
        if (world.rank() == 0) print(rpt, "  ", npot1pt, "  ", npot2pt);
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();


}
//*************************************************************************


//*************************************************************************
void test_G_R_vectors(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);

      std::vector<Vector<double,3> > gvecs = generate_G_vectors(world,50.0);
      std::vector<Vector<double,3> > rvecs = generate_R_vectors(world,150.0);

      for (unsigned int ir = 0; ir < rvecs.size(); ir++)
      {
        Vector<double,3> rvec = rvecs[ir];
        for (unsigned int ig = 0; ig < gvecs.size(); ig++)
        {
          Vector<double,3> gvec = gvecs[ig];
          double t1 = rvec[0]*gvec[0] + rvec[1]*gvec[1] + rvec[2]*gvec[2];
          double_complex t2 = std::exp(double_complex(0.0,t1));
          MADNESS_ASSERT(std::abs(std::imag(t2)) < 1e-10);
        }
      }

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();
}
//*************************************************************************

//*************************************************************************
void test_nuclear_energy(int argc, char** argv)
{
  initialize(argc, argv);

  World world(SafeMPI::COMM_WORLD);

  try {
      // Load info for MADNESS numerical routines
      startup(world,argc,argv);
      std::cout.precision(6);
      FunctionDefaults<3>::set_thresh(1e-6);
      FunctionDefaults<3>::set_k(8);
      FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
      FunctionDefaults<3>::set_cubic_cell(0,L);

      MolecularEntity mentity;
      mentity.add_atom(0.0,0.0,0.0,7,9);
      mentity.add_atom(L/2,L/2,L/2,1,3);
      mentity.center();

      if (world.rank() == 0) mentity.print();

//        int ntype = mentity.get_num_types();
//        int n1 = mentity.get_num_atoms_type(9);
//        int n2 = mentity.get_num_atoms_type(3);
//        printf("Number of type of atoms: %d\n", ntype);
//        printf("F: %d\n", n1);
//        printf("Li: %d\n", n2);

//        for (unsigned int i = 0; i < mentity.atoms.size(); i++)
//        {
//          printf("atom %d --- (%8.4f, %8.4f, %8.4f)\n");
//        }


//        for (unsigned int i = 0; i < 800; i++)
//        {
//          compute_madelung_energy2(world, mentity, 0.05 + i*0.005);
//        }

//        compute_madelung_energy_PWSCF(world,mentity,1.14017543,100.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.2,50.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.25,50.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.3,50.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.35,50.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.4,50.0,100.0);
      compute_madelung_energy_PWSCF(world,mentity,1.45,50.0,100.0);
//      compute_madelung_energy_PWSCF(world,mentity,0.3,50.0,100.0);

//        for (unsigned i = 0; i < 20; i++)
//        {
//          compute_madelung_energy2(world, mentity, 0.004, 0.0+i*50.0, 200.0);
//        }

//        compute_madelung_energy2(world, mentity, 1.2, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.3, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.4, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.5, 100.0, 100.0);

  } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
  } catch (const madness::MadnessException& e) {
      if (world.rank() == 0) print(e);
      error("caught a MADNESS exception");
  } catch (const madness::TensorException& e) {
      if (world.rank() == 0) print(e);
      error("caught a Tensor exception");
  } catch (const char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (char* s) {
      if (world.rank() == 0) print(s);
      error("caught a string exception");
  } catch (const std::string& s) {
      if (world.rank() == 0) print(s);
      error("caught a string (class) exception");
  } catch (const std::exception& e) {
      if (world.rank() == 0) print(e.what());
      error("caught an STL exception");
  } catch (...) {
      error("caught unhandled exception");
  }

  finalize();


}
//*************************************************************************

int main(int argc, char** argv)
{
//    test_nuclear_energy(argc,argv);
  //test_gaussian_num_coeffs(argc,argv);
  test_nuclear_potential3(argc,argv);
  //test_gence2(argc,argv);
  //test_G_R_vectors(argc,argv);
  return 0;
}

