#define WORLD_INSTANTIATE_STATIC_TEMPLATES

#include <mra/mra.h>
#include "mentity.h"

using namespace madness;

static double L = 6.5;

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
std::vector< Vector<double,3> > generate_R_vectors(World& world, double maxRlen = 50.0)
{
  const double t1 = L;

  std::vector< Vector<double,3> > rvecs;

  int rlo = -30;
  int rhi = 31;
  for (int ir1 = rlo; ir1 < rhi; ir1++)
  {
    for (int ir2 = rlo; ir2 < rhi; ir2++)
    {
      for (int ir3 = rlo; ir3 < rhi; ir3++)
      {
        double rlen = t1*std::sqrt(ir1*ir1 + ir2*ir2 + ir3*ir3);
        if (rlen <= maxRlen)
        {
          Vector<double,3> rvec = vec(t1*ir1, t1*ir2, t1*ir3);
          rvecs.push_back(rvec);
        }
      }
    }
  }
  std::sort(rvecs.begin(), rvecs.end(), vectorLengthFunctor());
//      if (_world.rank() == 0) printf("Size of vectors:  %d\n", rvecs.size());
//      if (_world.rank() == 0)
//        printf("\nR-vectors:\n");
  for (unsigned int ir = 0; ir < rvecs.size(); ir++)
  {
    Vector<double,3> rvec = rvecs[ir];
    double rlen = std::sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);
//        if (_world.rank() == 0)
//          printf("%10.5f %10.5f %10.5f   %10.5f\n",rvec[0],rvec[1],rvec[2],rlen);
  }

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

  int glo = -30;
  int ghi = 31;
  for (int ig1 = glo; ig1 < ghi; ig1++)
  {
    for (int ig2 = glo; ig2 < ghi; ig2++)
    {
      for (int ig3 = glo; ig3 < ghi; ig3++)
      {
        double glen = t1*std::sqrt(ig1*ig1 + ig2*ig2 + ig3*ig3);
        if (glen <= maxGlen)
        {
          Vector<double,3> gvec = vec(t1*ig1, t1*ig2, t1*ig3);
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

  return gvecs;

}
//*************************************************************************


void compute_madelung_energy(World& world, MolecularEntity mentity,
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

  // RECIPROCAL SPACE SUM
  int NG = 7;
  //      double_complex s1 = -12.0/alpha/4.0;
  double_complex s1 = -64.0/alpha/alpha/4.0;
  print("initial alpha = ", alpha);
  print("initial s1 = ", s1);
    // skip G=0
    //for (unsigned int ig = 1; ig < gvecs.size(); ig++)
    for (unsigned int ig = 1; ig < NG; ig++)
    {
      Vector<double,3> gvec = gvecs[ig];
      double G2 = gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2];
      double_complex rhon = double_complex(0.0,0.0);
      for (unsigned int ia = 0; ia < natoms; ia++)
      {
        Atom iatom = mentity.get_atom(ia);
        Vector<double,3> tvec = vec(1*(-iatom.x),
                                    1*(-iatom.y),
                                    1*(-iatom.z));
        double_complex t1 = std::exp(double_complex(0.0,gvec[0]*tvec[0] +
                                     gvec[1]*tvec[1] + gvec[2]*tvec[2]));
        rhon += iatom.q*t1;
      }
      print("rhon:   ", std::abs(rhon));
      s1 += std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2;
      print("extra amount:     ", std::abs(rhon)*std::abs(rhon)*std::exp(-G2/(4.0*alpha*alpha))/G2);
    }
  s1 *= 2.0*TWOPI/v;

//  // REAL SPACE SUM
//  double_complex s2 = 0.0;
//  // loop over R-lattice vectors
//  for (unsigned int ir = 0; ir < rvecs.size(); ir++)
//  {
//    Vector<double,3> rvec = rvecs[ir];
//    // loop through the atoms in the unit cell
//    for (unsigned int ia = 0; ia < natoms; ia++)
//    {
//      for (unsigned int ja = 0; ja < natoms; ja++)
//      {
//        // do not include term if R=(0.0,0.0,0.0) and i==j
//        if ((ir > 0) || ((ir==0) && (ia != ja)))
//        {
//          Atom iatom = mentity.get_atom(ia);
//          Atom jatom = mentity.get_atom(ja);
//          Vector<double,3> dvec = vec(L*(iatom.x-jatom.x),
//                                      L*(iatom.y-jatom.y),
//                                      L*(iatom.z-jatom.z));
//          Vector<double,3> tvec = vec(dvec[0]+rvec[0],
//                                      dvec[1]+rvec[1],
//                                      dvec[2]+rvec[2]);
//          double tnorm = std::sqrt(tvec[0]*tvec[0] + tvec[1]*tvec[1] +
//                                   tvec[2]*tvec[2]);
////              s2 += iatom.q*jatom.q*(1.0 - erf(alpha*tnorm))/tnorm;
//          s2 += iatom.q*jatom.q*(erfc(alpha*tnorm))/tnorm;
//        }
//      }
//    }
//  }
//  s2 *= 0.5;
//
//  double_complex s3 = 0.0;
//  double sqrtpi = std::sqrt(constants::pi);
//  for (unsigned int ia = 0; ia < natoms; ia++)
//  {
//    Atom iatom = mentity.get_atom(ia);
//    s3 += iatom.q*iatom.q*alpha/sqrtpi;
//  }
//
//  double energy = std::abs(s1) + std::abs(s2) - std::abs(s3);
//
//  if (world.rank()==0)
//  {
//    printf("alpha: %8.4f G-sum: %8.4f  R-sum: %8.4f  total energy:  %15.7f\n",
//        alpha, std::abs(s1)-std::abs(s3), std::abs(s2), energy);
//  }
}

int main(int argc, char** argv)
{
    initialize(argc, argv);

    World world(MPI::COMM_WORLD);

    try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        std::cout.precision(6);
        FunctionDefaults<3>::set_thresh(1e-6);
        FunctionDefaults<3>::set_k(8);
        FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
        FunctionDefaults<3>::set_cubic_cell(0,L);

        MolecularEntity mentity;
        mentity.add_atom(0.0,0.0,0.0,9,7);
        mentity.add_atom(L/2,L/2,L/2,3,3);
        mentity.center();

//        for (unsigned int i = 0; i < mentity.atoms.size(); i++)
//        {
//          printf("atom %d --- (%8.4f, %8.4f, %8.4f)\n");
//        }


//        for (unsigned int i = 0; i < 40; i++)
//        {
//          compute_madelung_energy(world, mentity, 1 + i*0.025);
//        }
        compute_madelung_energy(world,mentity,1.22474487,100.0,100.0);


    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (char* s) {
        print(s);
        error("caught a string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    finalize();

    return 0;
}

