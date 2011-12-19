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

  int rlo = -200;
  int rhi = 201;
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
//        if (world.rank() == 0)
//          printf("%10.5f %10.5f %10.5f   %10.5f\n",rvec[0],rvec[1],rvec[2],rlen);
  }

  unsigned int rsize = rvecs.size();
  Vector<double,3> rvec = rvecs[rsize-1];
  double maxRlen2 = std::sqrt(rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2]);
  print("R-max length requested:  ", maxRlen, "    R-max length:  ", maxRlen2);

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

  int glo = -200;
  int ghi = 201;
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

  unsigned int gsize = gvecs.size();
  Vector<double,3> gvec = gvecs[gsize-1];
  double maxGlen2 = std::sqrt(gvec[0]*gvec[0] + gvec[1]*gvec[1] + gvec[2]*gvec[2]);
  print("G-max length requested:  ", maxGlen, "    G-max length:  ", maxGlen2);
  return gvecs;

}
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

  mentity.print();

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
        Vector<double,3> tvec = vec(1*(-iatom.x),
                                    1*(-iatom.y),
                                    1*(-iatom.z));
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
          Vector<double,3> dvec = vec(L*(iatom.x-jatom.x),
                                      L*(iatom.y-jatom.y),
                                      L*(iatom.z-jatom.z));
          Vector<double,3> tvec = vec(dvec[0]+rvec[0],
                                      dvec[1]+rvec[1],
                                      dvec[2]+rvec[2]);
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
  printf("\n");
  double sqrtpi = std::sqrt(constants::pi);
  for (unsigned int ia = 0; ia < natoms; ia++)
  {
    Atom iatom = mentity.get_atom(ia);
//    printf("atom %d    charge: %8.4f\n",ia,iatom.q);
//    s3 += 2.0*iatom.q*iatom.q*alpha/sqrtpi;
    s3 += iatom.q*iatom.q*alpha*std::sqrt(8.0/TWOPI);
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
      Vector<double,3> tvec = vec(1*(iatom.x),
                                  1*(iatom.y),
                                  1*(iatom.z));
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
          Vector<double,3> dvec = vec(1*(iatom.x-jatom.x),
                                      1*(iatom.y-jatom.y),
                                      1*(iatom.z-jatom.z));
          Vector<double,3> tvec = vec(dvec[0]+rvec[0],
                                      dvec[1]+rvec[1],
                                      dvec[2]+rvec[2]);
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
        mentity.add_atom(0.0,0.0,0.0,7,9);
        mentity.add_atom(L/2,L/2,L/2,1,3);
        mentity.center();

        mentity.print();

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
        compute_madelung_energy_PWSCF(world,mentity,1.0,50.0,100.0);
        compute_madelung_energy_PWSCF(world,mentity,0.3,50.0,100.0);

//        for (unsigned i = 0; i < 20; i++)
//        {
//          compute_madelung_energy2(world, mentity, 0.004, 0.0+i*50.0, 200.0);
//        }

//        compute_madelung_energy2(world, mentity, 1.2, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.3, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.4, 100.0, 100.0);
//        compute_madelung_energy2(world, mentity, 1.5, 100.0, 100.0);

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

