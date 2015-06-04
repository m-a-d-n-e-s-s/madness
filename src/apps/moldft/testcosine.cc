//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>
using namespace madness;

static const double_complex I(0,1);
static const double pi = constants::pi;
static const double twopi = 2.0*constants::pi;

static const double L = 5.0; // Unit cell size in au for neon
static const double alpha = 2.5;

static const double thresh = 1e-2;
static const double kwavelet = 4;
static const int truncate_mode = 0;

static const double kx=0.5*twopi/L, ky=0.5*twopi/L, kz=0.5*twopi/L;
//static const double kx=0, ky=0, kz=0;

typedef SeparatedConvolution<double,3> operatorT;
typedef SeparatedConvolution<double_complex,3> coperatorT;
typedef std::shared_ptr<operatorT> poperatorT;

template <typename Q>
class ExpFunctor: public FunctionFunctorInterface<Q,3> {
private:
    Q qx;
    Q qy;
    Q qz;
public:
    ExpFunctor(Q qx, Q qy, Q qz) : qx(qx), qy(qy), qz(qz) {}
    Q operator()(const coord_3d& x) const {
      return std::exp(qx*x[0] + qy*x[1] + qz*x[2]);
    }
};

template <typename Q>
class ExpFunctor3d: public FunctionFunctorInterface<Q,3> {
private:
    Q q0;
    Q q1;
    Q q2;
public:
    ExpFunctor3d(Q q0, Q q1, Q q2) : q0(q0), q1(q1), q2(q2) {}
    Q operator()(const coord_3d& x) const {
      return std::exp(q0*x[0])*std::exp(q1*x[1])*std::exp(q2*x[2]);
    }
};

class PWFunctor: public FunctionFunctorInterface<double_complex,3> {
private:
   std::vector<coord_3d> gvecs; 
   real_tensor coeffs;
   double_complex factor;

public:
    PWFunctor(const std::vector<coord_3d>& gvecs, const real_tensor& coeffs, const double& L)
     : gvecs(gvecs), coeffs(coeffs), factor(double_complex(0.0,2*pi/L)) {};
    double_complex operator()(const coord_3d& x) const {
      double_complex s = 0.0;
      int ngvecs = gvecs.size();
      for (int ig = 0; ig < ngvecs; ig++) {
        coord_3d gv = gvecs[ig];
        s += coeffs(ig)*std::exp(factor*(gv[0]*x[0]+gv[1]*x[1]+gv[2]*x[2]));
      }
      MADNESS_ASSERT(false);
      return s;
    }
};

class CosPotentialFunctor : public FunctionFunctorInterface<double,3> {
public:
  double operator()(const coord_3d& x) const {
    return -alpha*(cos(2*pi*x[0]/L)*cos(2*pi*x[1]/L)*cos(2*pi*x[2]/L) + 1.0);
  }
};

std::vector<coord_3d> get_coeffs_pw(double maxKlen) {
    std::vector<coord_3d> coeffs;
    int maxI = 10;
    for (int ii0 = -maxI; ii0 <= maxI; ii0++) {
      for (int ii1 = -maxI; ii1 <= maxI; ii1++) {
        for (int ii2 = -maxI; ii2 <= maxI; ii2++) {
          double kvecLen = std::sqrt((double)(ii0*ii0 + ii1*ii1 + ii2*ii2))*2*pi/L;
          if (kvecLen < maxKlen) {
              coeffs.push_back({(double)ii0, (double)ii1, (double)ii2});
          }
        }
      }
    }
    sort(coeffs.begin(), coeffs.end(), [](const coord_3d& a, const coord_3d& b) -> bool {
      return a.normf() < b.normf();
    });
    return coeffs;
}

real_tensor make_pw_matrix(const std::vector<coord_3d>& gvecs) {
  int ngvecs = gvecs.size();
  real_tensor Hcos(ngvecs,ngvecs);
  std::vector<coord_3d> qvecs(8);
  qvecs[0] = { 1.0, 1.0, 1.0};
  qvecs[1] = { 1.0, 1.0,-1.0};
  qvecs[2] = { 1.0,-1.0, 1.0};
  qvecs[3] = { 1.0,-1.0,-1.0};
  qvecs[4] = {-1.0, 1.0, 1.0};
  qvecs[5] = {-1.0, 1.0,-1.0};
  qvecs[6] = {-1.0,-1.0, 1.0};
  qvecs[7] = {-1.0,-1.0,-1.0};
  for (int ig1 = 0; ig1 < ngvecs; ig1++) {
    double s = (2*pi/L)*gvecs[ig1].normf();
    Hcos(ig1,ig1) = 0.5*s*s;
    Hcos(ig1,ig1) += -alpha;
    //print("ig1: ", ig1, "    s: ", s, "    0.5*s*s: ", 0.5*s*s, "Hcos: ", Hcos(ig1,ig1));
    for (int ig2 = 0; ig2 < ngvecs; ig2++) {
      for (int iq = 0; iq < 8; iq++) {
        coord_3d gvec = gvecs[ig2]-gvecs[ig1]+qvecs[iq]; 
        if (abs(gvec.normf()) < 1e-8)
          Hcos(ig1,ig2) += -0.125*alpha; 
      }
    }
  }
  return Hcos;
}

vector_complex_function_3d make_basis(World& world, double maxKlen) {
    vector_complex_function_3d orbs;
    int maxI = 8;
    for (int ii0 = -maxI; ii0 <= maxI; ii0++) {
      for (int ii1 = -maxI; ii1 <= maxI; ii1++) {
        for (int ii2 = -maxI; ii2 <= maxI; ii2++) {
          double kvecLen = std::sqrt((double)(ii0*ii0 + ii1*ii1 + ii2*ii2))*2*pi/L;
          if (kvecLen < maxKlen) {
            complex_function_3d orb = complex_factory_3d(world).functor(complex_functor_3d(
              new ExpFunctor<double_complex>(I*(double)ii0*twopi/L, I*(double)ii1*twopi/L, I*(double)ii2*twopi/L))).truncate_on_project();
            orbs.push_back(orb);
          }
        }
      }
    }
    normalize(world,orbs);
    return orbs;
}

tensor_complex make_kinetic_matrix(World& world, const vector_complex_function_3d& v) {
    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    vector_complex_function_3d dvx = apply(world, Dx, v);
    vector_complex_function_3d dvy = apply(world, Dy, v);
    vector_complex_function_3d dvz = apply(world, Dz, v);

    // -1/2 (del + ik)^2 = -1/2 del^2 - i k.del + 1/2 k^2
    // -1/2 <p|del^2|q> = +1/2 <del p | del q>

    tensor_complex f1 = 0.5*(matrix_inner(world, dvx, dvx, true) +
                             matrix_inner(world, dvy, dvy, true) +
                             matrix_inner(world, dvz, dvz, true));

    tensor_complex f2 =
        (-I*kx)*matrix_inner(world, v, dvx, true) +
        (-I*ky)*matrix_inner(world, v, dvy, true) +
        (-I*kz)*matrix_inner(world, v, dvz, true);

    tensor_complex f3 = (0.5 * (kx*kx + ky*ky + kz*kz)) * matrix_inner(world, v, v, true);

    return f1 + f2 + f3;
}

vector_complex_function_3d apply_potential(const real_function_3d& potential, const vector_complex_function_3d& psi)
{
    vector_complex_function_3d vpsi;
    for (unsigned int i=0; i<psi.size(); i++)
        vpsi.push_back(potential*psi[i]);
    return vpsi;
}


void orthogonalize(World& world, vector_complex_function_3d& psi) {
    compress(world, psi);
    for (unsigned int i = 0; i<psi.size(); i++) {
        complex_function_3d& psi_i = psi[i];
        psi_i.scale(1.0/psi_i.norm2());
        for (unsigned int j = 0; j<i; j++) {
            complex_function_3d& psi_j = psi[j];
            double_complex s = inner(psi_j,psi_i);
            psi_i.gaxpy(1.0,psi_j,-s); // |i> = |i> - |j><j|i>
            psi_i.scale(1.0/psi_i.norm2());
        }
    }
}

// function to apply BSH with twisted PBC
// kx, ky, kz -- some k value in the 1BZ (e.g. 0.5*2.0*pi/L where L is the lattice constant)
// energy     -- bound state energy (should be negative)
// L          -- lattice constant
complex_function_3d apply_periodic_bsh(World& world, const complex_function_3d& f, 
                                       const double& kx, const double& ky, const double& kz,
                                       const double& energy, const double& L) {
  complex_function_3d phase_p = complex_factory_3d(world).functor(complex_functor_3d(
    new ExpFunctor3d<double_complex>(I*kx,I*ky,I*kz))).truncate_mode(0).truncate_on_project();
  complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(
    new ExpFunctor3d<double_complex>(-I*kx,-I*ky,-I*kz))).truncate_mode(0).truncate_on_project();
  SeparatedConvolution<double_complex,3> op = 
    PeriodicBSHOperator3D(world, {-kx*L, -ky*L, -kz*L}, sqrt(-2.0*(energy)),  1e-4, FunctionDefaults<3>::get_thresh());
  complex_function_3d g = phase_m*apply(op, phase_p*f);
  return g;
}

// DESTROYS VPSI
vector_complex_function_3d update(World& world,
                                  const vector_complex_function_3d& psi,
                                  vector_complex_function_3d& vpsi,
                                  const tensor_real& e,
                                  int iter)
{
    // psi = - 2 G(E+shift) * (V+shift) psi
    int nmo = psi.size();

    // Append additional terms for periodic case to the potential
    // -ik.del + 1/2 k^2
    double ksq = kx*kx + ky*ky + kz*kz;
    coord_3d k {kx, ky, kz};

    // determine shift to make homo <=-0.1
    double shift = 0.0;
    if (e(nmo-1) > -0.1) {
        shift = -0.1 - e(nmo-1);
        gaxpy(world, 1.0, vpsi, shift, psi);
    }

    // Do the BSH thing
    scale(world, vpsi, -2.0);
    truncate(world, vpsi);
    
    vector_complex_function_3d new_psi(nmo);
    for (int iorb = 0; iorb < nmo; iorb++) {
        new_psi[iorb] = apply_periodic_bsh(world, vpsi[iorb], kx, ky, kz, e[iorb]+shift, L);
    }


    // Step restriction
    double damp;
    if (iter < 10) damp = 0.95;
    else if (iter < 20) damp = 0.85;
    else damp = 0.75;
    damp = 0.25;
    if (world.rank() == 0) print("  shift", shift, "damp", damp, "\n");

    if (world.rank() == 0) printf("      eigenvalue    residual\n");
    for (int i=0; i<nmo; i++) {
        double rnorm = (psi[i]-new_psi[i]).norm2();
        if (world.rank() == 0) printf("%4d  %10.6f  %10.1e\n", i, e[i], rnorm);
        new_psi[i] = damp*psi[i] + (1.0-damp)*new_psi[i];
    }
    truncate(world,new_psi);
    normalize(world, new_psi);
    orthogonalize(world, new_psi);
    truncate(world,new_psi);
    normalize(world, new_psi);
    return new_psi;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_k(kwavelet);
    FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<3>::set_truncate_mode(truncate_mode);

    // Nuclear potential
    real_function_3d v = real_factory_3d(world).functor(real_functor_3d(new CosPotentialFunctor())).truncate_mode(0).truncate_on_project();
    v.truncate();

    // Make basis
    int nmo = 7;
    vector_complex_function_3d psi = make_basis(world, 7.0);
    print("initial size :", psi.size());
    print("reprojecting");
    for (int iter=0; iter<100; iter++) {
        if (world.rank() == 0) print("\n\n  Iteration",iter,"\n");
        vector_complex_function_3d vpsi = apply_potential(v, psi);

        tensor_complex ke_mat = make_kinetic_matrix(world, psi);
        tensor_complex pe_mat = matrix_inner(world, psi, vpsi, true);
        tensor_complex ov_mat = matrix_inner(world, psi, psi, true);

        //print("KE"); print(real(ke_mat));
        //print("PE"); print(real(pe_mat));
        //print("H"); print(real(ke_mat+pe_mat));
        //print("OV"); print(ov_mat);

        tensor_complex fock = ke_mat + pe_mat;
        // eliminate small off-diagonal elements and lift diagonal
        // degeneracies to reduce random mixing
        for (unsigned int i=0; i<psi.size(); i++) {
            fock(i,i) += i*thresh*1e-2;
            for (unsigned int j=0; j<i; j++) {
                if (std::abs(fock(i,j)) < thresh*1e-1 || std::abs(ov_mat(i,j)) < thresh*1e-1) {
                    fock(i,j) = fock(j,i) = 0.0;
                    ov_mat(i,j) = ov_mat(j,i) = 0.0;
                }
            }
        }

        tensor_complex c;
        tensor_real e;
        sygv(fock, ov_mat, 1, c, e);
        //print("eigenvectors"); print(c);
        //print("eigenvalues"); print(e);

        if (iter == 0) {
            c = copy(c(_,Slice(0,nmo-1))); // truncate to occupied states
            e = e(Slice(0,nmo-1));
        }

        psi = transform(world, psi, c);
        vpsi = transform(world, vpsi, c);

        if (iter == 0) {
          print("reprojecting ..");
            v = madness::project(v, kwavelet+2, thresh*1e-2, true); 
          for (int i = 0; i < nmo; i++) {
            FunctionDefaults<3>::set_k(kwavelet+2);
            FunctionDefaults<3>::set_thresh(thresh*1e-2);
            psi[i] = madness::project(psi[i], kwavelet+2, thresh*1e-2, true); 
            vpsi[i] = madness::project(vpsi[i], kwavelet+2, thresh*1e-2, true); 
          }
          print("done reprojecting ..");
        } else if (iter == 10) {
          print("reprojecting ..");
            v = madness::project(v, kwavelet+4, thresh*1e-4, true); 
          for (int i = 0; i < nmo; i++) {
            FunctionDefaults<3>::set_k(kwavelet+4);
            FunctionDefaults<3>::set_thresh(thresh*1e-4);
            psi[i] = madness::project(psi[i], kwavelet+4, thresh*1e-4, true); 
            vpsi[i] = madness::project(vpsi[i], kwavelet+4, thresh*1e-4, true); 
          }
          print("done reprojecting ..");
        } else if (iter == 20) {
          print("reprojecting ..");
            v = madness::project(v, kwavelet+6, thresh*1e-6, true); 
          for (int i = 0; i < nmo; i++) {
            FunctionDefaults<3>::set_k(kwavelet+6);
            FunctionDefaults<3>::set_thresh(thresh*1e-6);
            psi[i] = madness::project(psi[i], kwavelet+6, thresh*1e-6, true); 
            vpsi[i] = madness::project(vpsi[i], kwavelet+6, thresh*1e-6, true); 
          }
          print("done reprojecting ..");
        }

        psi = update(world, psi, vpsi, e, iter);
    }
    return 0;
}
