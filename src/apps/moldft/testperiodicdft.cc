/*
  free space evals (k=6,tol=1e-4) from moldft
  -3.0306e+01 -1.3228e+00 -4.9800e-01 -4.9800e-01 -4.9800e-01

  computed by testperiodic with gamma point L=30.0
  -3.0304e+01 -1.3213e+00 -4.9782e-01 -4.9782e-01 -4.9782e-01

 */


#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>
#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/potentialmanager.h>
#include<madness/chem/xcfunctional.h>

using namespace madness;

#include "subspace.h"

static const double_complex I(0,1);
static const double twopi = 2.0*constants::pi;

//static const double L = 5.0; // Unit cell size in au for neon
//static const double L = 8.37; // Unit cell size in au for neon
//static const double L = 7.65; // Unit cell size in au for LiF
//static const double L = 3.8; // Unit cell size in au for LiF
//static const double L = 8.0;
//static const double L = 10.26085381075144364474; // Unit cell size in au for Si
// static const double L = 10.3235; // Unit cell size in au for CaF2
static const double L = 10.6591; // Unit cell size in au for NaCl

static const int maxR = 3; // periodic sums from -R to +R inclusive
static const double thresh = 1e-5;
static const double kwavelet = 14;
static const int truncate_mode = 0;

static Molecule molecule;
static AtomicBasisSet aobasis;
static Subspace* subspace;

typedef SeparatedConvolution<double,3> operatorT;
typedef SeparatedConvolution<double_complex,3> coperatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef std::shared_ptr<coperatorT> pcoperatorT;
typedef std::pair<vector_complex_function_3d,vector_complex_function_3d> vcpairT;

static double ttt, sss;
static void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

static void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}

class SplitterFunctor: public FunctionFunctorInterface<double,3> {
public:
    double operator()(const coord_3d& r) const {
        return 2.0*r[2]*r[2] + r[0]*r[0];
    }
};

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

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
    const int maxR = 2;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis)
        : molecule(molecule), aobasis(aobasis) {}

    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        for (int i=-maxR; i<=+maxR; i++) {
            for (int j=-maxR; j<=+maxR; j++) {
                for (int k=-maxR; k<=+maxR; k++) {
                    sum += aobasis.eval_guess_density(molecule, x[0]+i*L, x[1]+j*L, x[2]+k*L);
                }
            }
        }
        return sum;
    }
};

class AtomicBasisFunctor : public FunctionFunctorInterface<double_complex,3> {
private:
    const AtomicBasisFunction aofunc;
    double kx, ky, kz;
    std::vector<coord_3d> specialpt;
public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc, double kx, double ky, double kz)
        : aofunc(aofunc), kx(kx), ky(ky), kz(kz)
    {
    	double x, y, z;
        aofunc.get_coords(x,y,z);
        coord_3d r;
        r[0]=x; r[1]=y; r[2]=z;
        specialpt=std::vector<coord_3d>(1,r);
    }

    double_complex operator()(const coord_3d& x) const {
        double_complex sum = 0.0;
        for (int i=-maxR; i<=+maxR; i++) {
            double xx = x[0]+i*L;
            double rsq = xx*xx;
            if (rsq < aofunc.rangesq()) {
                for (int j=-maxR; j<=+maxR; j++) {
                    double yy = x[1]+j*L;
                    rsq += yy*yy;
                    if (rsq < aofunc.rangesq()) {
                        for (int k=-maxR; k<=+maxR; k++) {
                            double zz = x[2]+k*L;
                            rsq += zz*zz;
                            if (rsq < aofunc.rangesq()) {
                                sum += exp(-I*(kx*xx+ky*yy+kz*zz))*aofunc(xx, yy, zz);
                            }
                        }
                    }
                }
            }
        }
        return sum;
    }

    std::vector<coord_3d> special_points() const {return specialpt;}
};

class AtomicOrbitalFunctor: public FunctionFunctorInterface<double,3> {
private:
    const AtomicBasisFunction aofunc;
    double R;
    int offx, offy, offz;
    std::vector<coord_3d> specialpt;
    double rangesq;
public:
    AtomicOrbitalFunctor(const AtomicBasisFunction& aofunc, double R, int offx, int offy, int offz)
        : aofunc(aofunc), R(R), offx(offx), offy(offy), offz(offz)
    {
        double x, y, z;
        aofunc.get_coords(x,y,z);
        coord_3d r;
        r[0]=x+offx*R; r[1]=y+offy*R; r[2]=z+offz*R;
        specialpt=std::vector<coord_3d>(1,r);
        rangesq = 3.0*aofunc.rangesq();
    }

    double operator()(const coord_3d& x) const {
        return aofunc(x[0] + offx*R, x[1] + offy*R, x[2] + offz*R);
    }

    std::vector<coord_3d> special_points() const {return specialpt;}
};

class KPeriodicBSHOperator {
private:
    double kx, ky, kz;
    double L;
    complex_function_3d phase_p;
    complex_function_3d phase_m;

public:
    KPeriodicBSHOperator(World& world, const double& kx, const double& ky, const double& kz, const double& L)
     : kx(kx), ky(ky), kz(kz), L(L) {
        phase_p = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(I*kx,I*ky,I*kz))).truncate_mode(0).truncate_on_project();
        phase_m = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(-I*kx,-I*ky,-I*kz))).truncate_mode(0).truncate_on_project();
    }

    KPeriodicBSHOperator(World& world, const Vector<double,3>& kpt, const double& L)
     : kx(kpt[0]), ky(kpt[1]), kz(kpt[2]), L(L) {
        phase_p = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(I*kx,I*ky,I*kz))).truncate_mode(0).truncate_on_project();
        phase_m = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(-I*kx,-I*ky,-I*kz))).truncate_mode(0).truncate_on_project();
    }

    vector_complex_function_3d apply(World& world, const vector_complex_function_3d& v, const tensor_real& evals, double shift = 0.0) {
        START_TIMER(world);
        int nmo = evals.dim(0);
        vector<pcoperatorT> ops(nmo);
        for(int i = 0;i < nmo; ++i){
            double eps = evals(i) + shift;
            ops[i] = pcoperatorT(PeriodicBSHOperatorPtr3D(world, vec(-kx*L, -ky*L, -kz*L), sqrt(-2.0*eps),  1e-4, FunctionDefaults<3>::get_thresh()));
        }
        vector_complex_function_3d t1 = mul(world,phase_p,v);
        truncate(world,t1);
        vector_complex_function_3d t2 = ::apply(world,ops,t1);
        vector_complex_function_3d t3 = mul(world,phase_m,t2);
        return t3;
        END_TIMER(world, "apply periodic bsh");
        throw "control reaches end of non-void function";
        return t2; // T2 is random choice ... WHAT SHOULD IT REALLY BE??????????????????????????????
    }
};

vector_complex_function_3d makeao_slow(World& world, const std::vector<Vector<double,3> >& kpoints) {
    int nkpt = kpoints.size();
    int nbf = aobasis.nbf(molecule);
    vector_complex_function_3d ao(nkpt*nbf);
    for (int ik = 0; ik < nkpt; ++ik) {
        Vector<double,3> kvec = kpoints[ik];
        double kx = kvec[0]; double ky = kvec[1]; double kz = kvec[2];
        for(int i = 0; i<aobasis.nbf(molecule); ++i) {
            aobasis.get_atomic_basis_function(molecule, i).print_me(std::cout);
            complex_functor_3d aofunc(
               new AtomicBasisFunctor(
               aobasis.get_atomic_basis_function(molecule, i), kx, ky, kz));
            ao[ik*nbf+i] = complex_factory_3d(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(0);
        }
    }
    world.gop.fence();
    normalize(world,ao);
    truncate(world,ao);
    return ao;
}

vector_complex_function_3d makeao(World& world, const std::vector<Vector<double,3> >& kpts, double R) {
    START_TIMER(world);
    int nkpt = kpts.size();
    int nbf = aobasis.nbf(molecule);
    double t1 = 1./std::sqrt((double)nbf);
    vector_complex_function_3d aos = zero_functions_compressed<double_complex,3>(world, nkpt*nbf);
    if (world.rank() == 0) print("nbf:  ", nbf);
    for (int ibf = 0; ibf < nbf; ibf++) {
        if (world.rank() == 0) print("\n\nibf: ", ibf);
        for (int i = -maxR; i <= maxR; i++) {
            for (int j = -maxR; j <= maxR; j++) {
                for (int k = -maxR; k <= maxR; k++) {
                    //AtomicBasisFunction abf = aobasis.get_atomic_basis_function(molecule, ibf); 
                    //if ((i*i+j*j+k*k)*L*L < 3.0*abf.rangesq()) {
                    if (true) {
                        real_functor_3d aofunc(
                           new AtomicOrbitalFunctor(
                           aobasis.get_atomic_basis_function(molecule, ibf), L, i, j, k));
                        real_function_3d ao = real_factory_3d(world).functor(aofunc).truncate_on_project().truncate_mode(0);
                        ao.compress();
                        for (int ik = 0; ik < nkpt; ik++) {
                            Vector<double,3> kpt = kpts[ik];
                            complex_function_3d t2 = t1*std::exp(double_complex(0.0,(i*kpt[0]+j*kpt[1]+k*kpt[2])*R))*ao;
                            aos[ik*nbf+ibf].gaxpy(1.0,t2,1.0);
                        }
                    }
                }
            }
        }
    }
    for (int ik = 0; ik < nkpt; ik++) {
        Vector<double,3> kpt = kpts[ik];
        complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(-I*kpt[0],-I*kpt[1],-I*kpt[2]))).truncate_mode(0).truncate_on_project();
        for (int ibf = 0; ibf < nbf; ibf++) {
            aos[ik*nbf+ibf] = aos[ik*nbf+ibf]*phase_m;
        }
    }
    END_TIMER(world, "makeao");
    return aos;
}

//vector_complex_function_3d makeao(World& world, const std::vector<Vector<double,3> >& kpts, double R) {
//    START_TIMER(world);
//    int nkpt = kpts.size();
//    int nbf = aobasis.nbf(molecule);
//    double t1 = 1./std::sqrt((double)nbf);
//    vector_complex_function_3d ao = zero_functions_compressed<double_complex,3>(world, std:pow(2*maxR+1, 3)*nbf);
//    vector_complex_function_3d aos = zero_functions_compressed<double_complex,3>(world, nkpt*nbf);
//    print("nbf:  ", nbf);
//    for (int ibf = 0; ibf < nbf; ibf++) {
//        print("\n\nibf: ", ibf);
//        for (int i = -maxR; i <= maxR; i++) {
//            for (int j = -maxR; j <= maxR; j++) {
//                for (int k = -maxR; k <= maxR; k++) {
//                    AtomicBasisFunction abf = aobasis.get_atomic_basis_function(molecule, ibf); 
//                    //if ((i*i+j*j+k*k)*L*L < 3.0*abf.rangesq()) {
//                    if (true) {
//                        real_functor_3d aofunc(
//                           new AtomicOrbitalFunctor(
//                           aobasis.get_atomic_basis_function(molecule, ibf), L, i, j, k));
//                        ao[i = real_factory_3d(world).functor(aofunc).truncate_on_project().truncate_mode(0);
//                    }
//                }
//            }
//        }
//    }
//    for (int ik = 0; ik < nkpt; ik++) {
//        Vector<double,3> kpt = kpts[ik];
//        complex_function_3d t2 = t1*std::exp(double_complex(0.0,(i*kpt[0]+j*kpt[1]+k*kpt[2])*R))*ao;
//        aos[ik*nbf+ibf].gaxpy(1.0,t2,1.0);
//    }
//    for (int ik = 0; ik < nkpt; ik++) {
//        Vector<double,3> kpt = kpts[ik];
//        complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(
//          new ExpFunctor3d<double_complex>(-I*kpt[0],-I*kpt[1],-I*kpt[2]))).truncate_mode(0).truncate_on_project();
//        for (int ibf = 0; ibf < nbf; ibf++) {
//            aos[ik*nbf+ibf] = aos[ik*nbf+ibf]*phase_m;
//        }
//    }
//    END_TIMER(world, "makeao");
//    return aos;
//}

tensor_complex make_kinetic_matrix(World& world, const vector_complex_function_3d& v, const Vector<double,3>& kpt) {
    START_TIMER(world);
    double kx = kpt[0]; double ky = kpt[1]; double kz = kpt[2];
    complex_derivative_3d Dx(world, 0);
    complex_derivative_3d Dy(world, 1);
    complex_derivative_3d Dz(world, 2);

    auto dvx = apply(world, Dx, v);
    auto dvy = apply(world, Dy, v);
    auto dvz = apply(world, Dz, v);

    // -1/2 (del + ik)^2 = -1/2 del^2 - i k.del + 1/2 k^2
    // -1/2 <p|del^2|q> = +1/2 <del p | del q>

    auto f1 = 0.5 * (matrix_inner(world, dvx, dvx, false) +
                     matrix_inner(world, dvy, dvy, false) +
                     matrix_inner(world, dvz, dvz, false));

    auto f2 =
        (-I*kx)*matrix_inner(world, v, dvx, false) +
        (-I*ky)*matrix_inner(world, v, dvy, false) +
        (-I*kz)*matrix_inner(world, v, dvz, false);

    auto f3 = (0.5 * (kx*kx + ky*ky + kz*kz)) * matrix_inner(world, v, v, true);
    END_TIMER(world, "kinetic energy");

    return f1 + f2 + f3;
}

vector_complex_function_3d apply_potential(World& world, const real_function_3d& potential, const vector_complex_function_3d& psi)
{
    START_TIMER(world);
    auto vpsi = mul(world, potential, psi, false);
    world.gop.fence();
    END_TIMER(world, "apply potential");
    return vpsi;
}

vector_complex_function_3d orth(World& world, const vector_complex_function_3d& v, double thresh = 1.e-10) {
    auto vsize = v.size();
    auto ov_mat = matrix_inner(world, v, v, true);
    for (unsigned int i=0; i<v.size(); i++) {
        for (unsigned int j=0; j<i; j++) {
            if (std::abs(ov_mat(i,j)) < thresh*1e-1) {
                ov_mat(i,j) = ov_mat(j,i) = 0.0;
            }
        }
    }
    tensor_complex U;
    tensor_real D;
    syev(ov_mat, U, D);
    print("D: ");
    print(D);
    auto indx = -1;
    for (unsigned int i = 0; i < vsize && indx < 0; i++) {
      if (std::abs(D(i)) > thresh) {
        indx = i;
      }
    }
    U = copy(U(_,Slice(indx,vsize-1)));
    auto R = transform(world, v, U);
    normalize(world, R, true);
    ov_mat = matrix_inner(world, R, R, true);
    print("new overlap: ");
    print(ov_mat);
    return R;
}

std::pair<vector_complex_function_3d, tensor_real> diag_and_transform(World& world, 
                                                                      const Vector<double,3> kpt, 
                                                                      const real_function_3d& v, 
                                                                      const vector_complex_function_3d& psik, 
                                                                      int nmo = 0) {
    auto vpsik = apply_potential(world, v, psik); 
    auto ke_mat = make_kinetic_matrix(world, psik, kpt);
    auto pe_mat = matrix_inner(world, psik, vpsik, true);
    auto ov_mat = matrix_inner(world, psik, psik, true);
    
    tensor_complex fock = ke_mat + pe_mat;
    // eliminate small off-diagonal elements and lift diagonal
    // degeneracies to reduce random mixing
    for (unsigned int i=0; i<psik.size(); i++) {
        fock(i,i) += i*thresh*1e-2;
        for (unsigned int j=0; j<i; j++) {
            if (std::abs(fock(i,j)) < thresh*1e-1 || std::abs(ov_mat(i,j)) < thresh*1e-1) {
                fock(i,j) = fock(j,i) = 0.0;
                ov_mat(i,j) = ov_mat(j,i) = 0.0;
            }
        }
    }

    // print("H:\n"); print(fock);
    // print("S:\n"); print(ov_mat);
    tensor_complex c;
    tensor_real e;
    sygv(fock, ov_mat, 1, c, e);

    if (nmo > 0) {
        c = copy(c(_,Slice(0,nmo-1))); // truncate to occupied states
        e = e(Slice(0,nmo-1));
    }

    auto new_psik = transform(world, psik, c);
    return std::pair<vector_complex_function_3d, tensor_real>(new_psik, e);
}
 
real_function_3d make_lda_potential(World& world, const real_function_3d &rho)
{
    START_TIMER(world);
    auto vlda = copy(rho);
    vlda.reconstruct();
    vlda.unaryop(xc_lda_potential());
    END_TIMER(world, "lda potential");
    return vlda;
}

real_function_3d make_coulomb_potential(World& world, const real_function_3d& rho)
{
    START_TIMER(world);
    real_convolution_3d op = CoulombOperator(world, 1e-4, thresh);
    END_TIMER(world, "hartree potential");
    return op(rho);
}

vector<poperatorT> make_bsh_operators(World & world, const tensor_real& evals, double shift)
{
    int nmo = evals.dim(0);
    vector<poperatorT> ops(nmo);
    for(int i = 0;i < nmo; ++i){
        double eps = evals(i) + shift;
        ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  1e-4, thresh));
    }
    return ops;
}

void orthogonalize(World& world, vector_complex_function_3d& psi) {
    START_TIMER(world);
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
    END_TIMER(world, "orthogonalize");
}

// function to apply BSH with twisted PBC
// kx, ky, kz -- some k value in the 1BZ (e.g. 0.5*2.0*pi/L where L is the lattice constant)
// energy     -- bound state energy (should be negative)
// L          -- lattice constant
//
// Obviously this is slow
complex_function_3d apply_periodic_bsh(World& world, const complex_function_3d& f, 
                                       const double& kx, const double& ky, const double& kz,
                                       const double& energy, const double& L) {
  complex_function_3d phase_p = complex_factory_3d(world).functor(complex_functor_3d(
    new ExpFunctor3d<double_complex>(I*kx,I*ky,I*kz))).truncate_mode(0).truncate_on_project();
  complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(
    new ExpFunctor3d<double_complex>(-I*kx,-I*ky,-I*kz))).truncate_mode(0).truncate_on_project();
  auto op = PeriodicBSHOperator3D(world, vec(-kx*L, -ky*L, -kz*L), sqrt(-2.0*(energy)),  1e-4, FunctionDefaults<3>::get_thresh());
  complex_function_3d g = phase_m*madness::apply(op, phase_p*f);
  return g;
}

// function to apply BSH with twisted PBC
complex_function_3d apply_periodic_bsh(World& world, const complex_function_3d& f, 
                                       const Vector<double,3>& kpt, 
                                       const double& energy, 
                                       const double& L) {
    return apply_periodic_bsh(world,f,kpt[0],kpt[1],kpt[2],energy,L);
}

// DESTROYS VPSI
vcpairT update(World& world,
               int ik,
               const vector_complex_function_3d& psi,
               vector_complex_function_3d& vpsi,
               const Vector<double,3>& kpt,
               const real_function_3d& v,
               const tensor_real& e)
{
    int nmo = psi.size();
    double kx = kpt[0]; double ky = kpt[1]; double kz = kpt[2];

    // determine shift to make homo <=-0.1
    double shift = 0.0;
    if (e(nmo-1) > -0.1) {
        shift = -0.1 - e(nmo-1);
        gaxpy(world, 1.0, vpsi, shift, psi);
    }

    // Do the BSH thing
    scale(world, vpsi, -2.0);
    //truncate(world, vpsi);
    
    //vector_complex_function_3d new_psi(nmo);
    //for (int iorb = 0; iorb < nmo; iorb++) {
    //  new_psi[iorb] = apply_periodic_bsh(world, vpsi[iorb], kx, ky, kz, e[iorb]+shift, L);
    //}
   
    KPeriodicBSHOperator kop(world, kx, ky, kz, L);
    vector_complex_function_3d new_psi = kop.apply(world, vpsi, e, shift);
    vector_complex_function_3d rm = sub(world, psi, new_psi);
    truncate(world,rm);

    if (world.rank() == 0) printf("kpoint:  %10.5f    %10.5f    %10.5f\n",kpt[0],kpt[1],kpt[2]);
    if (world.rank() == 0) printf("      eigenvalue    residual\n");
    for (int i=0; i<nmo; i++) {
        double rnorm = rm[i].norm2();
        if (world.rank() == 0) printf("%4d  %10.6f  %10.1e\n", i, e[i], rnorm);
    }
    return vcpairT(new_psi,rm);
}

real_function_3d make_density(World& world, const vector_complex_function_3d& v, double weight) {
    START_TIMER(world);
    real_function_3d rho(world);
    for (unsigned int i=0; i<v.size(); i++) {
        rho = rho + weight*abssq(v[i]);
    }
    rho.scale(2.0); // total closed-shell density
    rho.truncate();
    END_TIMER(world, "make density");
    return rho;
}

// Return all orbitals for all kpoints
// Modifies the density rho
vector_complex_function_3d initial_guess(World& world, const real_function_3d& vnuc, real_function_3d& rho, const std::vector<Vector<double,3> >& kpoints, int nst) {
    auto psi0 = makeao(world, kpoints, L);
    auto v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
    auto vpsi = apply_potential(world, v, psi0);

    int nkpt = kpoints.size();
    int nsize = psi0.size();
    int nst_initial = nsize/nkpt;
    MADNESS_CHECK(nst <= nst_initial);
    if (world.rank() == 0) print("nsize: ", nsize);
    if (world.rank() == 0) print("nst_initial: ", nst_initial);
    vector_complex_function_3d psi;
    for (int ik = 0; ik < nkpt; ++ik) {
        vector_complex_function_3d psik(psi0.begin()+ik*nst_initial, psi0.begin()+(ik+1)*nst_initial);
        vector_complex_function_3d vpsik(vpsi.begin()+ik*nst_initial, vpsi.begin()+(ik+1)*nst_initial);

        auto ke_mat = make_kinetic_matrix(world, psik, kpoints[ik]);
        auto pe_mat = matrix_inner(world, psik, vpsik, true);
        auto ov_mat = matrix_inner(world, psik, psik, true);
    
        auto fock = ke_mat + pe_mat;
        // eliminate small off-diagonal elements and lift diagonal
        // degeneracies to reduce random mixing
        for (unsigned int i=0; i<psik.size(); i++) {
            fock(i,i) += i*thresh*1e-2;
            for (unsigned int j=0; j<i; j++) {
                if (std::abs(fock(i,j)) < thresh*1e-1 && std::abs(ov_mat(i,j)) < thresh*1e-1) {
                    fock(i,j) = fock(j,i) = 0.0;
                    ov_mat(i,j) = ov_mat(j,i) = 0.0;
                }
            }
        }

        tensor_complex c;
        tensor_real e;
        fock = 0.5*(fock + transpose(fock));
        ov_mat = 0.5*(ov_mat + transpose(ov_mat));

        // print("Initial Fock: "); print(real(fock));
        // print("Initial Overlap: "); print(real(ov_mat));
        sygv(fock, ov_mat, 1, c, e);

        if (world.rank() == 0) print("initial_guess() ik = ", ik);
        if (world.rank() == 0) print(e);
    
        psik = transform(world, psik, c);
        vpsik = transform(world, vpsik, c);
        for (int ist = 0; ist < nst; ist++) {
            if (world.rank() == 0) print("pushing back ist = ", ist);
            psi.push_back(psik[ist]);    
        }
    }

    // compute new density
    double weights = 1.0/(double)kpoints.size();
    rho = make_density(world,psi,weights);
    return psi;
}

tensor_complex matrix_exponential(const tensor_complex& A) {
    const double tol = 1e-13;
    MADNESS_CHECK(A.dim(0) == A.dim(1));

    // Scale A by a power of 2 until it is "small"
    double anorm = A.normf();
    int n = 0;
    double scale = 1.0;
    while (anorm*scale > 0.1)
    {
        n++;
        scale *= 0.5;
    }
    tensor_complex B = scale*A;    // B = A*2^-n

    // Compute exp(B) using Taylor series
    tensor_complex expB = tensor_complex(2, B.dims());
    for (int i = 0; i < expB.dim(0); i++) expB(i,i) = double_complex(1.0,0.0);

    int k = 1;
    tensor_complex term = B;
    while (term.normf() > tol)
    {
        expB += term;
        term = inner(term,B);
        k++;
        term.scale(1.0/k);
    }

    // Repeatedly square to recover exp(A)
    while (n--)
    {
        expB = inner(expB,expB);
    }

    return expB;
}

void fixphases(World& world, tensor_real& e, tensor_complex& U) {
    int nmo = U.dim(0);
    long imax;
    for (long j = 0; j < nmo; j++)
    {
        // Get index of largest value in column
        U(_,j).absmax(&imax);
        double_complex ang = arg(U(imax,j));
        double_complex phase = std::exp(-ang*I);
        // Loop through the rest of the column and divide by the phase
        for (long i = 0; i < nmo; i++)
        {
            U(i,j) *= phase;
        }
    }

    // Within blocks with the same occupation number attempt to
    // keep orbitals in the same order (to avoid confusing the
    // non-linear solver).  Have to run the reordering multiple
    // times to handle multiple degeneracies.
    int maxpass = 15;
    for (int pass = 0; pass < maxpass; pass++)
    {
        long j;
        for (long i = 0; i < nmo; i++)
        {
            U(_, i).absmax(&j);
            if (i != j)
            {
              tensor_complex tmp = copy(U(_, i));
              U(_, i) = U(_, j);
              U(_, j) = tmp;
              //swap(e[i], e[j]);
              double ti = e[i];
              double tj = e[j];
              e[i] = tj; e[j] = ti;
            }
        }
    }

    // Rotations between effectively degenerate states confound
    // the non-linear equation solver ... undo these rotations
    long ilo = 0; // first element of cluster
    while (ilo < nmo-1) {
        long ihi = ilo;
        while (fabs(e[ilo]-e[ihi+1]) < thresh*700.0*std::max(fabs(e[ilo]),1.0)) {
            ihi++;
            if (ihi == nmo-1) break;
        }
        long nclus = ihi - ilo + 1;
        if (nclus > 1) {
            if (world.rank() == 0) print("   found cluster", ilo, ihi, e[ilo]);
            tensor_complex q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
            //print(q);
            // Special code just for nclus=2
            // double c = 0.5*(q(0,0) + q(1,1));
            // double s = 0.5*(q(0,1) - q(1,0));
            // double r = sqrt(c*c + s*s);
            // c /= r;
            // s /= r;
            // q(0,0) = q(1,1) = c;
            // q(0,1) = -s;
            // q(1,0) = s;

            // Iteratively construct unitary rotation by
            // exponentiating the antisymmetric part of the matrix
            // ... is quadratically convergent so just do 3
            // iterations
            tensor_complex rot = matrix_exponential(-0.5*(q - conj_transpose(q)));
            q = inner(q,rot);
            tensor_complex rot2 = matrix_exponential(-0.5*(q - conj_transpose(q)));
            q = inner(q,rot2);
            tensor_complex rot3 = matrix_exponential(-0.5*(q - conj_transpose(q)));
            q = inner(rot,inner(rot2,rot3));
            U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
        }
        ilo = ihi+1;
    }

}
void fixphases(World& world, tensor_real& e, tensor_complex& U, vector_complex_function_3d& psik) {
    int nmo = U.dim(0);
    long imax;
    for (long j = 0; j < nmo; j++)
    {
        // Get index of largest value in column
        U(_,j).absmax(&imax);
        double_complex ang = arg(U(imax,j));
        double_complex phase = std::exp(-ang*I);
        // Loop through the rest of the column and divide by the phase
        for (long i = 0; i < nmo; i++)
        {
            U(i,j) *= phase;
        }
    }

    //// Within blocks with the same occupation number attempt to
    //// keep orbitals in the same order (to avoid confusing the
    //// non-linear solver).  Have to run the reordering multiple
    //// times to handle multiple degeneracies.
    //int maxpass = 5;
    //for (int pass = 0; pass < maxpass; pass++)
    //{
    //    long j;
    //    for (long i = 0; i < nmo; i++)
    //    {
    //        U(_, i).absmax(&j);
    //        if (i != j)
    //        {
    //          tensor_complex tmp = copy(U(_, i));
    //          U(_, i) = U(_, j);
    //          U(_, j) = tmp;
    //          //swap(e[i], e[j]);
    //          double ti = e[i];
    //          double tj = e[j];
    //          e[i] = tj; e[j] = ti;
    //        }
    //    }
    //}

    // Rotations between effectively degenerate states confound
    // the non-linear equation solver ... undo these rotations
    long ilo = 0; // first element of cluster
    while (ilo < nmo-1) {
        long ihi = ilo;
        while (fabs(e[ilo]-e[ihi+1]) < thresh*10.0*std::max(fabs(e[ilo]),1.0)) {
            ihi++;
            if (ihi == nmo-1) break;
        }
        long nclus = ihi - ilo + 1;
//        if (nclus > 1) {
//            //if (world.rank() == 0) print("   found cluster", ilo, ihi);
//            tensor_complex q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
//            //print(q);
//            // Special code just for nclus=2
//            // double c = 0.5*(q(0,0) + q(1,1));
//            // double s = 0.5*(q(0,1) - q(1,0));
//            // double r = sqrt(c*c + s*s);
//            // c /= r;
//            // s /= r;
//            // q(0,0) = q(1,1) = c;
//            // q(0,1) = -s;
//            // q(1,0) = s;
//
//            // Iteratively construct unitary rotation by
//            // exponentiating the antisymmetric part of the matrix
//            // ... is quadratically convergent so just do 3
//            // iterations
//            tensor_complex rot = matrix_exponential(-0.5*(q - conj_transpose(q)));
//            q = inner(q,rot);
//            tensor_complex rot2 = matrix_exponential(-0.5*(q - conj_transpose(q)));
//            q = inner(q,rot2);
//            tensor_complex rot3 = matrix_exponential(-0.5*(q - conj_transpose(q)));
//            q = inner(rot,inner(rot2,rot3));
//            U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
//        }
        if (nclus > 1) {
            real_function_3d splitter = real_factory_3d(world).functor(real_functor_3d(
              new SplitterFunctor())).truncate_mode(0).truncate_on_project();
            vector_complex_function_3d psik_cluster(psik.begin()+ilo, psik.begin()+ihi+1);
            vector_complex_function_3d sfuncs = mul(world,splitter,psik_cluster);
            tensor_complex Mcluster = matrix_inner(world,psik_cluster,sfuncs);
            tensor_real eigs; tensor_complex uvecs;
            syev(Mcluster, uvecs, eigs);
            if (world.rank() == 0) {
                print("Found cluster of size: ", nclus, "with energy = ", e[ilo]);
                // print("cluster matrix:"); print(Mcluster);
                // print("cluster eigs:"); print(eigs);
                // print("cluster eigenvectors:"); print(uvecs);
            }
            psik_cluster = transform(world, psik_cluster, uvecs);
            for (int ist = 0; ist < nclus; ist++) {
                psik[ilo+ist] = psik_cluster[ist];
            }
        }
        ilo = ihi+1;
    }

}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_k(kwavelet);
    FunctionDefaults<3>::set_bc(BoundaryConditions<3>(BC_PERIODIC));
    FunctionDefaults<3>::set_cubic_cell(0,L);
    FunctionDefaults<3>::set_truncate_mode(truncate_mode);

    // kpoint list
    //int nkpt = 4;
    //std::vector<Vector<double,3> > kpoints(nkpt);
    //kpoints[0] = vec(0.0*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    //kpoints[1] = vec(0.0*twopi/L, 0.0*twopi/L, 0.5*twopi/L); 
    //kpoints[2] = vec(0.0*twopi/L, 0.5*twopi/L, 0.0*twopi/L); 
    //kpoints[3] = vec(0.0*twopi/L, 0.5*twopi/L, 0.5*twopi/L); 
    //double weight = 1.0/(double)nkpt;
    
    int nkpt = 8;
    std::vector<Vector<double,3> > kpoints(nkpt);
    kpoints[0] = vec(0.0*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    kpoints[1] = vec(0.0*twopi/L, 0.0*twopi/L, 0.5*twopi/L); 
    kpoints[2] = vec(0.0*twopi/L, 0.5*twopi/L, 0.0*twopi/L); 
    kpoints[3] = vec(0.0*twopi/L, 0.5*twopi/L, 0.5*twopi/L); 
    kpoints[4] = vec(0.5*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    kpoints[5] = vec(0.5*twopi/L, 0.0*twopi/L, 0.5*twopi/L); 
    kpoints[6] = vec(0.5*twopi/L, 0.5*twopi/L, 0.0*twopi/L); 
    kpoints[7] = vec(0.5*twopi/L, 0.5*twopi/L, 0.5*twopi/L); 
    double weight = 1.0/(double)nkpt;

    //int nkpt = 2;
    //std::vector<Vector<double,3> > kpoints(nkpt);
    //kpoints[0] = vec(0.0*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    //kpoints[1] = vec(0.0*twopi/L, 0.0*twopi/L, 0.5*twopi/L); 
    //double weight = 1.0/(double)nkpt;
    
    //int nkpt = 1;
    //std::vector<Vector<double,3> > kpoints(nkpt);
    ////kpoints[0] = vec(0.0*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    //kpoints[0] = vec(0.5*twopi/L, 0.5*twopi/L, 0.5*twopi/L); 
    //double weight = 1.0/(double)nkpt;
    
    // initialize subspace
    subspace = new Subspace[nkpt];

    // // FCC unit cell for ne
    molecule.add_atom(  0,  0,  0, 10.0, 10);
     molecule.add_atom(L/2,L/2,  0, 10.0, 10);
     molecule.add_atom(L/2,  0,L/2, 10.0, 10);
     molecule.add_atom(  0,L/2,L/2, 10.0, 10);

    // Cubic cell for LiF
    // molecule.add_atom(  0,  0,  0, 9.0, 9);
    // molecule.add_atom(L/2,L/2,  0, 9.0, 9);
    // molecule.add_atom(L/2,  0,L/2, 9.0, 9);
    // molecule.add_atom(  0,L/2,L/2, 9.0, 9);
    // molecule.add_atom(L/2,  0,  0, 3.0, 3);
    // molecule.add_atom(  0,L/2,  0, 3.0, 3);
    // molecule.add_atom(  0,  0,L/2, 3.0, 3);
    // molecule.add_atom(L/2,L/2,L/2, 3.0, 3);

    // Cubic cell for CaF2
    // molecule.add_atom(0.00*L, 0.00*L, 0.00*L, 20.0, 20);
    // molecule.add_atom(0.50*L, 0.50*L, 0.00*L, 20.0, 20);
    // molecule.add_atom(0.50*L, 0.00*L, 0.50*L, 20.0, 20);
    // molecule.add_atom(0.00*L, 0.50*L, 0.50*L, 20.0, 20);
    // molecule.add_atom(0.25*L, 0.25*L, 0.25*L, 9.0, 9); 
    // molecule.add_atom(0.75*L, 0.75*L, 0.75*L, 9.0, 9);
    // molecule.add_atom(0.75*L, 0.75*L, 0.25*L, 9.0, 9);
    // molecule.add_atom(0.25*L, 0.25*L, 0.75*L, 9.0, 9);
    // molecule.add_atom(0.75*L, 0.25*L, 0.75*L, 9.0, 9);
    // molecule.add_atom(0.25*L, 0.75*L, 0.25*L, 9.0, 9);
    // molecule.add_atom(0.25*L, 0.75*L, 0.75*L, 9.0, 9);
    // molecule.add_atom(0.75*L, 0.25*L, 0.25*L, 9.0, 9);    

    // Cubic cell for NaCl
//    molecule.add_atom(0.0*L, 0.0*L, 0.0*L, 11.0, 11);
//    molecule.add_atom(0.0*L, 0.5*L, 0.5*L, 11.0, 11);
//    molecule.add_atom(0.5*L, 0.0*L, 0.5*L, 11.0, 11);
//    molecule.add_atom(0.5*L, 0.5*L, 0.0*L, 11.0, 11);
//    molecule.add_atom(0.5*L, 0.5*L, 0.5*L, 17.0, 17);
//    molecule.add_atom(0.5*L, 0.0*L, 0.0*L, 17.0, 17);
//    molecule.add_atom(0.0*L, 0.5*L, 0.0*L, 17.0, 17);
//    molecule.add_atom(0.0*L, 0.0*L, 0.5*L, 17.0, 17);

    // Cubic cell for Si
    //molecule.add_atom(  0,     0,     0,     14.0, 14);
    //molecule.add_atom(  L/2,   L/2,   0,     14.0, 14);
    //molecule.add_atom(  L/2,   0,     L/2,   14.0, 14);
    //molecule.add_atom(  0,     L/2,   L/2,   14.0, 14);
    //molecule.add_atom(  L/4,   L/4,   L/4,   14.0, 14);
    //molecule.add_atom(  3*L/4, 3*L/4, L/4,   14.0, 14);
    //molecule.add_atom(  3*L/4, L/4,   3*L/4, 14.0, 14);
    //molecule.add_atom(  L/4,   3*L/4, 3*L/4, 14.0, 14);

    molecule.update_rcut_with_eprec(1e-3);

    // Load basis
    aobasis.read_file("sto-3g");

    // Nuclear potential
    real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearDensityFunctor(molecule, L))).truncate_mode(0).truncate_on_project();
    double nuclear_charge=vnuc.trace();
    if (world.rank() == 0) print("total nuclear charge", nuclear_charge);
    vnuc = -1.0*make_coulomb_potential(world, vnuc);
    vnuc.truncate();
    int nst = int(molecule.total_nuclear_charge() + 0.1)/2;

    // Guess density
    real_function_3d rho = real_factory_3d(world).functor(real_functor_3d(new MolecularGuessDensityFunctor(molecule,aobasis))).truncate_on_project();
    rho.truncate();
    double rhot = rho.trace();
    if (world.rank() == 0) print("total guess charge", rhot);
    rho.scale(molecule.total_nuclear_charge()/rhot);
    
    // Make AO basis functions
    auto psi = initial_guess(world, vnuc, rho, kpoints, nst);

    for (int iter=0; iter<100; iter++) {
        if (world.rank() == 0) print("\n\n  Iteration",iter,"\n");
        auto v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
        truncate(world, psi);
        auto vpsi = apply_potential(world, v, psi);

        vector_complex_function_3d new_psi(nst*nkpt);
        vector_complex_function_3d rm(nst*nkpt);
        for (int ik = 0; ik < nkpt; ++ik) {
            vector_complex_function_3d psik(psi.begin()+ik*nst, psi.begin()+(ik+1)*nst);
            vector_complex_function_3d vpsik(vpsi.begin()+ik*nst, vpsi.begin()+(ik+1)*nst);

            auto ke_mat = make_kinetic_matrix(world, psik, kpoints[ik]);
            auto pe_mat = matrix_inner(world, psik, vpsik, true);
            auto ov_mat = matrix_inner(world, psik, psik, true);
    
            auto fock = ke_mat + pe_mat;
            // eliminate small off-diagonal elements and lift diagonal
            // degeneracies to reduce random mixing
            for (unsigned int i=0; i<psik.size(); i++) {
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
            if (world.rank() == 0) print("main() ik = ", ik);
            if (world.rank() == 0) print(e);
  
            // if (world.rank() == 0) print("fock: "); 
            // if (world.rank() == 0) print(fock);
            // if (world.rank() == 0) print("overlap: "); 
            // if (world.rank() == 0) print(ov_mat);
            // if (world.rank() == 0) print("eigenvectors: "); 
            // if (world.rank() == 0) print(real(c));
            // if (world.rank() == 0) print(imag(c));
            fixphases(world, e, c); 

            psik = transform(world, psik, c);
            vpsik = transform(world, vpsik, c);
            vcpairT pair_update = update(world, ik, psik, vpsik, kpoints[ik], v, e);
            vector_complex_function_3d new_psik = pair_update.first; 
            vector_complex_function_3d rmk = pair_update.second;

            //subspace[ik].update_subspace(world, new_psik, psik, rmk);
            auto damp = 0.3;
            gaxpy(world,damp,psik,(1.0-damp),new_psik);

            //orthogonalize(world,psik);
            //truncate(world,psik);

            truncate(world,psik);
            orthogonalize(world,psik);

            for (int ist = 0; ist < nst; ist++) {
                psi[ik*nst+ist] = psik[ist];
            }
        }

        if (world.rank() == 0) print(nkpt,nst,psi.size());
        MADNESS_CHECK(nkpt*nst == (int) psi.size());

        if (iter == 20) {
          if (world.rank() == 0) print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+2, thresh*1e-2, true); 
            rho = madness::project(rho, kwavelet+2, thresh*1e-2, true); 
          for (unsigned int i = 0; i < psi.size(); i++) {
            FunctionDefaults<3>::set_k(kwavelet+2);
            FunctionDefaults<3>::set_thresh(thresh*1e-2);
            psi[i] = madness::project(psi[i], kwavelet+2, thresh*1e-2, true); 
          }
          if (world.rank() == 0) print("done reprojecting ..");
        }
        if (iter == 30) {
          if (world.rank() == 0) print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+4, thresh*1e-4, true); 
            rho = madness::project(rho, kwavelet+4, thresh*1e-4, true); 
          for (unsigned int i = 0; i < psi.size(); i++) {
            FunctionDefaults<3>::set_k(kwavelet+4);
            FunctionDefaults<3>::set_thresh(thresh*1e-4);
            psi[i] = madness::project(psi[i], kwavelet+4, thresh*1e-4, true); 
          }
          if (world.rank() == 0) print("done reprojecting ..");
        }
        if (iter == 50) {
          if (world.rank() == 0) print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+6, thresh*1e-6, true); 
            rho = madness::project(rho, kwavelet+6, thresh*1e-6, true); 
          for (unsigned int i = 0; i < psi.size(); i++) {
            FunctionDefaults<3>::set_k(kwavelet+6);
            FunctionDefaults<3>::set_thresh(thresh*1e-6);
            psi[i] = madness::project(psi[i], kwavelet+6, thresh*1e-6, true); 
          }
          if (world.rank() == 0) print("done reprojecting ..");
        }


        auto rho_new = make_density(world, psi, weight);
        double rdiff = (rho-rho_new).norm2();
        double s = 0.3;
        rho = s*rho + (1.0-s)*rho_new;
        if (world.rank() == 0) printf("electron density difference:  %15.8e\n\n", rdiff);
    }
    return 0;
}
