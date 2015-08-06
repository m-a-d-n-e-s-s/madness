/*
  free space evals (k=6,tol=1e-4) from moldft
  -3.0306e+01 -1.3228e+00 -4.9800e-01 -4.9800e-01 -4.9800e-01

  computed by testperiodic with gamma point L=30.0
  -3.0304e+01 -1.3213e+00 -4.9782e-01 -4.9782e-01 -4.9782e-01

 */


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>
using namespace madness;

#include <chem/molecule.h>
#include <chem/molecularbasis.h>
#include <chem/xcfunctional.h>

static const double_complex I(0,1);
static const double twopi = 2.0*constants::pi;

//static const double L = 5.0; // Unit cell size in au for neon
//static const double L = 8.37; // Unit cell size in au for neon
//static const double L = 7.65; // Unit cell size in au for LiF
static const double L = 3.8; // Unit cell size in au for LiF
//static const double L = 8.0;
//static const double L = 10.26085381075144364474; // Unit cell size in au for Si

static const int maxR = 2; // periodic sums from -R to +R inclusive
static const double thresh = 1e-6;
static const double kwavelet = 12;
static const int truncate_mode = 0;

static Molecule molecule;
static AtomicBasisSet aobasis;

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
public:
    AtomicOrbitalFunctor(const AtomicBasisFunction& aofunc, double R, int offx, int offy, int offz)
        : aofunc(aofunc), R(R), offx(offx), offy(offy), offz(offz)
    {
        double x, y, z;
        aofunc.get_coords(x,y,z);
        coord_3d r;
        r[0]=x+offx*R; r[1]=y+offy*R; r[2]=z+offz*R;
        specialpt=std::vector<coord_3d>(1,r);
    }

    double operator()(const coord_3d& x) const {
        return aofunc(x[0] + offx*R, x[1] + offy*R, x[2] + offz*R);
    }

    std::vector<coord_3d> special_points() const {return specialpt;}
};

//#define NTRANS 8
//class AtomicBasisFunctor : public FunctionFunctorInterface<double_complex,3> {
//private:
//    const AtomicBasisFunction aofunc;
//    const double kx;
//    const double ky;
//    const double kz;
//    const double R;
//    const double rangesq;
//    coord_3d r;
//    std::vector<coord_3d> specialpt;
//    Vector<std::complex<double>,2*NTRANS+1> tx; 
//    Vector<std::complex<double>,2*NTRANS+1> ty; 
//    Vector<std::complex<double>,2*NTRANS+1> tz; 
//
//public:
//    AtomicBasisFunctor(const AtomicBasisFunction& aofunc, const Vector<double,3>& kpt, 
//                       double R)
//        : aofunc(aofunc), kx(kpt[0]), ky(kpt[1]), kz(kpt[2]), R(R), rangesq(aofunc.rangesq()*2)
//    {
//    	double x, y, z;
//        aofunc.get_coords(x,y,z);
//        r[0]=x; r[1]=y; r[2]=z;
//        specialpt=std::vector<coord_3d>(1,r);
//        for (int ir = -NTRANS; ir <= NTRANS; ir += 1)
//        {
//          tx[ir+NTRANS] = exp(std::complex<double>(0.0, twopi*kx*ir * R));
//          ty[ir+NTRANS] = exp(std::complex<double>(0.0, twopi*ky*ir * R));
//          tz[ir+NTRANS] = exp(std::complex<double>(0.0, twopi*kz*ir * R));
//        }
//    }
//
//    double_complex operator()(const coord_3d& x) const {
//        double_complex sum = 0.0;
//        for (int i=-NTRANS; i<=+NTRANS; i++) {
//            const double xx = x[0]+i*R;
//            const double xxR = xx-r[0];
//            const double xxRsq = xxR*xxR;
//            if (xxRsq < rangesq) {
//                for (int j=-NTRANS; j<=+NTRANS; j++) {
//                    const double yy = x[1]+j*R;
//                    const double yyR = yy-r[1];
//                    const double yyRsq = yyR*yyR; 
//                    if (xxRsq+yyRsq < rangesq) { 
//                        for (int k=-NTRANS; k<=+NTRANS; k++) {
//                            const double zz = x[2]+k*R;
//                            const double zzR = zz-r[2];
//                            const double zzRsq = zzR*zzR;
//                            if (xxRsq+yyRsq+zzRsq < rangesq) {
//                                double ao = aofunc(xx, yy, zz);
//                                if (std::abs(ao) > 1e-8) {
//                                    sum += tx[xx+NTRANS]*ty[yy+NTRANS]*tz[zz+NTRANS]*ao;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        return sum*std::exp(-I*(kx*r[0]+ky*r[1]+kz*r[2]));
//    }
//
//    std::vector<coord_3d> special_points() const {return specialpt;}
//};

class NuclearDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const double R;
    std::vector<coord_3d> specialpt;
    const int maxR = 1;
public:
    NuclearDensityFunctor(const Molecule& molecule, double R)
        : molecule(molecule), R(R), specialpt(molecule.get_all_coords_vec())
    {}

    double operator()(const coord_3d& x) const {
        double big = 2*R + 6.0*molecule.smallest_length_scale();
        double sum = 0.0;
        for (int i=-maxR; i<=+maxR; i++) {
            double xx = x[0]+i*R;
            if (xx < big && xx > -big) {
                for (int j=-maxR; j<=+maxR; j++) {
                    double yy = x[1]+j*R;
                    if (yy < big && yy > -big) {
                        for (int k=-maxR; k<=+maxR; k++) {
                            double zz = x[2]+k*R;
                            if (zz < big && zz > -big)
                                sum += molecule.nuclear_charge_density(x[0]+i*R, x[1]+j*R, x[2]+k*R);
                        }
                    }
                }
            }
        }
        return sum;
    }

    std::vector<coord_3d> special_points() const {return specialpt;}

    Level special_level() {
        return 10;
    }

};

vector_complex_function_3d makeao(World& world, const std::vector<Vector<double,3> >& kpoints) {
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

//void makeao2(World& world, const std::vector<Vector<double,3> >& kpoints) {
//    int nkpt = kpoints.size();
//    int nbf = aobasis.nbf(molecule);
//    vector_real_function_3d aos(nbf*(2*maxR+1)*(2*maxR+1)*(2*maxR+1));
//    print("nbf:  ", nbf);
//    auto idx = 0;
//    for (int ibf = 0; ibf < nbf; ibf++) {
//        print("\n\nibf: ", ibf);
//        for (int i = -maxR; i <= maxR; i++) {
//            for (int j = -maxR; j <= maxR; j++) {
//                for (int k = -maxR; k <= maxR; k++) {
//                    print(idx,ibf,i,j,k);
//                    real_functor_3d aofunc(
//                       new AtomicOrbitalFunctor(
//                       aobasis.get_atomic_basis_function(molecule, ibf), L, i, j, k));
//                    real_function_3d ao = real_factory_3d(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(0);
//                    aos[idx++] = ao;
//                }
//            }
//        }
//    }
//}

vector_complex_function_3d makeao2(World& world, const std::vector<Vector<double,3> >& kpts, double R) {
    int nkpt = kpts.size();
    int nbf = aobasis.nbf(molecule);
    double t1 = 1./std::sqrt((double)nbf);
    vector_complex_function_3d aos = zero_functions_compressed<double_complex,3>(world, nkpt*nbf);
    print("nbf:  ", nbf);
    for (int ibf = 0; ibf < nbf; ibf++) {
        print("\n\nibf: ", ibf);
        for (int i = -maxR; i <= maxR; i++) {
            for (int j = -maxR; j <= maxR; j++) {
                for (int k = -maxR; k <= maxR; k++) {
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
    for (int ik = 0; ik < nkpt; ik++) {
        Vector<double,3> kpt = kpts[ik];
        complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(
          new ExpFunctor3d<double_complex>(-I*kpt[0],-I*kpt[1],-I*kpt[2]))).truncate_mode(0).truncate_on_project();
        for (int ibf = 0; ibf < nbf; ibf++) {
            aos[ik*nbf+ibf] = aos[ik*nbf+ibf]*phase_m;
        }
    }
    return aos;
}

tensor_complex make_kinetic_matrix(World& world, const vector_complex_function_3d& v, const Vector<double,3>& kpt) {
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

    return f1 + f2 + f3;
}

vector_complex_function_3d apply_potential(World& world, const real_function_3d& potential, const vector_complex_function_3d& psi)
{
    auto vpsi = mul(world, potential, psi, false);
    world.gop.fence();
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

    print("H:\n"); print(fock);
    print("S:\n"); print(ov_mat);
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
    auto vlda = copy(rho);
    vlda.reconstruct();
    vlda.unaryop(xc_lda_potential());
    return vlda;
}

real_function_3d make_coulomb_potential(World& world, const real_function_3d& rho)
{
    real_convolution_3d op = CoulombOperator(world, 1e-4, thresh);
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
  complex_function_3d g = phase_m*apply(op, phase_p*f);
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
vector_complex_function_3d update(World& world,
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
    truncate(world, vpsi);
    
    vector_complex_function_3d new_psi(nmo);
    
    for (int iorb = 0; iorb < nmo; iorb++) {
      new_psi[iorb] = apply_periodic_bsh(world, vpsi[iorb], kx, ky, kz, e[iorb]+shift, L);
    }

    // Step restriction
    double damp = 0.05;
    if (world.rank() == 0) print("  shift", shift, "damp", damp, "\n");

    if (world.rank() == 0) printf("kpoint:  %10.5f    %10.5f    %10.5f\n",kpt[0],kpt[1],kpt[2]);
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

//    // normalize(world, new_psi);
//    //new_psi.insert(new_psi.end(), psi.begin(), psi.end());    
//    new_psi = orth(world, new_psi);
//
//    auto result = diag_and_transform(world, kpt, v, new_psi, nmo);
//    tensor_real eigs = result.second;
//
//    for (int i=0; i<nmo; i++) {
//        double rnorm = (psi[i]-new_psi[i]).norm2();
//        if (world.rank() == 0) printf("%4d  %10.6f  %10.1e\n", i, eigs[i], rnorm);
//    }
    
    return new_psi;
}

real_function_3d make_density(World& world, const vector_complex_function_3d& v, double weight) {
    real_function_3d rho(world);
    for (unsigned int i=0; i<v.size(); i++) {
        rho = rho + weight*abssq(v[i]);
    }
    rho.scale(2.0); // total closed-shell density
    return rho;
}

// Return all orbitals for all kpoints
// Modifies the density rho
vector_complex_function_3d initial_guess(World& world, const real_function_3d& vnuc, real_function_3d& rho, const std::vector<Vector<double,3> >& kpoints, int nst) {
    print("gonna makeao2 now");
    auto psi0 = makeao2(world, kpoints, L);
    print("done makeao2");

    //print("gonna makeao now");
    //auto psi0 = makeao(world, kpoints);
    //print("done makeao");
    
    auto v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
    print("done making potential");
    auto vpsi = apply_potential(world, v, psi0);
    print("done applying potential");

    int nkpt = kpoints.size();
    int nsize = psi0.size();
    int nst_initial = nsize/nkpt;
    MADNESS_ASSERT(nst <= nst_initial);
    print("nsize: ", nsize);
    print("nst_initial: ", nst_initial);
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

        print("Initial Fock: "); print(real(fock));
        print("Initial Overlap: "); print(real(ov_mat));
        sygv(fock, ov_mat, 1, c, e);

        print("initial_guess() ik = ", ik);
        print(e);
    
        psik = transform(world, psik, c);
        vpsik = transform(world, vpsik, c);
        //psik = update(world, psik, vpsik, kpoints[ik], v, e);
        for (int ist = 0; ist < nst; ist++) {
            print("pushing back ist = ", ist);
            psi.push_back(psik[ist]);    
        }
    }
    print("HERE!!");

//    /////// BEGIN DEBUG CODE
//    vpsi = apply_potential(world, v, psi);
//    for (int ik = 0; ik < nkpt; ++ik) {
//        auto ke_mat = make_kinetic_matrix(world, psi, kpoints[ik]);
//        auto pe_mat = matrix_inner(world, psik, vpsi, true);
//        auto ov_mat = matrix_inner(world, psik, psi, true);
//        auto fock = ke_mat + pe_mat;
//        // eliminate small off-diagonal elements and lift diagonal
//        // degeneracies to reduce random mixing
//        for (unsigned int i=0; i<psik.size(); i++) {
//            fock(i,i) += i*thresh*1e-2;
//            for (unsigned int j=0; j<i; j++) {
//                if (std::abs(fock(i,j)) < thresh*1e-1 || std::abs(ov_mat(i,j)) < thresh*1e-1) {
//                    fock(i,j) = fock(j,i) = 0.0;
//                    ov_mat(i,j) = ov_mat(j,i) = 0.0;
//                }
//            }
//        }
//        tensor_complex c;
//        tensor_real e;
//        sygv(fock, ov_mat, 1, c, e);
//    }
//    /////// END DEBUG CODE

    // compute new density
    double weights = 1.0/(double)kpoints.size();
    rho = make_density(world,psi,weights);
    return psi;
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

    //// kpoint list
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
    //kpoints[0] = vec(0.0*twopi/L, 0.0*twopi/L, 0.0*twopi/L); 
    ////kpoints[0] = vec(0.5*twopi/L, 0.5*twopi/L, 0.5*twopi/L); 
    //double weight = 1.0/(double)nkpt;
    
    // // FCC unit cell for ne
    //molecule.add_atom(  0,  0,  0, 10.0, 10);
    // molecule.add_atom(L/2,L/2,  0, 10.0, 10);
    // molecule.add_atom(L/2,  0,L/2, 10.0, 10);
    // molecule.add_atom(  0,L/2,L/2, 10.0, 10);

    // Cubic cell for LiF
    molecule.add_atom(  0,  0,  0, 9.0, 9);
    //molecule.add_atom(L/2,L/2,  0, 9.0, 9);
    //molecule.add_atom(L/2,  0,L/2, 9.0, 9);
    //molecule.add_atom(  0,L/2,L/2, 9.0, 9);
    //molecule.add_atom(L/2,  0,  0, 3.0, 3);
    //molecule.add_atom(  0,L/2,  0, 3.0, 3);
    //molecule.add_atom(  0,  0,L/2, 3.0, 3);
    molecule.add_atom(L/2,L/2,L/2, 3.0, 3);

    // Cubic cell for Si
    //molecule.add_atom(  0,     0,     0,     14.0, 14);
    //molecule.add_atom(  L/2,   L/2,   0,     14.0, 14);
    //molecule.add_atom(  L/2,   0,     L/2,   14.0, 14);
    //molecule.add_atom(  0,     L/2,   L/2,   14.0, 14);
    //molecule.add_atom(  L/4,   L/4,   L/4,   14.0, 14);
    //molecule.add_atom(  3*L/4, 3*L/4, L/4,   14.0, 14);
    //molecule.add_atom(  3*L/4, L/4,   3*L/4, 14.0, 14);
    //molecule.add_atom(  L/4,   3*L/4, 3*L/4, 14.0, 14);

    molecule.set_eprec(1e-3);

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
        auto vpsi = apply_potential(world, v, psi);

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
            print("main() ik = ", ik);
            print(e);
    
            psik = transform(world, psik, c);
            vpsik = transform(world, vpsik, c);
            psik = update(world, psik, vpsik, kpoints[ik], v, e);

            for (int ist = 0; ist < nst; ist++) {
                psi[ik*nst+ist] = psik[ist];
                vpsi[ik*nst+ist] = vpsik[ist];
            }
        }

        print(nkpt,nst,psi.size());
        MADNESS_ASSERT(nkpt*nst == (int) psi.size());

        if (iter == 20) {
          print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+2, thresh*1e-2, true); 
          for (unsigned int i = 0; i < psi.size(); i++) {
            FunctionDefaults<3>::set_k(kwavelet+2);
            FunctionDefaults<3>::set_thresh(thresh*1e-2);
            psi[i] = madness::project(psi[i], kwavelet+2, thresh*1e-2, true); 
          }
          print("done reprojecting ..");
        }
        if (iter == 26) {
          print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+4, thresh*1e-4, true); 
          for (unsigned int i = 0; i < psi.size(); i++) {
            FunctionDefaults<3>::set_k(kwavelet+4);
            FunctionDefaults<3>::set_thresh(thresh*1e-4);
            psi[i] = madness::project(psi[i], kwavelet+4, thresh*1e-4, true); 
          }
          print("done reprojecting ..");
        }


        rho = make_density(world, psi, weight);
    }
    return 0;
}
