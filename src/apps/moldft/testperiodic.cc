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

static const int maxR = 3; // periodic sums from -R to +R inclusive
static const double thresh = 1e-4;
static const double kwavelet = 6;
static const int truncate_mode = 0;

//static const double kx=0.5*twopi/L, ky=0.5*twopi/L, kz=0.5*twopi/L;
static const double kx=0.5*twopi/L, ky=0.0, kz=0.0;
//static const double kx=0.0, ky=0.5*twopi/L, kz=0.0;
//static const double kx=0.0, ky=0.0, kz=0.0;

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
    std::vector<coord_3d> specialpt;
public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc)
        : aofunc(aofunc)
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
            for (int j=-maxR; j<=+maxR; j++) {
                double yy = x[1]+j*L;
                for (int k=-maxR; k<=+maxR; k++) {
                    double zz = x[2]+k*L;
                    sum += exp(-I*(kx*xx+ky*yy+kz*zz))*aofunc(xx, yy, zz);
                }
            }
        }
        return sum;
    }

    std::vector<coord_3d> special_points() const {return specialpt;}
};

class NuclearDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    std::vector<coord_3d> specialpt;
public:
    NuclearDensityFunctor(const Molecule& molecule)
        : molecule(molecule), specialpt(molecule.get_all_coords_vec())
    {}

    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        static const int R = std::max(::maxR,1);
        for (int i=-R; i<=+R; i++) {
            for (int j=-R; j<=+R; j++) {
                for (int k=-R; k<=+R; k++) {
                    sum += molecule.nuclear_charge_density(x[0]+i*L, x[1]+j*L, x[2]+k*L);
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

vector_complex_function_3d makeao(World& world) {
    vector_complex_function_3d ao(aobasis.nbf(molecule));
    for(int i = 0; i<aobasis.nbf(molecule); ++i) {
        complex_functor_3d aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule, i)));
        ao[i] = complex_factory_3d(world).functor(aofunc).truncate_on_project().truncate_mode(0);
    }
    return ao;
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

    tensor_complex f1 = 0.5 * (matrix_inner(world, dvx, dvx, false) +
                               matrix_inner(world, dvy, dvy, false) +
                               matrix_inner(world, dvz, dvz, false));

    tensor_complex f2 =
        (-I*kx)*matrix_inner(world, v, dvx, false) +
        (-I*ky)*matrix_inner(world, v, dvy, false) +
        (-I*kz)*matrix_inner(world, v, dvz, false);

    tensor_complex f3 = (0.5 * (kx*kx + ky*ky + kz*kz)) * matrix_inner(world, v, v, true);

    return f1 + f2 + f3;
}

vector_complex_function_3d apply_potential(World& world, const real_function_3d& potential, const vector_complex_function_3d& psi)
{
    vector_complex_function_3d vpsi = mul(world, potential, psi, false);
    world.gop.fence();
    return vpsi;
}


real_function_3d make_lda_potential(World& world, const real_function_3d &rho)
{
    real_function_3d vlda = copy(rho);
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
    //double ksq = kx*kx + ky*ky + kz*kz;
    //coord_3d k {kx, ky, kz};

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
//         //coperatorT op = PeriodicBSHOperator3D(world, {kx*L, ky*L, kz*L}, sqrt(-2.0*(e[iorb]+shift)),  1e-4, thresh);
//        //coperatorT op = PeriodicBSHOperator3D(world, {-kx*L, -ky*L, -kz*L}, sqrt(-2.0*(e[iorb]+shift)),  1e-4, thresh);
//        coperatorT op = PeriodicBSHOperator3D(world, {kx*L, ky*L, kz*L}, sqrt(-2.0*(e[iorb]+shift)),  1e-4, thresh);
//        operatorT op2 = BSHOperator3D(world, sqrt(-2.0*(e[iorb]+shift)),  1e-4, thresh);
//        complex_function_3d phase_p = complex_factory_3d(world).functor(complex_functor_3d(new ExpFunctor<double_complex>(I*kx,I*ky,I*kz))).truncate_mode(0).truncate_on_project();
//        complex_function_3d phase_m = complex_factory_3d(world).functor(complex_functor_3d(new ExpFunctor<double_complex>(-I*kx,-I*ky,-I*kz))).truncate_mode(0).truncate_on_project();
//        //new_psi[iorb] = phase_m*apply(op2, phase_p*vpsi[iorb]);
//        new_psi[iorb] = phase_p*apply(op2, phase_m*vpsi[iorb]);
      new_psi[iorb] = apply_periodic_bsh(world, vpsi[iorb], kx, ky, kz, e[iorb]+shift, L);
    }


    // Step restriction
    double damp;
    if (iter < 10) damp = 0.95;
    else if (iter < 20) damp = 0.85;
    else damp = 0.75;
    damp = 0.15;
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

real_function_3d make_density(World& world, const vector_complex_function_3d& v) {
    real_function_3d rho(world);
    for (unsigned int i=0; i<v.size(); i++) {
        rho = rho + abssq(v[i]);
    }
    rho.scale(2.0); // total closed-shell density
    return rho;
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

    // // FCC unit cell for ne
    //molecule.add_atom(  0,  0,  0, 10.0, 10);
    // molecule.add_atom(L/2,L/2,  0, 10.0, 10);
    // molecule.add_atom(L/2,  0,L/2, 10.0, 10);
    // molecule.add_atom(  0,L/2,L/2, 10.0, 10);

    // Cubic cell for LiF
    molecule.add_atom(  0,  0,  0, 9.0, 9);
    molecule.add_atom(L/2,L/2,  0, 9.0, 9);
    molecule.add_atom(L/2,  0,L/2, 9.0, 9);
    molecule.add_atom(  0,L/2,L/2, 9.0, 9);
    molecule.add_atom(L/2,  0,  0, 3.0, 3);
    molecule.add_atom(  0,L/2,  0, 3.0, 3);
    molecule.add_atom(  0,  0,L/2, 3.0, 3);
    molecule.add_atom(L/2,L/2,L/2, 3.0, 3);

    // Cubic cell for Si
    /*molecule.add_atom(  0,     0,     0,     14.0, 14);
    molecule.add_atom(  L/2,   L/2,   0,     14.0, 14);
    molecule.add_atom(  L/2,   0,     L/2,   14.0, 14);
    molecule.add_atom(  0,     L/2,   L/2,   14.0, 14);
    molecule.add_atom(  L/4,   L/4,   L/4,   14.0, 14);
    molecule.add_atom(  3*L/4, 3*L/4, L/4,   14.0, 14);
    molecule.add_atom(  3*L/4, L/4,   3*L/4, 14.0, 14);
    molecule.add_atom(  L/4,   3*L/4, 3*L/4, 14.0, 14);*/

   

    molecule.set_eprec(1e-3);

    // Load basis
    aobasis.read_file("sto-3g");

    // Nuclear potential
    real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearDensityFunctor(molecule))).truncate_mode(0).truncate_on_project();
    double nuclear_charge=vnuc.trace();
    if (world.rank() == 0) print("total nuclear charge", nuclear_charge);
    vnuc = -1.0*make_coulomb_potential(world, vnuc);
    vnuc.truncate();

    // Guess density
    real_function_3d rho = real_factory_3d(world).functor(real_functor_3d(new MolecularGuessDensityFunctor(molecule,aobasis))).truncate_on_project();
    rho.truncate();
    double rhot = rho.trace();
    if (world.rank() == 0) print("total guess charge", rhot);
    rho.scale(molecule.total_nuclear_charge()/rhot);

    int nmo = int(molecule.total_nuclear_charge() + 0.1)/2;

    // Make AO basis functions
    vector_complex_function_3d psi = makeao(world);
    vector_real norms = norm2s(world, psi);
    if (world.rank() == 0) print(norms);

    for (int iter=0; iter<100; iter++) {
        if (world.rank() == 0) print("\n\n  Iteration",iter,"\n");
        real_function_3d v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
        vector_complex_function_3d vpsi = apply_potential(world, v, psi);

        tensor_complex ke_mat = make_kinetic_matrix(world, psi);
        tensor_complex pe_mat = matrix_inner(world, psi, vpsi, true);
        tensor_complex ov_mat = matrix_inner(world, psi, psi, true);

        //print("KE"); print(ke_mat);
        //print("PE"); print(pe_mat);
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

        //for (unsigned int i = 0; i < psi.size(); i++) {
        //  for (unsigned int j = 0; j < psi.size(); j++) {
        //    printf("%15.8e    ", real(fock(i,j)));
        //  }
        //  printf("\n");
        //}
        //printf("\n\n");
        //for (unsigned int i = 0; i < psi.size(); i++) {
        //  for (unsigned int j = 0; j < psi.size(); j++) {
        //    printf("%15.8e    ", real(ov_mat(i,j)));
        //  }
        //  printf("\n");
        //}

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

        if (iter == 8) {
          print("reprojecting ..");
            vnuc = madness::project(vnuc, kwavelet+2, thresh*1e-2, true); 
          for (int i = 0; i < nmo; i++) {
            FunctionDefaults<3>::set_k(kwavelet+2);
            FunctionDefaults<3>::set_thresh(thresh*1e-2);
            psi[i] = madness::project(psi[i], kwavelet+2, thresh*1e-2, true); 
            vpsi[i] = madness::project(vpsi[i], kwavelet+2, thresh*1e-2, true); 
          }
          print("done reprojecting ..");
        }
        psi = update(world, psi, vpsi, e, iter);

        rho = make_density(world, psi);
    }
    return 0;
}
