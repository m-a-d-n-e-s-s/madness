/*
  free space evals (k=6,tol=1e-4) from moldft
  -3.0306e+01 -1.3228e+00 -4.9800e-01 -4.9800e-01 -4.9800e-01

  computed by testperiodic with gamma point L=30.0
  -3.0304e+01 -1.3213e+00 -4.9782e-01 -4.9782e-01 -4.9782e-01

 */


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>
using namespace madness;

#include <polar/molecule.h>
#include <polar/molecularbasis.h>
#include <polar/xcfunctional.h>

static const double_complex I(0,1);
static const double twopi = 2.0*constants::pi;

static const double L = 5.0; // Unit cell size in au for neon
//static const double L = 8.37; // Unit cell size in au for neon
//static const double L = 7.65; // Unit cell size in au for LiF
//static const double L = 3.8; // Unit cell size in au for LiF
//static const double L = 8.0;

static const int R = 1; // periodic sums from -R to +R inclusive
static const double thresh = 1e-4;
static const double kwavelet = 6;
static const int truncate_mode = 0;

static const double kx=0.5*twopi/L, ky=0.5*twopi/L, kz=0.5*twopi/L;
//static const double kx=0, ky=0, kz=0;

static Molecule molecule;
static AtomicBasisSet aobasis;

typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis)
        : molecule(molecule), aobasis(aobasis) {}

    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        for (int i=-R; i<=+R; i++) {
            for (int j=-R; j<=+R; j++) {
                for (int k=-R; k<=+R; k++) {
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
        for (int i=-R; i<=+R; i++) {
            double xx = x[0]+i*L;
            for (int j=-R; j<=+R; j++) {
                double yy = x[1]+j*L;
                for (int k=-R; k<=+R; k++) {
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
        static const int R = std::max(::R,1);
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

    tensor_complex f1 = 0.5 * (matrix_inner(world, dvx, dvx, true) +
                               matrix_inner(world, dvy, dvy, true) +
                               matrix_inner(world, dvz, dvz, true));

    tensor_complex f2 =
        (-I*kx)*matrix_inner(world, v, dvx, false) +
        (-I*ky)*matrix_inner(world, v, dvy, false) +
        (-I*kz)*matrix_inner(world, v, dvz, false);

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
    for (int i=0; i<nmo; i++) {
        for (int axis=0; axis<3; axis++) {
            complex_derivative_3d D(world, axis);
            vpsi[i] = vpsi[i] - (I*k[axis])*D(psi[i]);
        }
        vpsi[i] = vpsi[i] + (0.5*ksq)*psi[i];
    }

    // determine shift to make homo <=-0.1
    double shift = 0.0;
    if (e(nmo-1) > -0.1) {
        shift = -0.1 - e(nmo-1);
        gaxpy(world, 1.0, vpsi, shift, psi);
    }

    // Do the BSH thing
    scale(world, vpsi, -2.0);
    truncate(world, vpsi);
    vector<poperatorT> ops = make_bsh_operators(world, e, shift);
    vector_complex_function_3d new_psi = apply(world, ops, vpsi);

    // Step restriction
    double damp;
    if (iter < 10) damp = 0.95;
    else if (iter < 20) damp = 0.85;
    else damp = 0.75;
    print("  shift", shift, "damp", damp, "\n");

    printf("      eigenvalue    residual\n");
    for (int i=0; i<nmo; i++) {
        double rnorm = (psi[i]-new_psi[i]).norm2();
        printf("%4d  %10.6f  %10.1e\n", i, e[i], rnorm);
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
    // molecule.add_atom(L/2,L/2,  0, 9.0, 9);
    // molecule.add_atom(L/2,  0,L/2, 9.0, 9);
    // molecule.add_atom(  0,L/2,L/2, 9.0, 9);
    // molecule.add_atom(L/2,  0,  0, 3.0, 3);
    // molecule.add_atom(  0,L/2,  0, 3.0, 3);
    // molecule.add_atom(  0,  0,L/2, 3.0, 3);
    molecule.add_atom(L/2,L/2,L/2, 3.0, 3);

    molecule.set_eprec(1e-3);

    // Load basis
    aobasis.read_file("sto-3g");

    // Nuclear potential
    real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearDensityFunctor(molecule))).truncate_mode(0).truncate_on_project();
    print("total nuclear charge", vnuc.trace());
    vnuc = -1.0*make_coulomb_potential(world, vnuc);
    vnuc.truncate();

    // Guess density
    real_function_3d rho = real_factory_3d(world).functor(real_functor_3d(new MolecularGuessDensityFunctor(molecule,aobasis))).truncate_on_project();
    rho.truncate();
    double rhot = rho.trace();
    print("total guess charge", rhot);
    rho.scale(molecule.total_nuclear_charge()/rhot);

    int nmo = int(molecule.total_nuclear_charge() + 0.1)/2;

    // Make AO basis functions
    vector_complex_function_3d psi = makeao(world);
    vector_real norms = norm2s(world, psi);
    print(norms);

    for (int iter=0; iter<100; iter++) {
        print("\n\n  Iteration",iter,"\n");
        real_function_3d v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
        vector_complex_function_3d vpsi = apply_potential(v, psi);

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

        psi = update(world, psi, vpsi, e, iter);

        rho = make_density(world, psi);
    }
    return 0;
}
