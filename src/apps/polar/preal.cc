#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/solvers.h>
#include <tinyxml/tinyxml.h>
#include <world/worldmem.h>

using namespace madness;

#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>
#include <moldft/xcfunctional.h>

#include "gth_pseudopotential.h"

template <typename Q, int NDIM>
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
}


static const double twopi = 2.0*constants::pi;

static const double L = 51.143;
//static const double L = 50.0;

static double thresh = 1e-6;
static double kwavelet = 8;
static const int truncate_mode = 0;

static Molecule molecule;
static AtomicBasisSet aobasis;

typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

static std::vector<int> at_to_bf, at_nbf;

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis)
        : molecule(molecule), aobasis(aobasis) {}

    double operator()(const coord_3d& x) const {
        return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
    }
};

template <typename Q>
class AtomicBasisFunctor : public FunctionFunctorInterface<Q,3> {
private:
    const AtomicBasisFunction aofunc;

public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc)
        : aofunc(aofunc)
    {}

    Q operator()(const coord_3d& x) const {
        return aofunc(x[0], x[1], x[2]);
    }

    std::vector<coord_3d> special_points() const {
        return std::vector<coord_3d>(1,aofunc.get_coords_vec());
    }

    Level special_level() {
      return 8;
    }
};

class NuclearDensityFunctor : public FunctionFunctorInterface<double,3> {
  Molecule molecule;
  std::vector<coord_3d> specialpts;
public:
  NuclearDensityFunctor(const Molecule& molecule) : 
    molecule(molecule), specialpts(molecule.get_all_coords_vec()) {}
 
  double operator()(const Vector<double,3>& r) const {
    return molecule.mol_nuclear_charge_density(r[0], r[1], r[2]);
  }

  std::vector<coord_3d> special_points() const{
    return specialpts;
  }

  Level special_level() {
    return 15;
  }
};

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
public:
    MolecularPotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coord_3d& x) const {
        return molecule.nuclear_attraction_potential(x[0], x[1], x[2]);
    }

    std::vector<coord_3d> special_points() const {return molecule.get_all_coords_vec();}
};

vector_real_function_3d makeao(World& world) {
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
    vector_real_function_3d ao(aobasis.nbf(molecule));
    for(int i = 0; i<aobasis.nbf(molecule); ++i) {
        real_functor_3d aofunc(new AtomicBasisFunctor<double>(aobasis.get_atomic_basis_function(molecule, i)));
        ao[i] = real_factory_3d(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(0);
        AtomicBasisFunction atbf = aobasis.get_atomic_basis_function(molecule, i);
        //if (world.rank() == 0) {
        //  printf("%d:    ", i);
        //  atbf.print_me(std::cout);
        //}
    }
    world.gop.fence();
    truncate(world, ao);
    normalize(world, ao);
    print_meminfo(world.rank(), "makeao");

    return ao;
}

vector_real_function_3d makeao_real(World& world) {
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
    vector_real_function_3d ao(aobasis.nbf(molecule));
    for(int i = 0; i<aobasis.nbf(molecule); ++i) {
        print("basis function ", i);
        real_functor_3d aofunc(new AtomicBasisFunctor<double>(aobasis.get_atomic_basis_function(molecule, i)));
        ao[i] = real_factory_3d(world).functor(aofunc).truncate_on_project().truncate_mode(0);
    }
    print("world.gop.fence");
    world.gop.fence();
    //print("truncating");
    //truncate(world, ao);
    print("normalizing");
    normalize(world, ao);
    print_meminfo(world.rank(), "makeao");

    return ao;
}

real_tensor kinetic_energy_dir_matrix(World & world, int axis, const vector_real_function_3d& v)
{
        std::vector< std::shared_ptr<real_derivative_3d> > gradop = gradient_operator<double,3>(world);

        reconstruct(world, v);
        int n = v.size();
        real_tensor r(n, n);
        vector_real_function_3d dv = apply(world, *(gradop[axis]), v);
        r = matrix_inner(world, dv, dv, true);
        dv.clear();

        return r.scale(0.5);
}

real_tensor kinetic_energy_matrix(World& world, const vector_real_function_3d& v)
{
    std::vector< std::shared_ptr<real_derivative_3d> > gradop = gradient_operator<double,3>(world);
    reconstruct(world, v);
    int n = v.size();
    real_tensor r(n, n);
    for(int axis = 0;axis < 3;++axis){
        vector_real_function_3d dv = apply(world, *(gradop[axis]), v);
        r += matrix_inner(world, dv, dv, true);
        dv.clear();
    }
    return r.scale(0.5);
}

vector_real_function_3d apply_potential(World& world, const real_function_3d& potential, const vector_real_function_3d& psi)
{
    //vector_real_function_3d vpsi;
    //for (unsigned int i=0; i<psi.size(); i++)
    //    vpsi.push_back(potential*psi[i]);
    //return vpsi;
    return mul_sparse(world,potential,psi,FunctionDefaults<3>::get_thresh()*0.01);
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
    real_convolution_3d op = CoulombOperator(world, 1e-10, thresh);
    return op(rho);
}

vector<poperatorT> make_bsh_operators(World & world, const real_tensor& evals, double shift)
{
    int nmo = evals.dim(0);
    vector<poperatorT> ops(nmo);
    for(int i = 0;i < nmo; ++i){
        double eps = evals(i) + shift;
        ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  1e-10, thresh));
    }
    return ops;
}

void orthogonalize(World& world, vector_real_function_3d& psi) {
    compress(world, psi);
    for (unsigned int i = 0; i<psi.size(); i++) {
        real_function_3d& psi_i = psi[i];
        psi_i.scale(1.0/psi_i.norm2());
        for (unsigned int j = 0; j<i; j++) {
            real_function_3d& psi_j = psi[j];
            double s = inner(psi_j,psi_i);
            psi_i.gaxpy(1.0,psi_j,-s); // |i> = |i> - |j><j|i>
            psi_i.scale(1.0/psi_i.norm2());
        }
    }
}


// DESTROYS VPSI
vector_real_function_3d update(World& world,
                                  const vector_real_function_3d& psi,
                                  vector_real_function_3d& vpsi,
                                  const real_tensor& e,
                                  int iter,
                                  double& ravg)
{
    // psi = - 2 G(E+shift) * (V+shift) psi
    int nmo = psi.size();

    // determine shift to make homo <=-0.1
    double shift = 0.0;
    if (e(nmo-1) > -0.1) {
        shift = -0.1 - e(nmo-1);
        gaxpy(world, 1.0, vpsi, shift, psi);
        print("shifting vpsi by ", shift);
    }

    // Do the BSH thing
    scale(world, vpsi, -2.0);
    truncate(world, vpsi);
    vector<poperatorT> ops = make_bsh_operators(world, e, shift);
    vector_real_function_3d new_psi = apply(world, ops, vpsi);

    printf("\n\nBSH: traces of psi:\n");
    for (unsigned int i = 0; i < psi.size(); i++)
    {
      double tr = psi[i].trace();
      printf("%d      trace: %15.8e\n",i,tr);
    }
    printf("\n\nBSH: traces of vpsi:\n");
    for (unsigned int i = 0; i < vpsi.size(); i++)
    {
      double tr = vpsi[i].trace();
      printf("%d      trace: %15.8e\n",i,tr);
    }
    printf("\n\nBSH: traces of new_psi:\n");
    for (unsigned int i = 0; i < new_psi.size(); i++)
    {
      double tr = new_psi[i].trace();
      printf("%d      trace: %15.8e\n",i,tr);
    }

    // Step restriction
    double damp = 0.2;
    //if (iter < 10) damp = 0.95;
    //else if (iter < 20) damp = 0.85;
    //else damp = 0.75;
    //print("  shift", shift, "damp", damp, "\n");

    printf("      eigenvalue    residual\n");
    ravg = 0.0;
    for (int i=0; i<nmo; i++) {
        double rnorm = (psi[i]-new_psi[i]).norm2();
        ravg += rnorm/(double)nmo;
        printf("%4d  %15.10f  %10.1e\n", i, e[i], rnorm);
        new_psi[i] = damp*psi[i] + (1.0-damp)*new_psi[i];
    }

    truncate(world,new_psi);
    normalize(world, new_psi);
    orthogonalize(world, new_psi);
    truncate(world,new_psi);
    normalize(world, new_psi);
    return new_psi;
}

real_function_3d make_density(World& world, const vector_real_function_3d& v) {

   vector_real_function_3d vsq = square(world, v);
   compress(world, vsq);
   real_function_3d rho = real_factory_3d(world);
   rho.compress();
   for(unsigned int i = 0;i < vsq.size();++i){
     rho.gaxpy(1.0, vsq[i], 2.0, false);

   }
   world.gop.fence();
   vsq.clear();
   return rho;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(9);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_k(kwavelet);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_truncate_mode(truncate_mode);


    // Hydrogen atom
    // Also for hydrogen need to make changes to nmo (set to 1)
    // and make_density (multiply by 1.0 instead of 2.0)
    //molecule.add_atom(0.0, 0.0, 0.0, 1.0, 1); 

    // Hydrogen molecule 
    //molecule.add_atom(0.0, 0.0, 0.72372451, 1.0, 1); 
    //molecule.add_atom(0.0, 0.0, -0.72372451, 1.0, 1); 

    // Helium atom
    //molecule.add_atom(0.0, 0.0, 0.0, 2.0, 2); 
    //molecule.add_atom(0.0, 0.0, 5.0, 2.0, 2); 

    // Neon atom
    //molecule.add_atom(0, 0, 0, 8.0, 10);
    //molecule.add_atom(0, 0, 0, 10.0, 10);

    // Argon atom
    molecule.add_atom(0, 0, 0, 8.0, 18);

    // Calcium atom
    //molecule.add_atom(0, 0, 0, 2.0, 20);

    // Carbon atom
    //molecule.add_atom(0, 0, 0, 4.0, 6);

    // Oxygen atom
    //molecule.add_atom(0.0, 0.0, 2.286, 6.0, 8); 

    // Beryllium atom (at arbitrary location)
    //molecule.add_atom(1.23943, -0.3422, 5, 2.0, 4);

    // Water (doesn't work)
    //molecule.add_atom( 1.4375, 0, 1.15, 1.0, 1);
    //molecule.add_atom(-1.4375, 0, 1.15, 1.0, 1);
    //molecule.add_atom(0, 0, 0, 6.0, 8);
    
    // CO2 (doesn't work)
    //molecule.add_atom(0.0, 0.0, 0.0, 4.0, 6);
    //molecule.add_atom(0.0, 0.0, 2.24064, 6.0, 8);
    //molecule.add_atom(0.0, 0.0, -2.24064, 6.0, 8);

    // O2 (doesn't work)
    //molecule.add_atom(0.0, 0.0, 0.0  , 8.0, 8); 
    //molecule.add_atom(0.0, 0.0, 2.286, 8.0, 8); 

    molecule.orient();
    molecule.set_eprec(1e-5);

    molecule.print();

    // WSTHORNTON
    //real_function_3d f1 = real_factory_3d(world).f(proj_p_1_x_ne);
    //
    //GTHPseudopotential<double> ppotential2(world, molecule);
    //int maxL = ppotential2.radii[9].dim(0)-1;
    //RlmStore rlmstore(world, maxL, coord_3d(0.0));
    //ProjStore projstore(world, ppotential2.radii[9], coord_3d(0.0));
    ////real_function_3d f2a = real_factory_3d(world).functor(real_functor_3d(new ProjFunctor(0.214913, 1, 1, coord_3d(0.0))));
    ////real_function_3d f2b = real_factory_3d(world).functor(real_functor_3d(new RlmFunctor(1, 0, coord_3d(0.0))));
    //real_function_3d f2 = rlmstore.rlm(1,0)*projstore.nlproj(1,1);

    //coord_3d lo1 {-L,0.0,0.0}; coord_3d hi1 {L,0.0,0.0};
    //plot_line("error1.dat", 20001, lo1, hi1, f1);
    //plot_line("error2.dat", 20001, lo1, hi1, f2);
    //double err = (f1-f2).norm2();
    //print("err = ", err);
    //MADNESS_ASSERT(false);

    // Load basis
    //aobasis.read_file("sto-3g");
    aobasis.read_file("6-31g");

    // Nuclear potential (don't need if using pseudopotential)
    //real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new NuclearDensityFunctor(molecule))).truncate_mode(0).truncate_on_project();
    double safety = 0.1;
    double vtol = FunctionDefaults<3>::get_thresh() * safety;
    real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(new MolecularPotentialFunctor(molecule))).thresh(vtol).truncate_on_project();
    //vnuc.set_thresh(FunctionDefaults<3>::get_thresh());
    vnuc.reconstruct();
    print("total nuclear charge", vnuc.trace());
    //vnuc = -1.0*make_coulomb_potential(world, vnuc);

    // Create pseudopotential
    GTHPseudopotential<double> ppotential(world, molecule);

    // Guess density
    real_function_3d rho = real_factory_3d(world).functor(real_functor_3d(new MolecularGuessDensityFunctor(molecule,aobasis))).truncate_on_project();
    double nel = rho.trace();
    if(world.rank() == 0)
        print("guess dens trace", nel);
    int nmo = int(molecule.total_nuclear_charge() + 0.1)/2;
    rho.scale((2.0*nmo)/nel); 

    // Make AO basis functions
    vector_real_function_3d psi = makeao(world);
    vector_real norms = norm2s(world, psi);

    printf("\n\ntraces of atomic functions:\n");
    for (unsigned int i = 0; i < psi.size(); i++)
    {
      double tr = psi[i].trace();
      printf("%d      trace: %15.8e\n",i,tr);
    }

    coord_3d lo(-L); coord_3d hi(L);
    plot_line("initial_rho.dat", 20001, lo, hi, rho);
    plot_line("vnuc1.dat", 20001, lo, hi, vnuc);
    coord_3d lox(0.0); coord_3d hix(0.0);
    coord_3d loy(0.0); coord_3d hiy(0.0);
    coord_3d loz(0.0); coord_3d hiz(0.0);
    lox[0] = -L; hix[0] = L;
    loy[1] = -L; hiy[1] = L;
    loz[2] = -L; hiz[2] = L;
    for (unsigned int i = 0; i < psi.size(); i++)
    {
      real_function_3d rpsi = psi[i];
      char fnamex[25];
      sprintf(fnamex, "aox_%d.dat", i);
      plot_line(fnamex, 20001, lox, hix, rpsi);
      char fnamey[25];
      sprintf(fnamey, "aoy_%d.dat", i);
      plot_line(fnamey, 20001, loy, hiy, rpsi);
      char fnamez[25];
      sprintf(fnamez, "aoz_%d.dat", i);
      plot_line(fnamez, 20001, loz, hiz, rpsi);
    }

    double ravg = 10000.0;
    for (int iter=0; iter<100 && ravg > 1e-6; iter++) {
        print("\n\n  Iteration",iter,"\n");
        double rtr = rho.trace();
        print("rho trace:  ", rtr);
        //real_function_3d v = vnuc + make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
        //vector_real_function_3d vpsi = apply_potential(world, v, psi);

        real_function_3d v = make_coulomb_potential(world,rho) + make_lda_potential(world,rho);
        vector_real_function_3d vpsi = ppotential.apply_potential(world, v, psi);
        //vector_real_function_3d vpsi = ppotential.apply_potential_ne(world, v, psi);

        //vector_complex_function_3d cpsi(psi.size());
        //for (unsigned int i = 0; i < psi.size(); i++) cpsi[i] = function_real2complex(psi[i]);
        //vector_real_function_3d vpsi = ppotential.apply_potential_wsttiger(world, v, cpsi);

        printf("\n\ntraces of Vpsi:\n");
        for (unsigned int i = 0; i < psi.size(); i++)
        {
          double tr = vpsi[i].trace();
          printf("%d      trace: %15.8e\n",i,tr);
        }

        print("lo ", lo); print("hi ", hi);
        for (unsigned int i = 0; i < psi.size(); i++)
        {
          char fname[25];
          sprintf(fname, "vpsi_%d.dat", i);
          plot_line(fname, 20001, lo, hi, vpsi[i]);
        }

        real_tensor ke_mat = kinetic_energy_matrix(world, psi);
        real_tensor pe_mat = matrix_inner(world, psi, vpsi, true);
        real_tensor ov_mat = matrix_inner(world, psi, psi, true);

        //real_tensor ke0 = kinetic_energy_dir_matrix(world, 0, psi);
        //real_tensor ke1 = kinetic_energy_dir_matrix(world, 1, psi);
        //real_tensor ke2 = kinetic_energy_dir_matrix(world, 2, psi);

        real_tensor fock = ke_mat + pe_mat;
        fock = 0.5*(fock + transpose(fock));

        //real_tensor Utmp; real_tensor etmp;
        //double vtol = FunctionDefaults<3>::get_thresh()*safety;
        //vector_real_function_3d vpsi1 = mul_sparse(world, vnuc, psi, vtol);
        //real_tensor pe_mat1 = matrix_inner(world, psi, vpsi1, true);
        //vector_real_function_3d vpsi2 = mul_sparse(world, make_coulomb_potential(world,rho), psi, vtol);
        //real_tensor pe_mat2 = matrix_inner(world, psi, vpsi2, true);
        //vector_real_function_3d vpsi3 = mul_sparse(world, make_lda_potential(world,rho), psi, vtol);
        //real_tensor pe_mat3 = matrix_inner(world, psi, vpsi3, true);
        //vector_real_function_3d vrho = mul_sparse(world, rho, psi, vtol);
        //real_tensor re_mat = matrix_inner(world, psi, vrho, true);
        //syev(pe_mat, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs:");
        //  print(etmp);
        //}
        //syev(pe_mat1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vnuc)");
        //  print(etmp);
        //}
        //syev(pe_mat2, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vc)");
        //  print(etmp);
        //}
        //syev(pe_mat3, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vxc)");
        //  print(etmp);
        //}
        //syev(ke_mat, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("ke eigs:");
        //  print(etmp);
        //}
        //sygv(re_mat, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("re eigs: (overlap)");
        //  print(etmp);
        //}
        //sygv(pe_mat, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (overlap)");
        //  print(etmp);
        //}
        //sygv(pe_mat1, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vnuc)(overlap)");
        //  print(etmp);
        //}
        //sygv(pe_mat2, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vc)(overlap)");
        //  print(etmp);
        //}
        //sygv(pe_mat3, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("pe eigs: (vxc)(overlap)");
        //  print(etmp);
        //}
        //sygv(ke_mat, ov_mat, 1, Utmp, etmp);
        //if (world.rank() == 0) {
        //  printf("ke eigs: (overlap)");
        //  print(etmp);
        //}
        //if (world.rank() == 0) {
        //  coord_3d rr1(1.235);
        //  coord_3d rr2(0.142);
        //  print("rho(1) = ", rho(rr1));
        //  print("rho(2) = ", rho(rr2));
        //}
        //plot_line("rho.dat", 20001, lo, hi, rho);

        //print("KE"); print(ke_mat);
        //print("PE"); print(pe_mat);
        //print("OV"); print(ov_mat);
        //print("FOCK"); print(fock);

        // eliminate small off-diagonal elements and lift diagonal
        // degeneracies to reduce random mixing
        //for (unsigned int i=0; i<psi.size(); i++) {
        //    fock(i,i) += i*thresh*1e-2;
        //    for (unsigned int j=0; j<i; j++) {
        //        if (std::abs(fock(i,j)) < thresh*1e-1 || std::abs(ov_mat(i,j)) < thresh*1e-1) {
        //            fock(i,j) = fock(j,i) = 0.0;
        //            ov_mat(i,j) = ov_mat(j,i) = 0.0;
        //        }
        //    }
        //}

        real_tensor c;
        tensor_real e;
        sygv(fock, ov_mat, 1, c, e);
        double rfactor = 1e-3;
        for (int i = 0; i < fock.dim(0); i++)
        {
          double thresh = FunctionDefaults<3>::get_thresh();
          for (int j = 0; j < fock.dim(1); j++)
          {
            //if (std::abs(ke0(i,j)) < thresh*rfactor) ke0(i,j) = 0.0;
            //if (std::abs(ke1(i,j)) < thresh*rfactor) ke1(i,j) = 0.0;
            //if (std::abs(ke2(i,j)) < thresh*rfactor) ke2(i,j) = 0.0;

            if (std::abs(ov_mat(i,j)) < thresh*rfactor) ov_mat(i,j) = 0.0;
            if (std::abs(ke_mat(i,j)) < thresh*rfactor) ke_mat(i,j) = 0.0;
            if (std::abs(pe_mat(i,j)) < thresh*rfactor) pe_mat(i,j) = 0.0;
            if (std::abs(fock(i,j)) < thresh*0.1) fock(i,j) = 0.0;
          }
        }
        if (world.rank() == 0) {
          //print("ke0:");
          //print(real(ke0));
          //print("ke1:");
          //print(real(ke1));
          //print("ke2:");
          //print(real(ke2));
          //print("Overlap"); 
          //print(real(ov_mat));
          print("Kinetic (real part)"); 
          print(real(ke_mat));
          //print("Kinetic (imag part)"); 
          //print(imag(ke_mat));
          //print("Kinetic (real version)"); 
          //print(real(ke_mat_real));
          print("Potential (real part)"); 
          print(real(pe_mat));
          //print("Potential (imag part)"); 
          //print(imag(pe_mat));
          print("Fock (real part)"); 
          print(real(fock));
          //print("Fock (imag part)"); 
          //print(imag(fock));
          ////print("eigenvectors"); print(c);
          print("eigenvalues"); print(e);
        }

        if (iter == 0) {
            c = copy(c(_,Slice(0,nmo-1))); // truncate to occupied states
            e = e(Slice(0,nmo-1));
        }

        double trantol = vtol / std::min(30.0, double(psi.size()));
        psi = transform(world, psi, c, trantol, true);
        vpsi = transform(world, vpsi, c, trantol, true);

        psi = update(world, psi, vpsi, e, iter, ravg);
        truncate(world, psi);

        double rthresh = 50.0*thresh;
        if (ravg < rthresh) {
          print("reprojecting ...");
          kwavelet += 2;
          thresh /= 100.0;
          FunctionDefaults<3>::set_k(kwavelet);
          FunctionDefaults<3>::set_thresh(thresh);
          vnuc = madness::project(vnuc, kwavelet, thresh, true);
          ppotential.reproject(kwavelet, thresh);
          for (int i = 0; i < nmo; i++) psi[i] = madness::project(psi[i], kwavelet, thresh, true);  
        }

        //int md = psi[3].max_depth();
        //print("max depth:  ", md);

        rho = make_density(world, psi);
    }
    return 0;
}
