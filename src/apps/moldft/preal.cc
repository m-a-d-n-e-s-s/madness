#include <madness/mra/mra.h>
#include <madness/tensor/solvers.h>
#include <madness/external/tinyxml/tinyxml.h>
#include <madness/world/worldmem.h>

using namespace madness;

#include<madness/chem/molecule.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/xcfunctional.h>

#include<madness/chem/gth_pseudopotential.h>
#include "wst_functional.h"

static const double twopi = 2.0*constants::pi;

typedef std::shared_ptr<real_convolution_3d> poperatorT;

static std::vector<int> at_to_bf, at_nbf;

extern void drot(long n, double* MADNESS_RESTRICT a, double* MADNESS_RESTRICT b, double s, double c, long inc);

void drot3(long n, double* MADNESS_RESTRICT a, double* MADNESS_RESTRICT b, double s, double c, long inc) {
    if (inc == 1) {
        n*=3;
        for (long i=0; i<n; i+=3) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
    else {
        inc*=3;
        n*=inc;
        for (long i=0; i<n; i+=inc) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
}

// Functor for the initial guess atomic density
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

// Templated functor for the atomic Gaussian basis functions
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

    std::vector<coord_3d> special_points() const override final {
        return std::vector<coord_3d>(1,aofunc.get_coords_vec());
    }

    Level special_level() const override final {
      return 8;
    }
};

// Functor for the nuclear potential
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

/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coord_3d& x) const {
        return x[axis];
    }
};

struct BoysLocalization {

    std::vector<int> aset;

    BoysLocalization() {}

    BoysLocalization(unsigned int norbs) {
        aset = std::vector<int>(norbs, 0);
    }

    void set_size(unsigned int norbs) {
        aset = std::vector<int>(norbs, 0);
    }

    inline double DIP(const real_tensor& dip, int i, int j, int k, int l)
    {
        return dip(i, j, 0) * dip(k, l, 0) + dip(i, j, 1) * dip(k, l, 1) + dip(i, j, 2) * dip(k, l, 2);
    }

    real_tensor operator ()(World& world, 
                            const vector_real_function_3d& mo, 
                            const double thresh = 1e-9, 
                            const double thetamax = 0.5, 
                            const bool randomize = true)
    {
        //START_TIMER(world);
        const bool doprint = false;
        long nmo = mo.size();
        real_tensor dip(nmo, nmo, 3);
        double safety = 1e-1;
        double vtol = thresh*safety;
        for(int axis = 0;axis < 3;++axis){
            real_function_3d fdip = real_factory_3d(world).functor(real_functor_3d(new DipoleFunctor(axis))).initial_level(4);
            dip(_, _, axis) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
        }
        real_tensor U(nmo, nmo);
        if(world.rank() == 0){
            for(long i = 0;i < nmo;++i)
                U(i, i) = 1.0;

            double tol = thetamax;
            long ndone = 0;
            bool converged = false;
            for(long iter = 0;iter < 300;++iter){
                double sum = 0.0;
                for(long i = 0;i < nmo;++i){
                    sum += DIP(dip, i, i, i, i);
                }
                long ndone_iter = 0;
                double maxtheta = 0.0;
                if(doprint)
                    printf("iteration %ld sum=%.4f ndone=%ld tol=%.2e\n", iter, sum, ndone, tol);

                for(long i = 0;i < nmo;++i){
                    for(long j = 0;j < i;++j){
                        if (aset[i] == aset[j]) {
                            double g = DIP(dip, i, j, j, j) - DIP(dip, i, j, i, i);
                            double h = 4.0 * DIP(dip, i, j, i, j) + 2.0 * DIP(dip, i, i, j, j) - DIP(dip, i, i, i, i) - DIP(dip, j, j, j, j);
                            double sij = DIP(dip, i, j, i, j);
                            bool doit = false;
                            if(h >= 0.0){
                                doit = true;
                                if(doprint)
                                    print("             forcing negative h", i, j, h);

                                h = -1.0;
                            }
                            double theta = -g / h;
                            maxtheta = std::max<double>(std::abs(theta), maxtheta);
                            if(fabs(theta) > thetamax){
                                doit = true;
                                if(doprint)
                                    print("             restricting", i, j);

                                if(g < 0)
                                    theta = -thetamax;

                                else
                                    theta = thetamax * 0.8;

                            }
                            bool randomized = false;
                            if(randomize && iter == 0 && sij > 0.01 && fabs(theta) < 0.01){
                                randomized = true;
                                if(doprint)
                                    print("             randomizing", i, j);

                                theta += (RandomValue<double>() - 0.5);
                            }
                            if(fabs(theta) >= tol || randomized || doit){
                                ++ndone_iter;
                                if(doprint)
                                    print("     rotating", i, j, theta);

                                double c = cos(theta);
                                double s = sin(theta);
                                drot3(nmo, &dip(i, 0, 0), &dip(j, 0, 0), s, c, 1);
                                drot3(nmo, &dip(0, i, 0), &dip(0, j, 0), s, c, nmo);
                                drot(nmo, &U(i, 0), &U(j, 0), s, c, 1);
                            }
                        }
                    }
                }

                ndone += ndone_iter;
                if(ndone_iter == 0 && tol == thresh){
                    if(doprint)
                        print("Boys localization converged in", ndone,"steps");

                    converged = true;
                    break;
                }
                tol = std::max(0.1 * maxtheta, thresh);
            }

            if(!converged){
                print("warning: boys localization did not fully converge: ", ndone);
            }
            U = transpose(U);

	    bool switched = true;
	    while (switched) {
	      switched = false;
	      for (int i=0; i<nmo; i++) {
		for (int j=i+1; j<nmo; j++) {
		  if (aset[i] == aset[j]) {
		    double sold = U(i,i)*U(i,i) + U(j,j)*U(j,j);
		    double snew = U(i,j)*U(i,j) + U(j,i)*U(j,i);
		    if (snew > sold) {
		      real_tensor tmp = copy(U(_,i));
		      U(_,i) = U(_,j);
		      U(_,j) = tmp;
		      switched = true;
		    }
		  }
		}
	      }
	    }

        // Fix phases.
        for (long i=0; i<nmo; ++i) {
            if (U(i,i) < 0.0) U(_,i).scale(-1.0);
        }

        }

        world.gop.broadcast(U.ptr(), U.size(), 0);
        //END_TIMER(world, "Boys localize");
        return U;
    }
};

// Main mini MolDFT class
class MiniDFT {
private:
    Molecule molecule;
    AtomicBasisSet aobasis;
    BoysLocalization boys;
    XCfunctional xc;
    

public:
    MiniDFT(double              thresh,
            int                 kwavelet,
            int                 truncate_mode,
            double              boxsize,
            const Molecule&     molecule) : molecule(molecule) {
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_k(kwavelet);
        FunctionDefaults<3>::set_cubic_cell(-boxsize,boxsize);
        FunctionDefaults<3>::set_truncate_mode(truncate_mode);

	std::string xc_data;
	//xc_data="lda";
        xc_data = "GGA_X_PBE 1.0 GGA_C_PBE 1.0";
	//xc_data="GGA_X_PBE 1.";
	//xc_data="GGA_C_PBE 1.";
	//xc_data="GGA_X_B88 1.";
	xc.initialize(xc_data, false, world);
    }

    // Make the atomic basis functions
    vector_real_function_3d makeao(World& world, const Molecule& molecule) {
        aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);
        vector_real_function_3d ao(aobasis.nbf(molecule));
        for(int i = 0; i<aobasis.nbf(molecule); ++i) {
            real_functor_3d aofunc(new AtomicBasisFunctor<double>(aobasis.get_atomic_basis_function(molecule, i)));
            ao[i] = real_factory_3d(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(0);
            AtomicBasisFunction atbf = aobasis.get_atomic_basis_function(molecule, i);
        }
        world.gop.fence();
        truncate(world, ao);
        normalize(world, ao);
        print_meminfo(world.rank(), "makeao");
    
        return ao;
    }
   
    // Kinetic energy matrix (for the Fock matrix) 
    real_tensor kinetic_energy_matrix(World& world, const vector_real_function_3d& v)
    {
        std::vector< std::shared_ptr<real_derivative_3d> > gradop = gradient_operator<double,3>(world);
        reconstruct(world, v);
        auto n = v.size();
        real_tensor r(n, n);
        for(int axis = 0;axis < 3;++axis){
            vector_real_function_3d dv = apply(world, *(gradop[axis]), v);
            r += matrix_inner(world, dv, dv, true);
            dv.clear();
        }
        return r.scale(0.5);
    }
   
    // Apply a (local) potential to the orbitals 
    vector_real_function_3d apply_potential(World& world, const real_function_3d& potential, const vector_real_function_3d& psi)
    {
        return mul_sparse(world,potential,psi,FunctionDefaults<3>::get_thresh()*0.01);
    }
   
    // Construct simple LDA exchange-correlation potential 
    real_function_3d make_lda_potential(World& world, const real_function_3d &rho)
    {
        auto vlda = copy(rho);
        vlda.reconstruct();
        vlda.unaryop(xc_lda_potential());
        return vlda;
    }
   
    // Solve free-space Poisson's equation for a given charge density 
    real_function_3d make_coulomb_potential(World& world, const real_function_3d& rho)
    {
        auto thresh = FunctionDefaults<3>::get_thresh();
        auto op = CoulombOperator(world, 1e-10, thresh);
        return op(rho);
    }
   
    // BSH operators 
    vector<poperatorT> make_bsh_operators(World& world, const real_tensor& evals, double shift)
    {
        auto thresh = FunctionDefaults<3>::get_thresh();
        auto nmo = evals.dim(0);
        vector<poperatorT> ops(nmo);
        for(int i = 0;i < nmo; ++i){
            auto eps = evals(i) + shift;
            ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  1e-10, thresh));
        }
        return ops;
    }
   
    // orthogonalize orbitals using Gram-Schmidt 
    template<typename Q>
    void orthogonalize(World& world, std::vector<Function<Q,3> >& psi) {
        compress(world, psi);
        for (unsigned int i = 0; i < psi.size(); i++) {
            Function<Q,3>& psi_i = psi[i];
            psi_i.scale(1.0/psi_i.norm2());
            for (unsigned int j = 0; j < i; j++) {
                Function<Q,3>& psi_j = psi[j];
                Q s = inner(psi_j,psi_i);
                psi_i.gaxpy(1.0,psi_j,-s); // |i> = |i> - |j><j|i>
                psi_i.scale(1.0/psi_i.norm2());
            }
        }
    }
    
    
    // DESTROYS VPSI
    vector_real_function_3d update(World&                            world,
                                   const vector_real_function_3d&    psi,
                                   vector_real_function_3d&          vpsi,
                                   const real_tensor&                e,
                                   int                               iter,
                                   double&                           ravg)
    {
        // psi = - 2 G(E+shift) * (V+shift) psi
        int nmo = psi.size();
    
        // determine shift to make homo <=-0.1
        auto shift = 0.0;
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
    
        // Step restriction
        auto damp = 0.2;
        //if (iter < 10) damp = 0.95;
        //else if (iter < 20) damp = 0.85;
        //else damp = 0.75;
        //print("  shift", shift, "damp", damp, "\n");
    
        printf("      eigenvalue    residual\n");
        ravg = 0.0;
        for (int i=0; i<nmo; i++) {
            double rnorm = (psi[i]-new_psi[i]).norm2();
            ravg += rnorm/(double)nmo;
            printf("%4d  %15.10f  %18.6e\n", i, e[i], rnorm);
            new_psi[i] = damp*psi[i] + (1.0-damp)*new_psi[i];
        }
    
        truncate(world,new_psi);
        normalize(world, new_psi);
        orthogonalize<double>(world, new_psi);
        truncate(world,new_psi);
        normalize(world, new_psi);
        return new_psi;
    }
   
    // Construct charge density from the orbitals 
    real_function_3d make_density(World& world, const vector_real_function_3d& v) {
    
       auto vsq = square(world, v);
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

    // main routine of class
    void doit(World& world) {
    
        aobasis.read_file("6-31g");

        bool do_psp = true;

        auto kwavelet = FunctionDefaults<3>::get_k();
        auto thresh = FunctionDefaults<3>::get_thresh();
        if (world.rank() == 0) printf("wavelet order: %d\n", kwavelet); 
        if (world.rank() == 0) printf("thresh: %15.4e\n", thresh); 

        // Nuclear potential (don't need if using pseudopotential)
        //real_function_3d vnuc = real_factory_3d(world).functor(real_functor_3d(
        //  new NuclearDensityFunctor(molecule))).truncate_mode(0).truncate_on_project();
        auto safety = 0.1;
        auto vtol = FunctionDefaults<3>::get_thresh() * safety;

        real_function_3d vnuc;
        if (! do_psp){
          vnuc = real_factory_3d(world).functor(real_functor_3d(
            new MolecularPotentialFunctor(molecule))).thresh(vtol).truncate_on_project();
          vnuc.reconstruct();
          print("total nuclear charge", vnuc.trace());
        }
        //vnuc = -1.0*make_coulomb_potential(world, vnuc);
    
        // Create pseudopotential
        GTHPseudopotential<double> ppotential(world, molecule);
        if (do_psp){
          ppotential.make_pseudo_potential(world);
        }
    
        // Guess density
        real_function_3d rho = real_factory_3d(world).functor(real_functor_3d(
          new MolecularGuessDensityFunctor(molecule,aobasis))).truncate_on_project();
        double nel = rho.trace();
        if(world.rank() == 0)
            print("guess dens trace", nel);
        int nmo = int(molecule.total_nuclear_charge() + 0.1)/2;
        rho.scale((2.0*nmo)/nel); 
    
        // Make AO basis functions
        auto psi = makeao(world, molecule);
        auto norms = norm2s(world, psi);
    
        double ravg = 10000.0;
        for (int iter=0; iter<100 && ravg > 1e-6; iter++) {

            bool doboys = false;
            if (doboys) {
                if (iter==1) {
                    boys.set_size(psi.size());
                }            
                if (iter > 0) {
                    real_tensor U = boys(world, psi);
                    print(U);
                }
            }
 
            print("\n\n  Iteration",iter,"\n");
            auto rtr = rho.trace();
            print("rho trace:  ", rtr);

            WSTFunctional wstf;
            real_function_3d rho_half = 0.5 * rho;
            std::pair<real_function_3d, double> pair_xc = wstf.apply_xc(world,xc,rho_half);
            real_function_3d vxc = pair_xc.first;
            real_function_3d vxc2 = make_lda_potential(world,rho);
            double exc = pair_xc.second;

            vector_real_function_3d vpsi;
            if (!do_psp){
              real_function_3d v = vnuc + make_coulomb_potential(world,rho) + vxc;
              vpsi = apply_potential(world, v, psi);
            }
            else{
              double enl=0.0;
              tensorT occ = tensorT(nmo);
              for(int i = 0;i < nmo;++i)
                  occ[i] = 1.0;
              auto v = make_coulomb_potential(world,rho) + vxc + ppotential.vlocalp;
              vpsi = ppotential.apply_potential(world, v, psi, occ, enl);
            }

            bool doplots = false;
            /*if (doplots) {
                char rho_name[25]; sprintf(rho_name, "rho_%d.dat", iter);
                char vxc_name[25]; sprintf(vxc_name, "vxc_%d.dat", iter);
                char vxc2_name[25]; sprintf(vxc2_name, "vxc2_%d.dat", iter);
                char v_name[25]; sprintf(v_name, "v_%d.dat", iter);
                int ppnts = 5001; 
                plot_line(rho_name, ppnts, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, rho); 
                plot_line(vxc_name, ppnts, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vxc); 
                plot_line(vxc2_name, ppnts, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, vxc2); 
                plot_line(v_name, ppnts, {0.0,0.0,-50.0}, {0.0,0.0,50.0}, v); 
            }*/
   
            auto ke_mat = kinetic_energy_matrix(world, psi);
            auto pe_mat = matrix_inner(world, psi, vpsi, true);
            auto ov_mat = matrix_inner(world, psi, psi, true);
    
            auto fock = ke_mat + pe_mat;
            fock = 0.5*(fock + transpose(fock));
    
            real_tensor c;
            tensor_real e;
            sygv(fock, ov_mat, 1, c, e);
            double rfactor = 1e-3;
            for (int i = 0; i < fock.dim(0); i++)
            {
              auto thresh = FunctionDefaults<3>::get_thresh();
              for (int j = 0; j < fock.dim(1); j++)
              {
                if (std::abs(ov_mat(i,j)) < thresh*rfactor) ov_mat(i,j) = 0.0;
                if (std::abs(ke_mat(i,j)) < thresh*rfactor) ke_mat(i,j) = 0.0;
                if (std::abs(pe_mat(i,j)) < thresh*rfactor) pe_mat(i,j) = 0.0;
                if (std::abs(fock(i,j)) < thresh*0.1) fock(i,j) = 0.0;
              }
            }
            if (world.rank() == 0) {
              print("Kinetic (real part)"); 
              print(ke_mat);
              print("Potential (real part)"); 
              print(pe_mat);
              print("Fock (real part)"); 
              print(fock);
              print("eigenvalues"); print(e);
            }
    
            if (iter == 0) {
                c = copy(c(_,Slice(0,nmo-1))); // truncate to occupied states
                e = e(Slice(0,nmo-1));
            }
    
            auto trantol = vtol / std::min(30.0, double(psi.size()));
            psi = transform(world, psi, c, trantol, true);
            vpsi = transform(world, vpsi, c, trantol, true);
    
            psi = update(world, psi, vpsi, e, iter, ravg);
            truncate(world, psi);
    
            auto rthresh = 50.0*thresh;
            if (ravg < rthresh) {
              print("reprojecting ...");
              kwavelet += 2;
              thresh /= 100.0;
              FunctionDefaults<3>::set_k(kwavelet);
              FunctionDefaults<3>::set_thresh(thresh);
              if (!do_psp){
                vnuc = madness::project(vnuc, kwavelet, thresh, true);}
              else{
                ppotential.reproject(kwavelet, thresh);}
              for (int i = 0; i < nmo; i++) psi[i] = madness::project(psi[i], kwavelet, thresh, true);  
            }
    
            //int md = psi[3].max_depth();
            //print("max depth:  ", md);
    
            rho = make_density(world, psi);
        }
    }

};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(9);

    Molecule molecule;

    //print("Env: MADNESS_HAS_LIBXC =      ", MADNESS_HAS_LIBXC);

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
    
    // Argon atom
    //molecule.add_atom(0, 0, 0, 8.0, 18);
    
    // Calcium atom
    //molecule.add_atom(0, 0, 0, 2.0, 20);

    // Silicon atom
    //molecule.add_atom(0, 0, 0, 4.0, 14);
    
    // Carbon atom
    //molecule.add_atom(0, 0, 0, 4.0, 6);
    
    // Oxygen atom
    //molecule.add_atom(0.0, 0.0, 2.286, 6.0, 8); 
    
    // Beryllium atom (at arbitrary location)
    //molecule.add_atom(1.23943, -0.3422, 5, 2.0, 4);
    
    // Water (doesn't work)
    //molecule.add_atom( 1.4375, 0, 1.15, 1.0, 1);
    //molecule.add_atom(-1.4375, 0, 1.15, 1.0, 1);
    //molecule.add_atom(0, 0, 0, 8.0, 8);
   
    // H2 
    molecule.add_atom(0.0, 0.0,  0.7  , 1.0, 1); 
    molecule.add_atom(0.0, 0.0, -0.7  , 1.0, 1); 

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
   
    // MRA parameters
    double thresh = 1e-4;
    int kwavelet = 6;
    int truncate_mode = 0; 
    double L = 50.0;

    MiniDFT dft(thresh, kwavelet, truncate_mode, L, molecule);
    dft.doit(world);
    return 0;
}

