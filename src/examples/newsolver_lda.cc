
#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <type_traits>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

#include <madness/mra/lbdeux.h>
#include <madness/mra/qmprop.h>

#include <madness/misc/misc.h>
#include <madness/misc/ran.h>

#include <madness/tensor/systolic.h>
#include <madness/tensor/solvers.h>
#include <madness/tensor/elem.h>


#include <chem/xcfunctional.h>

#include <madness/mra/legendre.h>

#include <chem/xcfunctional.h>


using namespace madness;


typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Function<std::complex<double>,3> complex_functionT;
typedef std::vector<complex_functionT> cvecfuncT;
typedef Convolution1D<double_complex> complex_operatorT;


//static const double R = 1.4;    // bond length
//static const double L = 64.0*R; // box size
static const double L = 50.0; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-5; // precision
static const int maxiter=10;
static const int maxiterp=10;

std::vector< std::shared_ptr<real_derivative_3d> > gradop;
XCfunctional xc;

//prod/// simple projector class for 1- and 2-particle projectors
//prodtemplate<typename T, std::size_t NDIM>
//prodclass Projector {
//prod
//prod    int particle_;
//prod    std::vector<Function<T,NDIM> > p_;
//prod
//prodpublic:
//prod
//prod    Projector() : p_(std::vector<Function<T,NDIM> >()) {}
//prod
//prod    /// simple constructor with only one orbital to project out
//prod    Projector(const Function<T,NDIM>& p, const int particle=0)
//prod            : particle_(particle), p_(std::vector<Function<T,NDIM> >(1,p)) {
//prod        MADNESS_ASSERT(particle_==0 or particle_==1);
//prod        MADNESS_ASSERT(p_.size()>0);
//prod    }
//prod
//prod    /// constructor with a set of orbitals to project out
//prod    Projector(const std::vector<Function<T,NDIM> >& p, const int particle=0) : particle_(particle), p_(p) {
//prod        MADNESS_ASSERT(particle_==0 or particle_==1);
//prod        MADNESS_ASSERT(p_.size()>0);
//prod    }
//prod
//prod    int& particle() {return particle_;}
//prod    const int& particle() const {return particle_;}
//prod
//prod    /// get a const reference to the orbitals
//prod    const std::vector<Function<T,NDIM> >& p() const {return p_;}
//prod
//prod    /// project f on p: |result> =  | p><p | f>
//prod    template<std::size_t FDIM>
//prod    typename std::enable_if<NDIM==FDIM, Function<T,FDIM> >::type
//prod    operator()(const Function<T,FDIM>& f) const {
//prod
//prod        const double ovlp=inner(f,p_[0]);
//prod        Function<T,NDIM> sum=ovlp*p_[0];
//prod
//prod        for (unsigned int i=1; i<p_.size(); ++i) {
//prod            const double ovlp2=inner(f,p_[i]);
//prod            sum=(sum+ovlp2*p_[i]).truncate().reduce_rank();
//prod        }
//prod        return sum;
//prod    }
//prod
//prod    /// project p out of f: |result(1,2)> = sum_p | p(1)><p(1) | f(1,2)>
//prod    template<std::size_t FDIM>
//prod    typename std::enable_if<2*NDIM==FDIM, Function<T,FDIM> >::type
//prod    operator()(const Function<T,FDIM>& f) const {
//prod        real_function_6d sum=real_factory_6d(p_.begin()->world());
//prod        for (unsigned int i=0; i<p_.size(); ++i) {
//prod            const real_function_3d pf2=f.project_out(p_[i],particle_);
//prod            real_function_6d tmp;
//prod            MADNESS_EXCEPTION("Projector class: the hartree product is inaccurate -- don't use it",1);
//prod            if (particle_==0) tmp=hartree_product(p_[i],pf2);
//prod            else tmp=hartree_product(pf2,p_[i]);
//prod            sum=(sum+tmp);
//prod        }
//prod        sum.truncate();
//prod        return sum;
//prod    }
//prod};


double make_dft_energy(World & world, const vecfuncT& vf, int ispin)
{
	functionT vlda = multiop_values<double, xc_functional, 3>(xc_functional(xc, ispin), vf);
	return vlda.trace();
}

functionT make_dft_potential(World & world, const vecfuncT& vf, int ispin, int what)
{
	return multiop_values<double, xc_potential, 3>(xc_potential(xc, ispin, what), vf);
}

functionT make_dft_kernel(World & world, const vecfuncT& vf, int ispin, int what)
{
	return multiop_values<double, xc_kernel, 3>(xc_kernel(xc, ispin, what), vf);
}

//prodstatic double guess(const coord_3d& r) {
//prod    const double x=r[0], y=r[1], z=r[2];
//prod    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
//prod            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
//prod}
//prod
//prodstatic double V(const coord_3d& r) {
//prod    const double x=r[0], y=r[1], z=r[2];
//prod    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
//prod           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
//prod}

static double guess(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];

 	//return (x*x+y*y+z*z);
 	return (2.0*exp(-sqrt(x*x+y*y+z*z+1e-8)));
 	//return (2.0*exp(-sqrt(x*x+y*y+z*z+1e-8))+1e-8);
}

static double V(const coord_3d& r) {
	const double x=r[0], y=r[1], z=r[2];
	return  -2.0/sqrt(x*x+y*y+z*z+1e-8);
}

double rifunction(const coord_3d& r) {
    return r[2]; //z
}

double iterate_ground(World& world, NonlinearSolver& solver, 
                      functionT& V, functionT& psi, 
                      double& eps) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    functionT Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    functionT tmp = op(Vpsi).truncate();
    double norm = tmp.norm2();
    functionT r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi = solver.update(psi, r);
    psi.scale(1.0/psi.norm2());
    eps = eps_new;
    return rnorm;
}

double iterate_excite(World& world, NonlinearSolver& solver, 
                      functionT& V, functionT& psi, 
                      functionT& dpsi, double& eps, 
                      functionT& ri) {
    real_convolution_3d du = CoulombOperator(world, 0.001, 1e-6); 
    functionT Vdpsi = (V*dpsi);
// LDA
    functionT rho = psi*psi;
    vecfuncT vf;
    vf.push_back(rho);
    functionT fxc = make_dft_kernel(world, vf, 0, 0).truncate();

    functionT rhod = 2*dpsi*psi; //pert density

    //functionT rhs = Vdpsi + 2.*psi*du(psi*dpsi)  + ri*psi;
    functionT Gampsi = 2.*psi*du(rhod) + psi*fxc*rhod*2.*constants::pi + ri*psi;
    //Projector<double,3> rho0(psi);
    double ovlp = inner(psi,Gampsi);
    functionT rhs = Vdpsi + Gampsi - psi*ovlp;
//prod    functionT rhs = Vdpsi + Gampsi ;
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);

    rhs.scale(-2.0).truncate();
    //functionT r = op(rhs)*.7+dpsi*.3 - dpsi;
    functionT r = op(rhs) - dpsi;
    dpsi = solver.update(dpsi, r);
    dpsi.truncate();
    double dpsi_error = r.norm2();
    if (world.rank() == 0) {
        //printf("dpsi error = %f \t ovl = %f \n", dpsi_error, ovlp);
        printf("dpsi error = %f \n", dpsi_error);
    }
    return dpsi_error;
 }

//prod// This class is used to store information for the
//prod// non-linear solver used by the dynamic polarizability
//prodstruct F {
//prod    functionT x, y;
//prod
//prod    F(const functionT& x, const functionT& y) 
//prod        : x(x), y(y)
//prod    {}
//prod
//prod    F operator-(const F& b) {
//prod        return F(x-b.x, y-b.y);
//prod    }
//prod
//prod    F operator+=(const F& b) { // Operator+= necessary
//prod        x += b.x; y += b.y;
//prod        return *this;
//prod    }
//prod
//prod    F operator*(double a) { // Scale by a constant necessary
//prod        return F(x*a,y*a);
//prod    }
//prod};
//prod
//proddouble inner(const F& a, const F& b) {
//prod    return inner(a.x,b.x) + inner(a.y,b.y);
//prod}
//prod
//prod// The default constructor for functions does not initialize
//prod// them to any value, but the solver needs functions initialized
//prod// to zero for which we also need the world object.
//prodstruct allocator {
//prod    World& world;
//prod
//prod    allocator(World& world) : world(world) {}
//prod
//prod    F operator()() {
//prod        return F(functionT(world),functionT(world));
//prod    }
//prod};
//prod
//prod
//prodtemplate <class solverT>
//proddouble iterate_xy(World& world, solverT& solver, functionT& V, 
//prod                  functionT& psi, double& eps, 
//prod                  functionT& ri, functionT& x, 
//prod                  functionT& y, double& omega) {
//prod    real_convolution_3d gOp= BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
//prod    real_convolution_3d uop = CoulombOperator(world, 0.001, 1e-6);
//prod
//prod    functionT gp2 = uop(y*psi);
//prod    functionT gp3 = uop(x*psi);
//prod// LDA
//prod    vecfuncT vf;
//prod    functionT rho = psi*psi;
//prod    vf.push_back(rho);
//prod    functionT fxc = make_dft_kernel(world, vf, 0, 0);
//prod
//prod    functionT pp = (gp3 + gp2 + fxc*x + fxc*y + ri) * psi;
//prod
//prod    functionT xrhs = V*x + pp;
//prod    x.truncate();
//prod    xrhs.gaxpy(1.0, x, -omega,false);
//prod
//prod    xrhs.scale(-2.0).truncate();
//prod    functionT new_x = gOp(xrhs).truncate();
//prod    new_x = new_x - inner(psi,new_x)*psi;
//prod    double xerr = (x - new_x).norm2();
//prod    
//prod    functionT yrhs = V*y + pp;
//prod    y.truncate();
//prod    yrhs.gaxpy(1.0, y, omega,false);
//prod
//prod    yrhs.scale(-2.0).truncate();
//prod    functionT new_y = gOp(yrhs).truncate();
//prod    new_y = new_y - inner(psi,new_y)*psi;
//prod    double yerr = (y - new_y).norm2();
//prod    
//prod    F xy = solver.update(F(x,y), F(new_x-x, new_y-y));
//prod    
//prod    x = xy.x.truncate();
//prod    y = xy.y.truncate();
//prod    print("dynamic", xerr, yerr);
//prod    return xerr+yerr;
//prod}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    
    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);  
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);
    
    if (world.rank() == 0) print("\n  Solving for the DFT/LDA wave function\n");
    functionT Vnuc = real_factory_3d(world).f(V);
    functionT psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    std::string xc_data;
    xc_data="lda";
    xc.initialize(xc_data, false);

//    poperatorT coulop;
//    coulop = poperatorT(CoulombOperatorPtr(world, 1e-10, thresh));
    poperatorT coulop;
    coulop = poperatorT(CoulombOperatorPtr(world, 1e-10, thresh));
    gradop = gradient_operator<double,3>(world);
    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    double eps = -0.5;
    {
      NonlinearSolver solver;
      for (int iter=0; iter<maxiter; iter++) {
        functionT rho = psi*psi;
// LDA
        vecfuncT vf;
	vf.push_back(rho);
        //functionT fxc =  make_dft_kernel(world, vf, 0, 0);
        //functionT vxc =  fxc*rho*2;
        functionT vxc =  make_dft_potential(world, vf, 0, 0);

        functionT potential = Vnuc + 2.*op(rho).truncate() + vxc.truncate();
        double err  = iterate_ground(world, solver, potential, psi, eps);
	if (err < thresh) break;
      }
    }
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        functionT dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }

    functionT rho = square(psi);
    rho.reconstruct();
    vecfuncT vf;
    vf.push_back(rho);
    double exc=make_dft_energy(world, vf, 0); // Exc
    double two_electron_energy = 2.*inner(op(rho),rho); // <u|rho> = <phi phi | 1/r12 | phi phi>
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho); // <V|rho> = <phi|V|phi>
    double densii = inner(psi,psi); // <rho> 
   // double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy + 
            nuclear_attraction_energy +  exc;
   //     nuclear_attraction_energy + nuclear_repulsion_energy + exc;
    double virial = (nuclear_attraction_energy + two_electron_energy ) / kinetic_energy;
    //double virial = (nuclear_attraction_energy + two_electron_energy + nuclear_repulsion_energy) / kinetic_energy;

    if (world.rank() == 0) {
        print("");
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
//        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("                 XC energy ", exc);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
        print("                       dee ", densii);
    }

//    plot_line("psi.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, psi);


    world.gop.fence();

#if 1
    if (world.rank() == 0) print("\nSolving for the static response function\n");

    functionT ri  = real_factory_3d(world).f(rifunction);
//    functionT dpsi(world); //zero
//    dpsi = 1*psi;
    functionT dpsi = real_factory_3d(world); //zero
    rho = square(psi).truncate();
// LDA
    //vecfuncT vf;
    //vf.push_back(rho);
    functionT vxc =  make_dft_potential(world, vf, 0, 0);
    functionT potential = Vnuc + 2.*op(rho).truncate() + vxc.truncate();
    //functionT potential = Vnuc + 2.*op(rho).truncate() + vxc.truncate();
    
    {
        NonlinearSolver solver;
        for(int iter=1; iter<=maxiterp; iter++) {
            double err = iterate_excite(world, solver, potential, psi, dpsi, eps, ri);
            if (err < 10*thresh) break;
        }
    }

//prod    plot_line("dpsi.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, dpsi);
//prod    plot_line("rho.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, rho);

    double h1 = 0.0;
    for (int axis=0; axis<3; axis++) {
      real_derivative_3d D = free_space_derivative<double,3>(world, axis);
      functionT nadp = D(dpsi);
      h1 += 0.5*inner(nadp,nadp);
    }
    double h2 = inner(Vnuc*dpsi,dpsi);
    double ans1 = h1 + h2;

    functionT e1 = dpsi*psi;
    functionT e2 = psi*psi;
    functionT e3 = dpsi*dpsi;
    double ans2 = inner(e1,op(e1));
    double ans3 = inner(e3,op(e2));
    double ans4 = inner(dpsi,dpsi);
    double ans5 = inner(dpsi,ri*psi);
    double d2Ea = (4*ans1) + (8*ans2) + (4*ans3) - (4*(eps*ans4)) + (8*ans5);
    double d2Eb =  4*ans5;

    if (world.rank() == 0) {
        print("");
        print("< dpsi | h | dpsi >", ans1);
        print("< dpsi dpsi | psi psi >", ans2);
        print("< dpsi psi | dpsi psi >", ans3);
        print("< dpsi | dpsi >", ans4);
        print("< dpsi | r | psi > ",ans5);
        print("");
        print(" < 1 |  V - E1 | 0 > =", d2Ea);
        print("-< 0 | H0 - E0 | 0 > =", d2Eb);
        print("variational estimate =", 2*d2Ea - d2Eb);
    }
    
#endif
#if 0
    // For first frequency use zero as the initial guess but at subsequent
    // frequencies use the previous solution as the guess.
    functionT x = real_factory_3d(world); //zero
    functionT y = real_factory_3d(world); //zero
    for(int j=0; j<=1; j++) {
        //double omega = 0.365 + (j*0.005);
        //double omega = 0.1 + (j*0.01);
        double omega = 0.0;
        if (world.rank() == 0) print("\nSolving for the dynamic response function with omega =", omega,"\n");

        XNonlinearSolver<F,double,allocator> solver = XNonlinearSolver<F,double,allocator>(allocator(world));
        
        for(int iter=1; iter<=maxiterp; iter++) {
            double err = iterate_xy(world, solver, potential, psi, eps, ri, x, y, omega);
            if (err < 10*thresh) break;
        }

        functionT drho = (x*psi)+(psi*y);
        double alpha_dynamic = -2*inner(ri,drho);    
        if (world.rank() == 0) 
            print("\nalpha_dynamic omega = ",omega," alpha = ",alpha_dynamic);

        char fname[32];
        sprintf(fname,"x_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, x);
        sprintf(fname,"y_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, y);
        sprintf(fname,"drho_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, drho);
    }
#endif
    
    finalize();
    return 0;
}
