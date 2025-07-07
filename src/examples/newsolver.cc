#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-5; // precision

static double guess(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(-sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8)));
}

static double V(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/sqrt(x*x+y*y+(z-R/2)*(z-R/2)+1e-8)+
           -1.0/sqrt(x*x+y*y+(z+R/2)*(z+R/2)+1e-8);
}

double rifunction(const coord_3d& r) {
    return r[2]; //z
}

double iterate_ground(World& world, NonlinearSolver& solver, 
                      real_function_3d& V, real_function_3d& psi, 
                      double& eps) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_function_3d tmp = op(Vpsi).truncate();
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
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
                      real_function_3d& V, real_function_3d& psi, 
                      real_function_3d& dpsi, double& eps, 
                      real_function_3d& ri) {
    real_convolution_3d du = CoulombOperator(world, 0.001, 1e-6); 
    real_function_3d Vdpsi = (V*dpsi);
    real_function_3d rhs = Vdpsi + (psi*2*du(psi*dpsi)) + (ri*psi);
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);

    rhs.scale(-2.0).truncate();
    real_function_3d r = op(rhs) - dpsi;
    dpsi = solver.update(dpsi, r);
    dpsi.truncate();
    double dpsi_error = r.norm2();
    if (world.rank() == 0) {
        print("dpsi error = ", dpsi_error);
    }
    return dpsi_error;
 }

// This class is used to store information for the
// non-linear solver used by the dynamic polarizability
struct F {
    real_function_3d x, y;

    F(const real_function_3d& x, const real_function_3d& y) 
        : x(x), y(y)
    {}

    F operator-(const F& b) const {
        return F(x-b.x, y-b.y);
    }

    F operator+=(const F& b) { // Operator+= necessary
        x += b.x; y += b.y;
        return *this;
    }

    F operator*(double a) const { // Scale by a constant necessary
        return F(x*a,y*a);
    }
};

double inner(const F& a, const F& b) {
    return inner(a.x,b.x) + inner(a.y,b.y);
}

// The default constructor for functions does not initialize
// them to any value, but the solver needs functions initialized
// to zero for which we also need the world object.
struct allocator {
    World& world;

    allocator(World& world) : world(world) {}

    F operator()() {
        return F(real_function_3d(world),real_function_3d(world));
    }
};


/// solve the dynamic response equations

/// @param[in]	world	the world
/// @param[in]	solver	the KAIN solver
/// @param[in]	V		the local potential (here: V_nuc + 1/2 J )
/// @param[in]	psi		the orbital
/// @param[in]	eps		orbital energy for psi
/// @param[in]	ri		the external perturbation (here: z)
/// @param[in,out]	x	the x part of the response vector
/// @param[in,out]	y	the y part of the response vector
/// @param[in]	omega	the frequency of the external perturbation
/// @return		the current error in the residual of the response equations
template <class solverT>
double iterate_xy(World& world, solverT& solver, const real_function_3d& V,
                  const real_function_3d& psi,
                  double& eps, const real_function_3d& ri, real_function_3d& x,
                  real_function_3d& y, const double omega) {
    real_convolution_3d gOpx = BSHOperator3D(world, sqrt(-2*(eps+omega)), 0.001, 1e-6);
    real_convolution_3d gOpy = BSHOperator3D(world, sqrt(-2*(eps-omega)), 0.001, 1e-6);
    real_convolution_3d uop = CoulombOperator(world, 0.001, 1e-6);
    real_function_3d gp2 = uop(y*psi);
    real_function_3d gp3 = uop(x*psi);
    real_function_3d pp = (gp3 + gp2 + ri) * psi;
    
    real_function_3d xrhs = V*x + pp;
    xrhs.scale(-2.0).truncate();
    real_function_3d new_x = gOpx(xrhs).truncate();
    new_x = new_x - inner(psi,new_x)*psi;
    double xerr = (x - new_x).norm2();
    
    real_function_3d yrhs = V*y + pp;
    yrhs.scale(-2.0).truncate();
    real_function_3d new_y = gOpy(yrhs).truncate();
    new_y = new_y - inner(psi,new_y)*psi;
    double yerr = (y - new_y).norm2();
    
    F xy = solver.update(F(x,y), F(new_x-x, new_y-y));
    
    x = xy.x.truncate();
    y = xy.y.truncate();
    print("dynamic", xerr, yerr);
    return xerr+yerr;
}


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
    
    if (world.rank() == 0) print("\n  Solving for the HF wave function\n");
    real_function_3d Vnuc = real_factory_3d(world).f(V);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    double eps = -0.6;
    {
      NonlinearSolver solver;
      for (int iter=0; iter<10; iter++) {
        real_function_3d rho = square(psi).truncate();
        real_function_3d potential = Vnuc + op(rho).truncate();
        double err  = iterate_ground(world, solver, potential, psi, eps);
	if (err < thresh) break;
      }
    }
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }

    real_function_3d rho = square(psi);
    double two_electron_energy = inner(op(rho),rho); // <u|rho> = <phi phi | 1/r12 | phi phi>
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho); // <V|rho> = <phi|V|phi>
    double nuclear_repulsion_energy = 1.0/R;
    double total_energy = kinetic_energy + two_electron_energy + 
        nuclear_attraction_energy + nuclear_repulsion_energy;
    double virial = (nuclear_attraction_energy + two_electron_energy + nuclear_repulsion_energy) / kinetic_energy;

    if (world.rank() == 0) {
        print("");
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("       Two-electron energy ", two_electron_energy);
        print(" Nuclear  repulsion energy ", nuclear_repulsion_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", virial);
    }

    plot_line("psi.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, psi);

    world.gop.fence();

    if (world.rank() == 0) print("\nSolving for the static response function\n");

    real_function_3d ri  = real_factory_3d(world).f(rifunction);
    real_function_3d dpsi(world); //zero
    rho = square(psi).truncate();
    real_function_3d potential = Vnuc + op(rho).truncate();
    {
        NonlinearSolver solver;
        for(int iter=1; iter<=20; iter++) {
            double err = iterate_excite(world, solver, potential, psi, dpsi, eps, ri);
            if (err < 10*thresh) break;
        }
    }

    plot_line("dpsi.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, dpsi);
    plot_line("rho.dat", 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, rho);

    double h1 = 0.0;
    for (int axis=0; axis<3; axis++) {
      real_derivative_3d D = free_space_derivative<double,3>(world, axis);
      real_function_3d nadp = D(dpsi);
      h1 += 0.5*inner(nadp,nadp);
    }
    double h2 = inner(Vnuc*dpsi,dpsi);
    double ans1 = h1 + h2;

    real_function_3d e1 = dpsi*psi;
    real_function_3d e2 = psi*psi;
    real_function_3d e3 = dpsi*dpsi;
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
    
    // For first frequency use zero as the initial guess but at subsequent
    // frequencies use the previous solution as the guess.
    real_function_3d x = real_factory_3d(world); //zero
    real_function_3d y = real_factory_3d(world); //zero
    for(int j=0; j<=4; j++) {
        double omega = 0.365 + (j*0.005);
        if (world.rank() == 0) print("\nSolving for the dynamic response function with omega =", omega,"\n");

        XNonlinearSolver<F,double,allocator> solver = XNonlinearSolver<F,double,allocator>(allocator(world));
        
        for(int iter=1; iter<=20; iter++) {
            double err = iterate_xy(world, solver, potential, psi, eps, ri, x, y, omega);
            if (err < 10*thresh) break;
        }

        real_function_3d drho = (x*psi)+(psi*y);
        double alpha_dynamic = -2*inner(ri,drho);    
        if (world.rank() == 0) 
            print("\nalpha_dynamic omega = ",omega," alpha = ",alpha_dynamic);

        const std::size_t bufsize=32;
        char fname[bufsize];
        snprintf(fname, bufsize,"x_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, x);
        snprintf(fname, bufsize,"y_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, y);
        snprintf(fname, bufsize,"drho_%6.4f.dat", omega);
        plot_line(fname, 1001, {0.0,0.0,-20.0}, {0.0,0.0,20.0}, drho);
    }
    
    finalize();
    return 0;
}
