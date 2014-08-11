//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include "nonlinsol.h"

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 64.0*R; // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-4; // precision

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

double iterate_ground(World& world, NonlinearSolver& solver, real_function_3d& V, real_function_3d& psi, double& eps) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();
    real_function_3d tmp = op(Vpsi);
    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
    }
    psi = solver.update(psi, r);
    psi.truncate();
    psi.scale(1.0/psi.norm2());
    eps = eps_new;
    return rnorm;
}

double iterate_excite(World& world, NonlinearSolver& solver, real_function_3d& V, real_function_3d& psi, real_function_3d& dpsi, double& eps, real_function_3d& ri) {
    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    real_convolution_3d du = CoulombOperator(world, 0.001, 1e-6);
    real_function_3d rhs = V*dpsi + (du(psi*dpsi).scale(2.0) + ri)*psi;

    rhs.scale(-2.0).truncate();
    real_function_3d r = op(rhs) - dpsi;
    double rnorm = r.norm2();
    dpsi = solver.update(dpsi, r);
    dpsi.truncate();
    if (world.rank() == 0) {
        print("dpsi error = ", rnorm);
    }
    return rnorm;
 }

double iterate_xy(World& world, real_function_3d& V, real_function_3d& psi, real_function_3d& dpsi, double& eps, real_function_3d& ri, real_function_3d& x, real_function_3d& y, double omega) {
  real_convolution_3d uop = CoulombOperator(world, 0.001, 1e-6);
  real_convolution_3d gOpx = BSHOperator3D(world, sqrt(-2*(eps+omega)), 0.001, 1e-6);
  real_convolution_3d gOpy = BSHOperator3D(world, sqrt(-2*(eps-omega)), 0.001, 1e-6);
  real_function_3d gp2 = uop(y*psi);
  real_function_3d gp3 = uop(x*psi);

  real_function_3d pp = psi * (gp2 + gp3 + ri);

  real_function_3d xrhs = V*x + pp;
  xrhs.scale(-2.0).truncate();
  real_function_3d new_x = gOpx(xrhs);
  new_x = new_x - inner(psi,new_x)*psi;

  real_function_3d yrhs = V*y + pp;
  yrhs.scale(-2.0).truncate();
  real_function_3d new_y = gOpy(yrhs);
  new_y = new_y - inner(psi,new_y)*psi;

  double errx = (x-new_x).norm2();
  double erry = (y-new_y).norm2();

  if (world.rank() == 0)
      print("errx =", errx, "erry =", erry);

  x = new_x.truncate();
  y = new_y.truncate();

  return errx+erry;
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

    real_function_3d Vnuc = real_factory_3d(world).f(V);
    real_function_3d psi  = real_factory_3d(world).f(guess);
    psi.truncate();
    psi.scale(1.0/psi.norm2());

    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    double eps = -0.6;

    if (world.rank() == 0)
        print("\n  Solving for the ground state HF orbital\n");
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

    world.gop.fence();

    // for(int i=-100; i<=100; i++) {
    //   coord_3d r(0.0);
    //   r[2] = i*0.1;
    //   print(r[2],psi(r));
    // }

    if (world.rank() == 0)
        print("\n  Solving for the static response function (dphi/dFz)\n");

    real_function_3d ri  = real_factory_3d(world).f(rifunction);
    real_function_3d dpsi(world); //zero
    rho = square(psi).truncate();
    real_function_3d potential = Vnuc + op(rho).truncate();
    {
        NonlinearSolver solver;
        for(int iter=1; iter<=20; iter++) {
            double err = iterate_excite(world, solver, potential, psi, dpsi, eps, ri);
            if (err < thresh*10.0) break;
        }
    }
    // for (int i=-100; i<=100; i++) {
    //   coord_3d r(0.0);
    //   r[2] = i*0.1;
    //   print(r[2],dpsi(r));
    // }

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
    double d2Eb = 4*ans5;
    double d2E  = 2*d2Eb - d2Ea;
    double alpha_static = -d2E;

    if (world.rank() == 0) {
        print("< dpsi | h | dpsi >", ans1);
        print("< dpsi dpsi | psi psi >", ans2);
        print("< dpsi psi | dpsi psi >", ans3);
        print(" eps ",eps);
        print("< dpsi | dpsi >", ans4);
        print("< dpsi | r | psi > ",ans5);
        print("E2a = -<1|H0 - E0|1> =", d2Ea);
        print("E2b =  <1| V - E1|0> =", d2Eb);
        print("E2 variational form  =", d2E);
        print(" alpha_static ",alpha_static);
    }

    double omega = 0.0;
    for(int j=0; j<=4; j++) {
      omega = 0.365 + (j*0.005);

      if (world.rank() == 0)
          print("\n  Solving for the dynamic response function at omega =",omega,"\n");

      real_function_3d x = real_factory_3d(world); //zero
      real_function_3d y = real_factory_3d(world); //zero
      for(int iter=1; iter<=20; iter++) {
          double err = iterate_xy(world, potential, psi, dpsi, eps, ri, x, y, omega);
          if (err < 10*thresh) break;
      }
      //for (int i=-100; i<=100; i++) {
      //coord_3d r(0.0);
      //r[2] = i*0.1;
	//print(r[2],x(r));
	//print(r[2],y(r));
      //}

      //print("alpha_dynamic +w", 2.0*inner(ri,x*psi));
      //print("alpha_dynamic -w", 2.0*inner(ri,y*psi));

      real_function_3d drho = (x*psi)+(psi*y);
      double alpha_dynamic = -2*inner(ri,drho);
      print("alpha_dynamic omega = ",omega," alpha = ",alpha_dynamic);
    }

    finalize();
    return 0;
}
