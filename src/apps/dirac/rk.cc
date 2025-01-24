#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/world/vector.h>
#include <madness/mra/nonlinsol.h>
#include <cstdio>
#include "DKops.h"

void apbar_fit(double dx, double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta);
void pbar_fit(double dx, double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta);
void a_fit(double dx, double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta);
void tbar_fit(double dx, double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta);
void bshrel_fit(double epsilon, double dx, double thresh, double quadacc, std::vector<double>& coeffs, std::vector<double>& expnts, double& Cdelta);

using namespace madness;

// Needed for rebalancing
template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() < 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

static double Z = -1; 

// Nuclear half-charge radius of H (proton) is 1.3e-5 and the exponent
// of its finite-size Gaussian approximation is 2.12e9.  So this
// choice for dx resolves one-two orders of magnitude finer scale.
static const double op_dx=1e-16; 
static const double op_quadacc=1e-8; 
static const double op_thresh=1e-8;
    
static const long k = 10;        // wavelet order
static const double thresh = 1e-8; // precision
static const double L = 50.0;   // box size

static const bool FINITENUC = false;

// Pick ONE of the following three
static const bool NONREL = false;
static const bool RK = false;
static const bool DK1 = true;

class FiniteNucleusPotential : public FunctionFunctorInterface<double,3> {
    const double xi;

    static double expnt(int Z) {
        const Vector<double,109> xi
            {2.1248239171e+09,1.1671538870e+09,8.9266848806e+08,7.8788802914e+08,7.1178709563e+08,6.8077502929e+08,6.2865615725e+08,5.8631436655e+08,5.3546911034e+08,5.2105715255e+08,4.8349721509e+08,4.7254270882e+08,4.4335984491e+08,4.3467748823e+08,4.1117553148e+08,4.0407992047e+08,3.8463852873e+08,3.5722217300e+08,3.6228128110e+08,3.5722217300e+08,3.3451324570e+08,3.2263108827e+08,3.1181925878e+08,3.0842641793e+08,2.9881373610e+08,2.9578406371e+08,2.8716667270e+08,2.8996391416e+08,2.7665979354e+08,2.7419021043e+08,2.6267002737e+08,2.5235613399e+08,2.5042024280e+08,2.4130163719e+08,2.4305454351e+08,2.3461213272e+08,2.3301551109e+08,2.2839354730e+08,2.2690621893e+08,2.2544431039e+08,2.2120420724e+08,2.1458511597e+08,2.1458511597e+08,2.0965270287e+08,2.0846586999e+08,2.0500935221e+08,2.0389047621e+08,1.9648639618e+08,1.9548577691e+08,1.9067718154e+08,1.8975246242e+08,1.8193056289e+08,1.8444240538e+08,1.8030529331e+08,1.7950688281e+08,1.7565009043e+08,1.7490463170e+08,1.7416744147e+08,1.7343837120e+08,1.7129844956e+08,1.7060044589e+08,1.6591550422e+08,1.6527352089e+08,1.6215880671e+08,1.6155419421e+08,1.5977529080e+08,1.5977529080e+08,1.5636673634e+08,1.5581702004e+08,1.5314257850e+08,1.5262201512e+08,1.5008710340e+08,1.4959325643e+08,1.4813689532e+08,1.4671710337e+08,1.4442808782e+08,1.4398142103e+08,1.4309883584e+08,1.4223027307e+08,1.4011788914e+08,1.3888925203e+08,1.3768840081e+08,1.3729411599e+08,1.3729411599e+08,1.3690277000e+08,1.3242350205e+08,1.3206733609e+08,1.3101367628e+08,1.3066730974e+08,1.2897067480e+08,1.2930539512e+08,1.2700881714e+08,1.2733038109e+08,1.2512299012e+08,1.2543221826e+08,1.2420711085e+08,1.2420711085e+08,1.2301273547e+08,1.2271879740e+08,1.2127611477e+08,1.2099285491e+08,1.2071131346e+08,1.1987683191e+08,1.2015331850e+08,1.1987683191e+08,1.1960199758e+08,1.1987683191e+08,1.1905722195e+08,1.1878724932e+08};

        if (Z < 1) return 1.0/(op_dx*op_dx); // Use point charge if Z<1
        else return xi[Z];
    }
    
public:
    FiniteNucleusPotential(int Z) : xi(expnt(Z)) {}

    double operator()(const coord_3d& r) const {
        const double x=r[0], y=r[1], z=r[2];
        const double R=sqrt(x*x+y*y+z*z);
        const double arg=R*sqrt(xi);
        double q = -Z;
        if (arg<6) q *= std::erf(arg);
        return q/R;
    }
};

/// Given coefficients and exponents make the 3D operator adding delta part to last Gaussian
real_convolution_3d make_operator(World& world,
                                  const std::vector<double>& C,
                                  const std::vector<double>& T,
                                  double Cdelta) {
    Tensor<double> c(C.size()), t(T.size());

    for (unsigned int i=0; i<C.size(); i++) {
        c[i] = C[i];
        t[i] = T[i];
    }

    c[C.size()-1] += Cdelta * pow(t[C.size()-1]/constants::pi,1.5);

    return real_convolution_3d(world, c, t);
}

/// Factory function generating operator for convolution with grad(operator) in 3D

/// Returns a 3-vector containing the convolution operator for the
/// x, y, and z components of grad(operator) where the operator
/// is provided as a sum over Gaussians
static
std::vector<real_convolution_3d_ptr>
make_grad_operator(World& world,
		   std::vector<double> C, // TAKES A COPY SINCE WE MODIFY TI
		   const std::vector<double> T,
		   const double Cdelta,
		   int k=FunctionDefaults<3>::get_k())
{
  //typedef SeparatedConvolution<double,3> real_convolution_3d;
  //typedef std::shared_ptr<real_convolution_3d> real_convolution_3d_ptr;
  const double pi = constants::pi;
  const Tensor<double> width = FunctionDefaults<3>::get_cell_width();
  double hi = width.normf(); // Diagonal width of cell

  const int rank = C.size();

  C[rank-1] += Cdelta * pow(T[rank-1]/constants::pi,1.5);

  std::vector<real_convolution_3d_ptr> gradG(3);

  for (int dir=0; dir<3; dir++) {
    std::vector< ConvolutionND<double,3> > ops(rank);
    for (int mu=0; mu<rank; mu++) {
      // We cache the normalized operator so the factor is the value we must multiply
      // by to recover the coeff we want.
      double c = std::pow(T[mu]/pi,1.5); // Normalization coeff
      ops[mu].setfac(C[mu]/c/width[dir]);

      for (int d=0; d<3; d++) {
	if (d != dir)
	  ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, T[mu]*width[d]*width[d], 0, false));
      }
      ops[mu].setop(dir,GaussianConvolution1DCache<double>::get(k, T[mu]*width[dir]*width[dir], 1, false));
    }
    gradG[dir] = real_convolution_3d_ptr(new SeparatedConvolution<double,3>(world, ops));
  }

  return gradG;
}


/// Makes the relativistic equivalent of the BSH operator
//real_convolution_3d BSHrel(World& world, double epsilon) {
//    std::vector<double> C, T;
//    double Cdelta;
//    bshrel_fit(epsilon, op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_operator(world, C, T, Cdelta);
//}

/// Makes the Tbar operator (T_rel = E0-mc2 = Tbar T_nonrel = -1/2 Tbar del**2)
//real_convolution_3d Tbar(World& world) {
//    std::vector<double> C, T;
//    double Cdelta;
//    tbar_fit(op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_operator(world, C, T, Cdelta);
//}
    
/// Makes the Pbar operator
//real_convolution_3d Pbar(World& world) {
//    std::vector<double> C, T;
//    double Cdelta;
//    pbar_fit(op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_operator(world, C, T, Cdelta);
//}
    
/// Makes the APbar operator
//real_convolution_3d APbar(World& world) {
//    std::vector<double> C, T;
//    double Cdelta;
//    apbar_fit(op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_operator(world, C, T, Cdelta);
//}
    
/// Makes the gradient of the APbar operator
//std::vector<real_convolution_3d_ptr> gradAPbar(World& world) {
//    std::vector<double> C, T;
//    double Cdelta;
//    apbar_fit(op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_grad_operator(world, C, T, Cdelta);
//}
    

/// Makes the A operator (including the 1/sqrt(2) piece)
//real_convolution_3d A(World& world) {
//    std::vector<double> C, T;
//    double Cdelta;
//    a_fit(op_dx, op_thresh, op_quadacc, C, T, Cdelta);
//    return make_operator(world, C, T, Cdelta);
//}

/// Returns the exponent of the wavefunction cusp at the origin for RK
double rk_v(double Z) {
  if (NONREL) return 0.0;
  if (DK1) Z *= 0.5;

    static const double c = 137.0359895;
    double v = -2.0*Z/(constants::pi*c);
    while (true) {
        double vnew = -2.0*std::atan(Z/(c*(v+1.0)))/constants::pi;
        if (std::abs(v-vnew) < 1e-8) break;
        v = vnew;
    }
    return v;
}

class Guess : public FunctionFunctorInterface<double,3> {
    const double Z;
    const double fac;
    const double v;
    std::vector<coord_3d> specialpt;

public:
    Guess(double Z)
        : Z (Z)
        , fac(sqrt(Z*Z*Z/(8.0*constants::pi)))
        , v(rk_v(Z))
        , specialpt(1,coord_3d{0.0,0.0,0.0})
    {}

    double operator()(const coord_3d& r) const {
        const double x=r[0], y=r[1], z=r[2];
        const double R=sqrt(x*x+y*y+z*z+op_dx*op_dx);
        return fac*exp(-Z*R)*pow(R,-v);
        //return fac*exp(-Z*R);
    }

    std::vector<coord_3d> special_points() const override final {return specialpt;}

    Level special_level() const override final {return 5;}
};

// static double V(const coord_3d& r) {
//     const double x=r[0], y=r[1], z=r[2];
//     return -Z/(sqrt(x*x+y*y+z*z+1e-10));
// }

real_function_3d apply_potential(const real_function_3d& Vnuc, const real_function_3d& psi) {
    World& world = Vnuc.world();

    if (NONREL || RK) {
        return Vnuc*psi;
    }
    else if (DK1) {
        // A (V + Pp . V pP ) A psi = A V A psi + AP p . V p AP psi

        real_convolution_3d Aop = A(world);
        real_convolution_3d APbarop = PbarA(world);
        real_function_3d tempfunc(world);
        double tempdouble;
        double fac = pow(2.0*constants::pi,-3.0/2.0);
        
        tempfunc = (fac*fac)*Aop((Vnuc*Aop(psi)).truncate());
        //tempdouble = tempfunc.norm2();
        //tempdouble = (psi*tempfunc).trace();
        //if(world.rank()==0) print("expec AVApsi: ", tempdouble);
        real_function_3d Vpsi = tempfunc;

        tempfunc = fac*(1.0/sqrt(2.0))*Aop((Vnuc*psi).truncate());
        //tempdouble = tempfunc.norm2();
        //tempdouble = (psi*tempfunc).trace();
        //if(world.rank()==0) print("expec AVpsi: ", tempdouble);
        Vpsi += tempfunc;

        tempfunc = fac*(1.0/sqrt(2.0))*Vnuc*Aop(psi);
        //tempdouble = tempfunc.norm2();
        //tempdouble = (psi*tempfunc).trace();
        //if(world.rank()==0) print("expec VApsi: ", tempdouble);
        Vpsi += tempfunc;

        tempfunc = (1.0/2.0)*Vnuc*psi;
        //tempdouble = tempfunc.norm2();
        //tempdouble = (psi*tempfunc).trace();
        //if(world.rank()==0) print("expec Vpsi: ", tempdouble);
        Vpsi += tempfunc;


	    //// This works
         //real_function_3d APbarpsi = APbarop(psi);
         //for (int axis=0; axis<3; axis++) {
	    //   real_derivative_3d D = free_space_derivative<double,3>(world, axis);
         //   tempfunc = (fac*fac)*APbarop((D(Vnuc*D(APbarpsi))).truncate());
         //   //tempdouble = tempfunc.norm2();
         //   //if(world.rank()==0) print("expec PApVpPA[",axis,"]: ", tempdouble);
         //   Vpsi -= tempfunc;
	    //}
	
	// This combines derivative and integral operators and seems both faster and more accurate??
        std::vector<real_convolution_3d_ptr> gradAPbarop = gradPbarA(world);
        for (int axis=0; axis<3; axis++) {
	      real_convolution_3d& op = *(gradAPbarop[axis]);
	      tempfunc = (fac*fac)*op((Vnuc*op(psi)).truncate());
           Vpsi -= tempfunc; //should be minus?
        }

        return Vpsi.truncate();
    }
    else {
        throw "confusion";
    }
}

double compute_energy(World& world, const real_function_3d& psi, const real_function_3d& V, bool doprint=false) {
    real_convolution_3d Tbarop = Tbar(world);

    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);

        if (NONREL) {
            kinetic_energy += 0.5*inner(dpsi,dpsi);
        }
        else {
            kinetic_energy += 0.5*inner(dpsi,apply(Tbarop,dpsi));
        }
    }

    double nuclear_attraction_energy = inner(psi, apply_potential(V,psi));
    double total_energy = kinetic_energy + nuclear_attraction_energy;

    if (world.rank() == 0 && doprint) {
        if (NONREL) print("Non-relativistic");
        else if (RK) print("RK");
        else if (DK1) print("DK1");
        else throw "confused";
                         
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", nuclear_attraction_energy / kinetic_energy);
    }

    return total_energy;
}

double compute_energy_simple(World& world, const real_function_3d& psi, const real_function_3d& Vpsi, bool doprint=false) {
    real_convolution_3d Tbarop = Tbar(world);

    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);

        if (NONREL) {
            kinetic_energy += 0.5*inner(dpsi,dpsi);
        }
        else {
            kinetic_energy += 0.5*inner(dpsi,apply(Tbarop,dpsi));
        }
    }

    double nuclear_attraction_energy = inner(psi, Vpsi);
    double total_energy = kinetic_energy + nuclear_attraction_energy;

    if (world.rank() == 0 && doprint) {
        if (NONREL) print("Non-relativistic");
        else if (RK) print("RK");
        else if (DK1) print("DK1");
        else throw "confused";
                         
        print("            Kinetic energy ", kinetic_energy);
        print(" Nuclear attraction energy ", nuclear_attraction_energy);
        print("              Total energy ", total_energy);
        print("                    Virial ", nuclear_attraction_energy / kinetic_energy);
    }

    return total_energy;
}

real_function_3d apply_laplacian(World& world, const real_function_3d& psi) {
    real_function_3d delsqpsi(world);
    
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        delsqpsi += D(D(psi));
    }

    return delsqpsi;
}

real_convolution_3d make_bsh_operator(World& world, double eps, double& fac) {
    if (NONREL) {
        fac = -2.0;
        return BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
    }
    else {
        fac = -1.0;
        return Ebar(world, eps);
    }
}

real_function_3d iterate(World& world, const real_function_3d& V, const real_function_3d& psi, double& eps) {
    /*

      (Tbar T + V) psi = E psi

      (T + Tbar T - T + V) psi = E psi

      (T - E) psi = - (Tbar T - T + V) psi --> but this suffers from numerical noise

      instead

      psi = - (TbarT - E)^-1 V psi

     */

    real_function_3d Vpsi = apply_potential(V,psi);
    eps = compute_energy_simple(world, psi, Vpsi, false);

    double fac;
    real_convolution_3d op = make_bsh_operator(world, eps, fac);

    Vpsi.scale(fac).truncate();
    
    real_function_3d tmp = apply(op,Vpsi).truncate();
    double norm = tmp.norm2();
    real_function_3d r = psi-tmp;
    // double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    // eps = eps_new;

    return r;
}

void logplot(double Z, const real_function_3d& psi, double lo=1e-7, double hi=L*0.5) {
    char fname[256];

    if (NONREL) sprintf(fname, "nonrel-%3.3d.dat", int(Z));
    else if (RK) sprintf(fname, "rk-%3.3d.dat", int(Z));
    else if (DK1) sprintf(fname, "dk1-%3.3d.dat", int(Z));
    else throw "jdsalkfjalks";

    std::ofstream file;
    file.open (fname);
    file.precision(12);
    psi.reconstruct();
    const double fac = pow(10.0,1.0/16.0);
    double rr = lo;
    while (rr<hi) {
        coord_3d r{rr,0.0,0.0};
        file << rr << " " << psi(r) << std::endl;
        rr *= fac;
    }
    file.close();
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(16);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_truncate_mode(0); //See what happens to lighter elements with mode 1 (hyp more accuracy for lighter elements, but too long for higher elements)
    FunctionDefaults<3>::set_truncate_on_project(true);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

    double Zlist[] = {1.0,2.0,3.0,4.0,6.0,8.0,10.0,12.0,16.0,20.0,30.0,40.0,60.0,80.0};
    //double Zlist[] = {1.0,2.0,4.0,8.0,10.0,16.0,20.0,32.0,40.0,48.0,56.0,60.0,64.0,72.0,76.0,80.0};
    //double Zlist[] = {60.0,64.0,72.0,76.0,80.0};
    //double Zlist[] = {20.0,40.0,56.0,60.0,64.0,72.0,76.0,80.0};
    //double Zlist[] = {100.0};
    const int NumZs = sizeof(Zlist)/sizeof(Z);


    for (unsigned int i=0; i<NumZs; i++) {
        Z = Zlist[i];

	   if(world.rank()==0) print("!!!!!!!!!!!!!!!!!!!!!!!!!");
	   if(world.rank()==0) print("Z =", Z, "   v =", rk_v(Z), "   NONREL =", NONREL, "   RK =", RK, "   DK1 =", DK1, "   FINITENUC =", FINITENUC);

        real_function_3d psi  = real_factory_3d(world).functor(real_functor_3d(new Guess(Z)));
        psi.scale(1.0/psi.norm2());
        psi.truncate();

        real_function_3d Vnuc;


        if (FINITENUC) {
            Vnuc = real_factory_3d(world).functor(real_functor_3d(new FiniteNucleusPotential(Z))).truncate_mode(0);
        }
        else {
            Vnuc = real_factory_3d(world).functor(real_functor_3d(new FiniteNucleusPotential(-1))).truncate_mode(0);
        }
        

	   double eps;// = compute_energy(world, psi, Vnuc);
        //if(world.rank()==0) print("eps beforehand: ", eps);
	
	   NonlinearSolver solver(10);
	   for (int iter=0; iter<20; iter++) {
	      psi.scale(1.0/psi.norm2());
	      //double eps = compute_energy(world, psi, Vnuc);
	      real_function_3d residual = iterate(world, Vnuc, psi, eps);
	     double rnorm = residual.norm2();
	     real_function_3d psi_new = solver.update(psi, residual);

          if (rnorm > 0.1) {
	        psi = 0.5*psi_new + 0.5*psi;
          }
          else {
	        psi = psi_new;
          }

	     if (world.rank() == 0) {
             print(" eps=",eps," err(psi)=",rnorm);
	     }
	     //if (iter>3 && rnorm < pow(10,0.5*log10(thresh))) break;
	     if (iter>3 && rnorm < thresh*50.0) break;
          LoadBalanceDeux<3> lb(world);
          lb.add_tree(Vnuc, lbcost<double,3>(12.0,96.0),true);
          lb.add_tree(psi, lbcost<double,3>(12.0,96.0),true);
          FunctionDefaults<3>::redistribute(world,lb.load_balance(2),false);
	   }
	
        psi.scale(1.0/psi.norm2());
	   compute_energy(world, psi, Vnuc, true);

        logplot(Z,psi);
    }

    world.gop.fence();

    finalize();
    return 0;
}
