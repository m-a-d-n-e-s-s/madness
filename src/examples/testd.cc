#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <string>
#include <malloc.h>

// const char* BSP_DATA_PATH = "/gpfs/home/rharrison/madness/src/madness/mra/deriv-bsp.k=m+1.n=m+1";
// const char* PH_DATA_PATH  = "/gpfs/home/rharrison/madness/src/madness/mra/ph-spline-deriv.txt";
// const char* PH2_DATA_PATH = "/gpfs/home/rharrison/madness/src/madness/mra/ph-spline-deriv-2.txt";
// const char* PROL_DATA_PATH= "/gpfs/home/rharrison/madness/src/madness/mra/prolates-joel";

const char* BSP_DATA_PATH = "/home/rjh/Devel/mad-der/src/madness/mra/deriv-bsp.k=m+1.n=m+1";
const char* PH_DATA_PATH  = "/home/rjh/Devel/mad-der/src/madness/mra/ph-spline-deriv.txt";
const char* PH2_DATA_PATH = "/home/rjh/Devel/mad-der/src/madness/mra/ph-spline-deriv-2.txt";
const char* PROL_DATA_PATH= "/home/rjh/Devel/mad-der/src/madness/mra/prolates-joel";

using namespace madness;

const double L = 32.0;   // box size [-L,L]^3
const double a_gaussian = 1.0; 

double exponential(const coord_3d& r) {
    const double a = 0.5;
    const double small = 1e-5;
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x + y*y + z*z + small*small);
    return exp(-a*R);
}

double exponential1(const coord_3d& r) {
    const double a = 0.5;
    const double small = 1e-5;
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x + y*y + z*z + small*small);
    return -exp(-a*R)*a*x/R;
}

double exponential2(const coord_3d& r) {
    const double a = 0.5;
    const double small = 1e-5;
    const double x=r[0], y=r[1], z=r[2];
    const double R = std::sqrt(x*x + y*y + z*z + small*small);
    return exp(-a*R)*(a*a*x*x/(R*R) - a/R + a*x*x/(R*R*R));
}

double exponential3(const coord_3d& r) {
    const double a = 0.5;
    const double a2 = a*a;
    const double a3 = a2*a;
    const double small = 1e-5;
    const double x=r[0], y=r[1], z=r[2];
    const double x2 = x*x;
    const double x3 = x*x2;
    const double R = std::sqrt(x*x + y*y + z*z + small*small);
    const double R2 = R*R;
    const double R3 = R*R2;
    const double R4 = R*R3;
    const double R5 = R*R4;
    return exp(-a*R)*(-3*a*x3/R5 + 3*a*x/R3 - 3*a2*x3/R4 + 3*a2*x/R2 - a3*x3/R3);
}

double exponential_energy() {
    double a = 0.5;
    return constants::pi/a/3.0;
}


double flat(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2 = x*x + y*y + z*z;
    const double r4 = r2*r2;

    return erfc(r4);
}


double flat1(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2 = x*x + y*y + z*z;
    const double r4 = r2*r2;
    const double r8 = r4*r4;

    return -8*exp(-r8)* r2 * x/sqrt(constants::pi);
}

double flat2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2 = x*x + y*y + z*z;
    const double r4 = r2*r2;
    const double r8 = r4*r4;

    const double s2 = y*y + z*z;
    const double s4 = s2*s2;
    const double s6 = s4*s2;
    const double s8 = s6*s2;

    const double x2 = x*x;
    const double x4 = x2*x2;
    const double x6 = x4*x2;
    const double x8 = x6*x2;
    const double x10 = x8*x2;

    return 8*exp(-r8)*(8*s8*x2+32*s6*x4+48*s4*x6+32*s2*x8+8*x10-s2-3*x2)/sqrt(constants::pi);
}

double flat3(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2 = x*x + y*y + z*z;
    const double r4 = r2*r2;
    const double r8 = r4*r4;

    const double s2 = y*y + z*z;
    const double s4 = s2*s2;
    const double s6 = s4*s2;
    const double s8 = s6*s2;
    const double s10 = s8*s2;
    const double s12 = s10*s2;
    const double s14 = s12*s2;

    const double x2 = x*x;
    const double x4 = x2*x2;
    const double x6 = x4*x2;
    const double x8 = x6*x2;
    const double x10 = x8*x2;
    const double x12 = x10*x2;
    const double x14 = x12*x2;
    const double x16 = x14*x2;

    return -16*exp(-r8)*x*(32*s14*x2+224*s12*x4+672*s10*x6+1120*s8*x8+1120*s6*x10+672*s4*x12+224*s2*x14+32*x16-12*s8-88*s6*x2-192*s4*x4-168*s2*x6-52*x8+3)/sqrt(constants::pi);
}

double flat_energy() {
    return 4.60576993824350;
}


double gaussian(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return fac*exp(-a*(x*x+y*y+z*z));
}

double gaussian1(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return -2.0*a*x*fac*exp(-a*(x*x+y*y+z*z));
}

double gaussian2(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return (4.0*a*a*x*x - 2.0*a)*fac*exp(-a*(x*x+y*y+z*z));
}

double gaussian3(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return (12.0*a*a*x - 8.0*a*a*a*x*x*x)*fac*exp(-a*(x*x+y*y+z*z));
}

double gaussian_energy() {
    double a = a_gaussian;
    return a;
}

// const std::string funcname = "flat";
// double (*f )(const coord_3d&) = flat;
// double (*f1)(const coord_3d&) = flat1;
// double (*f2)(const coord_3d&) = flat2;
// double (*f3)(const coord_3d&) = flat3;
// double (*energy)() = flat_energy;

// const std::string funcname = "gaussian";
// double (*f )(const coord_3d&) = gaussian;
// double (*f1)(const coord_3d&) = gaussian1;
// double (*f2)(const coord_3d&) = gaussian2;
// double (*f3)(const coord_3d&) = gaussian3;
// double (*energy)() = gaussian_energy;

const std::string funcname = "exponential";
double (*f )(const coord_3d&) = exponential;
double (*f1)(const coord_3d&) = exponential1;
double (*f2)(const coord_3d&) = exponential2;
double (*f3)(const coord_3d&) = exponential3;
double (*energy)() = exponential_energy;

class F : public FunctionFunctorInterface<double,3> {
private:
    double(*f)(const coord_3d&);
    std::vector<coord_3d> pts;
public:
    F(double(*f)(const coord_3d&), const coord_3d& pt) : f(f), pts(1,pt) {}

    double operator()(const coord_3d& r) const {return f(r);}

    std::vector<coord_3d> special_points() const {return pts;}
};

char* p(char* buf, const char* name, int k, int initial_level, double thresh) {
    sprintf(buf, "%s-%02d-%02d-%.1e.dat", name, k, initial_level, thresh);
    return buf;
}

template <typename funcT>
Tensor<double> tabulate(funcT& f, const std::vector<coord_3d>& pts) {
  Tensor<double> v(pts.size());
  int i=0;
  for (auto pt : pts) {
    v[i++] = f(pt);
  }
  return v;
}

std::vector<coord_3d> make_pts(size_t n, double lo, double hi) {
  std::vector<coord_3d> pts;
  double h = (hi-lo)/(n-1);
  for (size_t i=0; i<n; i++) {
    double x = lo + h*i;
    pts.push_back(coord_3d({x,0.0,0.0}));
  }
  return pts;
}

void plotter(const std::string& fname, const std::vector<coord_3d>& pts, real_function_3d& g0, real_function_3d& g1, real_function_3d& g2, real_function_3d& g3) {
  Tensor<double> g0numer = tabulate(g0 , pts);
  Tensor<double> g1numer = tabulate(g1, pts);
  Tensor<double> g2numer = tabulate(g2, pts);
  Tensor<double> g3numer = tabulate(g3, pts);
  Tensor<double> g0exact = tabulate(f , pts);
  Tensor<double> g1exact = tabulate(f1, pts);
  Tensor<double> g2exact = tabulate(f2, pts);
  Tensor<double> g3exact = tabulate(f3, pts);

  FILE* file = fopen(fname.c_str(), "w");
  for (size_t i=0; i<pts.size(); i++) {
    fprintf(file,"%21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e \n", 
	    pts[i][0], pts[i][1], pts[i][2], g0numer[i], g0exact[i], g1numer[i], g1exact[i], g2numer[i], g2exact[i], g3numer[i], g3exact[i]);
  }

  fclose(file);
}

void test(World& world, int k, int initial_level, double thresh, int truncmode, bool refine) {

    char buf[256];

    print("\n\n\n");
    print("testing function:", funcname, "k:", k, "initial level:", initial_level, "thresh:", thresh, "truncmode:", truncmode, "refine:", refine);
    print("\n");

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_refine(refine);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_initial_level(initial_level);
    FunctionDefaults<3>::set_truncate_mode(truncmode); 
    FunctionDefaults<3>::set_truncate_on_project(true); // <<<<<<<<<

    const double lo=-L, hi=L;
    const size_t npts = 10001;
    const std::vector<coord_3d> pts = make_pts(npts, lo+1e-13, hi-1e-13); // Sampling points

    std::string model;

    FunctionDefaults<3>::set_cubic_cell(lo,hi);
    
    Derivative<double,3> D = free_space_derivative<double,3>(world, 0);
    
    double g0norm, g1norm, g2norm, g3norm;
    
    {
      real_functor_3d gfunctor(new F(f1,{0,0,0}));
      real_function_3d g1  = real_factory_3d(world).functor(gfunctor).thresh(1e-4);
      g1norm = g1.norm2();
    }
    {
      real_functor_3d gfunctor(new F(f2,{0,0,0}));
      real_function_3d g2  = real_factory_3d(world).functor(gfunctor).thresh(1e-4);
      g2norm = g2.norm2();
    }
    {
      real_functor_3d gfunctor(new F(f3,{0,0,0}));
      real_function_3d g3 = real_factory_3d(world).functor(gfunctor).thresh(1e-4);
      g3norm = g3.norm2();
    }

    real_functor_3d gfunctor(new F(f,{0,0,0}));
    real_function_3d g  = real_factory_3d(world).functor(gfunctor);
    g0norm = g.norm2();
    
    const double e = energy();
    
#define PRINT(what,value) print("@", initial_level, ",", k, ",", thresh, ",", model, ",", what, ",", value)
    
    model = "project";
    
    print("@ headings, initial_level, k, model, what, value");
    
    PRINT("g0norm", g0norm);
    PRINT("g1norm", g1norm);
    PRINT("g2norm", g2norm);
    PRINT("energyexact", e);
    
#define DOIT \
 	real_function_3d g1 = D(g);  \
 	real_function_3d g2 = D(g1); \
 	real_function_3d g3 = D(g2); \
 	PRINT("relerrg1", g1.err(f1)/g1norm); \
 	PRINT("relerrg2", g2.err(f2)/g2norm); \
 	PRINT("relerrg3", g3.err(f3)/g3norm); \
	PRINT("e-<g1|g1>", e-g1.inner(g1));				\
 	PRINT("e+<g|g2>", e+g.inner(g2));				\
        plotter(p(buf,model.c_str(),k,initial_level,thresh),  pts, g, g1, g2, g3);


    {
        model = "abgv";
	DOIT;
    }

    {
        model = "bsp";
        D.read_from_file(BSP_DATA_PATH);
	DOIT;
    }

    {
        model = "ph";
        D.read_from_file(PH_DATA_PATH);
	DOIT;
    }

    {
        model = "prol";
        D.read_from_file(PROL_DATA_PATH);
	DOIT;
    }

    // {
    //     model = "ph2";
    //     D.read_from_file(PH2_DATA_PATH);
    //     D.set_is_second(true);
    //     real_function_3d g2 = D(g);
    // 	PRINT("relerrg0norm", g.err(f)/g0norm); 
    // 	PRINT("relerrg2norm", g2.err(f2)/g2norm);
    // 	PRINT("e+<g|g2>", e+g.inner(g2));

    //     real_functor_3d g2functor(new F(f2,{0,0,0}));
    //     real_function_3d g2exact  = real_factory_3d(world).functor(g2functor);
    // 	PRINT("e+<g|g2exact>", e+g.inner(g2exact));
    //     plot_line(p(buf,"ph2d",k,initial_level,thresh),  10001, {lo,0,0}, {hi,0,0}, g, g2exact, g2);
    // }
}

int main(int argc, char** argv) {
    //if (!mallopt(M_MXFAST, 0)) return -1;
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    double threshes[] = {1e-4, 1e-6, 1e-8, 1e-10};

    for (double thresh : threshes) {
      for (int initial_level = 4; initial_level<=4; initial_level++) { // not used if adaptive refining
	for (int k=7; k<19; k+=2) {
	  test(world, k, initial_level, thresh, 1, true);
	}
      }
    }
            
    world.gop.fence();
    finalize();
    return 0;
}
