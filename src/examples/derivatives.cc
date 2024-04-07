
/* The following goes through each derivative model and runs the macro tests using that model.
 * To load the new derivative models, you create the derivative operator normally, then use the 
 * "set_deriv_type#()" function. See lines 413 - 437 for examples.
 */

#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <string>

using namespace madness;

// Constants for use in tests
const double L = 16.0;   // box size [-L,L]^3
const double a_gaussian = 1.0; 
const std::size_t bufsize=256;

// Analytic gaussian function
double gaussian(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return fac*exp(-a*(x*x+y*y+z*z));
}
// Analytic first derivative of gaussian function
double gaussian1(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return -2.0*a*x*fac*exp(-a*(x*x+y*y+z*z));
}
// Analytic second derivative of gaussian function
double gaussian2(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return (4.0*a*a*x*x - 2.0*a)*fac*exp(-a*(x*x+y*y+z*z));
}
// Analytic third derivative of gaussian function
double gaussian3(const coord_3d& r) {
    const double a = a_gaussian;
    const double fac = std::pow(2*a/constants::pi,0.75);
    const double x=r[0], y=r[1], z=r[2];
    return (12.0*a*a*x - 8.0*a*a*a*x*x*x)*fac*exp(-a*(x*x+y*y+z*z));
}
// Expectation value of gaussian
double gaussian_energy() {
    double a = a_gaussian;
    return a;
}

// Define functions
const std::string funcname = "gaussian";
double (*f )(const coord_3d&) = gaussian;
double (*f1)(const coord_3d&) = gaussian1;
double (*f2)(const coord_3d&) = gaussian2;
double (*f3)(const coord_3d&) = gaussian3;
double (*energy)() = gaussian_energy;

class F : public FunctionFunctorInterface<double,3> {
private:
    double(*f)(const coord_3d&);
    std::vector<coord_3d> pts;
public:
    F(double(*f)(const coord_3d&), const coord_3d& pt) : f(f), pts(1,pt) {}

    double operator()(const coord_3d& r) const {return f(r);}

    std::vector<coord_3d> special_points() const {return pts;}
};

char* p(char* buf, const char* name, int k, int initial_level, double thresh, int order) {
    snprintf(buf, bufsize, "%s-%02d-%02d-%.1e-%d.dat", name, k, initial_level, thresh, order);
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

void plotter1(const std::string& fname, const std::vector<coord_3d>& pts, real_function_3d& g0, real_function_3d& g1, real_function_3d& g2, real_function_3d& g3) {
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

void plotter2(const std::string& fname, const std::vector<coord_3d>& pts, real_function_3d& g0, real_function_3d& g2) {
  Tensor<double> g0numer = tabulate(g0 , pts);
  Tensor<double> g2numer = tabulate(g2, pts);
  Tensor<double> g0exact = tabulate(f , pts);
  Tensor<double> g2exact = tabulate(f2, pts);

  FILE* file = fopen(fname.c_str(), "w");
  for (size_t i=0; i<pts.size(); i++) {
    fprintf(file,"%21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e\n", 
	    pts[i][0], pts[i][1], pts[i][2], g0numer[i], g0exact[i], g2numer[i], g2exact[i]);
  }

  fclose(file);
}
 
void plotter3(const std::string& fname, const std::vector<coord_3d>& pts, real_function_3d& g0, real_function_3d& g3) {
  Tensor<double> g0numer = tabulate(g0 , pts);
  Tensor<double> g3numer = tabulate(g3, pts);
  Tensor<double> g0exact = tabulate(f , pts);
  Tensor<double> g3exact = tabulate(f3, pts);

  FILE* file = fopen(fname.c_str(), "w");
  for (size_t i=0; i<pts.size(); i++) {
    fprintf(file,"%21.12e %21.12e %21.12e %21.12e %21.12e %21.12e %21.12e \n", 
	    pts[i][0], pts[i][1], pts[i][2], g0numer[i], g0exact[i], g3numer[i], g3exact[i]);
}
 
  fclose(file);
}

void test(World& world, int k, int initial_level, double thresh, int truncmode, bool refine) {

    char buf[bufsize];

    print("\n\n\n");
    print("testing function:", funcname, "k:", k, "initial level:", initial_level, "thresh:", thresh, "truncmode:", truncmode, "refine:", refine);
    print("\n");

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_refine(refine);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_initial_level(initial_level);
    FunctionDefaults<3>::set_truncate_mode(truncmode); 
    FunctionDefaults<3>::set_truncate_on_project(true); 

    const double lo=-L, hi=L;
    const size_t npts = 10001;
    const std::vector<coord_3d> pts = make_pts(npts, lo+1e-13, hi-1e-13); // Sampling points

    std::string model;

    FunctionDefaults<3>::set_cubic_cell(lo,hi);
    
    Derivative<double,3> D = free_space_derivative<double,3>(world, 0);
    if (funcname == "sin"){
        FunctionDefaults<3>::set_bc(BC_PERIODIC);
    }
    
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

// Prints relevant values for a first derivative    
#define DOIT1 \
 	real_function_3d g1 = D(g);  \
 	real_function_3d g2 = D(g1); \
 	real_function_3d g3 = D(g2); \
 	PRINT("relerrg1", g1.err(f1)/g1norm); \
 	PRINT("relerrg2", g2.err(f2)/g2norm); \
 	PRINT("relerrg3", g3.err(f3)/g3norm); \
	PRINT("e-<g1|g1>", e-g1.inner(g1));   \
 	PRINT("e+<g|g2>", e+g.inner(g2));     \
        plotter1(p(buf,model.c_str(),k,initial_level,thresh,1), pts, g, g1, g2, g3);

// Prints relevant values for a second derivative
#define DOIT2 \
	real_function_3d g2 = D(g); \
 	PRINT("relerrg2", g2.err(f2)/g2norm); \
 	PRINT("e+<g|g2>", e+g.inner(g2));     \
        plotter2(p(buf,model.c_str(),k,initial_level,thresh,2), pts, g, g2);

// Prints relevant values for a third derivative
#define DOIT3 \
 	real_function_3d g3 = D(g); \
 	PRINT("relerrg3", g3.err(f3)/g3norm); \
        plotter3(p(buf,model.c_str(),k,initial_level,thresh,3),  pts, g, g3);

    {
        model = "abgv";
	DOIT1;
    }

    {
        model = "ble-1";
        D.set_ble1();
	DOIT1;
    }

    {
        model = "ble-2";
        D.set_ble2();
	DOIT2;
    }

    {
        model = "bspline-1";
        D.set_bspline1();
	DOIT1;
    }

    {
        model = "bspline-2";
        D.set_bspline2();
	DOIT2;
    }

    {
        model = "bspline-3";
        D.set_bspline3();
        DOIT3;
    }

}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(10);

    if(world.rank() == 0) print("\n\nBLE derivatives are available for k = 1 to k = 15.\nBspline derivatives are available for k = 1 to k = 18.\n");

    double thresh = 1e-4;
    int initial_level = 4; 
    for (int k=7; k<15; k+=2) {
       test(world, k, initial_level, thresh, 1, true);
    }

    world.gop.fence();
    finalize();
    return 0;
}
