
#include "madness.h"

using namespace madness;

// 1D sin(n*x)
class F : public FunctionFunctorInterface<double,1> 
{
   private:
      const int n;

   public:
      F(int n) : n(n) {}

      double operator()(const coord_1d& x) const 
      {
         return std::sin(n * x[0]);
      }    
};

// 1D gaussian, centered at point p
// exp(n * (x-p) * (x-p))
class G : public FunctionFunctorInterface<double,1> 
{
   private:
      const int n;
      const double p;

   public:
      G(int n, double p) : n(n), p(p) {}

      double operator()(const coord_1d& x) const 
      {
         return std::exp(-n * (x[0]-p) * (x[0]-p));
      }    
};


int main(int argc,char** argv)
{
   // MADNESS initialization boiler plate
   initialize(argc, argv);
   World world(SafeMPI::COMM_WORLD);
   startup(world, argc, argv);
   std::cout.precision(17);

   int k = 4;
   double thresh = 1e-6;

   // More MADNESS boiler plate
   FunctionDefaults<1>::set_k(k);
   FunctionDefaults<1>::set_thresh(thresh);
   FunctionDefaults<1>::set_cubic_cell(-6,6);

   // Make sure truncate on project is off
   FunctionDefaults<1>::set_truncate_on_project(false);

   // Our 1 parameter
   int n = 1;

   // Create the functions
   real_function_1d f = real_factory_1d(world).functor(real_functor_1d(new G(n, 0.5)));
   //real_function_1d f = real_factory_1d(world).functor(real_functor_1d(new F(n)));
   f = project(f,  k, thresh); 
   real_function_1d ftrunc = copy(f).truncate();
   real_function_1d eps = f - ftrunc;
   real_function_1d feps = f + eps;

   // Create an agbv derivative
   real_derivative_1d Dx(world, 0); 

   // Apply the derivative
   real_function_1d deps = apply(Dx, eps);
   real_function_1d dftrunc = apply(Dx, ftrunc);
   real_function_1d dfeps = apply(Dx, feps);
   real_function_1d df = apply(Dx, f);

   // Second derivative
   real_function_1d d2f = apply(Dx, df);
   double integral = -inner(f,d2f);
   print("integral =", integral);
   print("sqrt(pi/2) =", sqrt(constants::pi/2.0));

   // Test orthogonality between functions
   double eps_ftrunc = inner(eps, ftrunc);
   double deps_ftrunc = inner(deps, ftrunc); 
   double deps_dftrunc = inner(deps, dftrunc); 
   double dfeps_dfeps = inner(dfeps, dfeps);
   double dftrunc_dftrunc = inner(dftrunc, dftrunc);
   double dftrunc_deps = inner(dftrunc, deps);
   double deps_deps = inner(deps, deps);

   // Output 
   print("\n< eps | ftrunc > =", eps_ftrunc);
   print("\n< deps | ftrunc > =", deps_ftrunc);
   print("\n< deps |dftrunc > =", deps_dftrunc);
   print("");
   print("\n< dfeps | dfeps > =", dfeps_dfeps);
   print("\n< dftrunc | dftrunc > =", dftrunc_dftrunc);
   print("\n< dftrunc | deps > =", dftrunc_deps);
   print("\n< deps | dftrunc > =", deps_dftrunc);
   print("\n< deps | deps > =", deps_deps);

   // Line plots
   int npts = 1001;
   const Vector<double,1> lo(0.0);
   const Vector<double,1> hi( 1.0);
   plot_line("deriv_test_plots.txt", npts, lo, hi, f, ftrunc, eps, deps);

   // Printing trees in reconstructed form
   f.reconstruct();
   ftrunc.reconstruct();
   eps.reconstruct();
   deps.reconstruct();
 
   print("\n\nf (reconstructed) tree:\n");
   f.print_tree();
   print("\n\nftrunc (reconstructed) tree:\n");
   ftrunc.print_tree();
   print("\n\neps (reconstructed) tree:\n");
   eps.print_tree();
   print("\n\ndeps (reconstructed) tree:\n");
   deps.print_tree();

   // Printing trees in compressed form
   f.compress();
   ftrunc.compress();
   eps.compress();
   deps.compress();
 
   print("\n\nf (compressed) tree:\n");
   f.print_tree();
   print("\n\nftrunc (compressed) tree:\n");
   ftrunc.print_tree();
   print("\n\neps (compressed) tree:\n");
   eps.print_tree();
   print("\n\ndeps (compressed) tree:\n");
   deps.print_tree();

   // MADNESS finalization boiler plate
   world.gop.fence();
   finalize();
   return 0;
}
