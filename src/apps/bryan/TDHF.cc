/*
 *
 *   Written by: bsundahl
 *   Date: A long time ago...
 *
 */ 

#include "TDHF.h"
#include "Plot_VTK.h"
#include "TDHF_Basic_Operators.h"
#include "../chem/projector.h"             // For easy calculation of (1 - \hat{\rho}^0)
#include "../chem/potentialmanager.h"

using namespace madness;

// KAIN allocator for vectorfunctions
struct TDHF_allocator
{
   // Member variables
   World& world;
   const int num_vir;
   const int num_occ;

   // Constructor
   TDHF_allocator(World& world, const int num_vir, const int num_occ) : world(world), num_vir(num_vir), num_occ(num_occ) {}

   // Overloading () operator
   std::vector<std::vector<real_function_3d>> operator()()
   {
      std::vector<std::vector<real_function_3d>> f;

      for(int i = 0; i < num_vir; i++) f.push_back(zero_functions<double,3>(world,num_occ));

      return f;
   }

   // Copy constructor
   TDHF_allocator operator=(const TDHF_allocator &other)
   {
      TDHF_allocator tmp(world,other.num_occ,other.num_vir);
      return tmp;
   }
};

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

// Pulled from SCF.cc, starts a timer
void TDHF::start_timer(World& world)
{
   world.gop.fence();
   ttt.push_back(wall_time());
   sss.push_back(cpu_time());
}

// Needed for timers
double TDHF::pop(std::vector<double>& v)
{
   double x = v.back();
   v.pop_back();
   return x;
}

// Stops a timer
Tensor<double> TDHF::end_timer(World& world)
{
   Tensor<double> times(2);
   times[0] = wall_time() - pop(ttt);
   times[1] = cpu_time() - pop(sss);
   return times;
}

// Collective constructor
TDHF::TDHF(World & world,
         const char* filename) : TDHF(world, (world.rank() == 0 ? std::make_shared<std::ifstream>(filename) : nullptr))
{}

// Constructor that actually does stuff
TDHF::TDHF(World & world,
         std::shared_ptr<std::istream> input) 
{
   // Start the timer
   start_timer(world);

   // Try and open input file
   if(world.rank() == 0)
   {
      if (input->fail()) MADNESS_EXCEPTION("Response failed to open input stream", 0);
   
      // Welcome user (future ASCII art of Robert goes here) 
      print("\n   Preparing to solve the TDHF equations.\n"); 

      // Read input files
      Rparams.read(*input);

      // Print out what was read in
      Rparams.print_params();
   }

   // Broadcast to all other nodes
   world.gop.broadcast_serializable(Rparams, 0);

   // Read in archive
   Gparams.read(world, Rparams.archive);
   if(world.rank() == 0)
   {
      Gparams.print_params();
      print_molecule(world);
   }   

   // Set some function defaults   
   FunctionDefaults<3>::set_cubic_cell(-Gparams.L, Gparams.L);
   FunctionDefaults<3>::set_k(Gparams.order);
   FunctionDefaults<3>::set_thresh(Rparams.thresh);
   FunctionDefaults<3>::set_truncate_mode(1);   
}

// (Each state's norm should be 1, not the 
// individual functions norms)
void TDHF::normalize(World & world,
                    std::vector<std::vector<real_function_3d>> & f)
{
   // Run over rows
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Get the normalization constant
      // (Sum included inside inner) 
      double norm = inner(f[i], f[i]);
      norm = sqrt(norm);

      // And scale
      scale(world, f[i], 1.0/norm);
   }
};

// Prints norms of the given vector of vector of functions
void TDHF::print_norms(World & world,
                      std::vector<std::vector<real_function_3d>> f)
{
   // Container
   Tensor<double> norms(f.size(),f[0].size());

   // Calc the norms
   for(unsigned int i = 0; i < f.size(); i++)
   {
      for(unsigned int j = 0; j < f[0].size(); j++)
      {
         norms(i,j) = f[i][j].norm2();
      }
   }

   // Print em in a smart way
   if(world.rank() == 0) print(norms);
}

// Small function to print geometry of a molecule nicely
void TDHF::print_molecule(World &world)
{
   if(world.rank() == 0)
   {
      // Precision is set to 10 coming in, drop it to 5
      std::cout.precision(5);
      std::cout << std::fixed;

      // First get atoms
      const std::vector<Atom> atoms = Gparams.molecule.get_atoms();
      int num_atoms = atoms.size();

      // Now print
      print("\n   Geometry Information");
      print("   --------------------\n");
      print("   Units: a.u.\n");
      print(" Atom            x                 y                 z");
      print("----------------------------------------------------------------");
      for(int j = 0; j < num_atoms; j++)
      {
           Vector<double,3> coords = atoms[j].get_coords();
           std::cout << std::setw(3) << atomic_number_to_symbol(atoms[j].get_atomic_number());
           std::cout << std::setw(18) << std::right << coords[0] << std::setw(18) << coords[1] << std::setw(18) << coords[2] << endl;
      }
      print("");

      // Reset precision
      std::cout.precision(10);
      std::cout << std::scientific;
   }
}

// Returns a vector of vectors filled with zero functions 
// with proper sizes (an "m x n" matrix of zero functions)
std::vector<std::vector<real_function_3d>> TDHF::response_zero_functions(World & world,
                                                               int m,
                                                               int n)
{
   // Functions to be returned
   std::vector<std::vector<real_function_3d>> results(m);

   // Create vectors of zero functions
   for(int i = 0; i < m; i++) results[i] = zero_functions<double,3>(world, n);

   // Done
   return results;
}

// Radial function
static double radial(const coord_3d& r)
{
   return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

// Returns a list of symmetry related functions for the correct
// pointgroup of the provided molecule
std::vector<real_function_3d> TDHF::symmetry(World & world)
{
   // Container to return
   std::vector<real_function_3d> result;

   // Create the basic x, y, z
   real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{1,0,0})));
   real_function_3d y = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,1,0})));
   real_function_3d z = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,0,1})));
   real_function_3d r = real_factory_3d(world).f(radial);

   // Add in s function 
   result.push_back(r);

   // Add in p functions
   result.push_back(x);
   result.push_back(y);
   result.push_back(z);

   // Add in d functions
   result.push_back(x*y);
   result.push_back(x*z);
   result.push_back(y*z);
   result.push_back(x*x-y*y);
   result.push_back(z*z);

   // Done
   return result;
}

// Returns initial guess functions
// Probably needs a lot of work on how to do this intelligently
// Currently: Producing each combination of symmetry function on the orbitals
std::vector<std::vector<real_function_3d>> TDHF::create_trial_functions(World & world,
                                                                       int k,
                                                                       std::vector<real_function_3d> & orbitals,
                                                                       int print_level)
{
   // Get size
   int n = orbitals.size();

   // Create a vector of correct symmetry relate polynomials
   // Only going through the d symmetry functions
   std::vector<real_function_3d> symm = symmetry(world);

   // Determine how many functions will be created
   int size = (int)(n * symm.size()) >= k ? n * symm.size() : (k/n + 1) * n;

   // Container to return 
   std::vector<std::vector<real_function_3d>> trials = response_zero_functions(world, size, n);

   // Counter for number of trials created
   int count = 0;

   // Run over symmetry functions
   for(unsigned int i = 0; i < symm.size(); i++)
   {
      // Run over each occupied orbital
      for(int p = 0; p < n; p++)
      {
         trials[count][p] = symm[i] * orbitals[p];
         count++;
      }
   }

   // Make sure we have at least k functions by adding in powers of the symmetry 
   // functions times the orbitals
   int power = 1;
   while(count < size )
   {
      // Run over symmetry functions
      for(unsigned int i = 0; i < symm.size(); i++)
      {
         // Initial symmetry function
         real_function_3d x = symm[i];

         // Get the symmetry function to the right power
         for(int j = 0; j < power; j++) x = x * symm[i];

         // Run over each occupied orbital
         for(int p = 0; p < n; p++)
         {
            trials[count][p] = x * symm[i] * orbitals[p];
            count++;
         }
      }

      // Increase power of r
      power++;
   }

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) print("   Norms of guess functions:");
      print_norms(world, trials);
   }

   // Truncate 
   truncate(world, trials);

   // Done
   return trials;
}

// Returns the derivative of the coulomb operator, applied to ground state orbitals
std::vector<std::vector<real_function_3d>> TDHF::create_coulomb_derivative(World &world,
                                                                          std::vector<std::vector<real_function_3d>> & f,
                                                                          std::vector<real_function_3d> & orbitals,
                                                                          double small,
                                                                          double thresh)
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Zero function, to be returned
   std::vector<std::vector<real_function_3d>> deriv_j = response_zero_functions(world, m, n);

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Temperary storage
   real_function_3d rho = real_function_3d(world);

   // Need to run over virtual orbitals
   for(int k = 0; k < m; k++)
   {
      // Get transitionn density 
      rho = dot(world, f[k], orbitals);

      // Apply coulomb operator
      rho = apply(op, rho);

      // Need to run over all occupied orbitals
      for(int p = 0; p < n; p++)
      {
         // Multiply by ground state orbital p
         // and save the result
         deriv_j[k][p] = rho * orbitals[p];
      }
   }

   // Done
   return deriv_j;
}


// Does what it sounds like it does
std::vector<std::vector<real_function_3d>> TDHF::create_exchange_derivative(World &world,
                                                                            std::vector<std::vector<real_function_3d>> & f,
                                                                            std::vector<real_function_3d> & orbitals,
                                                                            double small,
                                                                            double thresh) 
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Zero function, to be returned
   std::vector<std::vector<real_function_3d>> deriv_k = response_zero_functions(world, m, n);

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);
  
   // Potential is not stored by default
   if(Rparams.store_potential)
   {
      // Need to run over occupied orbitals
      for(int p = 0; p < n; p++)
      {
         // Need to run over all virtual orbitals originating from orbital p
         for(int k = 0; k < m; k++)
         {
            // Need to sum over occupied orbitals
            for(int i = 0; i < n; i++)
            {
               // Multiply precalculated \int rho/r by response function (k,i)
               deriv_k[k][p] += stored_potential[i][p] * f[k][i];
            }
         }
      }
   }
   else // But the storage can be turned off...
   {
      // Need to run over occupied orbitals
      for(int p = 0; p < n; p++)
      {
         // Need to run over all virtual orbitals originating from orbital p
         for(int k = 0; k < m; k++)
         {
            // Need to sum over occupied orbitals
            for(int i = 0; i < n; i++)
            {
               // Get density (ground state orbitals)
               real_function_3d rho = orbitals[i] * orbitals[p];

               // Apply coulomb operator
               rho = apply(op, rho);

               // Multiply by response function (k,i)
               // and add to total
               deriv_k[k][p] += rho * f[k][i];
            }
         }
      }
   }
   
   // Done
   return deriv_k;
}

// Does what it sounds like it does
std::vector<std::vector<real_function_3d>> TDHF::create_exchange_derivative_tdhf(World &world,
                                                                                 std::vector<std::vector<real_function_3d>> & x,
                                                                                 std::vector<std::vector<real_function_3d>> & y,
                                                                                 std::vector<real_function_3d> & orbitals,
                                                                                 double small,
                                                                                 double thresh) 
{
   // Get sizes
   int m = x.size();
   int n = x[0].size();

   // Zero function, to be returned
   std::vector<std::vector<real_function_3d>> deriv_k = response_zero_functions(world, m, n);

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Two pieces: x states and y states
   // Calc. x states first
   // Need to run over occupied orbitals
   for(int p = 0; p < n; p++)
   {
      // Need to run over all virtual orbitals originating from orbital p
      for(int k = 0; k < m; k++)
      {
         // Need to sum over occupied orbitals
         for(int i = 0; i < n; i++)
         {
            // Get density (ground state orbitals)
            real_function_3d rho = orbitals[i] * orbitals[p];

            // Apply coulomb operator
            rho = apply(op, rho);

            // Multiply by response function (k,i)
            // and add to total
            deriv_k[k][p] += rho * x[k][i];
         }
      }
   }
   // Now calc. the y states 
   // Need to run over occupied orbitals
   for(int p = 0; p < n; p++)
   {
      // Need to run over all virtual orbitals originating from orbital p
      for(int k = 0; k < m; k++)
      {
         // Need to sum over occupied orbitals
         for(int i = 0; i < n; i++)
         {
            // Get density (ground state orbitals)
            real_function_3d rho = y[k][i] * orbitals[i]; // Changed index p to i 

            // Apply coulomb operator
            rho = apply(op, rho);

            // Multiply by response function (k,i)
            // and add to total
            deriv_k[k][p] += rho * orbitals[p]; // Changed index i to p
         }
      }
   }
   
   // Done
   return deriv_k;
}

// Computes gamma(r) given the ground state orbitals and response functions
// Only for TDA
std::vector<std::vector<real_function_3d>> TDHF::create_gamma(World &world,
                                                              std::vector<std::vector<real_function_3d>> & f,
                                                              std::vector<real_function_3d> & orbitals,
                                                              double small,
                                                              double thresh,
                                                              int print_level)
{
   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Creating Gamma\n");
   }

   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // The gamma function to be returned, intialized to zero
   std::vector<std::vector<real_function_3d>> gamma = response_zero_functions(world, m, n);

   // Gamma will have 2 terms for HF: dJ/drho[rho] and dK/drho[rho]
   // There is a different Gamma for each orbital-->virtual transition
   // Calculate both here
   std::vector<std::vector<real_function_3d>> deriv_J = create_coulomb_derivative(world, f, orbitals, small, thresh);
   std::vector<std::vector<real_function_3d>> deriv_K = create_exchange_derivative(world, f, orbitals, small, thresh);
   //std::vector<std::vector<real_function_3d>> deriv_XC = create_xc_derivative(world, f, orbitals, small, thresh);

   // Debugging output
   if (print_level >= 2)
   {
      if(world.rank() == 0) printf("   Coulomb Deriv matrix:\n");
      Tensor<double> temp =expectation(world, f, deriv_J);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Exchange Deriv matrix:\n");
      temp = expectation(world, f, deriv_K);
      if(world.rank() == 0) print(temp);
   }

   // Spin integration gives coefficients 
   // This is the spin restricted, singlet excitation coefficients
   gamma = scale(deriv_J, 2.0) - deriv_K;

   // Project out groundstate 
   QProjector<double, 3> projector(world, Gparams.orbitals);           // Projector to project out ground state
   for(int i = 0; i<m; i++) gamma[i] = projector(gamma[i]);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Gamma matrix:\n");
      Tensor<double> temp = expectation(world, f, gamma);
      if(world.rank() == 0) print(temp);
   }

   // Done
   return gamma;
}

// Computes gamma(r) given the ground state orbitals and response functions
// Only for TDHF
std::vector<std::vector<real_function_3d>> TDHF::create_gamma_tdhf(World &world,
                                                                   std::vector<std::vector<real_function_3d>> & x,
                                                                   std::vector<std::vector<real_function_3d>> & y,
                                                                   std::vector<real_function_3d> & orbitals,
                                                                   double small,
                                                                   double thresh,
                                                                   int print_level,
                                                                   std::string xy)
{
   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Creating Gamma for %s states\n", xy.c_str());
   }

   // Get sizes
   int m = x.size();
   int n = x[0].size();

   // The gamma function to be returned, intialized to zero
   std::vector<std::vector<real_function_3d>> gamma = response_zero_functions(world, m, n);

   // Gamma will have 2 terms for HF: dJ/drho[rho] and dK/drho[rho]
   // There is a different Gamma for each orbital-->virtual transition
   // Calculate both here
   std::vector<std::vector<real_function_3d>> deriv_Jx = create_coulomb_derivative(world, x, orbitals, small, thresh);
   std::vector<std::vector<real_function_3d>> deriv_Jy = create_coulomb_derivative(world, y, orbitals, small, thresh);
   std::vector<std::vector<real_function_3d>> deriv_J = deriv_Jx + deriv_Jy;

   std::vector<std::vector<real_function_3d>> deriv_K = create_exchange_derivative_tdhf(world, x, y, orbitals, small, thresh);

   // Spin integration gives coefficients 
   // This is the spin restricted, singlet excitation coefficients
   gamma = scale(deriv_J, 2.0) - deriv_K;

   // Project out groundstate 
   QProjector<double, 3> projector(world, Gparams.orbitals);           // Projector to project out ground state
   for(int i = 0; i<m; i++) gamma[i] = projector(gamma[i]);

   // Debugging output
   if (print_level >= 2 and xy == "x")
   {
      if(world.rank() == 0) printf("   Coulomb Deriv matrix for x states:\n");
      Tensor<double> temp = expectation(world, x, deriv_J);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Exchange Deriv matrix for x states:\n");
      temp = expectation(world, x, deriv_K);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Gamma matrix for x states:\n");
      temp = expectation(world, x, gamma);
      if(world.rank() == 0) print(temp);
   }
   else if (print_level >= 2 and xy == "y")
   {
      if(world.rank() == 0) printf("   Coulomb Deriv matrix for y states:\n");
      Tensor<double> temp = expectation(world, y, deriv_J);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Exchange Deriv matrix for y states:\n");
      temp = expectation(world, y, deriv_K);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Gamma matrix for y states:\n");
      temp = expectation(world, y, gamma);
      if(world.rank() == 0) print(temp);
   }

   // Done
   return gamma;
}


// Calculates ground state coulomb potential 
real_function_3d TDHF::coulomb(World& world)
{
   // Coulomb operator
   real_convolution_3d op = CoulombOperator(world, Rparams.small, Rparams.thresh);

   // Get density
   std::vector<real_function_3d> vsq = square(world, Gparams.orbitals);
   compress(world, vsq);
   real_function_3d rho = real_factory_3d(world);
   rho.compress();
   for (unsigned int i = 0; i < vsq.size(); ++i) {
       rho.gaxpy(1.0, vsq[i], 1.0, false);
   }
   world.gop.fence();
   vsq.clear();

   // Apply operator and truncate
   rho = apply(op, rho);
   rho.truncate();

   // Done 
   return rho;
}

// Calculates HF exchange between ground state orbitals and functions f 
std::vector<std::vector<real_function_3d>> TDHF::exchange(World& world,
                                                         std::vector<std::vector<real_function_3d>> & f)
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Adding this because localized orbitals need to run over 
   // all the ground state orbitals on inner loop below, but
   // wouldn't without this last size variable
   int q = Gparams.orbitals.size();

   // Coulomb operator
   real_convolution_3d op = CoulombOperator(world, Rparams.small, Rparams.thresh);

   // Container for results and others
   std::vector<std::vector<real_function_3d>> result = response_zero_functions(world, m, n);
   real_function_3d psif = real_function_3d(world);

   // Run over each excited state
   for(int k = 0; k < m; k++)
   {
      // And run over each occupied state
      for(int p = 0; p < n; p++)
      {
         for(int j = 0; j < q; j++)
         {
            // Get transition density
            psif = Gparams.orbitals[j] * f[k][p];

            // Apply coulomb operator 
            psif = apply(op, psif);

            // Final multiplication
            result[k][p] += Gparams.orbitals[j] * psif;
         }
      }
   }

   // Truncate
   truncate(world, result);

   // Done!
   return result;
}
// Returns the ground state potential applied to functions f
std::vector<std::vector<real_function_3d>> TDHF::create_potential(World & world,
                                                                 std::vector<std::vector<real_function_3d>> & f,
                                                                 int print_level,
                                                                 std::string xy)
{
   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Computing V0 * f(r) for %s states\n", xy.c_str());
   }

   // Computing \hat{V}^0 = v_nuc + v_coul + v_exch           
   // v_nuc first
   PotentialManager manager(Gparams.molecule, "a");
   manager.make_nuclear_potential(world);
   real_function_3d v_nuc = manager.vnuclear().truncate();

   // V_coul next
   // This does not include final multiplication of each orbital 
   // 2 is from integrating out spin
   real_function_3d v_coul = 2.0 * coulomb(world);

   // Sum coulomb (pre multiplied) and v_nuc
   // v_nuc comes out negative from potential manager, so add it
   real_function_3d v = v_coul + v_nuc;

   // Apply V to f functions
   std::vector<std::vector<real_function_3d>> V_x_resp = multiply(f, v);

   // V_exch last
   // Multiplication by f functions is included in construction
   std::vector<std::vector<real_function_3d>> v_exch = exchange(world, f);

   if(print_level >= 2)
   {
      // Print potential energy matrices
      if(world.rank() == 0) printf("   Nuclear potential matrix for %s states:\n", xy.c_str());
      std::vector<std::vector<real_function_3d>> temp1 = multiply(f,v_nuc);
      Tensor<double> temp = expectation(world, f, temp1);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Coulomb potential matrix for %s states:\n", xy.c_str());
      std::vector<std::vector<real_function_3d>> temp2 = multiply(f,v_coul);
      temp = expectation(world, f, temp2);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Exchange potential matrix for %s states:\n", xy.c_str());
      temp = expectation(world, f, v_exch);
      if(world.rank() == 0) print(temp);
   }

   // Subtract V_exch from V_x_resp
   V_x_resp = V_x_resp - v_exch; 
 
   // Print some basic output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Total Potential Energy matrix for %s states:\n", xy.c_str());
      Tensor<double> temp = expectation(world, f, V_x_resp);
      if(world.rank() == 0) print(temp);
   }

   truncate(world, V_x_resp);

   // Done
   return V_x_resp;
}


// Returns a Tensor of inner products, where
// result(i,j) = inner(a[i],b[j]).sum()
Tensor<double> TDHF::expectation(World &world,
                                std::vector<std::vector<real_function_3d>> & a,
                                std::vector<std::vector<real_function_3d>> & b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size());
   MADNESS_ASSERT(b[0].size() > 0);

   // Get sizes
   int dim_a = a.size();
   int dim_b = b.size();

   // Container for result
   Tensor<double> result(dim_a, dim_b);

   // Run over dimension one
   for(int p = 0; p < dim_a; p++)
   {
      // Run over dimension two
      for(int k = 0; k < dim_b; k++)
      {
         result(p,k) = inner(world, a[p],b[k]).sum();
      }
   }

   // Done
   return result;
}

// Creating overlap matrix for given function f
Tensor<double> TDHF::create_overlap(World & world,
                                   std::vector<std::vector<real_function_3d>> & f,
                                   int print_level,
                                   std::string xy)
{
   // Get sizes
   int m = f.size();

   // Container for answer to be returned
   Tensor<double> S(m,m);

   // Run over all virtuals i 
   for(int i = 0; i < m; i++)
   {
      // Run over all virtuals j
      for(int j = 0; j < m; j++)
      {
         // Using vmra.h line 627 function 
         // Sum included inside inner()
         S(i,j) = inner(f[i], f[j]);
      }
   }

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Overlap matrix for %s states:\n", xy.c_str());
      if(world.rank() == 0) print(S);
   }

   // Done
   return S;
}

// Returns the ground state fock operator applied to functions f
std::vector<std::vector<real_function_3d>> TDHF::create_fock(World & world,
                                                            std::vector<std::vector<real_function_3d>> & V,
                                                            std::vector<std::vector<real_function_3d>> & f,
                                                            int print_level,
                                                            std::string xy)
{
   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Creating perturbed fock matrix for %s states\n", xy.c_str());
   }

   // Size of fock matrix must match that of V
   int m = V.size();
   int n = V[0].size();

   // Container to return
   std::vector<std::vector<real_function_3d>> fock =  response_zero_functions(world, m, n);

   // Fock = (T + V) * orbitals
   // Already have V
   // Create T
   // Make the derivative operators in each direction
   real_derivative_3d Dx(world, 0);
   real_derivative_3d Dy(world, 1);
   real_derivative_3d Dz(world, 2);

   // Apply derivatives to orbitals
   std::vector<std::vector<real_function_3d>> dvx = apply(world, Dx, f);
   std::vector<std::vector<real_function_3d>> dvy = apply(world, Dy, f);
   std::vector<std::vector<real_function_3d>> dvz = apply(world, Dz, f);

   // Apply again for 2nd derivatives
   std::vector<std::vector<real_function_3d>> dvx2 = apply(world, Dx, dvx);
   std::vector<std::vector<real_function_3d>> dvy2 = apply(world, Dy, dvy);
   std::vector<std::vector<real_function_3d>> dvz2 = apply(world, Dz, dvz);

   // Add together derivatives
   fock = scale((dvx2 + dvy2 + dvz2), -0.5);

   // Debugging output   
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Kinetic energy matrix for %s states:\n", xy.c_str());
      Tensor<double> temp = expectation(world, f, fock);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Potential energy matrix for %s states:\n", xy.c_str());
      temp = expectation(world, f, V);
      if(world.rank() == 0) print(temp);
   }

   // Add in potential
   fock = fock + V;

   truncate(world, fock);

   truncate(world, fock);

   // Done
   return fock;
}

// Construct the Hamiltonian
Tensor<double> TDHF::create_hamiltonian(World & world,
                                       std::vector<std::vector<real_function_3d>> & gamma,
                                       std::vector<std::vector<real_function_3d>> & V,
                                       std::vector<std::vector<real_function_3d>> & f,
                                       std::vector<real_function_3d> & ground_orbitals,
                                       std::vector<real_function_3d> & full_ground_orbitals,
                                       Tensor<double> & energies, // 2D tensor that is ground state hamiltonian now
                                       int print_level,
                                       std::string xy)
{
   // Sizes inferred from V and gamma
   int m = V.size();

   // Container for A
   Tensor<double> A(m,m);

   // Create the ground-state fock operator on response orbitals
   std::vector<std::vector<real_function_3d>> fock_x_resp = create_fock(world, V, f, print_level, xy);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Ground Fock matrix for %s states:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world, f, fock_x_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Need to calculate energy * x_response
   // This now works for localized or canonical orbitals ???
   std::vector<std::vector<real_function_3d>> energy_x_resp = scale_2d(world, f, energies);

   // Verify this keeps orbitals in the virtual space
   // Verify this annihilates an occupied orbital (leave occupied orbital in occupied space at least)

   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Energy scaled response orbitals for %s states:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world,f,energy_x_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Construct intermediary
   std::vector<std::vector<real_function_3d>> temp = gamma + fock_x_resp - energy_x_resp;

   // Need to run over excited states
   for(int k = 0; k < m; k++)
   {
      // Need another run over excited states
      for(int j = 0; j < m; j++)
      {
         // Project out ground state densities
         // Shouldn't need this - taken care of elsewhere now
         //temp[j] = ground_state_density(temp[j]);

         // Run over all occupied orbitals to get their contribution
         // to the part of A we're calculating . Finally calculate 
         // \int dr x_p^k * temp (using vmra.h function, sum is included)
         A(k,j) = inner(f[k], temp[j]);
      }
   }

   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Hamiltonian matrix for %s states:\n", xy.c_str());
      if(world.rank() == 0) print(A);
   }

   // Done
   return A;
}
// Returns the shift needed to make sure that
// -2.0 * (ground_state_energy + excited_state_energy) 
// is negative. Please note: The same shift needs to 
// be applied to the potential.
Tensor<double> TDHF::create_shift(World & world,
                                 Tensor<double> & ground,
                                 Tensor<double> & omega,
                                 int print_level,
                                 std::string xy)
{
   // Get sizes
   int m = omega.size();
   int n = ground.size();

   // Container to hold shift
   Tensor<double> result(m,n);

   // Run over excited states
   for(int k = 0; k < m; k++)
   {
      // Run over ground states
      for(int p = 0; p < n; p++)
      {
         if(ground(p) + omega(k) > 0)
         {
            // Calculate the shift needed to get energy to -0.05,
            // which was arbitrary (same as moldft)
            result(k,p) = -(ground(p) + omega(k) + 0.05);

            // Basic output
            if(print_level >= 1)
            {
               if(world.rank() == 0) printf("   Shift needed for transition from ground orbital %d to state %s response orbital %d\n", p, xy.c_str(), k);
               if(world.rank() == 0) print("   Ground energy =", ground(p));
               if(world.rank() == 0) print("   Excited energy =", omega(k));
               if(world.rank() == 0) print("   Shifting by", result(k,p));
               if(world.rank() == 0) print("");
            }
         }
      }
   }

   // Done
   return result;
}
// Returns the given shift applied to the given potential
std::vector<std::vector<real_function_3d>> TDHF::apply_shift(World & world,
                                                            Tensor<double> & shifts,
                                                            std::vector<std::vector<real_function_3d>> & V,
                                                            std::vector<std::vector<real_function_3d>> & f)
{
   // Sizes inferred from V
   int n = V[0].size();
   int m = V.size();

   // Container to return
   std::vector<std::vector<real_function_3d>> shifted_V = response_zero_functions(world, m, n);

   // Run over occupied
   for(int k = 0; k < m; k++)
   {
      // Run over virtual
      for(int p = 0; p < n; p++)
      {
         shifted_V[k][p] =  V[k][p] + shifts(k,p) * f[k][p];
      }
   }

   truncate(world, shifted_V);

   // Done
   return shifted_V;
}

// Function to make a vector of BSH operators using ground and excited
// state energies
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> TDHF::create_bsh_operators(World & world,
                                                                                         Tensor<double> & shift,
                                                                                         Tensor<double> & ground,
                                                                                         Tensor<double> & omega,
                                                                                         double small,
                                                                                         double thresh)
{
   // Sizes inferred from ground and omega
   int n = ground.size();
   int m = omega.size();

   // Make the vector
   std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> operators;

   // Make a BSH operator for each response function
   // Run over excited states
   for(int k = 0; k < m; k++)
   {
      // Container for intermediary
      std::vector<std::shared_ptr<real_convolution_3d>> temp(n);

      // Run over occupied states
      for(int p = 0; p < n; p++)
      {
         temp[p] = std::shared_ptr<SeparatedConvolution<double,3>>(BSHOperatorPtr3D(world, sqrt(-2.0 * (ground(p) + omega(k) + shift(k,p))), small, thresh));
      }

      // Add intermediary to return container
      operators.push_back(temp);
   }

   // Done
   return operators;
}

// Returns the second order update to the energies of the excited states
// Not currently used.
Tensor<double> TDHF::calculate_energy_update(World & world,
                                            std::vector<std::vector<real_function_3d>> & rhs,
                                            std::vector<std::vector<real_function_3d>> & f_residuals,
                                            std::vector<std::vector<real_function_3d>> & new_f,
                                            int print_level,
                                            std::string xy)
{
   /*
    *  The correction is:
    *      \delta \omega^{(k)} = - \frac{ \sum_p\left< \hat{V}^0 x_p^{(k)}(r) + (1 - \hat{\rho}^0) \Gamma_p^{(k)}(r)\right|
    *                                         \left. x_p^{(k)} - \~{x}_p^{(k)} \right> }
    *                                   { \sum_p \left| \left| \~{x}_p^{(k)} \right| \right|^2 }
    */

   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Calculating energy residy residuals for %s states\n", xy.c_str());
   }

   // Size inferred
   int m = rhs.size();

   // Container for updates
   Tensor<double> updates(m);

   // Need to run over all functions in rhs and calculate inner products.
   // rhs contains the bra in the braket notation above, and f_residuals
   // is the ket.

   // Run over excited states
   for(int k = 0; k < m; k++)
   {
      // vmra.h function, line 627
      // Sum is included inside function call
      updates(k) = inner(f_residuals[k],rhs[k]);

      // Normalize update function
      // The -1.0 is the leading coefficient in the update formula
      // the 1/2 is to undo the scaling of V
      updates(k) = -1.0/2.0 * updates(k) / inner(new_f[k], new_f[k]);
   }

   if(print_level >= 1)
   {
      // Print energy deltas
      if(world.rank() == 0) printf("   Energy residuals for %s states:\n", xy.c_str());
      if(world.rank() == 0) print(updates);
   }

   // Done?
   return updates;
}

// Specialized for response calculations that returns orthonormalized
// functions
std::vector<std::vector<real_function_3d>> TDHF::gram_schmidt(World & world,
                                                              std::vector<std::vector<real_function_3d>> & f)
{
   // Sizes inferred
   int m = f.size();

   // Return container
   std::vector<std::vector<real_function_3d>> result = copy(world, f);

   // Orthogonalize
   for(int j = 0; j < m; j++)
   {
      // Need to normalize the row
      double norm = norm2(world, result[j]);

      // Now scale each entry      
      result[j] = result[j] * (1.0/norm);

      // Project out from the rest of the vectors
      for(int k = j+1; k < m; k++)
      {
         // Temp function to hold the sum
         // of inner products
         // vmra.h function, line 627
         double temp = inner(result[j], result[k]);

         // Now subtract 
         result[k] = result[k] - temp * result[j];
      }
   }

   truncate(world, result);

   // Done
   return result;
}

// Returns the max norm of the given vector of functions
double TDHF::calculate_max_residual(World & world,
                                   std::vector<std::vector<real_function_3d>> & f)
{
   // Container for max
   double max = 0.0;

   // Run over all functions in f
   for(unsigned int i = 0; i < f.size(); i++)
   {
      double temp = 0.0;

      for(unsigned int j = 0; j < f[0].size(); j++)
      {
         temp += pow(f[i][j].norm2(),2);
      }

      temp = sqrt(temp);

      if( temp > max) max = temp;
   }

   // Done
   return max;
}

// Selects the 'active' orbitals from ground state orbitals to be used in the calculation (based on energy distance
// from the HOMO). Function needs knowledge of Gparams.orbitals and Gparams.energies. Function sets act_orbitals and
// num_act_orbitals.
void TDHF::select_active_subspace(World & world)
{
   // Default output
   if( Rparams.print_level >= 0)
   {
      // Set print output to something reasonable
      std::cout.precision(2);
      std::cout << std::fixed;

      if(world.rank() == 0) print("   Selecting ground state subspace to excite from for states.");
      if(world.rank() == 0) print("   This is all orbitals between", Rparams.range_low, "and", Rparams.range_high, "\n");

      // Reset precision
      std::cout.precision(10);
      std::cout << std::scientific;
  }

   // Determine active orbitals based on energy differences 
   // from HOMO
   for(unsigned int i = 0; i < Gparams.num_orbitals; i++)
   {
      
      if(Rparams.range_low < Gparams.energies(i) and Gparams.energies(i) < Rparams.range_high)
      {
         // This orbital should be active, so add to list
         //active.push_back(i);
         active.push_back(i);
      }
   }

   // Make sure we have at least one ground state orbital to excite from
   MADNESS_ASSERT(active.size() > 0); 

   // Now that we know size, allocate act_ground_energies
   act_ground_energies = Tensor<double>(active.size());

   // Now to pull the functions and energies and add to act_orbitals and act_ground_energies
   for(unsigned int i = 0; i < active.size(); i++)
   {
      act_orbitals.push_back(Gparams.orbitals[active[i]]);
      act_ground_energies(i) = Gparams.energies(active[i]); // Put energies on diagonal
   }

   // Also set the active size
   act_num_orbitals = act_orbitals.size();

   print("Found", act_num_orbitals, "active orbitals.");
}

// Selects from a list of functions and energies the k functions with the lowest energy
std::vector<std::vector<real_function_3d>> TDHF::select_trial_functions(World & world,
                                                                       std::vector<std::vector<real_function_3d>> & f,
                                                                       Tensor<double> & energies,
                                                                       int k,
                                                                       int print_level)
{
   // Container for result
   std::vector<std::vector<real_function_3d>> answer;

   // Basic output
   if(print_level >= 0)
   {
      if(world.rank() == 0) print("\n   Guess functions created. Selecting the", k, "lowest energy states.\n");
   }

   // No energy updates or function differences, so create dummies for sort( ) function
   Tensor<double> dummy(energies.size());
   std::vector<std::vector<real_function_3d>> dummy2 = response_zero_functions(world, f.size(), f[0].size());

   // Sort by the energy
   // NOTE: sort() modifies in all its arguments 
   Tensor<int> selected = sort(world, energies, dummy, f, dummy2);

   // Pull out first k from selected.
   Tensor<int> k_selected(k);
   for(int i = 0; i < k; i++) k_selected(i) = selected(i);

   // Basic output
   if(print_level >= 0)
   {
      if(world.rank() == 0) print("   The selected states are:");
      if(world.rank() == 0) print(k_selected);
   }

   // Now just take the first k functions
   for(int i = 0; i < k; i++) answer.push_back(copy(world, f[i]));

   truncate(world, answer);

   // Done
   return answer;
}

// Calculate the exponentiation of a matrix through first order (I think)
Tensor<double> TDHF::matrix_exponential(const Tensor<double> & A)
{
   const double tol = 1e-13;
   MADNESS_ASSERT(A.dim((0) == A.dim(1)));

   // Scale A by a power of 2 until it is "small"
   double anorm = A.normf();
   int n = 0;
   double scale = 1.0;
   while (anorm * scale > 0.1)
   {
      ++n;
      scale *= 0.5;
   }
   Tensor<double> B = scale * A;    // B = A*2^-n

   // Compute exp(B) using Taylor series
   Tensor<double> expB = Tensor<double>(2, B.dims());
   for (int i = 0; i < expB.dim(0); ++i)
      expB(i, i) = 1.0;

   int k = 1;
   Tensor<double> term = B;
   while (term.normf() > tol)
   {
      expB += term;
      term = inner(term, B);
      ++k;
      term.scale(1.0 / k);
   }

   // Repeatedly square to recover exp(A)
   while (n--)
   {
      expB = inner(expB, expB);
   }

   return expB;
}

/// compute the unitary transformation that diagonalizes the fock matrix

/// @param[in]  world   the world
/// @param[in]  overlap the overlap matrix of the orbitals
/// @param[inout]       fock    the fock matrix; diagonal upon exit
/// @param[out] evals   the orbital energies
/// @param[in]  thresh_degenerate       threshold for orbitals being degenerate
/// @return             the unitary matrix U: U^T F U = evals
Tensor<double> TDHF::get_fock_transformation(World & world,
                                            const Tensor<double> & overlap,
                                            Tensor<double> & fock,
                                            Tensor<double> & evals,
                                            const double thresh_degenerate)
{
   // Diagonalize using lapack
   Tensor<double> U;
   sygvp(world, fock, overlap, 1, U, evals);

   long nmo = fock.dim(0);

   bool switched = true;
   while (switched)
   {
      switched = false;
      for (int i = 0; i < nmo; i++)
      {
         for (int j = i + 1; j < nmo; j++)
         {
            double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
            double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
            if (snew > sold)
            {
               Tensor<double> tmp = copy(U(_, i));
               U(_, i) = U(_, j);
               U(_, j) = tmp;
               std::swap(evals[i], evals[j]);
               switched = true;
            }
         }
      }
   }


   // Fix phases.
   for (long i = 0; i < nmo; ++i)
      if (U(i, i) < 0.0)
         U(_, i).scale(-1.0);

   // Rotations between effectively degenerate states confound
   // the non-linear equation solver ... undo these rotations
   long ilo = 0; // first element of cluster
   while (ilo < nmo - 1)
   {
      long ihi = ilo;
      while (fabs(evals[ilo] - evals[ihi + 1]) < thresh_degenerate * 10.0 * std::max(fabs(evals[ilo]), 1.0))
      {
         ++ihi;
         if (ihi == nmo - 1)
            break;
      }
      long nclus = ihi - ilo + 1;
      if (nclus > 1)
      {
         Tensor<double> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));

         // Polar Decomposition
         Tensor<double> VH(nclus, nclus);
         Tensor<double> W(nclus, nclus);
         Tensor<double> sigma(nclus);

         svd(q, W, sigma, VH);
         q = transpose(inner(W,VH));  // Should be conj. tranpose if complex
         U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);


         //  Robert's Rotations
         //          
         // Iteratively construct unitary rotation by
         // exponentiating the antisymmetric part of the matrix
         // ... is quadratically convergent so just do 3
         // iterations
         //Tensor<double> rot = matrix_exponential(-0.5 * (q - transpose(q)));
         //q = inner(q, rot);
         //Tensor<double> rot2 = matrix_exponential(-0.5 * (q - transpose(q)));
         //q = inner(q, rot2);
         //Tensor<double> rot3 = matrix_exponential(-0.5 * (q - transpose(q)));
         //q = inner(rot, inner(rot2, rot3));
         //U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
      }
      ilo = ihi + 1;
   }

   fock = 0;
   for (unsigned int i = 0; i < nmo; ++i)
      fock(i, i) = evals(i);
   return U;
}

/// diagonalize the fock matrix, taking care of degenerate states

/// Vpsi is passed in to make sure orbitals and Vpsi are in phase
/// @param[in]  world   the world
/// @param[inout]       fock    the fock matrix (diagonal upon exit)
/// @param[inout]       psi     the orbitals
/// @param[inout]       Vpsi    the orbital times the potential
/// @param[inout]       gamma   the orbital times the perturbed potential
/// @param[out] evals   the orbital energies
/// @param[in]  overlap the overlap matrix
/// @param[in]  thresh  threshold for rotation and truncation
/// @return             the unitary matrix U: U^T F U = evals
Tensor<double> TDHF::diag_fock_matrix(World & world,
                                     Tensor<double> & fock,
                                     std::vector<std::vector<real_function_3d>> & psi,
                                     std::vector<std::vector<real_function_3d>> & Vpsi,
                                     std::vector<std::vector<real_function_3d>> & gamma,
                                     Tensor<double> & evals,
                                     Tensor<double> & overlap,
                                     const double thresh)
{
    // compute the unitary transformation matrix U that diagonalizes
    // the fock matrix
    Tensor<double> U = get_fock_transformation(world, overlap, fock, evals, thresh);

    // transform the orbitals and the potential
    Vpsi = transform(world, Vpsi, U);
    gamma = transform(world, gamma, U);
    psi = transform(world, psi, U);

    // truncate all and normalize psi
    truncate(world, Vpsi);
    truncate(world, gamma);
    truncate(world, psi);
    normalize(world, psi);

    return U;
}

// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
std::vector<std::vector<real_function_3d>> TDHF::transform(World & world,
                                                          std::vector<std::vector<real_function_3d>> & f,
                                                          Tensor<double> & U)
{
   // Return container
   std::vector<std::vector<real_function_3d>> result;

   // Go element by element 
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Temp for the result of one row
      std::vector<real_function_3d> temp = zero_functions_compressed<double,3>(world, f[0].size());

      for(unsigned int j = 0; j < f.size(); j++)
      {
         gaxpy(world, 1.0, temp, U(j,i), f[j]);
      }

      // Add to temp to result
      result.push_back(temp);
   }

   truncate(world, result);

   // Done
   return result;
}

// Sorts the given Tensor of energies and vector of functions
// in place
Tensor<int> TDHF::sort(World & world,
                      Tensor<double> & vals,
                      Tensor<double> & val_residuals,
                      std::vector<std::vector<real_function_3d>> & f,
                      std::vector<std::vector<real_function_3d>> & f_diff)

{
   // Get relevant sizes
   int k = vals.size();

   // Tensor to hold selection order
   Tensor<int> selected(k);

   // Copy everything... 
   std::vector<std::vector<real_function_3d>> f_copy = copy(world, f);
   std::vector<std::vector<real_function_3d>> f_diff_copy = copy(world, f_diff);
   std::vector<double> vals_copy;
   for(int i = 0; i < k; i++) vals_copy.push_back(vals[i]);
   Tensor<double> vals_copy2 = copy(vals);
   Tensor<double> val_residuals_copy = copy(val_residuals);

   // Clear the vectors
   f.clear();
   f_diff.clear();

   // Now sort vals_copy
   std::sort(vals_copy.begin(), vals_copy.end());

   // Now sort the rest of the things, using the sorted energy list
   // to find the correct indices 
   for(int i = 0; i < k; i++)
   {
      // Find matching index in sorted vals_copy
      int j = 0;
      while(fabs(vals_copy[i] - vals_copy2[j]) > 1e-8 && j < k) j++;

      // Add in to list which one we're taking
      selected(i) = j;

      // Put corresponding function, difference function, value residual and value
      // in the correct place
      f.push_back(f_copy[j]);
      f_diff.push_back(f_diff_copy[j]);
      vals(i) =  vals_copy[i];
      val_residuals[i] = val_residuals_copy(j);

      // Change the value of vals_copy2[j] to help deal with duplicates?
      vals_copy2[j] = 10000.0;
   }

   // Done
   return selected;
}

// Iterates the response functions until converged or out of iterations
void TDHF::iterate(World & world)
{
   // Variables needed to iterate
   int iteration = 0;                                                // Iteration counter
   QProjector<double, 3> projector(world, Gparams.orbitals);         // Projector to project out ground state
   int n = Gparams.num_orbitals;                                     // Number of ground state orbitals
   int m = Rparams.states;                                           // Number of excited states
   bool converged = false;                                           // For convergence
   Tensor<double> old_x_energy(m);                                   // Holds previous iteration's energy
   Tensor<double> old_y_energy(m);                                   // Holds previous iteration's energy
   Tensor<double> energy_x_residuals;                                // Holds energy residuals 
   Tensor<double> energy_y_residuals;                                // Holds energy residuals 
   Tensor<double> x_shifts;                                          // Holds the shifted energy values
   Tensor<double> y_shifts;                                          // Holds the shifted energy values
   std::vector<std::vector<real_function_3d>> x_differences;         // Holds wave function corrections
   std::vector<std::vector<real_function_3d>> y_differences;         // Holds wave function corrections
   std::vector<std::vector<real_function_3d>> temp;                  // Used for step restriction 
   std::vector<std::vector<real_function_3d>> x_gamma;               // Holds the perturbed two electron piece
   std::vector<std::vector<real_function_3d>> y_gamma;               // Holds the perturbed two electron piece
   std::vector<std::vector<real_function_3d>> V_x_response;          // Holds V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> V_y_response;          // Holds V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> shifted_V_x_response;  // Holds the shifted V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> shifted_V_y_response;  // Holds the shifted V^0 applied to response functions

   // The KAIN solver
   XNonlinearSolver<std::vector<std::vector<real_function_3d>>, double, TDHF_allocator> kain(TDHF_allocator(world, m, n), Rparams.print_level);

   // Setting max sub size for KAIN solver
   kain.set_maxsub(5);

   // Get a start time
   Tensor<double> initial_time = end_timer(world);

   // Now to iterate
   while( iteration < Rparams.max_iter  && !converged)
   {
      // Start a timer
      Tensor<double> iter_time = end_timer(world);

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("\n   Iteration", iteration);
         if(world.rank() == 0) print("  --------------");
      }

      // Create gamma
      if(Rparams.tda) x_gamma = create_gamma(world, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh, Rparams.print_level);
      else
      {
         x_gamma = create_gamma_tdhf(world, x_response, y_response, Gparams.orbitals, Rparams.small, Rparams.thresh, Rparams.print_level, "x");
         y_gamma = create_gamma_tdhf(world, y_response, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh, Rparams.print_level, "y");
      }

      // Create \hat{V}^0 applied to x_response
      V_x_response = create_potential(world, x_response, Rparams.print_level, "x");

      if(not Rparams.tda) V_y_response = create_potential(world, y_response, Rparams.print_level, "y");

      // Load balance
      // Only balancing on x-states. Smart?
      if(Rparams.print_level >= 1)
      {
         if(world.rank()==0) print("\n   Load balancing using orbitals and the potential.");
      }
      LoadBalanceDeux<3> lb(world);
      for(int j = 0; j < n; j++)
      {
         for(int k = 0; k < Rparams.states; k++)
         {
            lb.add_tree(x_response[k][j], lbcost<double,3>(1.0,8.0),true);
            lb.add_tree(V_x_response[k][j], lbcost<double,3>(1.0,8.0), true);
            lb.add_tree(x_gamma[k][j], lbcost<double,3>(1.0,8.0), true);
         }
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));
      if(world.rank() == 0) print("");

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("   Solving Ax=Swx");
      }

      // Constructing S
      Tensor<double> S_x = create_overlap(world, x_response, Rparams.print_level, "x");
      
      Tensor<double> S_y;
      if(not Rparams.tda) S_y = create_overlap(world, y_response, Rparams.print_level, "y"); 

      // Constructing hamiltonian 
      Tensor<double> A_x = create_hamiltonian(world, x_gamma, V_x_response, x_response, Gparams.orbitals, Gparams.orbitals, hamiltonian, Rparams.print_level, "x");
      Tensor<double> A_y;
      if(not Rparams.tda) A_y = create_hamiltonian(world, y_gamma, V_y_response, y_response, Gparams.orbitals, Gparams.orbitals, hamiltonian, Rparams.print_level, "y");

      // Solve Ax = Sxw
      Tensor<double> U_x = diag_fock_matrix(world, A_x, x_response, V_x_response, x_gamma, x_omega, S_x, Rparams.thresh);
      Tensor<double> U_y;
      if(not Rparams.tda) U_y = diag_fock_matrix(world, A_y, y_response, V_y_response, y_gamma, y_omega, S_y, Rparams.thresh);

      // Debugging output
      if(Rparams.print_level >= 2 and world.rank() == 0)
      {
         print("   Eigenvector coefficients from diagonalization for x states:");
         print(U_x);

         if(not Rparams.tda)
         {
            print("   Eigenvector coefficients from diagonalization for y states:");
            print(U_y);
         }
      }

      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      x_shifts = create_shift(world, Gparams.energies, x_omega, Rparams.print_level, "x");

      if(not Rparams.tda) y_shifts = create_shift(world, Gparams.energies, y_omega, Rparams.print_level, "y");

      // Apply the shifts
      shifted_V_x_response = apply_shift(world, x_shifts, V_x_response, x_response);

      if(not Rparams.tda) shifted_V_y_response = apply_shift(world, y_shifts, V_y_response, y_response);

      // Construct RHS of equation
      std::vector<std::vector<real_function_3d>> rhs_x = x_gamma + shifted_V_x_response;
      std::vector<std::vector<real_function_3d>> rhs_y;

      if(not Rparams.tda)
      {        
         rhs_y = shifted_V_y_response + y_gamma;  

         // Scaling rhs_y by -1 to get the negative sign from eq. 33
         scale(rhs_y, -1.0); 
     }

      // Basic output
      if(Rparams.print_level >= 1 and world.rank() == 0)
      {
         print("   Response Orbital Energies:");
         print("   x states:"); 
         print(x_omega);
       
         if(not Rparams.tda)
         {
            print("   y states:");
            print(y_omega);
         }
      }

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms of RHS of main equation:");
         if(world.rank() == 0) print("   x states:");
         print_norms(world, rhs_x);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y states:");
            print_norms(world, rhs_y);
         }
      }

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_operators = create_bsh_operators(world, x_shifts, Gparams.energies, x_omega, Rparams.small, Rparams.thresh);
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_y_operators;

      if(not Rparams.tda) bsh_y_operators = create_bsh_operators(world, y_shifts, Gparams.energies, y_omega, Rparams.small, Rparams.thresh);
     
      // TESTING
      for(int i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
      if(not Rparams.tda) for(int i = 0; i < m; i++) y_response[i] = projector(y_response[i]);
 
      // Apply BSH operators to RHS of equation
      if(Rparams.print_level >= 1 and world.rank() == 0)
      {
         print("   Applying BSH operators\n");
      }
      std::vector<std::vector<real_function_3d>> new_x_response = apply(world, bsh_x_operators, rhs_x);
      std::vector<std::vector<real_function_3d>> new_y_response;

      if(not Rparams.tda) new_y_response = apply(world, bsh_y_operators, rhs_y);

      // Scale by -2.0 (coefficient in eq. 37 of reference paper)
      new_x_response = scale(new_x_response, -2.0);
      if(not Rparams.tda) new_y_response = scale(new_y_response, -2.0);

      // Project out ground state
      for(int i = 0; i < m; i++) new_x_response[i] = projector(new_x_response[i]);
      if(not Rparams.tda) for(int i = 0; i < m; i++) new_y_response[i] = projector(new_y_response[i]);

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms after application of BSH");
         if(world.rank() == 0) print("   x-states:");
         print_norms(world, new_x_response);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y-states:");
            print_norms(world, new_y_response);
         }
      }

      // Get the difference between old and new
      x_differences = x_response - new_x_response;

      if(not Rparams.tda) y_differences = y_response - new_y_response;

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("   Response function residuals:");
         if(world.rank() == 0) print("   x states:");
         print_norms(world, x_differences);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y states:");
            print_norms(world, y_differences);
         }
      }

      // KAIN solver update 
      // Returns next set of orbitals
      // If not kain, save the new orbitals
      if(Rparams.kain)
      {
         x_response = new_x_response; //kain.update(x_response, x_differences); 
         y_response = new_y_response; //kain.update(y_response, y_differences);
      }
      else
      {
         x_response = new_x_response;
         y_response = new_y_response;
      }

      // Calculate energy residual and update old_energy 
      energy_x_residuals = abs(x_omega - old_x_energy);
      old_x_energy = x_omega;

      if(not Rparams.tda)
      {
         energy_y_residuals = abs(y_omega - old_y_energy);
         old_y_energy = y_omega;
      }

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("   Energy residuals:");
         if(world.rank() == 0) print("   x states:");
         if(world.rank() == 0) print(energy_x_residuals);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y states:");
            if(world.rank() == 0) print(energy_y_residuals);
         }
      }

      // Check convergence
      if(not Rparams.tda)
      {
         if(iteration >=1 && energy_x_residuals.absmax() < Rparams.econv && energy_y_residuals.absmax() < Rparams.econv)
            converged = true;
      }
      else
      {
         if(iteration >= 1 && energy_x_residuals.absmax() < Rparams.econv) 
            converged = true;
      }

      // Update counter
      iteration += 1;

      // TESTING
      // Orthogonalize
      //x_response = gram_schmidt(world, x_response);
      //if(not Rparams.tda) y_response = gram_schmidt(world, y_response);

      // Done with the iteration.. normalize and truncate
      normalize(world, x_response);
      truncate(world, x_response);

      if(not Rparams.tda)
      {
         normalize(world, y_response);
         truncate(world, y_response);
      }

      // Basic output
      if(Rparams.print_level >= 1)
      {
         // Precision is set to 10 coming in, drop it to 2
         std::cout.precision(2);
         std::cout << std::fixed;

         Tensor<double> this_time = end_timer(world);
         Tensor<double> current_time = this_time - iter_time;
         Tensor<double> total_time = this_time - initial_time;
         if(world.rank() == 0) print("   Time this iteration:", current_time[0],"s");
         if(world.rank() == 0) print("   Total time in iterations:", total_time[0],"s\n");

         // Reset precision
         std::cout.precision(10);
         std::cout << std::scientific;
      }
   }

   if(world.rank() == 0) print("\n");
   if(world.rank() == 0) print("   Finished TDHF Calculation ");
   if(world.rank() == 0) print("   ------------------------");
   if(world.rank() == 0) print("\n");

   // Did we converge?
   if(iteration == Rparams.max_iter && not converged)
   {
      if(world.rank() == 0) print("   Failed to converge. Reason:");
      if(world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
      if(world.rank() == 0) print("    Running analysis on current values.\n");
   }

   // Sort values and functions into ascending order based on values
   sort(world, x_omega, energy_x_residuals, x_response, x_differences);
   if(not Rparams.tda) sort(world, y_omega, energy_y_residuals, y_response, y_differences);

   // Print final things 
   if(world.rank() == 0) print(" Final x-state energies:");
   if(world.rank() == 0) print(x_omega);
   if(world.rank() == 0) print(" Final x-state energy residuals:");
   if(world.rank() == 0) print(energy_x_residuals);
   if(world.rank() == 0) print(" Final x-state response function residuals:");
   print_norms(world, x_differences);
   
   if(not Rparams.tda)
   {
      if(world.rank() == 0) print(" Final y-state energies:");
      if(world.rank() == 0) print(y_omega);
      if(world.rank() == 0) print(" Final y-state energy residuals:");
      if(world.rank() == 0) print(energy_y_residuals);
      if(world.rank() == 0) print(" Final y-state response function residuals:");
      print_norms(world, y_differences);
   }

   // A little more detailed analysis
   if(Rparams.tda) analysis_tda(world, x_response, x_omega);
   else analysis_tdhf(world, x_response, y_response, x_omega);

}   // Done with iterate. 

// More detailed analysis of the response functions
void TDHF::analysis_tda(World & world,
                        std::vector<std::vector<real_function_3d>> f,
                        Tensor<double> energies)
{
   // Sizes get used a lot here, so lets get a local copy
   int n = f[0].size();
   int m = f.size();

   // Per response function, want to print the contributions from each ground state
   // So print the norm of each function?
   Tensor<double> norms(m, n); 

   // Calculate the inner products
   for(int i = 0; i < m; i++)
   {
      for(int j = 0; j < n; j++)
      {
         norms(i,j) = f[i][j].norm2();
      }
   }

   // Need these to calculate dipole/quadrapole
   real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{1,0,0})));
   real_function_3d y = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,1,0})));
   real_function_3d z = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,0,1})));

   // Calculate transition dipole moments for each response function
   Tensor<double> dipoles(m, 3);

   // Run over each excited state
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         dipoles(i,0) += inner(Gparams.orbitals[j], x * f[i][j]);
         dipoles(i,1) += inner(Gparams.orbitals[j], y * f[i][j]);
         dipoles(i,2) += inner(Gparams.orbitals[j], z * f[i][j]);
      }
   }

   // Calculate transition quadrapole moments
   Tensor<double> quadrapoles(m,3,3);

   // Run over each excited state 
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         quadrapoles(i,0,0) += inner(Gparams.orbitals[j], x * x * f[i][j]);
         quadrapoles(i,0,1) += inner(Gparams.orbitals[j], x * y * f[i][j]);
         quadrapoles(i,0,2) += inner(Gparams.orbitals[j], x * z * f[i][j]);
         quadrapoles(i,1,0) += inner(Gparams.orbitals[j], y * x * f[i][j]);
         quadrapoles(i,1,1) += inner(Gparams.orbitals[j], y * y * f[i][j]);
         quadrapoles(i,1,2) += inner(Gparams.orbitals[j], y * z * f[i][j]);
         quadrapoles(i,2,0) += inner(Gparams.orbitals[j], z * x * f[i][j]);
         quadrapoles(i,2,1) += inner(Gparams.orbitals[j], z * y * f[i][j]);
         quadrapoles(i,2,2) += inner(Gparams.orbitals[j], z * z * f[i][j]);
      }
   }

   // Now print?
   if(world.rank() == 0)
   {
      for(int i = 0; i < m; i++)
      {
         printf("   Response Function %d\t\t%7.8f a.u.", i, energies[i]);
         print ("\n   --------------------------------------------");

         print("\n   Transition Dipole Moments");
         printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", dipoles(i,0), dipoles(i,1), dipoles(i,2));

         print("\n   Transition Quadrapole Moments");
         printf("   %16s %16s %16s\n", "X", "Y", "Z");
         printf(" X %16.8f %16.8f %16.8f\n", quadrapoles(i,0,0), quadrapoles(i,0,1), quadrapoles(i,0,2));
         printf(" Y %16.8f %16.8f %16.8f\n", quadrapoles(i,1,0), quadrapoles(i,1,1), quadrapoles(i,1,2));
         printf(" Z %16.8f %16.8f %16.8f\n", quadrapoles(i,2,0), quadrapoles(i,2,1), quadrapoles(i,2,2));

         // Print contributions
         print("\n   Norms of the Components:");
         for(int j = 0; j < n; j++)
         {
            printf("   Occupied %d  --->  Virtual %d   %7.8f\n", j, i, norms(i,j));
         }

         print("\n");
      }
   }
}

// More detailed analysis of the response functions
void TDHF::analysis_tdhf(World & world,
                         std::vector<std::vector<real_function_3d>> f,
                         std::vector<std::vector<real_function_3d>> g,
                         Tensor<double> energies)
{
   // Sizes get used a lot here, so lets get a local copy
   int n = f[0].size();
   int m = f.size();

   // Per response function, want to print the contributions from each ground state
   // So print the norm of each function?
   Tensor<double> norms(m, n); 

   // Calculate the inner products
   for(int i = 0; i < m; i++)
   {
      for(int j = 0; j < n; j++)
      {
         norms(i,j) = f[i][j].norm2();
      }
   }

   // Need these to calculate dipole/quadrapole
   real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{1,0,0})));
   real_function_3d y = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,1,0})));
   real_function_3d z = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,0,1})));

   // Calculate transition dipole moments for each response function
   Tensor<double> dipoles(m, 3);

   // Run over each excited state
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         dipoles(i,0) += inner(Gparams.orbitals[j], x * f[i][j]) + inner(Gparams.orbitals[j], x * g[i][j]);
         dipoles(i,1) += inner(Gparams.orbitals[j], y * f[i][j]) + inner(Gparams.orbitals[j], y * g[i][j]);
         dipoles(i,2) += inner(Gparams.orbitals[j], z * f[i][j]) + inner(Gparams.orbitals[j], z * g[i][j]);
      }
   }

   // Calculate transition quadrapole moments
   Tensor<double> quadrapoles(m,3,3);

   // Run over each excited state 
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         quadrapoles(i,0,0) += inner(Gparams.orbitals[j], x * x * f[i][j]) + inner(Gparams.orbitals[j], x * x * g[i][j]);
         quadrapoles(i,0,1) += inner(Gparams.orbitals[j], x * y * f[i][j]) + inner(Gparams.orbitals[j], x * y * g[i][j]);
         quadrapoles(i,0,2) += inner(Gparams.orbitals[j], x * z * f[i][j]) + inner(Gparams.orbitals[j], x * z * g[i][j]);
         quadrapoles(i,1,0) += inner(Gparams.orbitals[j], y * x * f[i][j]) + inner(Gparams.orbitals[j], y * x * g[i][j]);
         quadrapoles(i,1,1) += inner(Gparams.orbitals[j], y * y * f[i][j]) + inner(Gparams.orbitals[j], y * y * g[i][j]);
         quadrapoles(i,1,2) += inner(Gparams.orbitals[j], y * z * f[i][j]) + inner(Gparams.orbitals[j], y * z * g[i][j]);
         quadrapoles(i,2,0) += inner(Gparams.orbitals[j], z * x * f[i][j]) + inner(Gparams.orbitals[j], z * x * g[i][j]);
         quadrapoles(i,2,1) += inner(Gparams.orbitals[j], z * y * f[i][j]) + inner(Gparams.orbitals[j], z * y * g[i][j]);
         quadrapoles(i,2,2) += inner(Gparams.orbitals[j], z * z * f[i][j]) + inner(Gparams.orbitals[j], z * z * g[i][j]);
      }
   }

   // Now print?
   if(world.rank() == 0)
   {
      for(int i = 0; i < m; i++)
      {
         printf("   Response Function %d\t\t%7.8f a.u.", i, energies[i]);
         print ("\n   --------------------------------------------");

         print("\n   Transition Dipole Moments");
         printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", dipoles(i,0), dipoles(i,1), dipoles(i,2));

         print("\n   Transition Quadrapole Moments");
         printf("   %16s %16s %16s\n", "X", "Y", "Z");
         printf(" X %16.8f %16.8f %16.8f\n", quadrapoles(i,0,0), quadrapoles(i,0,1), quadrapoles(i,0,2));
         printf(" Y %16.8f %16.8f %16.8f\n", quadrapoles(i,1,0), quadrapoles(i,1,1), quadrapoles(i,1,2));
         printf(" Z %16.8f %16.8f %16.8f\n", quadrapoles(i,2,0), quadrapoles(i,2,1), quadrapoles(i,2,2));

         // Print contributions
         print("\n   Norms of the Components:");
         for(int j = 0; j < n; j++)
         {
            printf("   Occupied %d  --->  Virtual %d   %7.8f\n", j, i, norms(i,j));
         }

         print("\n");
      }
   }
}

// Diagonalizes the given functions
void TDHF::diagonalize_guess(World & world,
                            std::vector<std::vector<real_function_3d>> & f,
                            Tensor<double> & omega,
                            std::vector<real_function_3d> & orbitals,
                            std::vector<real_function_3d> & full_orbitals,
                            Tensor<double> & energies,
                            double thresh,
                            double small,
                            int print_level,
                            std::string xy)
{
   // Create gamma 
   std::vector<std::vector<real_function_3d>> gamma = create_gamma(world, f, orbitals, small, thresh, print_level);

   // Create \hat{V}^0 applied to guess functions 
   std::vector<std::vector<real_function_3d>> V_response = create_potential(world, f, print_level, xy);

   // Constructing S
   Tensor<double> S = create_overlap(world, f, print_level, xy);

   // Constructing hamiltonian
   Tensor<double> A = create_hamiltonian(world, gamma, V_response, f, orbitals, full_orbitals, energies, print_level, xy);

   // Solve Ax = Sxw  
   diag_fock_matrix(world, A, f, V_response, gamma, omega, S, thresh);
}

// Adds in random noise to a vector of vector of functions
std::vector<std::vector<real_function_3d>> TDHF::add_randomness(World & world,
                                                               std::vector<std::vector<real_function_3d>> & f)
{
   // Copy input functions
   std::vector<std::vector<real_function_3d>> f_copy = copy(world, f);

   // Lambda function to add in noise
   auto lambda = [](const Key<3> & key, Tensor<double> & x) mutable
   {
      Tensor<double> y(x.size());
      y.fillrandom();
      y.scale(1e2);
      x = x + y;
   };

   // Go through each function in f_copy and add in random noise
   for(unsigned int i = 0; i < f_copy.size(); i++)
   {
      for(unsigned int j = 0; j < f_copy[0].size(); j++)
      {
         // Add in random noise using rng and a the defined lambda function
         f_copy[i][j].unaryop(lambda);
      }
   }

   // Done
   return f_copy;
}

// Creates the ground state hamiltonian from given functions f
void TDHF::create_ground_hamiltonian(World & world,
                                    std::vector<real_function_3d> f,
                                    int print_level)

{
   // Basic output
   if(print_level > 0)
   {
      if(world.rank() == 0) print("   Creating the ground state hamiltonian.");
   }

   // Get sizes
   int m = f.size();

   // Calculate T
   // Make the derivative operators in each direction
   real_derivative_3d Dx(world, 0);
   real_derivative_3d Dy(world, 1);
   real_derivative_3d Dz(world, 2);

   // Apply derivatives once, and take inner products
   // according to this formula (faster / less noise):
   //  < f | \nabla^2 | f > = - < \nabla f | \nabla f >
   std::vector<real_function_3d> fx = apply(world, Dx, f);
   std::vector<real_function_3d> fy = apply(world, Dy, f);
   std::vector<real_function_3d> fz = apply(world, Dz, f);

   // Construct T according to above formula
   // Note: No negative as the formula above 
   // has one as well, so they cancel
   Tensor<double> T = 1.0/2.0 * (matrix_inner(world, fx, fx) +
                                 matrix_inner(world, fy, fy) +
                                 matrix_inner(world, fz, fz));

   // Construct V 
   // v_nuc first
   PotentialManager manager(Gparams.molecule, "a");
   manager.make_nuclear_potential(world);
   real_function_3d v_nuc = manager.vnuclear().truncate();

   // V_coul next
   // This does not include final multiplication of each orbital 
   // 2 is from integrating out spin
   real_function_3d v_coul = 2.0 * coulomb(world);

   // Sum coulomb (pre multiplied) and v_nuc
   // v_nuc comes out negative from potential manager, so add it
   real_function_3d v = v_coul + v_nuc;

   // Apply V to f functions
   std::vector<real_function_3d> vf = v * f;

   // exchange last
   // 'small memory' algorithm from SCF.cc 
   real_convolution_3d op = CoulombOperator(world, Rparams.small, Rparams.thresh);
   std::vector<real_function_3d> Kf = zero_functions_compressed<double,3>(world, m);
   for(int i=0; i<m; ++i)
   {
      std::vector<real_function_3d> psif = mul_sparse(world, f[i], f, Rparams.thresh);
      truncate(world, psif);
      psif = apply(world, op, psif);
      truncate(world, psif);

      // Save the potential here if we are saving it
      if(Rparams.store_potential)
      {
         stored_potential.push_back(psif);
      }

      psif = mul_sparse(world, f[i], psif, Rparams.thresh);
      gaxpy(world, 1.0, Kf, 1.0, psif);
   }

   // Construct V
   Tensor<double> V = matrix_inner(world, f, vf) - matrix_inner(world, f, Kf);

   // Now create the hamiltonian
   hamiltonian = T + V;

   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) print(hamiltonian);
   }

}

// Main function, makes sure everything happens in correcct order
void TDHF::solve(World & world)
{
   // Get start time
   Tensor<double> start_time = end_timer(world);

   // Plotting input orbitals
   if(Rparams.plot_initial)
   {
      if(world.rank() == 0) print("\n   Plotting ground state densities.\n");
      do_vtk_plots(world, 202, Gparams.L/2.0, 0, Gparams.num_orbitals, Gparams.molecule, square(world, Gparams.orbitals), "ground");
   }

   // Create initial guesses
   if(world.rank() == 0)
   {
      print("\n\n   TDHF Response Calculation");
      print("   ------------------------");
   }

   // Create the active subspace (select which ground state orbitals to calculate excitations from)
   //if(Rparams.e_window) select_active_subspace(world);

   // Create hamiltonian from ground state orbitals (need the matrix for both local and canonical orbitals)
   create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level);

   // Create trial functions by
   // creating a large number of symmetry included guesses 
   // or use random guesses
   if(world.rank() == 0)
   {
      print("   Creating trial functions.\n");
   }
   std::vector<std::vector<real_function_3d>> x_guesses, y_guesses;
   if(Rparams.random)
   {
      x_guesses = response_zero_functions(world, Rparams.states, Gparams.num_orbitals);
      x_guesses = add_randomness(world, x_guesses);
   }
   else
   {
      x_guesses = create_trial_functions(world, Rparams.states, Gparams.orbitals, Rparams.print_level);
   } 
  
   // Project out groundstate from guesses
   QProjector<double, 3> projector(world, Gparams.orbitals);
   for(unsigned int i = 0; i < x_guesses.size(); i++) x_guesses[i] = projector(x_guesses[i]);

   // Normalize
   normalize(world, x_guesses);

   // Only do this if using the constructed guess (not random)
   if(not Rparams.random)
   {
      // Basic output
      if(world.rank() == 0) print("\n   Diagonalizing trial functions for an improved initial guess.\n");

      Tensor<double> guess_x_omega(x_guesses.size());
      Tensor<double> guess_y_omega(y_guesses.size());

      // Diagonalize
      // Inplace modificaiton of guesses and guess_omega
      // Using the Tamm-Danchof approximation in this, should still be good enough for a first guess
      diagonalize_guess(world, x_guesses, guess_x_omega, Gparams.orbitals, Gparams.orbitals, hamiltonian, Rparams.thresh, Rparams.small, Rparams.print_level, "x");

      // Basic output
      if(Rparams.print_level >= 0)
      {
         if(world.rank() == 0)
         {
            print("   Initial response energies:");
            print(guess_x_omega);
         }
      }

      // Now we need to choose the Rparam.states lowest energy states
      x_response = select_trial_functions(world, x_guesses, guess_x_omega, Rparams.states, Rparams.print_level);
   }
   else
   {
      if(world.rank() == 0) print("   Using a random guess for initial response functions.");
      x_response = copy(world, x_guesses);
   }

   // Create y states as a copy of x states
   // Probably need to do something smarter?
   if( not Rparams.tda)
   {
      y_response = copy(world, x_response); 
   }

   // Ready to iterate! 
   // Protocol things go here?
   iterate(world);

   // Plot the response function if desired
   if(Rparams.plot)
   {
      // Need to sum contributions to get densities first
      std::vector<real_function_3d> densities = zero_functions<double, 3>(world, Rparams.plot_data.size());
      for(int i : Rparams.plot_data)
      {
         for(unsigned int j = 0; j < Gparams.num_orbitals; j++)
         {
            densities[i] = densities[i] + Gparams.orbitals[j] * x_response[i][j];

            if(not Rparams.tda)
            {
               densities[i] = densities[i] + Gparams.orbitals[j] * y_response[i][j];
            }
         }
      }

      // Now plot
      if(world.rank() == 0) print("\n   Plotting response state densities.\n");
      do_vtk_plots(world, 150, Gparams.L, 0, Rparams.plot_data.size(), Gparams.molecule, densities, "response-state");   
   }

   // Print total time
   // Precision is set to 10 coming in, drop it to 2
   std::cout.precision(2);
   std::cout << std::fixed;

   // Get start time
   Tensor<double> current_time = end_timer(world);
   current_time = current_time - start_time;

   //Tensor<double> current_time = end_timer(world);
   if(world.rank() == 0) print("   Total time:", current_time[0],"\n");
}

// Deuces



