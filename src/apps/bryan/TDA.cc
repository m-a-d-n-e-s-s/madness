#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <vector>
#include <math.h> 
#include <stdio.h> 
#include <iomanip>
#include <complex>
#include <cmath>
#include "../chem/molecule.h"
#include "../chem/potentialmanager.h"
#include "../chem/projector.h"          // For easy calculation of (1 - \hat{\rho}^0)
#include "../../examples/nonlinsol.h"   // For the KAIN solver
#include "TDA.h"
#include "TDA_Basic_Operators.h"
#include "Plot_VTK.h"                   // The plotting function I wrote

using namespace madness;

// KAIN allocator for vectorfunctions
struct TDA_allocator
{
   // Member variables
   World& world;
   const int num_vir;
   const int num_occ;
   
   // Constructor
   TDA_allocator(World& world, const int num_vir, const int num_occ) : world(world), num_vir(num_vir), num_occ(num_occ) {}
   
   // Overloading () operator
   std::vector<std::vector<real_function_3d>> operator()()
   {
      std::vector<std::vector<real_function_3d>> f;

      for(int i = 0; i < num_vir; i++) f.push_back(zero_functions<double,3>(world,num_occ));
   
      return f;
   }
   
   // Copy constructor
   TDA_allocator operator=(const TDA_allocator &other)
   {
      TDA_allocator tmp(world,other.num_occ,other.num_vir);
      return tmp;
   }
};

// Pulled from SCF.cc, starts a timer
void TDA::start_timer(World& world)
{
     world.gop.fence();
     tda_ttt.push_back(wall_time());
     tda_sss.push_back(cpu_time());
}

// Needed for timers
double TDA:: pop(std::vector<double>& v)
{
   double x = v.back();
   v.pop_back();
   return x;
}

// Stops a timer
Tensor<double> TDA::end_timer(World& world)
{
     Tensor<double> times(2);
     times[0] = wall_time() - pop(tda_ttt);
     times[1] = cpu_time() - pop(tda_sss);
     return times;
}


// Constructor
TDA::TDA(World & world,
         char* input_file, // Note: Leave off the .00000 for this to work. Not sure why.
         int k = 1,
         int print_level = 0)
{
      // Start the timers
      start_timer(world);

      // Dummy variable I don't care about
      double dummy1;
      std::vector<int> dummy2;

      // Everything else gets pulled out of the input archive
      archive::ParallelInputArchive input(world,input_file);
      input & dummy1;              // double
      input & tda_spinrestricted;  // bool
      input & tda_num_orbitals;    // int
      input & tda_ground_energies; // Tensor<double>    orbital energies
      input & tda_occ;             // Tensor<double>    orbital occupations
      input & dummy2;              // std::vector<int>  sets of orbitals(?)
      input & tda_L;               // double            box size
      input & tda_order;           // int               wavelet order
      input & tda_molecule;        // Molecule   

      // Check that order is positive and less than 30
      if (tda_order < 1 or tda_order > 30)
      {
         if(world.rank() == 0) print("\n ***PLEASE NOTE***\n   Invalid wavelet order read from archive, setting to 8.\n   This seems to happen when the default wavelet order is used in moldft."); 
         tda_order = 8;
      }

      // Derived parameters
      tda_thresh = pow(10, 2 - tda_order);
      tda_small = 1e-4; // Chosen arbitrarily

      // Set certain MADNESS parameters
      FunctionDefaults<3>::set_cubic_cell(-tda_L,tda_L);
      FunctionDefaults<3>::set_thresh(tda_thresh); 
      FunctionDefaults<3>::set_truncate_mode(1);
      FunctionDefaults<3>::set_k(tda_order);
 
      // Read in ground state orbitals
      for(int i = 0; i < tda_num_orbitals; i++)
      {
           real_function_3d reader;
           input & reader;
           tda_orbitals.push_back(reader); 
      }

      // User input things
      // Currently set at defaults
      tda_num_excited = k;
      tda_print_level = print_level;
      tda_max_iterations = 300;
      tda_energy_threshold = tda_thresh; 
}

// Prints norms of the given vector
void TDA::print_norms(World & world,
                      std::vector<std::vector<real_function_3d>> f)
{
   // Container
   Tensor<double> norms(f.size(),f[0].size());
   
   // Calc the norms
   for(int i = 0; i < f.size(); i++)
   {
      for(int j = 0; j < f[0].size(); j++)
      {
         norms(i,j) = f[i][j].norm2();
      }
   }

   // Print em in a smart way
   if(world.rank() == 0) print(norms);

}

// Small function to print MADNESS parameters
void TDA::print_madness_params(World &world)
{
   if(world.rank() == 0)
   {
      print("\n   MADNESS Calculation Parameters");
      print("   ------------------------------");
      print("\n            Box Width: ", 2.0 * tda_L); 
      print("        Wavelet Order: ", FunctionDefaults<3>::get_k());
      print(" Refinement Threshold: ", pow(10,-FunctionDefaults<3>::get_k()+2));
      print("Convergence Tolerance: ", FunctionDefaults<3>::get_thresh(), "\n");     
    }
}

// Small function to print geometry of a molecule nicely
void TDA::print_molecule(World &world) 
{
   if(world.rank() == 0)
   {
      // Precision is set to 10 coming in, drop it to 5
      std::cout.precision(5);
      std::cout << std::fixed;

      // First get atoms
      const std::vector<Atom> atoms = tda_molecule.get_atoms(); 
      int num_atoms = atoms.size();
 
      // Now print
      print("\n   Geometry Information");
      print("   --------------------\n");
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

// Small function to print response parameters
void TDA::print_response_params(World & world)
{
   if(world.rank() == 0)
   {
      print("\n   Response Information");
      print("   --------------------");
      print("           Method:", "Hartree-Fock");
      print("       Functional:", "None");
      print(" States Requested:", tda_num_excited);
   }
}


// Returns a vector of vectors filled with zero functions 
// with proper sizes (an "m x n" matrix of zero functions)
std::vector<std::vector<real_function_3d>> TDA::tda_zero_functions(World & world,
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


// Returns the given function f raised to power n
real_function_3d TDA::power_function(World & world,
                                     real_function_3d f,
                                     int n)
{
   MADNESS_ASSERT(n > 0);

   // Container to return
   real_function_3d result = copy(f);

   // Now multiply for each n-1 times needed
   for(int i = 0; i < n-1; i++) result = result * f;

   // Done
   return result;
}

// Returns a list of symmetry related functions for the correct
// pointgroup of the provided molecule
std::vector<real_function_3d> TDA::symmetry(World & world,
                                            Molecule molecule)
{
   // Container to return
   std::vector<real_function_3d> result;

   // Ensure the pointgroup is stored
   // This prints the pointgroup
   molecule.identify_point_group();

   // Get the pointgroup
   std::string pg = molecule.pointgroup_;

   // Create the basic x, y, and z functions
   std::vector<int> f(3,0); 
   f[0] = 1;
   real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
   f[0] = 0; f[1] = 1;
   real_function_3d y = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
   f[1] = 0; f[2] = 1;
   real_function_3d z = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));

   // Based on pointgroup, append functions to result
   // Functions pulled from character table of Cotton's "Chemical
   // Applications of Group Theory", 3rd Edition
   if(pg == "Ci" or pg == "C2h" or pg == "Cs" or pg == "C2" or pg == "C2v" or pg == "D2")
   {
      result.push_back(x*x);
      result.push_back(y*y);
      result.push_back(z*z);
      result.push_back(x*y);
      result.push_back(x*z);
      result.push_back(y*z);
   }
   else if(pg == "D2h")
   {
      result.push_back(x*x-y*y);
      result.push_back(z*z);
      result.push_back(x*z);
      result.push_back(y*z);
      result.push_back(x*x-y*y);
      result.push_back(x*y);
   }
   else // Assuming C1 if none of the above
   {
      // No symmetry relations in C1   
   }

   // Done
   return result;
}

// Base case for generating chebyshev guesses
std::vector<real_function_3d> TDA::chebyshev(World & world,
                                             int k)
{
   // Basic output
   if(tda_print_level >= 0)
   {
      if(world.rank() == 0) print("   Creating Chebyshev polynomials.");
   }

   // Create an empty vector of functions to return
   std::vector<real_function_3d> T;

   // A zero function
   real_function_3d zero = real_factory_3d(world);

   // If k >= 2, add in first 2
   if(k >= 1) T.push_back(zero.add_scalar(1.0));
   if(k >= 2) T.push_back(real_factory_3d(world).functor(GaussianGuess<3>(Vector<double,3>(0.0), 0.1, std::vector<int>{1,1,1} )));
 
   while(T.size() < k)
   {
      // Add in next polynomial
      // Formula is T_{n+1} = x * T_n - T_{n-1}
      // Remember: T[1] = x, so use that to get
      // T_{n+1} = T[1] * T_n - T_{n-1}
      // Also note: should only get here if T[1]
      // already exists
      T.push_back(T[1] * T[T.size() - 1] - T[T.size() - 2]);
   }

   // Done
   return T;
}

// Returns initial guess functions
// Probably needs a lot of work on how to do this intelligently
std::vector<std::vector<real_function_3d>> TDA::create_trial_functions1(World & world,
                                                                        int k,
                                                                        std::vector<real_function_3d> orbitals,
                                                                        int axis,
                                                                        int print_level)
{
   // Number of ground state orbitals
   int n = orbitals.size();

   // Holder for r*orbitals
   std::vector<std::vector<real_function_3d>> trials;

   // Tells functor which axis to use
   std::vector<int> f(3,0);
   f[axis] = true;

   // Creating the functor  
   real_function_3d r = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));

   // Create k functions
   for(int i = 0; i < k; i++)
   {
      trials.push_back(power_function(world, r, i+1) * tda_orbitals);      
   }

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("Norms of guess functions:");
      print_norms(world, trials); 
   }

   // Truncate 
   truncate(world, trials);

   // Done?
   return trials;    
}

// Returns initial guess functions
// Probably needs a lot of work on how to do this intelligently
// Currently: Producing each combination of symmetry function on the orbitals
//            times two (produces two orbitals for each symmetry function * ground
//            state orbital)
std::vector<std::vector<real_function_3d>> TDA::create_trial_functions2(World & world,
                                                                        int k,
                                                                        std::vector<real_function_3d> orbitals,
                                                                        int print_level)
{
   // Get size
   int n = orbitals.size();

   // Container to return 
   std::vector<std::vector<real_function_3d>> trials;

   // Create a vector of correct symmetry relate polynomials
   std::vector<real_function_3d> symm = symmetry(world, tda_molecule);
   int sym_size = symm.size();

   // Default case
   if(sym_size > 0 && sym_size * 2 > k)
   {
      if(tda_print_level >= 0)
      {
         if(world.rank() == 0) print("   Found", sym_size, " symmetry functions for this molecule.");
         if(world.rank() == 0) print("   Constructing 2 *", sym_size, " guesses.");
      }

      // Run over symmetry functions
      for(int i = 0; i < sym_size; i++)
      {
         std::vector<real_function_3d> temp(n);

         // Run over each occupied orbital
         for(int p = 0; p < n; p++)
         {
            temp[p] = symm[i] * orbitals[p];
         } 

         // Run over ground state functions
         trials.push_back(temp);      
      } 
   }
   else if(sym_size == 0) // C1 symmetry case
   {
      if(tda_print_level >= 0)
      {
         if(world.rank() == 0) print("   Found", sym_size, " symmetry functions for this molecule.");
         if(world.rank() == 0) print("   Constructing 2 *", k, " guesses.");
      }
   }
   else // Case when symmetry doesn't provide enough guesses
   {
      if(tda_print_level >= 0)
      {
         if(world.rank() == 0) print("   Found", sym_size, " symmetry functions for this molecule,");
         if(world.rank() == 0) print("   but", k, " functions were requested. Creating additional");
         if(world.rank() == 0) print("   functions.");
         if(world.rank() == 0) print("   Constructing 2 *", k, " guesses.");
      }
   } 

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("Norms of guess functions:");
      print_norms(world, trials); 
   }

   // Truncate 
   truncate(world, trials);

   // Done?
   return trials;    
}

// Constructing guesses from Chebyshev polynomials
std::vector<std::vector<real_function_3d>> TDA::create_trial_functions3(World & world,
                                                                        int k,
                                                                        std::vector<real_function_3d> orbitals,
                                                                        int print_level)
{
   // Get size from input orbitals
   int n = orbitals.size();
   
   // Container to return 
   std::vector<std::vector<real_function_3d>> trials;

   // Create a vector of the first k chebyshev polynomials   
   // The +1 is so we don't use the first chebyshev polynomial,
   // which is just 1 
   std::vector<real_function_3d> cheby = chebyshev(world, k+1); 

   // Make an outer product between chebyshev and orbitals 
   for(int i = 1; i < cheby.size(); i++)
   {
      // Container for the row
      std::vector<real_function_3d> temp(n);

      for(int p = 0; p < n; p++)
      {
         // Create the row
         temp[p] = cheby[i] * orbitals[p];
      }

      // Save the row
      trials.push_back(temp);
   }

   // Normalize the guesses
   normalize(world, trials);

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("Norms of guess functions:");
      print_norms(world, trials); 
   }

   // Truncate 
   truncate(world, trials);

   // Done?
   return trials;    
}

// Returns the derivative of the coulomb operator, applied to ground state orbitals
std::vector<std::vector<real_function_3d>> TDA::create_coulomb_derivative(World &world,
                                                                          int m,
                                                                          int n,
                                                                          std::vector<std::vector<real_function_3d>> f,
                                                                          std::vector<real_function_3d> orbitals,
                                                                          double small,
                                                                          double thresh)
{
   // Zero function, to be returned
   std::vector<std::vector<real_function_3d>> deriv_j = tda_zero_functions(world, m, n);

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Temperary storage
   real_function_3d rho = real_function_3d(world);

   // Need to run over virtual orbitals
   for(int k = 0; k < m; k++)
   {
      // Need to run over all occupied orbitals
      for(int p = 0; p < n; p++)
      {
         // Get transitionn density 
         rho = dot(world, f[k], orbitals);

         // Apply coulomb operator
         rho = apply(op, rho);
         
         // Multiply by ground state orbital p
         rho = rho * orbitals[p];

         // Save the result
         deriv_j[k][p] = rho;
      }
   }

   // Done
   return deriv_j;
}

// Returns the derivative of the coulomb operator, applied to ground state orbitals
// Overloaded to get defaults correct
std::vector<std::vector<real_function_3d>> TDA::create_coulomb_derivative(World &world)
{
   return create_coulomb_derivative(world, tda_num_excited, tda_num_orbitals, tda_x_response, tda_orbitals, tda_small, tda_thresh);
}


// Does what it sounds like it does
std::vector<std::vector<real_function_3d>> TDA::create_exchange_derivative(World &world,
                                                                           int m,
                                                                           int n,
                                                                           std::vector<std::vector<real_function_3d>> f,
                                                                           std::vector<real_function_3d> orbitals,
                                                                           double small,
                                                                           double thresh)
{
   // Zero function, to be returned
   std::vector<std::vector<real_function_3d>> deriv_k = tda_zero_functions(world, m, n); 

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);
 
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
            rho = rho * f[k][i];

            // Add to total
            deriv_k[k][p] += rho;
         }
      }
   }

   // Done
   return deriv_k;
}

// Does what it sounds like it does
// Overloaded to get default values right
std::vector<std::vector<real_function_3d>> TDA::create_exchange_derivative(World &world)
{
  return create_exchange_derivative(world, tda_num_excited, tda_num_orbitals, tda_x_response, tda_orbitals, tda_small, tda_thresh); 
}

// Computes gamma(r) given the ground state orbitals and response functions
std::vector<std::vector<real_function_3d>> TDA::create_gamma(World &world,
                                                             int m,
                                                             int n,
                                                             std::vector<std::vector<real_function_3d>> f,
                                                             std::vector<real_function_3d> orbitals,
                                                             double small,
                                                             double thresh,
                                                             int print_level)
{
   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) print("   Creating Gamma");
   }

   // The gamma function to be returned, intialized to zero
   std::vector<std::vector<real_function_3d>> gamma = tda_zero_functions(world, m, n); 

   // Gamma will have 2 terms for HF: dJ/drho[rho] and dK/drho[rho]
   // There is a different Gamma for each orbital-->virtual transition
   // Calculate both here
   std::vector<std::vector<real_function_3d>> deriv_J = create_coulomb_derivative(world, m, n, f, orbitals, small, thresh);
   std::vector<std::vector<real_function_3d>> deriv_K = create_exchange_derivative(world, m, n, f, orbitals, small, thresh);

   // Debugging output
   if (print_level >= 2)
   {
      if(world.rank() == 0) print("Norms of coulomb deriv:");
      print_norms(world, scale(deriv_J, 2.0)); 
      if(world.rank() == 0) print("Norms of exchange deriv:");
      print_norms(world, deriv_K);
   }

   // Spin integration gives coefficients 
   // This is the spin restricted, singlet excitation coefficients
   gamma = scale(deriv_J, 2.0) - deriv_K; 
   
   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   Gamma function norms:");
      print_norms(world, gamma);
   }

   // Done
   world.gop.fence();
   return gamma;
}

// Computes gamma(r) given the ground state orbitals and response functions
// Overloaded to get defaults right
std::vector<std::vector<real_function_3d>> TDA::create_gamma(World &world)
{
   return create_gamma(world, tda_num_excited, tda_num_orbitals, tda_x_response, tda_orbitals, tda_small, tda_thresh, tda_print_level);
}

// Calculates ground state coulomb potential 
real_function_3d TDA::coulomb(World& world)
{
   // Coulomb operator
   real_convolution_3d op = CoulombOperator(world, tda_small, tda_thresh);

   // Create the density
   real_function_3d rho = dot(world, tda_orbitals, tda_orbitals); 
  
   // Apply operator to density
   rho = apply(op,rho);
   rho.truncate();

   // Done 
   return rho;
}

// Calculates HF exchange between ground state orbitals and response orbitals 
std::vector<std::vector<real_function_3d>> TDA::exchange(World& world,
                                                         int m,
                                                         int n,
                                                         std::vector<std::vector<real_function_3d>> f)
{
     // Coulomb operator
     real_convolution_3d op = CoulombOperator(world, tda_small, tda_thresh);

     // Container for results and others
     std::vector<std::vector<real_function_3d>> result = tda_zero_functions(world, m, n);
     real_function_3d psif = real_function_3d(world);

     // Run over each excited state
     for(int k = 0; k < m; k++)
     {
        // And run over each occupied state
        for(int p = 0; p < n; p++)
        {
           // Get transition density
           psif = tda_orbitals[p] * f[k][p];
           
           // Apply coulomb operator 
           psif = apply(op, psif);

           // Final multiplication
           result[k][p] = tda_orbitals[p] * psif;
        }
     }      
   
     // Done!
     return result;
}


// Calculates HF exchange between ground state orbitals and response orbitals 
// Overloaded to get defaults right
std::vector<std::vector<real_function_3d>> TDA::exchange(World& world)
{
   return exchange(world, tda_num_excited, tda_num_orbitals, tda_x_response);
}

// Returns the ground state potential applied to response functions
std::vector<std::vector<real_function_3d>> TDA::create_potential(World & world,
                                                                 int m,
                                                                 int n,
                                                                 std::vector<std::vector<real_function_3d>> f,
                                                                 int print_level)
{
   // Basic output
   if(print_level >= 1) 
   {
      if(world.rank() == 0) print("   Computing V0 * x(r)");
   }   

   // Computing \hat{V}^0 = v_nuc + v_coul + v_exch           
   // v_nuc first
   PotentialManager manager(tda_molecule, "a");
   manager.make_nuclear_potential(world);
   real_function_3d v_nuc = manager.vnuclear().truncate(); 
   
   // V_coul next
   // This does not include final multiplication of each orbital 
   real_function_3d v_coul = coulomb(world); 
   
   // Sum coulomb (pre multiplied) and v_nuc
   // 2 is from integrating out spin
   // v_nuc comes out negative from potential manager, so add it
   real_function_3d v = 2.0 * v_coul + v_nuc;
   
   // Apply V to response functions
   std::vector<std::vector<real_function_3d>> V_x_resp = multiply(f, v); 
  
   // V_exch last
   // Multiplication by x_response functions is included in construction
   std::vector<std::vector<real_function_3d>> v_exch = exchange(world, m, n, f);
   
   if(print_level >= 2) 
   {
      // Print potential energy matrices
      if(world.rank() == 0) print("Nuclear:");
      print_norms(world, multiply(f, multiply(f, v_nuc)));
      if(world.rank() == 0) print("Coulomb:");
      print_norms(world, scale(multiply(f, multiply(f, v_coul)), 2.0));
      if(world.rank() == 0) print("Exchange:");
      print_norms(world, multiply(f, v_exch));
   }   

   // Subtract V_exch and V_x_resp
   V_x_resp = V_x_resp - v_exch;
   
   // Print some basic output
   if(print_level >= 1) 
   {
      if(world.rank() == 0) print("   V*x(r) norms:");
      print_norms(world, V_x_resp);
   }

   // Done?
   return V_x_resp;
}

// Returns the ground state potential applied to response functions
// Overloaded for default parameters
std::vector<std::vector<real_function_3d>> TDA::create_potential(World & world)
{
   return create_potential(world, tda_num_excited, tda_num_orbitals, tda_x_response, tda_print_level);
}

// Returns a Tensor of inner products, where
// result(i,j) = inner(a[i][j],b[i][j])
Tensor<double> TDA::m_m_inner(World &world,
                              std::vector<std::vector<real_function_3d>> a,
                              std::vector<std::vector<real_function_3d>> b)
{
   MADNESS_ASSERT(a.size() > 0);
   MADNESS_ASSERT(a.size() == b.size());
   MADNESS_ASSERT(b[0].size() > 0);

   // Get sizes
   int dim_a = a.size();
   int dim_b = b[0].size();

   // Container for result
   Tensor<double> result(dim_a, dim_b);

   // Run over dimension one
   for(int p = 0; p < dim_a; p++)
   {
      // Run over dimension two
      for(int k = 0; k < dim_b; k++)
      {
         result(p,k) = inner(a[p][k],b[p][k]);
      }
   }

   // Done
   return result;
}

// Creating matrix S for first guess at omega
Tensor<double> TDA::create_S(World & world,
                             int m,
                             int n,
                             std::vector<std::vector<real_function_3d>> f,
                             int print_level)
{
   // Container for answer to be returned
   Tensor<double> S(m,m);
   
   // Temp variable for summing
   double temp;

   // Run over all virtuals i 
   for(int i = 0; i < m; i++)
   {
      // Run over all virtuals j
      for(int j = 0; j < m; j++)
      {        
         // Using vmra.h line 627 function 
         S(i,j) = inner(f[i], f[j]); 
      }
   }

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   S matrix:");
      print(S);
   }
 
   // Done
   return S;
}

// Creating matrix S for first guess at omega
// Overloaded for default parameters
Tensor<double> TDA::create_S(World & world)
{
   return create_S(world, tda_num_excited, tda_num_orbitals, tda_x_response, tda_print_level);
}

// Returns the ground state fock operator applied to response functions
std::vector<std::vector<real_function_3d>> TDA::create_fock(World & world,
                                                            std::vector<std::vector<real_function_3d>> V,
                                                            std::vector<std::vector<real_function_3d>> f,
                                                            int print_level)
{
   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   Preparing to create perturbed fock matrix");
   }

   // Size of fock matrix must match that of V
   int m = V.size();
   int n = V[0].size();

   // Container to return
   std::vector<std::vector<real_function_3d>> fock =  tda_zero_functions(world, m, n);

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

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   Norms of d^2/dx (response functions):");
      print_norms(world, dvx2);
      if(world.rank() == 0) print("   Norms of d^2/dy (response functions):");
      print_norms(world, dvy2);
      if(world.rank() == 0) print("   Norms of d^2/dz (response functions):");
      print_norms(world, dvz2);
   }

   // Add together derivatives
   fock = scale((dvx2 + dvy2 + dvz2), -0.5);

   // Debugging output   
   if(tda_print_level >= 2) 
   {
      if(world.rank() == 0) print("   Kinetic matrix:");
      print(m_m_inner(world, f, fock));
      if(world.rank() == 0) print("   Potential matrix:"); 
      print(m_m_inner(world, f, V));
   }

   // Add in potential
   fock = fock + V; 

   // Done
   return fock;
}


// Creating matrix A for first guess at omega
Tensor<double> TDA::create_A(World & world,
                             std::vector<std::vector<real_function_3d>> gamma,
                             std::vector<std::vector<real_function_3d>> V,
                             std::vector<std::vector<real_function_3d>> f,
                             std::vector<real_function_3d> ground_orbitals,
                             Tensor<double> energies,
                             int print_level)
{
   // Sizes inferred from V and gamma
   int n = V[0].size();
   int m = V.size();

   // Container for A
   Tensor<double> A(m,m); 

   // Projector on the unperturbed density
   QProjector<double,3> ground_state_density(world, ground_orbitals);
   
   // Create the ground-state fock operator on response orbitals
   std::vector<std::vector<real_function_3d>> fock_x_resp = create_fock(world, V, f, print_level);

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   Perturbed Fock:");
      print(m_m_inner(world, f, fock_x_resp));
   }

   // Need to calculate energy * x_response
   std::vector<std::vector<real_function_3d>> energy_x_resp = scale(f, energies); 

   // Construct intermediary
   std::vector<std::vector<real_function_3d>> temp = gamma + fock_x_resp - energy_x_resp;

   // Need to run over excited states
   for(int k = 0; k < m; k++)
   {
      // Need another run over excited states
      for(int j = 0; j < m; j++)
      {
         // Project out all overlap with ground state orbitals
         temp[j] = ground_state_density(temp[j]);

         // Run over all occupied orbitals to get their contribution
         // to the part of A we're calculating 
         for(int p = 0; p < n; p++)
         {            
            // Finally calculate \int dr x_p^k * temp
            A(j,k) += inner(f[k][p], temp[j][p]);
         }
      }
   }

   // Debugging output
   if(print_level >= 2) 
   {
      if(world.rank() == 0) print("   A matrix:");
      if(world.rank() == 0) print(A);
   }

   // Done
   return A;
}

// Overloaded to get default parameters right
Tensor<double> TDA::create_A(World & world,
                             std::vector<std::vector<real_function_3d>> gamma,
                             std::vector<std::vector<real_function_3d>> V)
{
   return create_A(world, gamma, V, tda_x_response, tda_orbitals, tda_ground_energies, tda_print_level);
}


// Returns the shift needed to make sure that
// -2.0 * (ground_state_energy + excited_state_energy) 
// is negative. Please note: The same shift needs to 
// be applied to the potential.
Tensor<double> TDA::create_shift(World & world,
                                 int m,
                                 int n,
                                 Tensor<double> ground,
                                 Tensor<double> omega,
                                 int print_level)
{
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

            if(print_level >= 1) 
            {
               // Let user know whats up
               if(world.rank() == 0) print("   ground energy =", ground(p));
               if(world.rank() == 0) print("   excited energy =", omega(k));
               if(world.rank() == 0) print("   Potential shift needed. Shifting by", result(k,p));
               if(world.rank() == 0) print("");
            }
         }
      }
   }

   // Done
   return result;
}

// Overloaded for default parameters
Tensor<double> TDA::create_shift(World & world)
{
   return create_shift(world, tda_num_excited, tda_num_orbitals, tda_ground_energies, tda_omega, tda_print_level);
}

// Returns the given shift applied to the given potential
std::vector<std::vector<real_function_3d>> TDA::apply_shift(World & world,
                                                            Tensor<double> shifts,
                                                            std::vector<std::vector<real_function_3d>> V,
                                                            std::vector<std::vector<real_function_3d>> f)
{
   // Sizes inferred from V
   int n = V[0].size();
   int m = V.size();

   // Container to return
   std::vector<std::vector<real_function_3d>> shifted_V = tda_zero_functions(world, m, n); 

   // Run over occupied
   for(int k = 0; k < m; k++)
   {
      // Run over virtual
      for(int p = 0; p < n; p++)
      {
         shifted_V[k][p] =  V[k][p] + shifts(k,p) * f[k][p];
      }
   }

   // Done
   return shifted_V;
}

// Returns the given shift applied to the given potential
std::vector<std::vector<real_function_3d>> TDA::apply_shift(World & world,
                                                            Tensor<double> shifts,
                                                            std::vector<std::vector<real_function_3d>> V)
{
   return apply_shift(world, shifts, V, tda_x_response);
}

// Function to make a vector of BSH operators using ground and excited
// state energies
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> TDA::make_bsh_operators(World & world,
                                                                                        Tensor<double> shift,
                                                                                        Tensor<double> ground,
                                                                                        Tensor<double> omega,
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

// Overloaded for function defaults 
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> TDA::make_bsh_operators(World & world,
                                                                                       Tensor<double> shift)
{
   return make_bsh_operators(world, shift, tda_ground_energies, tda_omega, tda_small, tda_thresh);
}


// Returns the second order update to the energies of the excited states
Tensor<double> TDA::calculate_energy_update(World & world, 
                                            std::vector<std::vector<real_function_3d>> gamma,
                                            std::vector<std::vector<real_function_3d>> f_residuals,
                                            std::vector<std::vector<real_function_3d>> new_f,
                                            int print_level) 
{
   /*
    *  The correction is:
    *      \delta \omega^{(k)} = - \frac{ \sum_p\left< \hat{V}^0 x_p^{(k)}(r) + (1 - \hat{\rho}^0) \Gamma_p^{(k)}(r)\right|
    *                                         \left. x_p^{(k)} - \~{x}_p^{(k)} \right> }
    *                                   { \sum_p \left| \left| \~{x}_p^{(k)} \right| \right|^2 }
    */
  
   if(print_level >= 1) 
   { 
      // Basic output
      if(world.rank() == 0) print("   Updating energies");
   }   

   // Size inferred
   int m = gamma.size();

   // Container for updates
   Tensor<double> updates(m);
   double normalization_constant;
     
   // Need to run over all functions in gamma and calculate innter products
   // Run over excited states
   for(int k = 0; k < m; k++)
   {
      // Zero the norm constant
      normalization_constant = 0.0;

      // vmra.h function, line 627
      updates(k) = inner(gamma[k],f_residuals[k]);
     
      // Normalize update function
      // The -1.0 is the leading coefficient in the update formula
      updates(k) = -1.0 * updates(k) / pow(norm2(world, new_f[k]), 2);
   }
   
   if(print_level >= 1) 
   {
      // Print energy deltas
      if(world.rank() == 0) print("   Energy changes:");
      if(world.rank() == 0) print(updates);
   }

   // Done?
   return updates;
}


// Specialized for response calculations that returns orthonormalized
// functions
std::vector<std::vector<real_function_3d>> TDA::gram_schmidt(World & world,
                                                              std::vector<std::vector<real_function_3d>> f)
{
   // Sizes inferred
   int n = f[0].size();
   int m = f.size();

   // Return container
   std::vector<std::vector<real_function_3d>> result(m);
   for(int i = 0; i < m; i++) result[i] = copy(world, f[i]);

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

   // Done
   return result;
}

// Returns the max norm of the given vector of functions
double TDA::calculate_max_residual(World & world,
                                   std::vector<std::vector<real_function_3d>> f)
{
   // Container for max
   double max = 0.0;

   // Run over all functions in f
   for(int i = 0; i < f.size(); i++)
   {
      for(int j = 0; j < f[0].size(); j++)
      {
         double temp = f[i][j].norm2();
         if( temp > max) max = temp;
      }
   }

   // Done
   return max;
}

// Applies step restriction if neccessary
std::vector<std::vector<real_function_3d>> TDA::step_restriction(World & world,
                                                                std::vector<std::vector<real_function_3d>> old_f,
                                                                std::vector<std::vector<real_function_3d>> new_f)
{
   // Determine if step restriction is necessary first
   // Calc max norm of difference between old and new
   std::vector<std::vector<real_function_3d>> diff = new_f - old_f;
   Tensor<double> norms(diff.size(), diff[0].size());

   // Container to return
   std::vector<std::vector<real_function_3d>> result(diff.size());
   for(int i = 0; i < diff.size(); i++) result[i] = zero_functions<double,3>(world,diff[0].size());

   // Calculate norms
   for(int i = 0; i < diff.size(); i++) 
   {
      for(int j = 0; j < diff[0].size(); j++)
      {
         norms(i,j) = diff[i][j].norm2();
      }
   }

   // Check all to see if step is neccessary
   for(int i = 0; i < result.size(); i++)
   { 
      for(int j = 0; j < result[0].size(); j++)
      {
         if(norms(i,j) > 0.5) // Arbitrary choice here
         {
            // Step restriction necessary
            // Take average
            result[i][j] = (new_f[i][j] + old_f[i][j]) * 0.5;
         }
         else
         {
            result[i][j] = copy(new_f[i][j]);
         }
      }
   }

   // Done
   return result;
} 
 
// Selects from a list of functions and energies the k functions with the lowest energy
std::vector<std::vector<real_function_3d>> TDA::select_trial_functions(World & world,
                                                                       std::vector<std::vector<real_function_3d>> f,
                                                                       Tensor<double> energies,
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

   // Container to hold the ordering of energies
   Tensor<int> ordering(k);   

   // Simple sorting thing
   // Copy first tda_num_excited 
   for(int i = 0; i < k; i++)
   {
      ordering(i) = i;
   }

   // Now sort 
   // Run over first bloc of requested states
   for(int i = 0; i < k; i++)
   {
      // Compare that energy to rest of bloc
      // to sort
      for(int j = i + 1; j < k; j++)
      {
         // Check if switch is needed
         if(energies(j) < energies(i))
         {
            // Switch 
            int temp = ordering(i);
            ordering(i) = ordering(j);
            ordering(j) = temp;
         }
      }
   }

   // Now go through the rest of guesses and put them in the right 
   // spot if lower in energy
   for(int i = k; i < energies.size(); i++)
   {
      // Run over the guesses already selected
      for(int j = 0; j < k; j++)
      {
         // See if new guess should be kept
         if(energies(i) < energies(j))
         {
            // It should be kept. Need to move all 
            // guesses after this one up one and
            // remove last one
            int temp = ordering(j);
            ordering(j) = i;             

            for(int a = j+1; a < k; a++)
            {
               int temp2 = ordering(a);
               ordering(a) = temp;
               temp = temp2; 
            }
         }
        
         // Need to break the j indexed for loop
         break;
      }
   }

   if(print_level >= 2)
   {
      if(world.rank() == 0) print("   These functions were chosen:");
      if(world.rank() == 0) print(ordering);
      if(world.rank() == 0) print("   And they have energies (use the above indices in the following array):");
      if(world.rank() == 0) print(energies);
   }

   // Now pull out the states listed in ordering and put into tda_x_response
   for(int i = 0; i < ordering.size(); i++) answer.push_back(f[ordering(i)]);

   if(print_level >= 2)
   {
      if(world.rank() == 0) print("   Initial guess norms:");
      print_norms(world, tda_x_response); 
   }
   
   // Done
   return answer;
}

// Calculate the exponentiation of a matrix through first order (I think)
Tensor<double> TDA::matrix_exponential(const Tensor<double> & A)
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
    while (n--) {
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
Tensor<double> TDA::get_fock_transformation(World & world, 
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
    while (switched) {
        switched = false;
        for (int i = 0; i < nmo; i++) {
            for (int j = i + 1; j < nmo; j++) {
               
                    double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                    double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                    if (snew > sold) {
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
           
            // Iteratively construct unitary rotation by
            // exponentiating the antisymmetric part of the matrix
            // ... is quadratically convergent so just do 3
            // iterations
            Tensor<double> rot = matrix_exponential(-0.5 * (q - transpose(q)));
            q = inner(q, rot);
            Tensor<double> rot2 = matrix_exponential(-0.5 * (q - transpose(q)));
            q = inner(q, rot2);
            Tensor<double> rot3 = matrix_exponential(-0.5 * (q - transpose(q)));
            q = inner(rot, inner(rot2, rot3));
            U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
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
Tensor<double> TDA::diag_fock_matrix(World & world,
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
   
    // transform the orbitals and the orbitals times the potential
    Vpsi = transform(world, Vpsi, U);
    gamma = transform(world, gamma, U);
    psi = transform(world, psi, U);

    // truncate
    truncate(world, Vpsi);
    truncate(world, gamma);
    truncate(world, psi);
    normalize(world, psi);

    return U;
}


// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
std::vector<std::vector<real_function_3d>> TDA::transform(World & world,
                                                          std::vector<std::vector<real_function_3d>> f,
                                                          Tensor<double> U)
{
   // Return container
   std::vector<std::vector<real_function_3d>> result;

   // Go element by element 
   for(int i = 0; i < f.size(); i++)
   {
      // Temp for the result of one row
      std::vector<real_function_3d> temp = zero_functions_compressed<double,3>(world, f[0].size());

      for(int j = 0; j < f.size(); j++)
      {
         gaxpy(world, 1.0, temp, U(j,i), f[j]);
      }
      
      // Add to temp to result
      result.push_back(temp);
   }

   // Done?
   return result;
} 

// Sorts the given Tensor of energies and vector of functions
// in place
void TDA::sort(World & world,
               Tensor<double> vals,
               std::vector<std::vector<real_function_3d>> f)
{
   // Get relevant sizes
   int k = vals.size();

   // Copy f and vals 
   std::vector<std::vector<real_function_3d>> f_copy;
   for(int i = 0; i < k; i++) f_copy.push_back(f[i]);
   Tensor<double> vals_copy = copy(vals);

   // Container to hold the ordering of energies
   Tensor<int> ordering(k);   

   // Initial order
   for(int i = 0; i < k; i++) ordering(i) = i;

   // Now sort 
   for(int i = 0; i < k; i++)
   {
      // Compare that energy to rest of bloc
      // to sort
      for(int j = i + 1; j < k; j++)
      {
         // Check if switch is needed
         if(vals(j) < vals(i))
         {
            // Switch 
            int temp = ordering(i);
            ordering(i) = ordering(j);
            ordering(j) = temp;
         }
      }
   }

   // Now pull out the states listed in ordering and put into tda_x_response
   for(int i = 0; i < k; i++)
   {
      f[i] = f_copy[ordering(i)];
      vals[i] = vals_copy[ordering(i)];
   }
}

// Prints headers for standard output of iterate()
void TDA::print_iterate_headers(World & world)
{
   if(world.rank() == 0)
   {
      print("");
      print("                 Max            Max   ");
      print("                Energy        Function");
      print(" Iteration     Residual       Residual        Time");
      print("---------------------------------------------------");
   }
}

// Prints headers for standard output of iterate()
void TDA::print_guess_iterate_headers(World & world)
{
   if(world.rank() == 0)
   {
      print("");
      print("                 Max            Max   ");
      print("   Guess        Energy        Function");
      print(" Iteration     Residual       Residual        Time");
      print("---------------------------------------------------");
   }
}

// Iterates the response functions until converged or out of iterations
void TDA::iterate(World & world)
{
   // Variables needed to iterate
   int iteration = 0;                                                            // Iteration counter
   QProjector<double, 3> projector(world, tda_orbitals);                         // Projector to project out ground state
   int n = tda_num_orbitals;                                                     // Number of ground state orbitals
   int m = tda_num_excited;                                                      // Number of excited states
   bool converged = false;                                                       // For convergence
   Tensor<double> energy_updates;                                                // Holds energy corrections
   std::vector<std::vector<real_function_3d>> differences;                       // Holds wave function corrections
   std::vector<std::vector<real_function_3d>> temp;                              // Used for step restriction 
 
   // The KAIN solver
   XNonlinearSolver<std::vector<std::vector<real_function_3d>>, double, TDA_allocator> kain(TDA_allocator(world, m, n), tda_print_level);  
   bool kain_on = true; // controls when to turn KAIN on

   // Setting max sub size for KAIN solver
   kain.set_maxsub(5);

   // Create gamma 
   //std::vector<std::vector<real_function_3d>> gamma = create_gamma(world);

   // Create \hat{V}^0 applied to tda_x_response
   //std::vector<std::vector<real_function_3d>> V_response = create_potential(world); 

   // Print headers
   if(tda_print_level == 0) print_iterate_headers(world);

   // Get a start time
   Tensor<double> initial_time = end_timer(world);

   // Now to iterate
   while( iteration < tda_max_iterations  && !converged)
   {
      // Basic output
      if(tda_print_level >= 1)
      {
         if(world.rank() == 0) print("\n   Iteration", iteration);
         if(world.rank() == 0) print("  --------------");
      }

      // Create gamma 
      std::vector<std::vector<real_function_3d>> gamma = create_gamma(world);

      // Project out ground state
      for(int i = 0; i < m; i++) gamma[i] = projector(gamma[i]);

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms of gamma after projecting out ground state:");
         print_norms(world, gamma);
      }

      // Create \hat{V}^0 applied to tda_x_response
      std::vector<std::vector<real_function_3d>> V_response = create_potential(world);

      // Basic output
      if(tda_print_level >= 1) 
      {
         if(world.rank() == 0) print("   Solving Ax=Swx for initial energies");
      }
      
      if(iteration == 0)
      {
         // Constructing S
         Tensor<double> S = create_S(world);
         
         // Constructing A
         Tensor<double> A = create_A(world, gamma, V_response);

         // Symmeterize A
         A = 0.5 * (A + transpose(A));

         // Solve Ax = Sxw
         Tensor<double> e_vecs;
         sygv(A, S, 1, e_vecs, tda_omega);
         //diag_fock_matrix(world, A, tda_x_response, V_response, gamma, tda_omega, S, tda_thresh);             
      }

      // Basic output
      if(tda_print_level >= 1) 
      {
         if(world.rank() == 0) print("   Response Orbital Energies:");
         if(world.rank() == 0) print(tda_omega);
      }

      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      Tensor<double> shifts = create_shift(world);  

      // Apply the shifts
      V_response = apply_shift(world, shifts, V_response);

      // Construct RHS of equation
      std::vector<std::vector<real_function_3d>> rhs = gamma + V_response;

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms of RHS of main equation:");
         print_norms(world, rhs);
      }

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_operators = make_bsh_operators(world, shifts);

      // Apply BSH operators to RHS of equation
      if(tda_print_level >= 1) 
      {
         if(world.rank() == 0) print("   Applying BSH operators");
      }
      std::vector<std::vector<real_function_3d>> new_x_response = apply(world, bsh_operators, rhs);

      // Scale by -2.0 (coefficient in eq. 37 of reference paper)
      new_x_response = scale(new_x_response, -2.0);

      // Project out ground state
      for(int i = 0; i < m; i++) new_x_response[i] = projector(new_x_response[i]);

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms after application of BSH");
         print_norms(world, new_x_response);
      }

      // Get the difference between old and new
      differences = tda_x_response - new_x_response;

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Response function residuals:");
         print_norms(world, differences);
      }

      // Calculate energy correction
      // NOTE: full formula is calculated, so just add to omega to update 
      energy_updates = calculate_energy_update(world, rhs, differences, new_x_response, tda_print_level);

      // Update energies now
      tda_omega += energy_updates;

      // Basic output
      if(tda_print_level >= 1) 
      {
         // Energy update
         if(world.rank() == 0) print("   Updated energies:"); 
         if(world.rank() == 0) print(tda_omega);
         if(world.rank() == 0) print("");
      }
 
      // KAIN solver update
      if(kain_on)
      { 
         //new_x_response = kain.update(new_x_response, differences); 
         //temp = kain.update(new_x_response, differences); 

         // Step restriction
         //new_x_response = step_restriction(world, new_x_response, temp);
      }      

      // Basic output
      if(tda_print_level >= 1) 
      {   
         // Orthonormalize functions
         if(world.rank() == 0 ) print("   Orthogonalizing response functions");
      }

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms before orthogonalization:");
         print_norms(world, new_x_response);
      }

      // Modified gram-schmidt
      tda_x_response = gram_schmidt(world, new_x_response);

      // Debugging output
      if(tda_print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms after orthogonalization:");
         if(world.rank() == 0) print("      Note: The row will be normalized,\n      not the individual functions.");
         print_norms(world, tda_x_response);
      }

      // Check convergence
      if(energy_updates.absmax() < tda_energy_threshold) converged = true;

      // Update counter
      iteration += 1;

      // Default print
      if(tda_print_level == 0) 
      {
         Tensor<double> current_time = end_timer(world);
         if(world.rank() == 0) printf("%6d      % 10.6e   % 10.6e    %5.2f\n", iteration, energy_updates.absmax(), calculate_max_residual(world, differences), current_time[0] - initial_time[0]);
      }

      // Done with the iteration.. truncate?
      truncate(world, tda_x_response);
   }
   
   if(world.rank() == 0) print("\n");

   // Did we converge?
   if(iteration == tda_max_iterations && energy_updates.absmax() > tda_energy_threshold) print("\n  ***  Ran out of iterations  ***\n");

   // Sort values and functions into ascending order based on values
   sort(world, tda_omega, tda_x_response);

   // Print final things
   if(world.rank() == 0) print(" Final energies:\n"); 
   if(world.rank() == 0) print(tda_omega);
   if(world.rank() == 0) print(" Final energy residuals:");
   if(world.rank() == 0) print(energy_updates);
   if(world.rank() == 0) print(" Final response function residuals:");
   print_norms(world, differences);   

   // Done with iterate
}

// Iterates the response functions until converged or out of iterations
Tensor<double> TDA::guess_iterate(World & world,
                                  std::vector<std::vector<real_function_3d>> guess,
                                  std::vector<real_function_3d> orbitals,
                                  Tensor<double> energies,
                                  double thresh,
                                  double small,
                                  int print_level)
{
   // Variables needed to iterate
   int iteration = 0;                                                            // Iteration counter
   int max_iterations = 50;                                                      // Chosen arbitrarily
   bool converged = false;                                                       // For convergence
   QProjector<double, 3> projector(world, orbitals);                             // Projector to project out ground state
   int m = guess.size();                                                         // Number of guess orbitals (symmetry functions * ground state) 
   int n = guess[0].size();                                                      // Number of ground state orbitals
   Tensor<double> guess_omega(m);                                                // Holds energy corrections
   std::vector<std::vector<real_function_3d>> differences;                       // Holds wave function corrections
   std::vector<std::vector<real_function_3d>> temp;                              // Used for step restriction 

   // Print headers
   if(print_level == 0) print_guess_iterate_headers(world);

   // Create gamma 
   //std::vector<std::vector<real_function_3d>> gamma = create_gamma(world, m, n, guess, orbitals, small, thresh, print_level);

   // Create \hat{V}^0 applied to guess functions 
   //std::vector<std::vector<real_function_3d>> V_response = create_potential(world, m, n, guess, print_level);

   // Get a start time
   Tensor<double> initial_time = end_timer(world);

   // Now to iterate
   while( iteration < max_iterations && !converged)
   {
      // Basic output
      if(print_level >= 1)
      {
         if(world.rank() == 0) print("\n   Guess Iteration", iteration);
         if(world.rank() == 0) print("  --------------------");
      }

      // Create gamma 
      std::vector<std::vector<real_function_3d>> gamma = create_gamma(world, m, n, guess, orbitals, small, thresh, print_level);

      // Project out ground state
      for(int i = 0; i < m; i++) gamma[i] = projector(gamma[i]);

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms of gamma after projecting out ground state:");
         print_norms(world, gamma);
      }

      // Create \hat{V}^0 applied to guess functions 
      std::vector<std::vector<real_function_3d>> V_response = create_potential(world, m, n, guess, print_level);

      // Basic output
      if(print_level >= 1) 
      {
         if(world.rank() == 0) print("   Solving Ax=Swx for initial energies");
      }

      // Constructing S
      Tensor<double> S = create_S(world, m, n, guess, print_level);
      
      // Constructing A
      Tensor<double> A = create_A(world, gamma, V_response, guess, orbitals, energies, print_level);

      // Symmeterize A
      A = 0.5 * (A + transpose(A));

      // Solve Ax = Sxw
      Tensor<double> e_vecs;
      sygv(A, S, 1, e_vecs, guess_omega);
      //diag_fock_matrix(world, A, guess, V_response, gamma, guess_omega, S, thresh);             

      // Basic output
      if(print_level >= 1) 
      {
         if(world.rank() == 0) print("   Guess Orbital Energies:");
         if(world.rank() == 0) print(guess_omega);
      }

      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      Tensor<double> shifts = create_shift(world, m, n, energies, guess_omega, print_level);  
     
      // Apply the shifts
      V_response = apply_shift(world, shifts, V_response, guess);

      // Construct RHS of equation
      std::vector<std::vector<real_function_3d>> rhs  = gamma + V_response;

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms of RHS of main equation:");
         print_norms(world, rhs);
      }

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_operators = make_bsh_operators(world, shifts, energies, guess_omega, small, thresh);

      // Apply BSH operators to RHS of equation
      if(print_level >= 1)
      {
         if(world.rank() == 0) print("   Applying BSH operators");
      }
      std::vector<std::vector<real_function_3d>> new_guess = apply(world, bsh_operators, rhs);

      // Scale by -2.0 (coefficient in eq. 37 of reference paper)
      new_guess = scale(new_guess, -2.0);

      // Project out ground state
      for(int i = 0; i < m; i++) new_guess[i] = projector(new_guess[i]);

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms after application of BSH");
         print_norms(world, new_guess);
      }

      // Get the difference between old and new
      differences = guess - new_guess;

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Guess function residuals:");
         print_norms(world, differences);
      }

      // Calculate energy correction
      // NOTE: full formula is calculated, so just add to omega to update 
      Tensor<double> energy_updates = calculate_energy_update(world, rhs, differences, guess, print_level);

      // Update energies now
      guess_omega += energy_updates;

      // Basic output
      if(print_level >= 1) 
      {
         // Energy update
         if(world.rank() == 0) print("   Updated guess energies:"); 
         if(world.rank() == 0) print(guess_omega);
         if(world.rank() == 0) print("");
      }

      // Basic output
      if(print_level >= 1) 
      {   
         // Orthonormalize functions
         if(world.rank() == 0 ) print("   Orthogonalizing guess functions");
      }

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms before orthogonalization:");
         print_norms(world, new_guess);
      }

      // Modified gram-schmidt
      guess = gram_schmidt(world, new_guess);

      // Debugging output
      if(print_level >= 2) 
      {
         if(world.rank() == 0) print("   Norms after orthogonalization:");
         print_norms(world, guess);
      }

      // Update counter
      iteration += 1;

      // Default print
      if(print_level == 0) 
      {
         Tensor<double> current_time = end_timer(world);
         if(world.rank() == 0) printf("%6d      % 10.6e   % 10.6e    %5.2f\n", iteration, energy_updates.absmax(), calculate_max_residual(world, differences), current_time[0] - initial_time[0]);
      }

      // Check for convergence
      if(energy_updates.absmax() < thresh) converged = true;

      // Done with the iteration.. truncate?
      truncate(world, guess);
   }

   // Done with guess_iterate
   return guess_omega;
}

// Calculates the full polarizability tensor given the 
Tensor<double> TDA::polarizability(World & world,
                                   std::vector<std::vector<std::vector<real_function_3d>>> f)
{
   // Return container
   Tensor<double> polarizability(3,3);

   // Run over all axes
   //for(int axis1 = 0; axis1 < 3; axis1++)
   //{
   //   // Now get the density of each perturbation wave function
   //   real_function_3d density = sum(world, f[axis1]);

   //   // Now run over each dipole operator
   //   for(int axis2 = 0; axis2 < 3; axis2++)
   //   {
   //      // Tells functor which axis to use
   //      std::vector<int> f(3,0);
   //      f[axis1] = true;

   //      // Creating the functor  
   //      real_function_3d r = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));

   //      // Calculate the polarizability tensor values
   //      polarizability(axis1, axis2) = -inner(density, r*density); 
   //   }
   //}
   
   // Done
   return polarizability;
}


// Main function, makes sure everything happens in correcct order
void TDA::solve(World & world)
{
   // Welcome user (future ASCII art of Robert goes here) 
   if(world.rank() == 0) print("\n   Preparing to solve the TDA equations.\n");
   print_madness_params(world);
   print_molecule(world);
   if(world.rank() == 0) print("   Ground state energies:\n", tda_ground_energies);
   print_response_params(world);

   // Plotting input orbitals
   //print("\n   Plotting ground state densities.\n");
   //do_vtk_plots(world, 200, tda_L, 0, tda_num_orbitals, tda_molecule, square(world, tda_orbitals), "ground");

   // Create initial guesses
   if(world.rank() == 0)
   {
      print("");
      print("   Creating trial functions.");
   }      

   // Create large number of symmetry included guesses 
   // Doubling the requested states to help deal with root skipping
   //std::vector<std::vector<real_function_3d>> guesses = create_trial_functions1(world, tda_num_excited, tda_orbitals, 2, tda_print_level);
   std::vector<std::vector<real_function_3d>> guesses = create_trial_functions3(world, 2*tda_num_excited, tda_orbitals, tda_print_level);

   // Normalize
   normalize(world, guesses);

   // Basic output
   if(world.rank() == 0) print("\n   Iterating trial functions for an improved initial guess.");

   // Now iterate a few times to let the functions differentiate themselves
   Tensor<double> guess_omega = guess_iterate(world, guesses, tda_orbitals, tda_ground_energies, 1e-2, tda_small, tda_print_level);

   // Now we need to choose the tda_num_excited lowest energy states
   tda_x_response = select_trial_functions(world, guesses, guess_omega, tda_num_excited, tda_print_level);

   // Ready to iterate!
   iterate(world);

   // Any post processing would go here? 
   // Plot the response function
   // Need to sum contributions to get actual orbital first
   //std::vector<real_function_3d> x = zero_functions<double, 3>(world, tda_num_excited);
   //for(int i = 0; i < tda_num_excited; i++)
   //{
   //   for(int j = 0; j < tda_num_orbitals; j++)
   //   {
   //      x[i] = x[i] + tda_x_response[j] * tda_x_response[j];
   //   }
   //}

   //// Now plot
   //print("\n   Plotting excited state densities.\n");
   //do_vtk_plots(world, 200, tda_L, 0, tda_num_excited, tda_molecule, x, "excited");

}


// Dueces



