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

// Masking function to switch from 0 to 1 smoothly at boundary
// Pulled from SCF.h
inline double mask1(double x) {
    /* Iterated first beta function to switch smoothly
       from 0->1 in [0,1].  n iterations produce 2*n-1
       zero derivatives at the end points. Order of polyn
       is 3^n.
       
       Currently use one iteration so that first deriv.
       is zero at interior boundary and is exactly representable
       by low order multiwavelet without refinement */
    
    x = (x*x*(3.-2.*x));
    return x;
}

static double mask3(const coord_3d& ruser) {
    coord_3d rsim;
    user_to_sim(ruser, rsim);
    double x= rsim[0], y=rsim[1], z=rsim[2];
    double lo = 0.0625, hi = 1.0-lo, result = 1.0;
    double rlo = 1.0/lo;
    
    if (x<lo)
        result *= mask1(x*rlo);
    else if (x>hi)
        result *= mask1((1.0-x)*rlo);
    if (y<lo)
        result *= mask1(y*rlo);
    else if (y>hi)
        result *= mask1((1.0-y)*rlo);
    if (z<lo)
        result *= mask1(z*rlo);
    else if (z>hi)
        result *= mask1((1.0-z)*rlo);
    
    return result;
}


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

   // Create the mask function
   mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());
}

// Save the current response calculation
void TDHF::save(World & world)
{
   // Archive to write everything to
   archive::ParallelOutputArchive ar(world, "restart_response", 1); // Just going to enforce 1 io server

   // Saving, in this order;
   //  string           ground-state archive name (garch_name)
   //  bool             TDA flag
   //  int              number of ground state orbitals (n)
   //  int              number of excited state orbitals (m)
   //  Tensor<double>   energies of m x-states
   //  for i from 0 to m-1
   //     for j from 0 to n-1
   //        Function<double,3> x_response[i][j]
   //  (If TDA flag == True)
   //  (Tensor<double>  energies of m y-states    )
   //  (for i from 0 to m-1                       ) 
   //  (   for j from 0 to n-1                    )
   //  (      Function<double,3> y_response[i][j] )
   ar & Gparams.inFile;
   ar & Rparams.tda;
   ar & Gparams.num_orbitals;
   ar & Rparams.states;
   ar & x_omega; 

   //for(int i=0; i<Rparams.states; i++)
   //   ar & x_omega(i);
   for(int i=0; i<Rparams.states; i++)
      for(unsigned int j=0; j<Gparams.num_orbitals; j++)
         ar & x_response[i][j];
   if(Rparams.tda)
   {
      //for(int i=0; i<Rparams.states; i++)
      //   ar & y_omega(i);
      ar & y_omega;
      for(int i=0; i<Rparams.states; i++)
         for(unsigned int j=0; j<Gparams.num_orbitals; j++)
            ar & y_response[i][j]; 
   }
}

// Load a response calculation
//void TDHF::load(World& world,
//                std::string archive)
//{
//   // The archive to read from
//   archive::ParallelInputArchive ar(world, archive.c_str());
//
//   // Reading in, in this order;
//   //  string           ground-state archive name (garch_name)
//   //  bool             TDA flag
//   //  int              number of ground state orbitals (n)
//   //  int              number of excited state orbitals (m)
//   //  Tensor<double>   energies of m x-states
//   //  for i from 0 to m-1
//   //     for j from 0 to n-1
//   //        Function<double,3> x_response[i][j]
//   //  (If TDA flag == True)
//   //  (Tensor<double>  energies of m y-states    )
//   //  (for i from 0 to m-1                       ) 
//   //  (   for j from 0 to n-1                    )
//   //  (      Function<double,3> y_response[i][j] )
//   
//   ar & Rparams.archive;
//   ar & Rparams.tda;
//   ar & Gparams.num_orbitals;
//   ar & Rparams.states;
//   ar & x_omega;   
//
//   for(int i=0; i<Rparams.states; i++)
//      for(unsigned int j=0; j<Gparams.num_orbitals; j++)
//         ar & x_response[i][j];
//
//   if(Rparams.tda)
//   {
//      ar & y_omega;
//      
//      for(int i=0; i<Rparams.states; i++)
//         for(unsigned int j=0; j<Gparams.num_orbitals; j++)
//            ar & y_response[i][j];
//   }
//}

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
      // Get transition density 
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

// Creates the off diagonal (letter B) portions of response matrix
std::vector<std::vector<real_function_3d>> TDHF::create_B(World &world,
                                                          std::vector<std::vector<real_function_3d>> & f, 
                                                          std::vector<real_function_3d> & orbitals,
                                                          double small,
                                                          double thresh) 
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Zero function
   std::vector<std::vector<real_function_3d>> deriv_j = response_zero_functions(world, m, n);
   std::vector<std::vector<real_function_3d>> deriv_k = response_zero_functions(world, m, n);
   real_function_3d rho;

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Two pieces: coulomb and exchange 
   // Exchange first 
   // Need to run over virtual orbitals
   for(int k = 0; k < m; k++)
   {
      // Get transition density 
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

   // Coulomb 
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
            real_function_3d rho = f[k][i] * orbitals[p];

            // Apply coulomb operator
            rho = apply(op, rho);

            // Multiply by response function (k,i)
            // and add to total
            deriv_k[k][p] += rho * orbitals[i];
         }
      }
   }
   
   // Take care of coeficients
   deriv_j = scale(deriv_j, 2.0) - deriv_k; 
 
   // Project out the ground state
   QProjector<double, 3> projector(world, orbitals);          
   for(int i = 0; i<m; i++) deriv_j[i] = projector(deriv_j[i]);

   // Done
   return deriv_j; 
}

// Computes gamma(r) given the ground state orbitals and response functions
// Only for TDA
std::vector<std::vector<real_function_3d>> TDHF::create_gamma(World &world,
                                                              std::vector<std::vector<real_function_3d>> & f,
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
      Tensor<double> temp = expectation(world, f, deriv_J);
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
   //int m = V.size();
   //int n = V[0].size();

   // Container to return
   //std::vector<std::vector<real_function_3d>> fock =  response_zero_functions(world, m, n);
   std::vector<std::vector<real_function_3d>> fock;

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

   // Done
   return fock;
}

// Construct the Hamiltonian
Tensor<double> TDHF::create_response_matrix(World & world,
                                            std::vector<std::vector<real_function_3d>> & gamma,
                                            std::vector<std::vector<real_function_3d>> & V,
                                            std::vector<std::vector<real_function_3d>> & f,
                                            std::vector<real_function_3d> & ground_orbitals,
                                            Tensor<double> & hamiltonian, // Ground state 
                                            int print_level,
                                            std::string xy)
{
   // Sizes inferred from V and gamma
   int m = V.size();

   // Container for A
   Tensor<double> A(m,m);

   // Create the ground-state fock operator on response orbitals
   std::vector<std::vector<real_function_3d>> fock_resp = create_fock(world, V, f, print_level, xy);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Ground Fock matrix for %s states:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world, f, fock_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Need to calculate hamiltonian * x_response
   // Name of function sounds strange, I know...
   std::vector<std::vector<real_function_3d>> energy_resp = scale_2d(world, f, hamiltonian);    

   // Verify this keeps orbitals in the virtual space
   // Verify this annihilates an occupied orbital (leave occupied orbital in occupied space at least)

   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Energy scaled response orbitals for %s states:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world, f, energy_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Construct intermediary
   std::vector<std::vector<real_function_3d>> temp = gamma + fock_resp - energy_resp;

   // Need to run over excited states
   for(int k = 0; k < m; k++)
   {
      // Need another run over excited states
      for(int j = 0; j < m; j++)
      {
         // Run over all occupied orbitals to get their contribution
         // to the part of A we're calculating . Finally calculate 
         // \int dr f_p^k * temp (using vmra.h function, sum is included)
         A(k,j) = inner(f[k], temp[j]);
      }
   }

   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) printf("   Response matrix for %s states:\n", xy.c_str());
      if(world.rank() == 0) print(A);
   }

   // Done
   return A;
}

// Constructs full response matrix of
// [ A  B ] [ X ] = w [ X ]
// [-B -A ] [ Y ]     [ Y ]
Tensor<double> TDHF::create_full_response_matrix(World & world, 
                                                 std::vector<std::vector<real_function_3d>> x_b, // x perturbed two electron piece 
                                                 std::vector<std::vector<real_function_3d>> Vx,  // potential * x
                                                 std::vector<std::vector<real_function_3d>> x,   // x response functions
                                                 std::vector<std::vector<real_function_3d>> y_b, // y perturbed two electron piece
                                                 std::vector<std::vector<real_function_3d>> Vy,  // potential * y
                                                 std::vector<std::vector<real_function_3d>> y,   // y response functions
                                                 std::vector<real_function_3d> ground_orbitals,  // ground state orbitals
                                                 Tensor<double> ground_ham,                      // full ground state hamiltonian
                                                 double small,
                                                 double thresh,
                                                 int print_level)
{
   // Get size
   int m = x.size();

   // Create the A pieces (A_x is top left, A_y is bottom right)
   // The -1 suppresses output
   Tensor<double> A_x = create_response_matrix(world, x_b, Vx, x, ground_orbitals, ground_ham, -1, "x"); 
   Tensor<double> A_y = create_response_matrix(world, y_b, Vy, y, ground_orbitals, ground_ham, -1, "y");
   
   // Construct matrix rep. of B (y is first row, x is second)
   std::vector<std::vector<real_function_3d>> tmp1 = create_B(world, x, ground_orbitals, small, thresh);
   std::vector<std::vector<real_function_3d>> tmp2 = create_B(world, y, ground_orbitals, small, thresh);
   Tensor<double> B_x = expectation(world, x, tmp1);
   Tensor<double> B_y = expectation(world, y, tmp2);  

   // Construct the large, 2*m x 2*m size matrix to be returned
   // (m is number of states requested)
   // Use madness slicing
   Tensor<double> response_matrix(2*m, 2*m); 

   // Place top left A
   response_matrix(Slice(0, m-1, 1), Slice(0, m-1, 1)) = A_x;
   
   // Place top right B
   response_matrix(Slice(0, m-1, 1), Slice(m, 2*m-1, 1)) = B_y;
   
   // Place bot left B (adjoint here if complex)
   response_matrix(Slice(m, 2*m-1, 1), Slice(0, m-1, 1)) = -B_x;
   
   // Place bot right A (adjoint here if complex)
   response_matrix(Slice(m, 2*m-1, 1), Slice(m, 2*m-1, 1)) = -A_y;

   // Print matrix if user requests
   if(world.rank() == 0 and print_level >= 1)
   {
      print("   Full Coupled Response Matrix:");
      print(response_matrix);
   }

   // Done
   return response_matrix;
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

// Diagonalize the full response matrix, taking care of degenerate states
Tensor<double> TDHF::diag_full_response(World & world,
                                        Tensor<double> & full_response,
                                        std::vector<std::vector<real_function_3d>> & x,
                                        std::vector<std::vector<real_function_3d>> & Vx,
                                        std::vector<std::vector<real_function_3d>> & x_g,
                                        std::vector<std::vector<real_function_3d>> & y,
                                        std::vector<std::vector<real_function_3d>> & Vy,
                                        std::vector<std::vector<real_function_3d>> & y_g,
                                        Tensor<double> & x_evals, 
                                        Tensor<double> & y_evals,
                                        const double thresh)
{
    // Get sizes
    int m = x.size();
    int n = x[0].size();

    // Place the response functions into a 2*m vector of vectors (X first, then Y)
    // Also need to do this for V_x, V_y, x_gamma, y_gamma (for rotations)
    std::vector<std::vector<real_function_3d>> response_vector(n, std::vector<real_function_3d>(2*m)); 
    std::vector<std::vector<real_function_3d>> V_response(n, std::vector<real_function_3d>(2*m));
    std::vector<std::vector<real_function_3d>> gamma(n, std::vector<real_function_3d>(2*m));
     
    // To make this work with current transform methods,
    // actually need to make this a transpose 
    for(int i = 0; i < n; i++) 
    {
       for(int j = 0; j < m; j++) 
       {
          // Response functions
          response_vector[i][j]   = x[j][i];
          response_vector[i][j+m] = y[j][i];

          // Gamma
          gamma[i][j] = x_g[j][i];
          gamma[i][j+m] = y_g[j][i];

          // Potentials
          V_response[i][j] = Vx[j][i];
          V_response[i][j+m] = Vy[j][i];
       }
    }

    // Create overlap matrix of everything
    std::vector<std::vector<real_function_3d>> resp_vec_trans = transpose(response_vector); 
    Tensor<double> overlap = create_overlap(world, resp_vec_trans, Rparams.print_level, "x and y");
    
    // compute the unitary transformation matrix U that diagonalizes
    // the fock matrix
    Tensor<double> evals(2*m);
    Tensor<double> vecs = get_full_response_transformation(world, overlap, full_response, evals, thresh);

    // Copy energies into the correct tensors
    for(int i = 0; i < m; i++) 
    {
       x_evals(i) = evals(i+m); 
       y_evals(i) = evals(i);
    }

    // Transform the vectors of functions
    V_response = scale_2d(world, V_response, vecs);
    gamma = scale_2d(world, gamma, vecs);
    response_vector = scale_2d(world, response_vector, vecs);

    // Now put everything back where it belongs 
    for(int i = 0; i < n; i++) 
    {
       for(int j = 0; j < m; j++) 
       {
          // Response functions
          x[j][i] = response_vector[i][j+m]; 
          y[j][i] = response_vector[i][j];

          // Gamma
          x_g[j][i] = gamma[i][j+m];
          y_g[j][i] = gamma[i][j];

          // Potentials
          Vx[j][i] = V_response[i][j+m];
          Vy[j][i] = V_response[i][j];
       }
    }

    // Normalize (x and y only) and truncate all the new functions
    truncate(world, Vx);
    truncate(world, Vy);
    truncate(world, x_g);
    truncate(world, y_g);
    truncate(world, x);
    truncate(world, y);
    normalize(world, x);
    normalize(world, y);

    // Return eigenvector tensor
    return vecs;
}

// Similar to what robert did above in "get_fock_transformation"
Tensor<double> TDHF::get_full_response_transformation(World& world,
                                                      Tensor<double>& overlap,
                                                      Tensor<double>& full_response,
                                                      Tensor<double>& evals,
                                                      const double thresh_degenerate)    
{
    // Get size
    int m = overlap.dim(0);

    // Diagonalize (NOT A SYMMETRIC DIAGONALIZATION!!!!)
    // Potentially complex eigenvalues come out of this
    Tensor<std::complex<double>> omega(m);
    Tensor<double> U(m,m);
    ggevp(world, full_response, overlap, U, omega);

    // Eigenvectors come out oddly packaged if there are 
    // complex eigenvalues.
    // Currently only supporting real valued eigenvalues
    // so throw an error if any imaginary components are 
    // not zero enough
    double max_imag = abs(imag(omega)).max();
    print("   Max imaginary component of eigenvalues:", max_imag, "\n");
    MADNESS_ASSERT(max_imag == 0); // MUST BE REAL!
    evals = real(omega);

    bool switched = true;
    while (switched)
    {
       switched = false;
       for (int i = 0; i < m; i++)
       {
          for (int j = i + 1; j < m; j++)
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
    for (long i = 0; i < m; ++i)
       if (U(i, i) < 0.0)
          U(_, i).scale(-1.0);

    // Rotations between effectively degenerate states confound
    // the non-linear equation solver ... undo these rotations
    long ilo = 0; // first element of cluster
    while (ilo < m - 1)
    {
       long ihi = ilo;
       while (fabs(evals[ilo] - evals[ihi + 1]) < thresh_degenerate * 10.0 * std::max(fabs(evals[ilo]), 1.0))
       {
          ++ihi;
          if (ihi == m - 1)
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
       }
       ilo = ihi + 1;
    }

    full_response = 0;
    for (int i = 0; i < m; ++i)
       full_response(i, i) = evals(i);

    // Finally, lets sort the eigenvalues
    // and eigenvectors
    sort_eigenvalues(world, evals, U);

    return U;
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

// Sorts the given Tensor of energies and vector of functions
// in place
Tensor<int> TDHF::sort_eigenvalues(World & world,
                                   Tensor<double> & vals,
                                   Tensor<double> & vecs)
{
   // Get relevant sizes
   int k = vals.size();

   // Tensor to hold selection order
   Tensor<int> selected(k);

   // Copy everything...    
   std::vector<double> vals_copy;
   for(int i = 0; i < k; i++) vals_copy.push_back(vals[i]);
   Tensor<double> vals_copy2 = copy(vals);
   Tensor<double> vecs_copy = copy(vecs);

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

      // Put corresponding things in the correct place
      vals(i) =  vals_copy[i];
      vecs(_,i) = vecs_copy(_,j);

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
   std::vector<std::vector<real_function_3d>> step;                  // Used for step restriction 
   std::vector<std::vector<real_function_3d>> x_gamma;               // Holds the perturbed two electron piece
   std::vector<std::vector<real_function_3d>> y_gamma;               // Holds the perturbed two electron piece
   std::vector<std::vector<real_function_3d>> V_x_response;          // Holds V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> V_y_response;          // Holds V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> shifted_V_x_response;  // Holds the shifted V^0 applied to response functions
   std::vector<std::vector<real_function_3d>> shifted_V_y_response;  // Holds the shifted V^0 applied to response functions
   Tensor<double> S_x;                                               // Overlap matrix of response orbitals for x states

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
      x_gamma = create_gamma(world, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh, Rparams.print_level, "x");
      if(!Rparams.tda)
      {
         y_gamma = create_gamma(world, y_response, Gparams.orbitals, Rparams.small, Rparams.thresh, Rparams.print_level, "y");
      }

      // Create \hat{V}^0 applied to response functions
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

      // TDA approximation
      if(Rparams.tda)
      {
         // Basic output
         if(Rparams.print_level >= 1 and world.rank() == 0) print("   Solving Ax=Swx");

         // Constructing S
         S_x = create_overlap(world, x_response, Rparams.print_level, "x"); 
      
         // Constructing response matrix
         Tensor<double> A_x = create_response_matrix(world, x_gamma, V_x_response, x_response, Gparams.orbitals, hamiltonian, Rparams.print_level, "x");

         // Solve Ax = Sxw
         Tensor<double> U_x = diag_fock_matrix(world, A_x, x_response, V_x_response, x_gamma, x_omega, S_x, Rparams.thresh);

         // Debugging output
         if(world.rank() == 0 and Rparams.print_level >= 2)
         {
            print("   Eigenvector coefficients from diagonalization for x states:");
            print(U_x);
         }
      }
      // Full TDHF
      else 
      {
         // Basic output
         if(Rparams.print_level >= 1 and world.rank() == 0) print("   Solving\n   [ A  B ][ X ] = S w [ X ]\n   [-B -A ][ Y ]       [ Y ]\n");

         // Construct full response matrix 
         Tensor<double> full_response = create_full_response_matrix(world, x_gamma, V_x_response, x_response,
                                                                    y_gamma, V_y_response, y_response, Gparams.orbitals,
                                                                    hamiltonian, Rparams.small, Rparams.thresh, Rparams.print_level);

         // Diagonalize         
         // Overlap matrix is constructed inside here
         Tensor<double> vecs = diag_full_response(world, full_response, x_response, V_x_response, x_gamma,
                                                  y_response, V_y_response, y_gamma, x_omega, y_omega, Rparams.thresh);

         // Debugging output
         if(world.rank() == 0 and Rparams.print_level >= 1)
         {
            print("   Eigenvector coefficients from diagonalization for x and y states:");
            print(vecs);
         } 
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

      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      x_shifts = create_shift(world, Gparams.energies, x_omega, Rparams.print_level, "x");

      // Negative here is to ensure we are looking at eps - omega 
      // (function is written explicity for eps + omega)
      if(not Rparams.tda) 
      {
         y_omega = -y_omega;
         y_shifts = create_shift(world, Gparams.energies, y_omega, Rparams.print_level, "y");
      }

      // Apply the shifts
      shifted_V_x_response = apply_shift(world, x_shifts, V_x_response, x_response);

      if(not Rparams.tda) shifted_V_y_response = apply_shift(world, y_shifts, V_y_response, y_response);

      // Construct RHS of equation
      std::vector<std::vector<real_function_3d>> rhs_x = x_gamma + shifted_V_x_response;
      std::vector<std::vector<real_function_3d>> rhs_y;
      if(not Rparams.tda)
      {
         // Add in coupling
         rhs_x = rhs_x +  create_B(world, y_response, Gparams.orbitals, Rparams.small, Rparams.thresh);

         // And construct y 
         rhs_y = shifted_V_y_response + y_gamma + create_B(world, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh);
      }

      // Add in localized orbital piece if using localized orbitals
      // This should be all off diagonal elements of ground state Fock
      // matrix 
      if(Rparams.localized)
      {
         std::vector<std::vector<real_function_3d>> temp = scale_2d(world, x_response, ham_no_diag); 
         rhs_x = rhs_x - temp;

         // Debugging output
         if(Rparams.print_level >= 2)
         {
            if(world.rank() == 0) print("   Norms of localized orbital correction for x states:");
            print_norms(world, temp);
         }

         if(not Rparams.tda)
         {
            temp = scale_2d(world, y_response, ham_no_diag);
            rhs_y = rhs_y - temp; 

            // Debugging output
            if(Rparams.print_level >= 2)
            {
               if(world.rank() == 0) print("   Norms of localized orbital correction for y states:");
               print_norms(world, temp);
            }
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
 
      if(not Rparams.tda) 
      {  
         bsh_y_operators = create_bsh_operators(world, y_shifts, Gparams.energies, y_omega, Rparams.small, Rparams.thresh);
      }     

      // Project out ground state 
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
         x_response = kain.update(x_response, x_differences); 
         if(not Rparams.tda) y_response = kain.update(y_response, y_differences);
      }
      else
      {
         x_response = new_x_response;
         if(not Rparams.tda) y_response = new_y_response;
      }

      // Apply mask
      for(int i = 0; i < m; i++) x_response[i] = mask * x_response[i];
      if(not Rparams.tda)
      { 
         for(int i = 0; i < m; i++) y_response[i] = mask * y_response[i];
      }

      // Calculate energy residual and update old_energy 
      energy_x_residuals = abs(x_omega - old_x_energy);
      old_x_energy = copy(x_omega);

      if(not Rparams.tda)
      {
         energy_y_residuals = abs(y_omega - old_y_energy);
         old_y_energy = copy(y_omega);
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
      truncate(world, x_response);
      normalize(world, x_response);

      if(not Rparams.tda)
      {
         truncate(world, y_response);
         normalize(world, y_response);
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

// TESTING
//      if(world.rank() == 0) print("Plotting orbitals");
//      // Need to sum contributions to get densities first
//      std::vector<real_function_3d> densities = zero_functions<double, 3>(world, Rparams.plot_data.size());
//      for(int i : Rparams.plot_data)
//      {
//         for(unsigned int j = 0; j < Gparams.num_orbitals; j++)
//         {
//            densities[i] = densities[i] + Gparams.orbitals[j] * x_response[i][j];
//         }
//      }
//      std::string iter_count = "reg-local-iter" + std::to_string(iteration) + "-";
//      densities = square(world, densities);
//      // Now plot
//      do_vtk_plots(world, 151, 5.0, 0, Rparams.plot_data.size(), Gparams.molecule, densities, iter_count);   
//
// END TESTING
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
   //sort(world, x_omega, energy_x_residuals, x_response, x_differences);
   //if(not Rparams.tda) sort(world, y_omega, energy_y_residuals, y_response, y_differences);

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
   analysis(world);

}   // Done with iterate. 

// More detailed analysis of the response functions
// Uses member variables
void TDHF::analysis(World & world)
{
   // Sizes get used a lot here, so lets get a local copy
   int n = x_response[0].size();
   int m = x_response.size();

   // Per response function, want to print the contributions from each ground state
   // So print the norm of each function?
   Tensor<double> x_norms(m, n);
   Tensor<double> y_norms(m, n); 

   // Calculate the inner products
   for(int i = 0; i < m; i++)
   {
      for(int j = 0; j < n; j++)
      {
         x_norms(i,j) = x_response[i][j].norm2();

         if(not Rparams.tda) y_norms(i,j) = y_response[i][j].norm2();
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
         dipoles(i,0) += inner(Gparams.orbitals[j], x * x_response[i][j]);
         dipoles(i,1) += inner(Gparams.orbitals[j], y * x_response[i][j]);
         dipoles(i,2) += inner(Gparams.orbitals[j], z * x_response[i][j]);

         if(not Rparams.tda) 
         {
            dipoles(i,0) += inner(Gparams.orbitals[j], x * y_response[i][j]);
            dipoles(i,1) += inner(Gparams.orbitals[j], y * y_response[i][j]);
            dipoles(i,2) += inner(Gparams.orbitals[j], z * y_response[i][j]);
         }
      }
   }

   // Calculate oscillator strength
   Tensor<double> oscillator(m);
   for(int i = 0; i < m; i++)
   {
      oscillator(i) = 2.0/3.0 * (dipoles(i,0)*dipoles(i,0) + dipoles(i,1)*dipoles(i,1) + dipoles(i,2)*dipoles(i,2)) * Gparams.energies(i);
   }

   // Calculate transition quadrapole moments
   Tensor<double> quadrapoles(m,3,3);

   // Run over each excited state 
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         quadrapoles(i,0,0) += inner(Gparams.orbitals[j], x * x * x_response[i][j]);
         quadrapoles(i,0,1) += inner(Gparams.orbitals[j], x * y * x_response[i][j]);
         quadrapoles(i,0,2) += inner(Gparams.orbitals[j], x * z * x_response[i][j]);
         quadrapoles(i,1,0) += inner(Gparams.orbitals[j], y * x * x_response[i][j]);
         quadrapoles(i,1,1) += inner(Gparams.orbitals[j], y * y * x_response[i][j]);
         quadrapoles(i,1,2) += inner(Gparams.orbitals[j], y * z * x_response[i][j]);
         quadrapoles(i,2,0) += inner(Gparams.orbitals[j], z * x * x_response[i][j]);
         quadrapoles(i,2,1) += inner(Gparams.orbitals[j], z * y * x_response[i][j]);
         quadrapoles(i,2,2) += inner(Gparams.orbitals[j], z * z * x_response[i][j]);

         if(not Rparams.tda)
         {
            quadrapoles(i,0,0) += inner(Gparams.orbitals[j], x * x * y_response[i][j]);
            quadrapoles(i,0,1) += inner(Gparams.orbitals[j], x * y * y_response[i][j]);
            quadrapoles(i,0,2) += inner(Gparams.orbitals[j], x * z * y_response[i][j]);
            quadrapoles(i,1,0) += inner(Gparams.orbitals[j], y * x * y_response[i][j]);
            quadrapoles(i,1,1) += inner(Gparams.orbitals[j], y * y * y_response[i][j]);
            quadrapoles(i,1,2) += inner(Gparams.orbitals[j], y * z * y_response[i][j]);
            quadrapoles(i,2,0) += inner(Gparams.orbitals[j], z * x * y_response[i][j]);
            quadrapoles(i,2,1) += inner(Gparams.orbitals[j], z * y * y_response[i][j]);
            quadrapoles(i,2,2) += inner(Gparams.orbitals[j], z * z * y_response[i][j]);
         }
      }
   }

   // Now print?
   if(world.rank() == 0)
   {
      for(int i = 0; i < m; i++)
      {
         printf("   Response Function %d\t\t%7.8f a.u.", i, x_omega(i));
         print ("\n   --------------------------------------------");

         print("\n   Transition Dipole Moments");
         printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", dipoles(i,0), dipoles(i,1), dipoles(i,2));

         printf("\n   Dipole Oscillator Strength: %7.8f\n", oscillator(i));

         print("\n   Transition Quadrapole Moments");
         printf("   %16s %16s %16s\n", "X", "Y", "Z");
         printf("  X %16.8f %16.8f %16.8f\n", quadrapoles(i,0,0), quadrapoles(i,0,1), quadrapoles(i,0,2));
         printf("  Y %16.8f %16.8f %16.8f\n", quadrapoles(i,1,0), quadrapoles(i,1,1), quadrapoles(i,1,2));
         printf("  Z %16.8f %16.8f %16.8f\n", quadrapoles(i,2,0), quadrapoles(i,2,1), quadrapoles(i,2,2));

         // Print contributions
         // Unique for tda/tdhf
         if(Rparams.tda)
         {
            print("\n   Norms of the Components:");
            for(int j = 0; j < n; j++)
            {
               printf("   Occupied %d  --->  Virtual %d   %7.8f\n", j, i, x_norms(i,j));
            }

            print("\n");
         }
         else
         {
            print("\n   Norms of the Components:");
            print("                                          x          y");
            for(int j = 0; j < n; j++)
            {
               printf("   Occupied %d  --->  Virtual %d   %7.8f %7.8f\n", j, i, x_norms(i,j), y_norms(i,j));
            }

            print("\n");

         }
      }
   }
}

// Diagonalizes the given functions
void TDHF::diagonalize_guess(World & world,
                            std::vector<std::vector<real_function_3d>> & f,
                            Tensor<double> & omega,
                            std::vector<real_function_3d> & orbitals,
                            Tensor<double> & energies,
                            double thresh,
                            double small,
                            int print_level,
                            std::string xy)
{
   // Create gamma 
   std::vector<std::vector<real_function_3d>> gamma = create_gamma(world, f, orbitals, small, thresh, print_level, xy);

   // Create \hat{V}^0 applied to guess functions 
   std::vector<std::vector<real_function_3d>> V_response = create_potential(world, f, print_level, xy);

   // Constructing S
   Tensor<double> S = create_overlap(world, f, print_level, xy);

   // Constructing response matrix
   Tensor<double> A = create_response_matrix(world, gamma, V_response, f, orbitals, energies, print_level, xy);

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
      y.scale(1e3);
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
      
      // Apply mask to get boundary condition right
      f_copy[i] = mask * f_copy[i];
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

   // Clear stored_potential
   stored_potential.clear();

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

   // If using localized orbitals, just save a matrix that is
   // (T+V) - Lambda * eye (so we can multiply this for RHS)
   if(Rparams.localized)
   {
      // Copy hamiltonian and zero the diagonal 
      ham_no_diag = copy(hamiltonian); 
      for(int i = 0; i < m; i++) ham_no_diag(i,i) = 0.0; 
   }
    
   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) print(hamiltonian);
   }

}

// Creates the transition density
// Uses member variables, not input parameters
std::vector<real_function_3d> TDHF::transition_density(World& world)
{
   // Get sizes
   int m = x_response.size();
   int n = Gparams.orbitals.size();

   // Return container 
   std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

   // Run over virtual...
   for(int i =0; i < m; i++)
      {
         // Run over occupied...
         for(int j = 0; j < n; j++)
         {
            densities[i] = densities[i] + Gparams.orbitals[j] * x_response[i][j];

            // Add in de-excitation if applicable
            if(not Rparams.tda)
            {
               densities[i] = densities[i] + Gparams.orbitals[j] * y_response[i][j];
            }
         }
      }

   // Done!
   return densities;
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
      if(Rparams.plot_L > 0.0) do_vtk_plots(world, Rparams.plot_pts, Rparams.plot_L, 0, Gparams.num_orbitals, Gparams.molecule, square(world, Gparams.orbitals), "ground");
      else do_vtk_plots(world, Rparams.plot_pts, Gparams.L/2.0, 0, Gparams.num_orbitals, Gparams.molecule, square(world, Gparams.orbitals), "ground");
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
   // Class variable "hamiltonian" is set here
   create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level);

   // Create trial functions by
   // creating a large number of symmetry included guesses 
   // or use random guesses
   if(Rparams.random)
   {
      if(world.rank() == 0) print("   Using a random guess for initial response functions.");
      x_response = response_zero_functions(world, Rparams.states, Gparams.num_orbitals);
      x_response = add_randomness(world, x_response);

      if(not Rparams.tda)
      {
         y_response = response_zero_functions(world, Rparams.states, Gparams.num_orbitals);
         y_response = add_randomness(world, y_response);
      }

      // Project out groundstate from guesses
      QProjector<double, 3> projector(world, Gparams.orbitals);
      for(unsigned int i = 0; i < x_response.size(); i++) x_response[i] = projector(x_response[i]);
      if(not Rparams.tda) for(unsigned int i = 0; i < y_response.size(); i++) y_response[i] = projector(y_response[i]);
  
      // Normalize
      normalize(world, x_response);
      if(not Rparams.tda) normalize(world, y_response);
   }
   else
   {
      if(world.rank() == 0) print("   Creating trial functions.\n");

      std::vector<std::vector<real_function_3d>> x_guesses;
      x_guesses = create_trial_functions(world, Rparams.states, Gparams.orbitals, Rparams.print_level);
     
      // Project out groundstate from guesses
      QProjector<double, 3> projector(world, Gparams.orbitals);
      for(unsigned int i = 0; i < x_guesses.size(); i++) x_guesses[i] = projector(x_guesses[i]);
  
      // Normalize
      normalize(world, x_guesses);

      // Basic output
      if(world.rank() == 0) print("\n   Diagonalizing trial functions for an improved initial guess.\n");

      Tensor<double> guess_x_omega(x_guesses.size());

      // Diagonalize
      // Inplace modificaiton of guesses and guess_omega
      // Using the Tamm-Danchof approximation in this, should still be good enough for a first guess
      diagonalize_guess(world, x_guesses, guess_x_omega, Gparams.orbitals, hamiltonian, Rparams.thresh, Rparams.small, Rparams.print_level, "x");

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

      // Create y states as a copy of x states
      // Probably need to do something smarter?
      if( not Rparams.tda)
      {
         y_response = copy(world, x_response); 
      }
   }

   // Initialize x and y omega
   x_omega = Tensor<double>(x_response.size());
   if(!Rparams.tda) y_omega = Tensor<double>(y_response.size());

   // Ready to iterate! 
   // Protocol things go here?
   iterate(world);

// NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// DONE SPECIFICALLY FOR BE-LOCAL ONLY!!!!!!!!
//
//   // Original ground state hamiltonian
//   if(world.rank() == 0) print("Original ground state hamiltonian");
//   create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level);
//   if(world.rank() == 0) print("Eigenvalues:\n", Gparams.energies);
//
//   std::vector<std::vector<real_function_3d>> x_gamma = create_gamma(world, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh, -2, "x");   
//   Tensor<double> gamma_matrix = expectation(world, x_response, x_gamma);
//   Tensor<double> gamma_vecs, gamma_vals;
//   syev(gamma_matrix, gamma_vecs, gamma_vals);
//   if(world.rank() == 0) print("Gamma eigenvalues:\n", gamma_vals);
//
//   std::vector<std::vector<real_function_3d>> V_x_response = create_potential(world, x_response, -2, "x");
//   Tensor<double> Vx_matrix = expectation(world, x_response, V_x_response);
//   Tensor<double> Vx_vecs, Vx_vals;
//   syev(Vx_matrix, Vx_vecs, Vx_vals);
//   if(world.rank() == 0) print("Vx eigenvalues:\n", Vx_vals);
//
//   std::vector<std::vector<real_function_3d>> fock_resp = create_fock(world, V_x_response, x_response, 0, "x");
//   fock_resp = fock_resp - V_x_response; // Leaves just kinetic
//   Tensor<double> fock_matrix = expectation(world, x_response, fock_resp);
//   Tensor<double> fock_vecs, fock_vals;
//   syev(fock_matrix, fock_vecs, fock_vals);
//   if(world.rank() == 0) print("kinetic eigenvalues:\n", fock_vals);
//
//   Tensor<double> A_x = create_response_matrix(world, x_gamma, V_x_response, x_response, Gparams.orbitals, hamiltonian, -2, "x");
//   if(world.rank() == 0) print("Original response matrix:");
//   if(world.rank() == 0) print(A_x);
//
//   // Look at eigenvalues
//   Tensor<double> overlap = create_overlap(world, x_response, -2, "x");
//   Tensor<double> vals;
//   Tensor<double> vecs = diag_fock_matrix(world, A_x, x_response, V_x_response, x_gamma, vals, overlap, Rparams.thresh);
//   if(world.rank() == 0) print("Response Eigenvalues:\n", vals);
//
//
//
//   // Rotate here
//   // Rotating ground state by pi/4
//   if(world.rank() == 0) print("\nRotating by pi/4.");
//   Tensor<double> rotation(2,2);
//   rotation(0,0) = cos(constants::pi/4.);  rotation(0,1) = -sin(constants::pi/4.);
//   rotation(1,0) = sin(constants::pi/4.);  rotation(1,1) =  cos(constants::pi/4.);
//   Gparams.orbitals = madness::transform(world, Gparams.orbitals, rotation, false);
//
//   // Rotating response orbitals by -pi/4
//   x_response = scale_2d(world, x_response, rotation);
// 
//   // Turn on off diagonal elements
//   Rparams.localized = true;
// 
//   // Remake groundstate hamiltonian 
//   create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level);
//
//   // Make diagonal the new eigenvalues
//   for(int i = 0; i < 2; i++) Gparams.energies(i) = hamiltonian(i,i);
//
//   // Diagonalize hamiltonian
//   syev(hamiltonian, vecs, Gparams.energies);
//   if(world.rank() == 0) print("Eigenvalues:\n", Gparams.energies);
//   
//   // Remake response matrix
//   x_gamma = create_gamma(world, x_response, Gparams.orbitals, Rparams.small, Rparams.thresh, 0, "x");
//   gamma_matrix = expectation(world, x_response, x_gamma); 
//   syev(gamma_matrix, gamma_vecs, gamma_vals);
//   if(world.rank() == 0) print("Gamma eigenvalues:\n", gamma_vals);
//
//   V_x_response = create_potential(world, x_response, 0, "x");
//   Vx_matrix = expectation(world, x_response, V_x_response);
//   syev(Vx_matrix, Vx_vecs, Vx_vals);
//   if(world.rank() == 0) print("Vx eigenvalues:\n", Vx_vals);
//
//   fock_resp = create_fock(world, V_x_response, x_response, 0, "x");
//   fock_resp = fock_resp - V_x_response; // Leaves just kinetic
//   fock_matrix = expectation(world, x_response, fock_resp); 
//   syev(fock_matrix, fock_vecs, fock_vals);
//   if(world.rank() == 0) print("kinetic eigenvalues:\n", fock_vals);
//
//   A_x = create_response_matrix(world, x_gamma, V_x_response, x_response, Gparams.orbitals, hamiltonian, 0, "x");
//   if(world.rank() == 0) print("Rotated response matrix:");
//   if(world.rank() == 0) print(A_x);
//
//   // Make a copy so originals aren't changed in diagonalization
//   // V_x and x_gamma will be recalculated so they can change here and it won't matter
//   std::vector<std::vector<real_function_3d>> b = copy(world, x_response);
//
//   // Look at eigenvalues
//   overlap = create_overlap(world, x_response, -2, "x");
//   vecs = diag_fock_matrix(world, A_x, b, V_x_response, x_gamma, vals, overlap, Rparams.thresh);
//   if(world.rank() == 0) print("Response Eigenvalues:\n", vals);
//
//
//   // Iterate again 
//   iterate(world);
//
//DONE TESTING


// TESTING WITH WATER

// Clear old everything 
//Gparams.orbitals.clear();
//Gparams.energies.clear();
//
//// Change ground state to be localized orbitals
//Gparams.read(world, "archives/water-local");
//
//// Just to be sure
//if(world.rank() == 0) Gparams.print_params();
//
//// Turn on localization functions
//Rparams.localized = true;
//
//// Recreate ground state hamiltonian
//create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level);
//
//// Diagonalize and get rotation matrix
//Tensor<double> vecs, vals;
//syev(hamiltonian, vecs, vals);
//
//// Transpose?
//vecs = transpose(vecs);
//
//// Rotate response orbitals
//x_response = scale_2d(world, x_response, vecs);
//
//// iterate again?
//iterate(world);
//
// END TESTING

   // Plot the response function if desired
   if(Rparams.plot)
   {
      // Need to get densities first
      std::vector<real_function_3d> densities = transition_density(world); 

      // For the instance where we don't plot all the orbitals
      std::vector<real_function_3d> plot_densities;

      for(int i : Rparams.plot_data)
      {
         plot_densities.push_back(densities[i]);
      }

      // Densities are the square of the wave function
      densities = square(world, plot_densities);

      // Now plot
      if(world.rank() == 0) print("\n   Plotting response state densities.\n");
      if(Rparams.plot_L > 0.0) do_vtk_plots(world, Rparams.plot_pts, Rparams.plot_L, 0, Rparams.plot_data.size(), Gparams.molecule, plot_densities, "response-state");   
      else do_vtk_plots(world, Rparams.plot_pts, Gparams.L, 0, Rparams.plot_data.size(), Gparams.molecule, plot_densities, "response-state");   
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



