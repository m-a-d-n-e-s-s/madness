/*
 *
 *   Written by: bsundahl
 *   Date: A long time ago...
 *
 */ 

#include "TDHF2.h"
#include "Plot_VTK.h"
#include "TDHF_Basic_Operators2.h"
#include "NWChem.h"                // For nwchem interface
#include "projector.h"             // For easy calculation of (1 - \hat{\rho}^0)
#include "potentialmanager.h"
#include "../chem/SCFOperators.h"

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
   ResponseFunction operator()()
   {
      ResponseFunction f(world, num_vir, num_occ);;

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
static std::vector<double> ttt, sss;
static void start_timer(World& world)
{
   world.gop.fence();
   ttt.push_back(wall_time());
   sss.push_back(cpu_time());
}

// Needed for timers
static double pop(std::vector<double>& v)
{
   double x = v.back();
   v.pop_back();
   return x;
}

// Stops a timer
static void end_timer(World& world,
                      const char* msg)
{ 
   double wall = wall_time() - pop(ttt);
   double cpu  = cpu_time() - pop(sss);
   if(world.rank() == 0) printf("   timer: %20.20s %8.2fs %8.2fs\n", msg, cpu, wall); 
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
   FunctionDefaults<3>::set_truncate_mode(1);   
   FunctionDefaults<3>::set_truncate_on_project(true);

   // Initialize response state xcfuntion object
   xcf.initialize(Rparams.xc, false, world, true);

   // Create the masking function
   mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());

   // Load balance on orbitals
   if(world.size() > 1)
   {
      // Start a timer
      if(Rparams.print_level >= 1) start_timer(world); 
      if(world.rank() == 0) print(""); // Makes it more legible
 
      LoadBalanceDeux<3> lb(world);
      for(unsigned int j = 0; j < Gparams.num_orbitals; j++)
      {
         lb.add_tree(Gparams.orbitals[j], lbcost<double,3>(1.0,8.0),true);
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

      if(Rparams.print_level >= 1) end_timer(world, "Load balancing:");
   }
}

// Save the current response calculation
void TDHF::save(World & world)
{
   // Archive to write everything to
   archive::ParallelOutputArchive ar(world, "resp_restart", 1); // Just going to enforce 1 io server

   // Saving, in this order;
   //  string           ground-state archive name (garch_name)
   //  bool             TDA flag
   //  int              number of ground state orbitals (n)
   //  int              number of excited state orbitals (m)
   //  Tensor<double>   energies of m x-components
   //  for i from 0 to m-1
   //     for j from 0 to n-1
   //        Function<double,3> x_response[i][j]
   //  (If TDA flag == True)
   //  (Tensor<double>  energies of m y-components    )
   //  (for i from 0 to m-1                       ) 
   //  (   for j from 0 to n-1                    )
   //  (      Function<double,3> y_response[i][j] )
   ar & Gparams.inFile;
   ar & Rparams.tda;
   ar & Gparams.num_orbitals;
   ar & Rparams.states;
   ar & omega; 

   for(int i=0; i<Rparams.states; i++)
      for(unsigned int j=0; j<Gparams.num_orbitals; j++)
         ar & x_response[i][j];
   if(not Rparams.tda)
   {
      for(int i=0; i<Rparams.states; i++)
         for(unsigned int j=0; j<Gparams.num_orbitals; j++)
            ar & y_response[i][j]; 
   }
}

// Load a response calculation
void TDHF::load(World& world,
                std::string name)
{
   // The archive to read from
   archive::ParallelInputArchive ar(world, name.c_str());

   // Reading in, in this order;
   //  string           ground-state archive name (garch_name)
   //  bool             TDA flag
   //  int              number of ground state orbitals (n)
   //  int              number of excited state orbitals (m)
   //  Tensor<double>   energies of m x-components
   //  for i from 0 to m-1
   //     for j from 0 to n-1
   //        Function<double,3> x_response[i][j]
   //  (If TDA flag == True)
   //  (Tensor<double>  energies of m y-components    )
   //  (for i from 0 to m-1                       ) 
   //  (   for j from 0 to n-1                    )
   //  (      Function<double,3> y_response[i][j] )
   
   ar & Rparams.archive;
   ar & Rparams.tda;
   ar & Gparams.num_orbitals;
   ar & Rparams.states;
   ar & omega;   

   x_response = ResponseFunction(world, Rparams.states, Gparams.num_orbitals);

   for(int i=0; i<Rparams.states; i++) 
      for(unsigned int j=0; j<Gparams.num_orbitals; j++)
         ar & x_response[i][j];
   world.gop.fence();

   y_response = ResponseFunction(world, Rparams.states, Gparams.num_orbitals);
   if(not Rparams.tda)
   {      
      for(int i=0; i<Rparams.states; i++)
         for(unsigned int j=0; j<Gparams.num_orbitals; j++)
            ar & y_response[i][j];
      world.gop.fence();
   }
}

// (Each state's norm should be 1, not the 
// individual functions norms)
void TDHF::normalize(World & world,
                     ResponseFunction & f)
{
   // Run over rows
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Get the normalization constant
      // (Sum included inside inner) 
      double norm = inner(f[i], f[i]);
      norm = sqrt(norm);

      // Doing this to deal with zero functions. 
      // Maybe not smrt.
      if( norm == 0) continue;

      // And scale
      scale(world, f[i], 1.0/norm);
   }
};

// (Each state's norm should be 1, not the 
// individual functions norms)
void TDHF::normalize(World & world,
                     ResponseFunction & f,
                     ResponseFunction & g)
{
   // Run over rows
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Get the normalization constant
      // (Sum included inside inner) 
      double normf = inner(f[i], f[i]);
      double normg = inner(g[i], g[i]);
      double norm = sqrt(normf - normg);

      // Doing this to deal with zero functions. 
      // Maybe not smrt.
      if( norm == 0) continue;

      // And scale
      scale(world, f[i], 1.0/norm);
      scale(world, g[i], 1.0/norm);
   }
};

// Prints norms of the given vector of vector of functions
void TDHF::print_norms(World & world,
                       ResponseFunction f)
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


// Radial function
static double kronecker(int l, int n)
{
   if(l == n) return 1.0;
   return 0.0;
}



// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
std::map<std::vector<int>,real_function_3d> TDHF::solid_harmonics(World & world,
                                                             int n)
{
   // Container to return
   std::map<std::vector<int>, real_function_3d> result;

   // Create the basic x, y, z, constant and zero
   real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{1,0,0})));
   real_function_3d y = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,1,0})));
   real_function_3d z = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,0,1})));
   real_function_3d c = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(std::vector<int>{0,0,0})));
   real_function_3d zero = real_factory_3d(world);

   // Add in first few, since they're simple
   // Assuming n >= 1
   result[std::vector<int>{0,0}] = copy(c);
   result[std::vector<int>{0,-1}] = zero;
   result[std::vector<int>{0, 1}] = zero;
   result[std::vector<int>{-1,0}] = zero;

   // Generate the solid harmonics recursively from here
   for(int l = 0; l < n; l++)
   {
      // Calculate ends of this row first
      result[std::vector<int>{l+1, l+1}] = sqrt(pow(2,kronecker(l,0)*(2*l)/(2*l+1)))*(x*result[std::vector<int>{l,l}]-(1-kronecker(l,0)*y*result[std::vector<int>{l,-l}]));
      result[std::vector<int>{l+1,-l-1}] = sqrt(pow(2,kronecker(l,0)*(2*l)/(2*l+1)))*(y*result[std::vector<int>{l,l}]+(1-kronecker(l,0)*x*result[std::vector<int>{l,-l}]));

      // Formula below calls for some functions that don't exist.
      // Need zeroes where that would occur
      result[std::vector<int>{l+1,  l+2}] = zero;
      result[std::vector<int>{l+1, -l-2}] = zero; 

      // Run over quantum number m
      for(int m = -l; m < l+1; m++)
      {
         // Calculate remaining terms
         result[std::vector<int>{l+1,m}] = 1.0/std::sqrt((l+m+1)*(l-m+1)) * ((2*l+1) * z * result[std::vector<int>{l,m}] - sqrt((l+m)*(l-m)) * (x*x + y*y + z*z) * result[std::vector<int>{l-1,m}]);
      }
   }

   // Get rid of any zero functions we added
   for(auto it = result.begin(); it != result.end(); )
   {
      if (it->second.norm2() == 0) 
         it = result.erase(it);
      else 
         ++it;
   }

   // Also get rid of the constant
   result.erase(std::vector<int>{0,0});

   // Done
   return result;
}

// Returns initial guess functions as
// ground MO * solid harmonics
ResponseFunction TDHF::create_trial_functions(World & world,
                                              int k,
                                              std::vector<real_function_3d> & orbitals,
                                              int print_level)
{
   // Get size
   int n = orbitals.size();

   // Create solid harmonics such that num. solids * num. orbitals > k.
   // The total number of solid harmonics that exist up to level n is 
   // (n+1)^2 (because we count from zero)
   // Always do at least 8 (through the d orbital angular momentum functions,
   // minus )
   std::map<std::vector<int>, real_function_3d> solids = solid_harmonics(world, std::max(2.0,ceil(sqrt(k/n)-1)));

   // Useful info.
   if(world.rank() == 0) print("   Created", solids.size(), "solid harmonics.\n");

   // Container to return 
   ResponseFunction trials;

   // Counter for number of trials created 
   int count = 0;

   // Multiply each solid harmonic onto a ground state orbital
   for(int i = 0; i < n; i++)
   {
      // For each solid harmonic 
      for(auto key : solids)
      { 
         // Temp zero functions
         std::vector<real_function_3d> temp = zero_functions_compressed<double,3>(world, n);

         // Create one non-zero function and add to trials
         temp[count % n] = key.second * orbitals[n - count % n - 1];
         trials.push_back(temp);
         count++;
      }

      // Stop when we first get beyond k components
      if(count >= k) break;
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


// Returns dipole operator * molecular orbitals 
ResponseFunction TDHF::dipole_guess(World &world,
                                    std::vector<real_function_3d> orbitals)
{
   // Return container
   ResponseFunction dipole_guesses(world, 3, orbitals.size());

   for(int axis = 0; axis < 3; axis++)
   {
      // Create dipole operator in the 'axis' direction
      std::vector<int> f(3,0);
      f[axis] = true;
      real_function_3d dip = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));

      reconstruct(world, orbitals);

      // Create guesses
      for(unsigned int i = 0; i < dipole_guesses[0].size(); i++)
      {
         dipole_guesses[axis][i] = mul_sparse(dip, orbitals[i], false); 
      }
      world.gop.fence();
   }
   //dipole_guesses.truncate_rf();

   // Done
   return dipole_guesses;
}


// Returns the derivative of the coulomb operator, applied to ground state orbitals
ResponseFunction TDHF::create_coulomb_derivative(World &world,
                                                 ResponseFunction & f,
                                                 std::vector<real_function_3d> & orbitals,
                                                 double small,
                                                 double thresh)
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Zero function, to be returned
   ResponseFunction deriv_j(world, m, n);

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Temperary storage
   real_function_3d rho = real_function_3d(world);

   // Need to run over each state 
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
ResponseFunction TDHF::create_exchange_derivative(World &world,
                                                  ResponseFunction & f,
                                                  std::vector<real_function_3d> & orbitals,
                                                  double small,
                                                  double thresh) 
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Zero function, to be returned
   ResponseFunction deriv_k(world, m, n);

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

// Creates diagonal (letter A) portions of response matrix
ResponseFunction TDHF::create_A(World &world,
                                ResponseFunction & fe,
                                ResponseFunction & gamma,
                                ResponseFunction & V,
                                ResponseFunction & f,
                                std::vector<real_function_3d> & ground_orbitals,
                                Tensor<double> & hamiltonian, // Ground state 
                                int print_level,
                                std::string xy)
{
   // Sizes inferred from V and gamma
   int m = V.size();

   // Create A
   Tensor<double> A(m,m);

   // Create the ground-state fock operator on response components
   ResponseFunction fock_resp = create_fock(world, V, f, print_level, xy);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Ground Fock matrix for %s components:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world, f, fock_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Need to calculate hamiltonian * x_response
   // Name of function sounds strange, I know...
   ResponseFunction energy_resp = scale_2d(world, f, hamiltonian);    

   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Energy scaled response %s components:\n", xy.c_str());
      Tensor<double> temp2 = expectation(world, f, energy_resp);
      if(world.rank() == 0) print(temp2);
   }

   // Saving this here for larger subspace calculations
   fe = fock_resp - energy_resp;

   // And return the sum 
   return gamma + fe; 
}

// Creates the off diagonal (letter B) portions of response matrix
// Very similiar to create_gamma, but the order of ground state and
// response components are different inside the integrals
ResponseFunction TDHF::create_B(World &world,
                                ResponseFunction & f, 
                                ResponseFunction & g, 
                                std::vector<real_function_3d> & orbitals,
                                double small,
                                double thresh,
                                int print_level) 
{
   // Start a timer
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Initialize function
   ResponseFunction deriv_j(world, m, n);
   ResponseFunction deriv_k(world, m, n);
   ResponseFunction deriv_xc(world, m, n);
   ResponseFunction gamma;
   real_function_3d rho;

   // Need the coulomb operator 
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   // Two pieces: coulomb and exchange 
   // Coulomb first 
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

   // Exchange 
   // Determine if including HF exchange
   if(xcf.hf_exchange_coefficient())
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
               real_function_3d rho = f[k][i] * orbitals[p];

               // Apply coulomb operator
               rho = apply(op, rho);

               // Multiply by response function (k,i)
               // and add to total
               deriv_k[k][p] += rho * orbitals[i];
            }
         }
      }
   }   
   // Determine if DFT potential is needed
   if(xcf.hf_exchange_coefficient() != 1.0)
   {
      // Get v_xc
      std::vector<real_function_3d> vxc = create_fxc(world, orbitals, f, g);

      // Apply xc kernel to ground state orbitals
      for(int i = 0; i < m; i++)
      {
         deriv_xc[i] = mul_sparse(world, vxc[i], orbitals, thresh, false);
      }
      world.gop.fence();
   }

   // Take care of coeficients
   gamma = deriv_j * 2.0 + deriv_xc - deriv_k * xcf.hf_exchange_coefficient();  
 
   // Project out the ground state
   QProjector<double, 3> projector(world, orbitals);          
   for(int i = 0; i<m; i++) gamma[i] = projector(gamma[i]);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   B coulomb deriv matrix:\n");
      ResponseFunction t = deriv_j * 2.0;
      Tensor<double> temp = expectation(world, f, t);
      if(world.rank() == 0) print(temp);
      if(Rparams.xc == "hf")
      {
         if(world.rank() == 0) printf("   B exchange deriv matrix:\n");
         temp = expectation(world, f, deriv_k);
      }
      else
      {
         if(world.rank() == 0) printf("   B XC deriv matrix:\n");
         temp = expectation(world, f, deriv_xc);
      }
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   B gamma matrix:\n");
      temp = expectation(world, f, gamma);
      if(world.rank() == 0) print(temp);
   }

   // End timer
   if(print_level >= 1) end_timer(world, "   Creating B:");

   // Done
   return gamma; 
}

// Computes gamma(r) given the ground state orbitals and response functions
ResponseFunction TDHF::create_gamma(World &world,
                                    ResponseFunction & f,
                                    ResponseFunction & g,
                                    std::vector<real_function_3d> & orbitals,
                                    double small,
                                    double thresh,
                                    int print_level,
                                    std::string xy)
{
   // Start timer
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // The gamma function to be returned, intialized to zero
   ResponseFunction gamma(world, m, n);

   // Intermediaries, initialized to zero
   ResponseFunction deriv_K(world, m, n);
   ResponseFunction deriv_XC(world, m, n);

   // Perturbed coulomb piece 
   ResponseFunction deriv_J = create_coulomb_derivative(world, f, orbitals, small, thresh);

   // If including any HF exchange:
   if(xcf.hf_exchange_coefficient())
   {
      deriv_K = create_exchange_derivative(world, f, orbitals, small, thresh); 
   }
   // Get the DFT contribution
   if(xcf.hf_exchange_coefficient() != 1.0) 
   {
      // Get v_xc
      std::vector<real_function_3d> vxc = create_fxc(world, orbitals, f, g);

      // Apply xc kernel to ground state orbitals      
      for(int i = 0; i < m; i++)
      {
          deriv_XC[i] = mul_sparse(world, vxc[i], orbitals, thresh, false);
      }
      world.gop.fence();
   }

   // Now assemble pieces together to get gamma
   // Spin integration gives 2.0
   gamma = deriv_J * 2.0 + deriv_XC - deriv_K * xcf.hf_exchange_coefficient();   

   // Project out groundstate 
   QProjector<double, 3> projector(world, Gparams.orbitals);
   for(int i = 0; i<m; i++) gamma[i] = projector(gamma[i]);

   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Coulomb Deriv matrix for %s components:\n", xy.c_str());
      ResponseFunction t = deriv_J * 2.0;
      Tensor<double> temp = expectation(world, f, t);
      if(world.rank() == 0) print(temp);
      if(Rparams.xc == "hf")
      {
         if(world.rank() == 0) printf("   Exchange Deriv matrix for %s components:\n", xy.c_str());
         temp = expectation(world, f, deriv_K);
      }
      else
      {
         if(world.rank() == 0) printf("   XC Deriv matrix for %s components:\n", xy.c_str());
         temp = expectation(world, f, deriv_XC);
      }
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Gamma matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, gamma);
      if(world.rank() == 0) print(temp);
   }

   // Basic output
   if(print_level >= 1) end_timer(world, "Creating gamma:");

   truncate(world, gamma);

   // Done
   return gamma;
}


// Calculates ground state coulomb potential 
real_function_3d TDHF::coulomb(World& world)
{
   // Coulomb operator
   real_convolution_3d op = CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());

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
ResponseFunction TDHF::exchange(World& world,
                                ResponseFunction & f)
{
   // Get sizes
   int m = f.size();
   int n = f[0].size();

   // Coulomb operator
   real_convolution_3d op = CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());

   // Container for results and others
   ResponseFunction result(world, m, n);
   real_function_3d psif = real_function_3d(world);

   // Modified 'small memory' algorithm from SCF.cc
   f.reconstruct_rf();

   // Run over each excited state
   for(int k = 0; k < m; k++)
   {
      // And run over each occupied state
      for(int j = 0; j < n; j++)
      {
         // Get a vector of transition densities     
         auto phix = mul_sparse(world, Gparams.orbitals[j], f[k], FunctionDefaults<3>::get_thresh());

         // Clean up
         truncate(world, phix);         

         // Apply operator to each member of vector
         phix = apply(world, op, phix);

         // Clean up
         truncate(world, phix);

         // Final multiplication of each member of vector by a single function
         phix = mul_sparse(world, Gparams.orbitals[j], phix, FunctionDefaults<3>::get_thresh());

         // Add the vector to result
         gaxpy(world, 1.0, result[k], 1.0, phix);
      }
   }

   // Truncate
   truncate(world, result);

   // Done!
   return result;
}

// Returns the ground state potential applied to functions f
ResponseFunction TDHF::create_potential(World & world,
                                        ResponseFunction & f,
                                        XCOperator xc,
                                        int print_level,
                                        std::string xy)
{
   // Start a timer
   if(print_level >= 1) start_timer(world);

   // Return container
   ResponseFunction V_x_resp;

   // Computing \hat{V}^0 = v_nuc + v_coul + v_exch           
   // v_nuc first
   real_function_3d v_nuc, v_coul;
   if(!Rparams.store_potential)
   {
      PotentialManager manager(Gparams.molecule, "a");
      manager.make_nuclear_potential(world);
      //v_nuc = manager.vnuclear().truncate();
      v_nuc = manager.vnuclear();
      v_nuc.truncate();

      // v_coul next
      // This does not include final multiplication of each orbital 
      // 2.0 scale is from spin integration
      v_coul = coulomb(world);
      v_coul.scale(2.0);
   }
   else // Already pre-computed
   {
      v_nuc = stored_v_nuc;
      v_coul = stored_v_coul;
   }

   // Intermediaries
   ResponseFunction v_exch(world, f.size(), f[0].size());
   real_function_3d v_xc = real_factory_3d(world);

   // If including any exact HF exchange
   if(xcf.hf_exchange_coefficient())
   {
      // V_exch last
      // Multiplication by f functions is included in construction
      v_exch = exchange(world, f);
   }
   if(xcf.hf_exchange_coefficient() != 1.0) 
   {
      // Calculate DFT potential
      v_xc = xc.make_xc_potential();
   }

   // Assemble all the pieces for V_x
   V_x_resp = f * (v_coul + v_nuc + v_xc) - v_exch * xcf.hf_exchange_coefficient();

   // Debugging output
   if(print_level >= 2)
   {
      // Print potential energy matrices
      if(world.rank() == 0) printf("   Nuclear potential matrix for %s components:\n", xy.c_str());
      ResponseFunction temp1 = f * v_nuc;
      Tensor<double> temp = expectation(world, f, temp1);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Coulomb potential matrix for %s components:\n", xy.c_str());
      ResponseFunction temp2 = f * v_coul;
      temp = expectation(world, f, temp2);
      if(world.rank() == 0) print(temp);
      if(xcf.hf_exchange_coefficient())
      {
         if(world.rank() == 0) printf("   Exchange potential matrix for %s components:\n", xy.c_str());
         temp = expectation(world, f, v_exch);
      }
      else if(xcf.hf_exchange_coefficient() != 1.0)
      {
         if(world.rank() == 0) printf("   XC potential matrix for %s components:\n", xy.c_str());
         v_exch = f * v_xc;
         temp = expectation(world, f, v_exch); 
      }
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Total Potential Energy matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, V_x_resp);
      if(world.rank() == 0) print(temp);
   }
   truncate(world, V_x_resp);

   // Basic output
   if(print_level >= 1) end_timer(world, "Creating V0 * x:");

   truncate(world, V_x_resp);

   // Done
   return V_x_resp;
}


// Returns a Tensor of inner products, where
// result(i,j) = inner(a[i],b[j]).sum()
Tensor<double> TDHF::expectation(World &world,
                                 ResponseFunction & a,
                                 ResponseFunction & b)
{
   // Get sizes
   int dim_1 = a.size();
   int dim_2 = a[0].size();

   // Need to take transpose of each input ResponseFunction
   ResponseFunction aT(world, dim_2, dim_1);
   ResponseFunction bT(world, dim_2, dim_1);
   for(int i = 0; i < dim_1; i++)
   {
      for(int j = 0; j < dim_2; j++)
      {
         aT[j][i] = a[i][j];
         bT[j][i] = b[i][j];
      }
   }

   // Container for result
   Tensor<double> result(dim_1, dim_1);

   // Run over dimension two 
   for(int p = 0; p < dim_2; p++)
   {
      result += matrix_inner(world, aT[p], bT[p]);
   }

   // Done
   return result;
}

// Returns the ground state fock operator applied to functions f
ResponseFunction TDHF::create_fock(World & world,
                                   ResponseFunction & V,
                                   ResponseFunction & f,
                                   int print_level,
                                   std::string xy)
{
   // Debugging output
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Creating perturbed fock matrix for %s components\n", xy.c_str());
   }

   // Container to return
   ResponseFunction fock;

   // Fock = (T + V) * orbitals
   // Already have V (input parameter)
   // Create T
   // Make the derivative operators in each direction
   real_derivative_3d Dx(world, 0);
   real_derivative_3d Dy(world, 1);
   real_derivative_3d Dz(world, 2);

   // Apply derivatives to orbitals
   f.reconstruct_rf();
   ResponseFunction dvx = apply(world, Dx, f);
   ResponseFunction dvy = apply(world, Dy, f);
   ResponseFunction dvz = apply(world, Dz, f);

   // Apply again for 2nd derivatives
   ResponseFunction dvx2 = apply(world, Dx, dvx);
   ResponseFunction dvy2 = apply(world, Dy, dvy);
   ResponseFunction dvz2 = apply(world, Dz, dvz);

   // Add together derivatives
   fock = (dvx2 + dvy2 + dvz2) * (-0.5);

   // Debugging output   
   if(print_level >= 2)
   {
      if(world.rank() == 0) printf("   Kinetic energy matrix for %s components:\n", xy.c_str());
      Tensor<double> temp = expectation(world, f, fock);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Potential energy matrix for %s components:\n", xy.c_str());
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
					    ResponseFunction & fe,
                                            ResponseFunction & gamma,
                                            ResponseFunction & V,
                                            ResponseFunction & f,
                                            std::vector<real_function_3d> & ground_orbitals,
                                            Tensor<double> & hamiltonian, // Ground state 
                                            int print_level,
                                            std::string xy)
{
   // Start a timer
   if(print_level >= 1) start_timer(world);

   // Construct intermediary
   // Sets fe to be (\hat{fock} - eps)*f
   ResponseFunction temp = create_A(world, fe, gamma, V, f, ground_orbitals, hamiltonian, print_level, xy); 

   // Make the matrix
   Tensor<double> A = expectation(world, f, temp);

   // Debugging output
   if(print_level >= 2 and world.rank() == 0)
   {
      printf("\n   Response matrix for %s components:\n", xy.c_str());
      print(A);
   }

   // End timer
   if(print_level >= 1) end_timer(world, "Create resp. matrix:");

   // Done
   return A;
}

// Constructs the matrix, so really it does
// [ X Y ] [ A  B ] [ X ]
//         [ B  A ] [ Y ] 
Tensor<double> TDHF::create_full_response_matrix(World & world, 
                                                 ResponseFunction & x_gamma, // x perturbed two electron piece 
                                                 ResponseFunction & Vx,      // potential * x
                                                 ResponseFunction & B_x,     // mixing terms for y
                                                 ResponseFunction & fe_x,    // eps * x 
                                                 ResponseFunction & x,       // x response functions
                                                 ResponseFunction & y_gamma, // y perturbed two electron piece
                                                 ResponseFunction & Vy,      // potential * y
                                                 ResponseFunction & B_y,     // mixing terms for x
                                                 ResponseFunction & fe_y,    // eps * y 
                                                 ResponseFunction & y,       // y response functions
                                                 std::vector<real_function_3d> & ground_orbitals, // ground state orbitals
                                                 Tensor<double> & ground_ham,                     // full ground state hamiltonian
                                                 double small,
                                                 double thresh,
                                                 int print_level)
{
   // Start timer
   if(print_level >= 1) start_timer(world);

   // Create the A pieces
   // (Sets fe_x and fe_y to be (\hat{F}-eps) * resp. funcs.
   ResponseFunction A_x = create_A(world, fe_x, x_gamma, Vx, x, ground_orbitals, ground_ham, print_level, "x"); 
   ResponseFunction A_y = create_A(world, fe_y, y_gamma, Vy, y, ground_orbitals, ground_ham, print_level, "y");
   
   // Create the B pieces 
   B_x = create_B(world, x, y, ground_orbitals, small, thresh, print_level);
   B_y = create_B(world, y, x, ground_orbitals, small, thresh, print_level);

   // Debugging output 
   if(print_level >= 2)
   {
      Tensor<double> t = expectation(world, x, A_x);
      if(world.rank() == 0) 
      {
         print("\n   A_x:");
         print(t);
      }

      t = expectation(world, x, B_y);
      if(world.rank() == 0) 
      {
         print("\n   B_y:");
         print(t);
      }

      t = expectation(world, y, A_y);
      if(world.rank() == 0)
      {
         print("\n   A_y:");
         print(t);
      }

      t = expectation(world, y, B_x);
      if(world.rank() == 0)
      {
         print("\n   B_x:");
         print(t);
      }
   }


   // Add corresponding pieces together
   A_x = A_x + B_y;
   A_y = B_x + A_y; // Needs adjoint if complex

   // Construct matrix to be returned
   Tensor<double> response_matrix = expectation(world, x, A_x) + expectation(world, y, A_y); 

   // Debugging output 
   if(print_level >= 2 and world.rank() == 0)
   {
      print("\n   Full Coupled Response Matrix:");
      print(response_matrix);
   }

   // End timer 
   if(print_level >= 1) end_timer(world, "Create resp. matrix:");
   
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
   // Start a timer
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = omega.size();
   int n = ground.size();

   // Container to hold shift
   Tensor<double> result(m,n);

   // Run over excited components
   for(int k = 0; k < m; k++)
   {
      // Run over ground components
      for(int p = 0; p < n; p++)
      {
         if(ground(p) + omega(k) > 0)
         {
            // Calculate the shift needed to get energy to -0.05,
            // which was arbitrary (same as moldft)
            result(k,p) = -(ground(p) + omega(k) + 0.05);

            // Basic output
            if(print_level >= 2)
            {
               if(world.rank() == 0) printf("   Shift needed for transition from ground orbital %d to response %s state %d\n", p, xy.c_str(), k);
               if(world.rank() == 0) print("   Ground energy =", ground(p));
               if(world.rank() == 0) print("   Excited energy =", omega(k));
               if(world.rank() == 0) print("   Shifting by", result(k,p));
               if(world.rank() == 0) print("");
            }
         }
      }
   }

   // End timer
   if(print_level >= 1) end_timer(world, "Create shift:");

   // Done
   return result;
}

// Returns the shift needed to make sure that
// (ground_state_energy + excited_state_energy + shift) = target 
// Please note: The same shift needs to be applied to the potential.
Tensor<double> TDHF::create_shift_target(World & world,
                                         Tensor<double> & ground,
                                         Tensor<double> & omega,
                                         double target,
                                         int print_level,
                                         std::string xy)
{
   // Start a timer
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = omega.size();
   int n = ground.size();

   // Container to hold shift
   Tensor<double> result(m,n);

   // Run over excited components
   for(int k = 0; k < m; k++)
   {
      // Run over ground components
      for(int p = 0; p < n; p++)
      {
         // Calculate the shift needed to get energy to target 
         result(k,p) = -(ground(p) + omega(k) - target);

         // Basic output
         if(print_level >= 2)
         {
            if(world.rank() == 0) printf("   Shift needed for transition from ground orbital %d to response %s state %d\n", p, xy.c_str(), k);
            if(world.rank() == 0) print("   Ground energy =", ground(p));
            if(world.rank() == 0) print("   Excited energy =", omega(k));
            if(world.rank() == 0) print("   Shifting by", result(k,p));
            if(world.rank() == 0) print("");
         }         
      }
   }

   // End timer
   if(print_level >= 1) end_timer(world, "Create shift:");

   // Done
   return result;
}

// Returns the given shift applied to the given potential
ResponseFunction TDHF::apply_shift(World & world,
                                   Tensor<double> & shifts,
                                   ResponseFunction & V,
                                   ResponseFunction & f)
{
   // Start timer
   if(Rparams.print_level >= 1) start_timer(world);

   // Sizes inferred from V
   int n = V[0].size();
   int m = V.size();

   // Container to return
   ResponseFunction shifted_V(world, m, n);

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

   // End timer
   if(Rparams.print_level >= 1) end_timer(world, "Apply shift:");
   
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
   // Start timer
   if(Rparams.print_level >= 1) start_timer(world);

   // Sizes inferred from ground and omega
   int n = ground.size();
   int m = omega.size();

   // Make the vector
   std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> operators;

   // Make a BSH operator for each response function
   // Run over excited components
   for(int k = 0; k < m; k++)
   {
      // Container for intermediary
      std::vector<std::shared_ptr<real_convolution_3d>> temp(n);

      // Run over occupied components
      for(int p = 0; p < n; p++)
      {
         temp[p] = std::shared_ptr<SeparatedConvolution<double,3>>(BSHOperatorPtr3D(world, sqrt(-2.0 * (ground(p) + omega(k) + shift(k,p))), small, thresh));
      }

      // Add intermediary to return container
      operators.push_back(temp);
   }

   // End timer
   if(Rparams.print_level >= 1) end_timer(world, "Creating BSH ops:");

   // Done
   return operators;
}


// Returns the second order update to the energies of the excited components
// Not currently used.
Tensor<double> TDHF::calculate_energy_update(World & world,
                                             ResponseFunction & rhs,
                                             ResponseFunction & f_residuals,
                                             ResponseFunction & new_f,
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
      if(world.rank() == 0) printf("   Calculating energy residuals for %s components\n", xy.c_str());
   }

   // Size inferred
   int m = rhs.size();

   // Container for updates
   Tensor<double> updates(m);

   // Need to run over all functions in rhs and calculate inner products.
   // rhs contains the bra in the braket notation above, and f_residuals
   // is the ket.

   // Run over excited components
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
      if(world.rank() == 0) printf("   Energy residuals for %s components:\n", xy.c_str());
      if(world.rank() == 0) print(abs(updates));
   }

   // Done?
   return updates;
}

// Specialized for response calculations that returns orthonormalized
// functions
ResponseFunction TDHF::gram_schmidt(World & world,
                                    ResponseFunction & f)
{
   // Sizes inferred
   int m = f.size();

   // Return container
   ResponseFunction result = f.copy();

   // Orthogonalize
   for(int j = 0; j < m; j++)
   {
      // Need to normalize the row
      double norm = norm2(world, result[j]);

      // Now scale each entry      
      scale(world, result[j], 1.0/norm);

      // Project out from the rest of the vectors
      for(int k = j+1; k < m; k++)
      {
         // Temp function to hold the sum
         // of inner products
         // vmra.h function, line 627
         double temp = inner(result[j], result[k]);

         // Now subtract 
         gaxpy(world, 1.0, result[k], -temp, result[j]);
      }
   }

   truncate(world, result);

   // Done
   return result;
}

// Specialized for response calculations that returns orthonormalized
// functions
void TDHF::gram_schmidt(World & world,
                        ResponseFunction & f,
                        ResponseFunction & g)
{
   // Sizes inferred
   int m = f.size();

   // Orthogonalize
   for(int j = 0; j < m; j++)
   {
      // Need to normalize the row
      double norm = inner(f[j],f[j]) - inner(g[j],g[j]);

      // Now scale each entry      
      scale(world, f[j], 1.0/sqrt(norm));
      scale(world, g[j], 1.0/sqrt(norm));

      // Project out from the rest of the vectors
      for(int k = j+1; k < m; k++)
      {
         // Temp function to hold the sum
         // of inner products
         // vmra.h function, line 627
         double temp = inner(f[j], f[k]) - inner(g[j], g[k]);

         // Now subtract 
         gaxpy(world, 1.0, f[k], -temp, f[j]);
         gaxpy(world, 1.0, g[k], -temp, g[j]);
      }
   }

   truncate(world, f);
   truncate(world, g);
}

// Returns the max norm of the given vector of functions
double TDHF::calculate_max_residual(World & world,
                                    ResponseFunction & f)
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

      if(world.rank() == 0) print("   Selecting ground state subspace to excite from for components.");
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
ResponseFunction TDHF::select_functions(World & world,
                                        ResponseFunction & f,
                                        Tensor<double> & energies,
                                        int k,
                                        int print_level)
{
   // Container for result
   ResponseFunction answer;

   // Debugging output
   if(print_level >= 1)
   {
      if(world.rank() == 0) print("\n   Selecting the", k, "lowest excitation energy components.\n");
   }

   // Get rid of extra functions and save
   // the first k
   while(int(f.size()) > k) f.pop_back();
   answer = f;
   truncate(world, answer);

   // Get rid of extra energies and save 
   // the first k
   energies = energies(Slice(0,k-1)); 

   // Basic output
   if(print_level >= 1)
   {
      if(world.rank() == 0) print("   The selected components have excitation energies:");
      if(world.rank() == 0) print(energies);
   }

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
                                             Tensor<double> & overlap,
                                             Tensor<double> & fock,
                                             Tensor<double> & evals,
                                             const double thresh_degenerate)
{
   // Run an SVD on the overlap matrix and ignore values
   // less than thresh_degenerate
   Tensor<double> r_vecs;
   Tensor<double> s_vals;
   Tensor<double> l_vecs;
   Tensor<double> overlap_copy = copy(overlap);
   svd(overlap_copy, l_vecs, s_vals, r_vecs);

   // Debugging output
   if(Rparams.print_level >= 2 and world.rank() == 0)
   {
      print("\n   Singular values of overlap matrix:");
      print(s_vals);
      print("   Left singular vectors of overlap matrix:");
      print(l_vecs);
   }

   // Check how many singular values are less than 10*thresh_degen
   int num_sv = 0;
   for(int i = 0; i < s_vals.dim(0); i++)
   {
      if(s_vals(i) < 10 * thresh_degenerate)
      {
         if(world.rank() == 0 and num_sv == 0) print("");
         if(world.rank() == 0) printf("   Detected singular value (%.8f) below threshold (%.8f). Reducing subspace size.\n", 
                                      s_vals(i), 10*thresh_degenerate);
         num_sv++;
      }
      if(world.rank() == 0 and i == s_vals.dim(0) - 1 and num_sv > 0) print("");
   }

   // Going to use these a lot here, so just calculate them
   int size_l = s_vals.dim(0);
   int size_s = size_l - num_sv;
   Tensor<double> l_vecs_s(size_l, num_sv);

   // Transform into this smaller space if necessary
   if(num_sv > 0)
   {
      // Cut out the singular values that are small
      // (singular values come out in descending order)
      overlap = Tensor<double>(size_s, size_s);
      for(int i = 0; i < size_s; i++) overlap(i,i) = s_vals(i); 

      // Copy the active vectors to a smaller container
      l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

      // Debugging output
      if(Rparams.print_level >= 2 and world.rank() == 0)
      {
         print("   Reduced size left singular vectors of overlap matrix:");
         print(l_vecs_s);
      }

      // Transform
      Tensor<double> work(size_l, size_s);
      mxm(size_l, size_s, size_l, work.ptr(), fock.ptr(), l_vecs_s.ptr());
      fock = Tensor<double>(size_s, size_s);
      Tensor<double> l_vecs_t = transpose(l_vecs);
      mxm(size_s, size_s, size_l, fock.ptr(), l_vecs_t.ptr(), work.ptr());
   }

   // Diagonalize using lapack
   Tensor<double> U;
   sygv(fock, overlap, 1, U, evals);

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

   // Rotations between effectively degenerate components confound
   // the non-linear equation solver ... undo these rotations
   long ilo = 0; // first element of cluster
   while (ilo < nmo - 1)
   {
      long ihi = ilo;
      while (fabs(evals[ilo] - evals[ihi + 1]) < thresh_degenerate * 100.0 * std::max(fabs(evals[ilo]), 1.0))
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
      }
      ilo = ihi + 1;
   }

   fock = 0;
   for (unsigned int i = 0; i < nmo; ++i)
      fock(i, i) = evals(i);

   // If we transformed into the smaller subspace, time to transform back
   if(num_sv > 0)
   {
      // Temp. storage
      Tensor<double> temp_U(size_l, size_l);
      Tensor<double> U2(size_l, size_l);

      // Copy U back to larger size
      temp_U(Slice(0, size_s-1), Slice(0, size_s-1)) = copy(U);
      for(int i = size_s; i < size_l; i++) temp_U(i,i) = 1.0;

      // Transform back
      mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());

      U = copy(U2);
   }

   return U;
}

// Sorts the given Tensor of energies
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



/// diagonalize the fock matrix, taking care of degenerate components

/// Vpsi is passed in to make sure orbitals and Vpsi are in phase
/// @param[in]  world    the world
/// @param[in]           fock    
/// @param[inout]        psi     the orbitals
/// @param[inout]        Vpsi    the orbital times the potential
/// @param[inout]        gamma   the orbital times the perturbed potential
/// @param[out] evals    the orbital energies
/// @param[in]  overlap  the overlap matrix
/// @param[in]  thresh   threshold for rotation and truncation
Tensor<double> TDHF::diag_fock_matrix(World & world,
                                      Tensor<double> & fock,
                                      ResponseFunction & psi,
                                      ResponseFunction & Vpsi,
                                      ResponseFunction & gamma,
                                      ResponseFunction & fe,
                                      Tensor<double> & evals,
                                      Tensor<double> & overlap,
                                      const double thresh)
{
    // compute the unitary transformation matrix U that diagonalizes
    // the fock matrix
    Tensor<double> U = get_fock_transformation(world, overlap, fock, evals, thresh);

    // Sort into ascending order
    Tensor<int> selected = sort_eigenvalues(world, evals, U);

    // Debugging output
    if(Rparams.print_level >= 2 and world.rank() == 0)
    {
       print("   U:");
       print(U);
    }

    // Start timer
    if(Rparams.print_level >= 1) start_timer(world);

    // transform the orbitals and the potential
    // Truncate happens inside here
    Vpsi = transform(world, Vpsi, U);
    gamma = transform(world, gamma, U);
    fe = transform(world, fe, U);
    psi = transform(world, psi, U);

    // End timer
    if(Rparams.print_level >= 1) end_timer(world, "Transform orbs.:");

    // Normalize x
    normalize(world, psi);

    // Debugging output
    if(Rparams.print_level >= 2 and world.rank() == 0)
    {
       print("   Eigenvector coefficients from diagonalization:");
       print(U);
    }

    return U;
}

// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
ResponseFunction TDHF::transform(World & world,
                                 ResponseFunction & f,
                                 Tensor<double> & U)
{
   // Return container
   ResponseFunction result;

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

// If using a larger subspace to diagonalize in, this will put everything in the right spot
void TDHF::augment(World & world,
                   Tensor<double> & S_x,     
                   Tensor<double> & A_x,     
                   ResponseFunction & x_gamma,
                   ResponseFunction & x_response,
                   ResponseFunction & V_x_response,
                   ResponseFunction & x_fe,         // Contains fock and energy scaled orbitals
                   Tensor<double> & old_S, 
                   Tensor<double> & old_A, 
                   ResponseFunction & old_x_gamma, 
                   ResponseFunction & old_x_response, 
                   ResponseFunction & old_V_x_response, 
                   ResponseFunction & old_x_fe,
                   int print_level)
{
   // Basic output
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = x_gamma.size();

   // Create work space, will overwrite S and A in the end
   Tensor<double> temp_S(2*m, 2*m);
   Tensor<double> temp_A(2*m, 2*m);

   // Need to create off diagonal blocks of A, so
   // create temps that are the sums of current and
   // old components respectively
   ResponseFunction temp_cur = x_gamma + x_fe;
   ResponseFunction temp_old = old_x_gamma + old_x_fe;

   // Calculate correct inner products of upper off diagonal
   Tensor<double> off = expectation(world, x_response, temp_old);

   // Use slicing to put in correct spot 
   temp_A(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);

   // Now for lower off diagonal 
   off = expectation(world, old_x_response, temp_cur); 

   // Use slicing to put in correct spot 
   temp_A(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);

   // Put together the rest of A
   temp_A(Slice(0, m-1), Slice(0, m-1)) = copy(A_x);
   temp_A(Slice(m, 2*m-1), Slice(m, 2*m-1)) = copy(old_A);

   // Debugging output
   if(print_level >= 2 and world.rank() == 0)
   {
      print("   Before symmeterizing A:");
      print(temp_A);
   }

   // Save temp_A as A_x
   // Need to symmeterize A as well (?)
   A_x = 0.5 * (temp_A + transpose(temp_A)); 

   // Now create upper off diagonal block of S
   off = expectation(world, x_response, old_x_response); 

   // Use slicing to put in correct spot 
   temp_S(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);

   // Now the lower off diagonal block
   // (Go ahead and cheat and use the transpose...)
   off = transpose(off);   

   // Use slicing to put in correct spot 
   temp_S(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);

   // Put together the rest of S
   temp_S(Slice(0, m-1), Slice(0, m-1)) = copy(S_x);
   temp_S(Slice(m, 2*m-1), Slice(m, 2*m-1)) = copy(old_S);

   // Save temp_S as S_x
   S_x = copy(temp_S);

   // Add in old vectors to current vectors for the appropriate ones 
   for(int i = 0; i < m; i++)
   {
      x_response.push_back(old_x_response[i]);
      x_gamma.push_back(old_x_gamma[i]);
      V_x_response.push_back(old_V_x_response[i]);
      x_fe.push_back(old_x_fe[i]);
   }

   // End the timer
   if(print_level >= 1) end_timer(world, "Aug. resp. matrix:");

   // Debugging output
   if(print_level >= 2 and world.rank() == 0) 
   {
      print("\n   Augmented response matrix:");
      print(A_x);
   }

   // Debugging output
   if(print_level >= 2 and world.rank() == 0) 
   {
      print("   Augmented overlap matrix:");
      print(S_x);
   }

   // SUPER debugging
   if(print_level >= 3)
   {
      if(world.rank() == 0) print("   Calculating condition number of aug. response matrix");

      
   }
}

// If using a larger subspace to diagonalize in, this will put everything in the right spot
void TDHF::augment_full(World & world,
                        Tensor<double> & S,
                        Tensor<double> & A,
                        ResponseFunction & B_x, 
                        ResponseFunction & x_gamma,
                        ResponseFunction & x_response,
                        ResponseFunction & V_x_response,
                        ResponseFunction & x_fe,  // Contains V_x_response
                        ResponseFunction & B_y, 
                        ResponseFunction & y_gamma,
                        ResponseFunction & y_response,
                        ResponseFunction & V_y_response,
                        ResponseFunction & y_fe, // Contains V_y_response
                        Tensor<double> & old_S, 
                        Tensor<double> & old_A, 
                        ResponseFunction & old_B_x,
                        ResponseFunction & old_x_gamma, 
                        ResponseFunction & old_x_response, 
                        ResponseFunction & old_V_x_response, 
                        ResponseFunction & old_x_fe,
                        ResponseFunction & old_B_y,
                        ResponseFunction & old_y_gamma, 
                        ResponseFunction & old_y_response, 
                        ResponseFunction & old_V_y_response, 
                        ResponseFunction & old_y_fe,
                        int print_level)
{
   // Basic output
   if(print_level >= 1) start_timer(world);

   // Get sizes
   int m = x_gamma.size();

   // Create work space, will overwrite S and A in the end
   Tensor<double> temp_S(2*m, 2*m);
   Tensor<double> temp_A(2*m, 2*m);

   // Need to create off diagonal blocks of A, so
   // create temps that are the sums of current and
   // old components respectively
   ResponseFunction temp_cur_x = x_gamma + x_fe + B_y;
   ResponseFunction temp_cur_y = y_gamma + y_fe + B_x;
   ResponseFunction temp_old_x = old_x_gamma + old_x_fe + old_B_y;
   ResponseFunction temp_old_y = old_y_gamma + old_y_fe + old_B_x;

   // Calculate correct inner products of upper off diagonal
   Tensor<double> off = expectation(world, x_response, temp_old_x) +
                        expectation(world, y_response, temp_old_y);

   // Use slicing to put in correct spot 
   temp_A(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);
   //temp_A(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);

   // Now for lower off diagonal 
   off = expectation(world, old_x_response, temp_cur_x) + 
         expectation(world, old_y_response, temp_cur_y);

   // Use slicing to put in correct spot 
   temp_A(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);
   //temp_A(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);

   // Put together the rest of A 
   temp_A(Slice(0, m-1), Slice(0, m-1)) = copy(A);
   temp_A(Slice(m, 2*m-1), Slice(m, 2*m-1)) = copy(old_A);

   // Save temp_A into A 
   A = copy(temp_A);

   // Now create upper off diagonal block of S
   off = expectation(world, x_response, old_x_response) -
         expectation(world, y_response, old_y_response); 

   // Use slicing to put in correct spot 
   temp_S(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);

   // Now the lower off diagonal block
   // (Cheating by using the transpose...)
   off = transpose(off);   

   // Use slicing to put in correct spot 
   temp_S(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);

   // Put together the rest of S
   temp_S(Slice(0, m-1), Slice(0, m-1)) = copy(S);
   temp_S(Slice(m, 2*m-1), Slice(m, 2*m-1)) = copy(old_S);

   // Save temp_S as S
   S = copy(temp_S);

   // Finally, add in old vectors to current vectors for the appropriate ones 
   for(int i = 0; i < m; i++)
   {
      x_response.push_back(old_x_response[i]);
      x_gamma.push_back(old_x_gamma[i]);
      V_x_response.push_back(old_V_x_response[i]);
      x_fe.push_back(old_x_fe[i]);
      B_x.push_back(old_B_x[i]);

      y_response.push_back(old_y_response[i]);
      y_gamma.push_back(old_y_gamma[i]);
      V_y_response.push_back(old_V_y_response[i]);
      y_fe.push_back(old_y_fe[i]);
      B_y.push_back(old_B_y[i]);
   }

   // End the timer
   if(print_level >= 1) end_timer(world, "Aug. resp. matrix:");

   // Debugging output 
   if(print_level >= 2 and world.rank() == 0) 
   {
      print("\n   Augmented response matrix:");
      print(A);
   }

   // Debugging output
   if(print_level >= 2 and world.rank() == 0) 
   {
      print("   Augmented overlap matrix:");
      print(S);
   }
}


// If using a larger subspace to diagonalize in, after diagonalization this will put everything in the right spot
void TDHF::unaugment(World & world,
                     int m,
                     int iter,
                     Tensor<double> & omega,
                     Tensor<double> & S_x,     
                     Tensor<double> & A_x,     
                     ResponseFunction & x_gamma,
                     ResponseFunction & x_response,
                     ResponseFunction & V_x_response,
                     ResponseFunction & x_fe,         // Contains fock and energy scaled orbitals
                     Tensor<double> & old_S, 
                     Tensor<double> & old_A, 
                     ResponseFunction & old_x_gamma, 
                     ResponseFunction & old_x_response, 
                     ResponseFunction & old_V_x_response, 
                     ResponseFunction & old_x_fe,
                     int print_level)        
{
   // Basic output
   if(print_level >= 1) start_timer(world);

   // Note: the eigenvalues and vectors were sorted after diagonalization
   // and hence all the functions are sorted in ascending order of energy

   // Quick copy of m lowest eigenvalues        
   omega = omega(Slice(0,m-1));

   // Pop off the "m" vectors off the back end of appropriate vectors
   // (only after first iteration)
   if(iter > 0)
   {
      for(int i = 0; i < m; i++)
      {
         x_fe.pop_back();
         V_x_response.pop_back();
         x_gamma.pop_back();
         x_response.pop_back();
      }
   }

   // Save the "current" into "old"
   old_x_fe = x_fe.copy();
   old_x_gamma = x_gamma.copy();
   old_V_x_response = V_x_response.copy();

   // Project out ground state  
   //QProjector<double, 3> projector(world, Gparams.orbitals);
   //for(int i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
   old_x_response = x_response.copy();

   // Copy S into old_S
   // S is identity
   old_S = expectation(world, x_response, x_response);

   // Copy A into old_A
   // A is (nearly?) diagonal
   old_A = Tensor<double>(m,m);
   for(int i = 0; i < m; i++) old_A(i,i) = omega(i);

   // End the timer
   if(print_level >= 1) end_timer(world, "Unaug. resp. mat.:");
}

// If using a larger subspace to diagonalize in, after diagonalization this will put everything in the right spot
void TDHF::unaugment_full(World & world,
                          int m,
                          int iter,
                          Tensor<double> & U,
                          Tensor<double> & omega,
                          Tensor<double> & S,     
                          Tensor<double> & A, 
                          ResponseFunction & x_gamma,
                          ResponseFunction & x_response,
                          ResponseFunction & V_x_response,
                          ResponseFunction & x_fe,
                          ResponseFunction & B_x, 
                          ResponseFunction & y_gamma,
                          ResponseFunction & y_response,
                          ResponseFunction & V_y_response,
                          ResponseFunction & y_fe,
                          ResponseFunction & B_y,  
                          Tensor<double> & old_S, 
                          Tensor<double> & old_A, 
                          ResponseFunction & old_x_gamma, 
                          ResponseFunction & old_x_response, 
                          ResponseFunction & old_V_x_response, 
                          ResponseFunction & old_x_fe,
                          ResponseFunction & old_B_x,
                          ResponseFunction & old_y_gamma, 
                          ResponseFunction & old_y_response, 
                          ResponseFunction & old_V_y_response, 
                          ResponseFunction & old_y_fe,
                          ResponseFunction & old_B_y,
                          int print_level)        
{
   // Basic output
   if(print_level >= 1) start_timer(world);

   // Note: the eigenvalues and vectors were sorted after diagonalization
   // and hence all the functions are sorted in ascending order of energy

   // Quick copy of m lowest eigenvalues        
   omega = omega(Slice(0,m-1));

   // Pop off the "m" vectors off the back end of appropriate vectors
   // (only after first iteration)
   if(iter > 0)
   {
      for(int i = 0; i < m; i++)
      {
         x_fe.pop_back();
         V_x_response.pop_back();
         x_gamma.pop_back();
         x_response.pop_back();
         B_x.pop_back();

         y_fe.pop_back();
         V_y_response.pop_back();
         y_gamma.pop_back();
         y_response.pop_back();
         B_y.pop_back();
      }
   }

   // Save the "current" into the "old"
   old_x_fe = x_fe.copy();
   old_x_gamma = x_gamma.copy();
   old_V_x_response = V_x_response.copy();
   old_B_x = B_x.copy();

   old_y_fe = y_fe.copy();
   old_y_gamma = y_gamma.copy();
   old_V_y_response = V_y_response.copy();
   old_B_y = B_y.copy();

   // Project out ground components
   //QProjector<double, 3> projector(world, Gparams.orbitals);
   //for(int i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
   //for(int i = 0; i < m; i++) y_response[i] = projector(y_response[i]);
   old_x_response = x_response.copy();
   old_y_response = y_response.copy();

   // Now to pull out correct values from S and A (both are size 2*m by 2*m, 
   // and only want m by m values)
   old_S = expectation(world, x_response, x_response) - expectation(world, y_response, y_response);  

   // And construct old_A
   ResponseFunction t1 = old_x_fe + old_x_gamma + old_B_y;
   ResponseFunction t2 = old_y_fe + old_y_gamma + old_B_x;
   old_A = expectation(world, old_x_response, t1) + expectation(world, old_y_response, t2); 
  
   // End the timer
   if(print_level >= 1) end_timer(world, "Unaug. resp. mat.:");
}


// Diagonalize the full response matrix, taking care of degenerate components
Tensor<double> TDHF::diag_full_response(World & world,
                                        Tensor<double> & S,
                                        Tensor<double> & A,
                                        ResponseFunction & x,
                                        ResponseFunction & Vx,
                                        ResponseFunction & x_g,
                                        ResponseFunction & x_fe,
                                        ResponseFunction & B_x,
                                        ResponseFunction & y,
                                        ResponseFunction & Vy,
                                        ResponseFunction & y_g,
                                        ResponseFunction & y_fe,
                                        ResponseFunction & B_y,
                                        Tensor<double> & omega, 
                                        const double thresh,
                                        int print_level)
{  
    // compute the unitary transformation matrix U that diagonalizes
    // the response matrix
    Tensor<double> vecs = get_full_response_transformation(world, S, A, omega, thresh);    

    // Start timer
    if(Rparams.print_level >= 1) start_timer(world);    

    // Transform the vectors of functions
    // Truncate happens in here
    Vx = transform(world, Vx, vecs);
    x_g = transform(world, x_g, vecs);
    x_fe = transform(world, x_fe, vecs);
    x = transform(world, x, vecs);
    B_x = transform(world, B_x, vecs);
    Vy = transform(world, Vy, vecs);
    y_g = transform(world, y_g, vecs);
    y_fe = transform(world, y_fe, vecs);
    y = transform(world, y, vecs);
    B_y = transform(world, B_y, vecs);

    // End timer
    if(Rparams.print_level >= 1) end_timer(world, "Transform orbs.:");

    // Normalize x and y 
    normalize(world, x, y);

    // Debugging output
    if(world.rank() == 0 and print_level >= 2)
    {
       print("   Eigenvector coefficients from diagonalization:");
       print(vecs);
    } 

    // Return the selected functions 
    return vecs;
}

// Similar to what robert did above in "get_fock_transformation"
Tensor<double> TDHF::get_full_response_transformation(World& world,
                                                      Tensor<double>& S,
                                                      Tensor<double>& A,
                                                      Tensor<double>& evals,
                                                      const double thresh_degenerate)    
{
    // Start timer
    if(Rparams.print_level >= 1) start_timer(world);
 
    // Get size
    int m = S.dim(0);

    // Run an SVD on the overlap matrix and ignore values
    // less than thresh_degenerate
    Tensor<double> r_vecs, s_vals, l_vecs;
    Tensor<double> S_copy = copy(S);
    svd(S_copy, l_vecs, s_vals, r_vecs);

    // Debugging output
    if(Rparams.print_level >= 2 and world.rank() == 0)
    {
       print("\n   Singular values of overlap matrix:");
       print(s_vals);
       print("   Left singular vectors of overlap matrix:");
       print(l_vecs);
    }

    // Check how many singular values are less than 10*thresh_degen
    int num_sv = 0;
    for(int i = 0; i < s_vals.dim(0); i++)
    {
       if(s_vals(i) < 10 * thresh_degenerate)
       {
          if(world.rank() == 0 and num_sv == 0) print("");
          if(world.rank() == 0) printf("   Detected singular value (%.8f) below threshold (%.8f). Reducing subspace size.\n", 
                                       s_vals(i), 10*thresh_degenerate);
          num_sv++;
       }
       if(world.rank() == 0 and i == s_vals.dim(0) - 1 and num_sv > 0) print("");
    }

    // Going to use these a lot here, so just calculate them
    int size_l = s_vals.dim(0);
    int size_s = size_l - num_sv;
    Tensor<double> l_vecs_s(size_l, num_sv);
    Tensor<double> copyA = copy(A);

    // Transform into this smaller space if necessary
    if(num_sv > 0)
    {
       // Cut out the singular values that are small
       // (singular values come out in descending order)
       S = Tensor<double>(size_s, size_s);
       for(int i = 0; i < size_s; i++) S(i,i) = s_vals(i); 

       // Copy the active vectors to a smaller container
       l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

       // Debugging output
       if(Rparams.print_level >= 2 and world.rank() == 0)
       {
          print("   Reduced size left singular vectors of overlap matrix:");
          print(l_vecs_s);
       }

       // Transform
       Tensor<double> work(size_l, size_s);
       mxm(size_l, size_s, size_l, work.ptr(), A.ptr(), l_vecs_s.ptr());
       copyA = Tensor<double>(size_s, size_s);
       Tensor<double> l_vecs_t = transpose(l_vecs);
       mxm(size_s, size_s, size_l, copyA.ptr(), l_vecs_t.ptr(), work.ptr());
    
       // Debugging output
       if(Rparams.print_level >= 2 and world.rank() == 0) 
       {
          print("   Reduced response matrix:");
          print(copyA);
          print("   Reduced overlap matrix:");
          print(S);
       }
    }
    // Diagonalize (NOT A SYMMETRIC DIAGONALIZATION!!!!)
    // Potentially complex eigenvalues come out of this
    Tensor<std::complex<double>> omega(size_s);
    Tensor<double> U(size_s, size_s);
    ggevp(world, copyA, S, U, omega);

    // Eigenvectors come out oddly packaged if there are 
    // complex eigenvalues.
    // Currently only supporting real valued eigenvalues
    // so throw an error if any imaginary components are 
    // not zero enough
    double max_imag = abs(imag(omega)).max();
    if(world.rank() == 0 and Rparams.print_level >= 2) print("\n   Max imaginary component of eigenvalues:", max_imag, "\n");
    MADNESS_ASSERT(max_imag == 0); // MUST BE REAL!
    evals = real(omega);

    // Easier to just resize here
    m = evals.dim(0);

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

    // Rotations between effectively degenerate components confound
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

    // If we transformed into the smaller subspace, time to transform back
    if(num_sv > 0)
    {
       // Temp. storage
       Tensor<double> temp_U(size_l, size_l);
       Tensor<double> U2(size_l, size_l);

       // Copy U back to larger size
       temp_U(Slice(0, size_s-1), Slice(0, size_s-1)) = copy(U);
       for(int i = size_s; i < size_l; i++) temp_U(i,i) = 1.0;

       // Transform U back
       mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());
       U = copy(U2);
    }

    // Sort into ascending order
    Tensor<int> selected = sort_eigenvalues(world, evals, U);

    // End timer
    if(Rparams.print_level >= 1) end_timer(world, "Diag. resp. mat.");

    return U;
}



// Sorts the given tensor of eigenvalues and 
// response functions 
void TDHF::sort(World & world,
                Tensor<double> & vals,
                ResponseFunction& f)
{
   // Get relevant sizes
   int k = vals.size();

   // Copy everything... 
   ResponseFunction f_copy(f); 
   Tensor<double> vals_copy = copy(vals);
   Tensor<double> vals_copy2 = copy(vals);

   // Now sort vals_copy
   std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

   // Now sort the rest of the things, using the sorted energy list
   // to find the correct indices 
   for(int i = 0; i < k; i++)
   {
      // Find matching index in sorted vals_copy
      int j = 0;
      while(fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

      // Put corresponding function, difference function, value residual and value
      // in the correct place
      f[i] = f_copy[j]; 
      vals(i) =  vals_copy(i);

      // Change the value of vals_copy2[j] to help deal with duplicates?
      vals_copy2(j) = 10000.0;
   }
}


// Creates the XCOperator object and initializes it with correct parameters
XCOperator TDHF::create_xcoperator(World& world,
                                   std::vector<real_function_3d> orbitals,
                                   std::string xc)
{
   // First calculate the ground state density
   std::vector<real_function_3d> vsq = square(world, Gparams.orbitals);
   compress(world, vsq);
   real_function_3d rho = real_factory_3d(world);
   rho.compress();
   for (unsigned int i = 0; i < vsq.size(); ++i) {
       rho.gaxpy(1.0, vsq[i], 1.0, false);
   }
   world.gop.fence();

   // And create the object using Gparams.xc
   XCOperator xcop(world, xc, false, rho, rho);

   return xcop;
} 

// Uses an XCOperator to construct v_xc for the ground state density 
// Returns d^2/d rho^2 E_xc[rho]
std::vector<real_function_3d> TDHF::create_fxc(World& world,
                                               std::vector<real_function_3d>& orbitals,
                                               ResponseFunction& f,
                                               ResponseFunction& g)
{
   // Create the xcop
   XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

   // Next need the perturbed density
   std::vector<real_function_3d> drho = transition_density(world, orbitals, f, g); 

   // Return container
   std::vector<real_function_3d> vxc;

   // Finally create the functions we want, one per response state 
   // (xc_args_prep_response happens inside this call)
   for(unsigned int i = 0; i < f.size(); i++)
   {
      vxc.push_back(xc.apply_xc_kernel(drho[i]));
   }

   return vxc;
}

// Iterates the response functions until converged or out of iterations
void TDHF::iterate(World & world)
{
   // Variables needed to iterate
   int iteration = 0;                                         // Iteration counter
   QProjector<double, 3> projector(world, Gparams.orbitals);  // Projector to project out ground state
   int n = Gparams.num_orbitals;                              // Number of ground state orbitals
   int m = Rparams.states;                                    // Number of excited states
   bool all_converged = false;                                // For convergence
   bool relax = false;                                        // For convergence
   int relax_start=Rparams.max_iter+1;                        // For convergence
   int num_conv=0;                                            // For convergence
   std::vector<bool>converged(m, false);                      // For convergence
   Tensor<double> old_energy(m);                              // Holds previous iteration's energy
   Tensor<double> energy_residuals;                           // Holds energy residuals 
   Tensor<double> x_norms(m);                                 // Holds the norms of x function residuals (for convergence)
   Tensor<double> y_norms(m);                                 // Holds the norms of y function residuals (for convergence)
   Tensor<double> x_shifts;                                   // Holds the shifted energy values
   Tensor<double> y_shifts;                                   // Holds the shifted energy values
   ResponseFunction bsh_x_resp;                               // Holds wave function corrections
   ResponseFunction bsh_y_resp;                               // Holds wave function corrections
   ResponseFunction x_differences;                            // Holds wave function corrections
   ResponseFunction y_differences;                            // Holds wave function corrections
   ResponseFunction x_gamma;                                  // Holds the perturbed two electron piece
   ResponseFunction y_gamma;                                  // Holds the perturbed two electron piece
   ResponseFunction x_fe;                                     // Holds the ground state-fock and energy scaled x response components 
   ResponseFunction y_fe;                                     // Holds the ground state-fock and energy scaled y response components 
   ResponseFunction V_x_response;                             // Holds V^0 applied to response functions
   ResponseFunction V_y_response;                             // Holds V^0 applied to response functions
   ResponseFunction B_x;                                      // Holds the off diagonal perturbed piece of y equation
   ResponseFunction B_y;                                      // Holds the off diagonal perturbed piece of x equation
   ResponseFunction shifted_V_x_response;                     // Holds the shifted V^0 applied to response functions
   ResponseFunction shifted_V_y_response;                     // Holds the shifted V^0 applied to response functions
   ResponseFunction old_x_response(world, m, n);              // Holds the old x_response vector of vectors
   ResponseFunction old_y_response(world, m, n);              // Holds the old y_response vector of vectors
   Tensor<double> S;                                          // Overlap matrix of response components for x states
   real_function_3d v_xc;                                     // For TDDFT

   // Versions from previous iteration that need to be stored
   // in order to diagonalize in a larger subspace
   ResponseFunction old_x_gamma; 
   ResponseFunction old_V_x_response;
   ResponseFunction old_x_fe;
   ResponseFunction old_B_x;
   ResponseFunction old_y_gamma;   
   ResponseFunction old_V_y_response;
   ResponseFunction old_y_fe;
   ResponseFunction old_B_y;
   Tensor<double> old_A;
   Tensor<double> old_S;

   // If DFT, initialize the XCOperator
   XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

   // The KAIN solver
   XNonlinearSolver<ResponseFunction, double, TDHF_allocator> kain(TDHF_allocator(world, Rparams.tda ? m : 2*m, n), false); 

   // Setting max sub size for KAIN solver
   if(Rparams.kain) kain.set_maxsub(Rparams.maxsub);

   // Set y things if not doing TDA
   if(not Rparams.tda) old_y_response = ResponseFunction(world, m, n);

   // Now to iterate
   while( iteration < Rparams.max_iter  && !all_converged)
   {
      // Start a timer for this iteration
      start_timer(world); 

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) printf("\n   Iteration %d at time %.1fs\n", iteration, wall_time());
         if(world.rank() == 0) print(" -------------------------------");
      }

      // Project out ground state 
      for(int i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
      if(not Rparams.tda) for(int i = 0; i < m; i++) y_response[i] = projector(y_response[i]);

      // Truncate before doing expensive things
      truncate(world, x_response);
      if(not Rparams.tda) truncate(world, y_response);

      // Normalize after projection
      if(Rparams.tda) normalize(world, x_response);
      else normalize(world, x_response, y_response);

      // Create gamma
      x_gamma = create_gamma(world, x_response, y_response, Gparams.orbitals, Rparams.small, 
                             FunctionDefaults<3>::get_thresh(), Rparams.print_level, "x");
      if(!Rparams.tda) y_gamma = create_gamma(world, y_response, x_response, Gparams.orbitals, Rparams.small, 
                                              FunctionDefaults<3>::get_thresh(), Rparams.print_level, "y");

      // Create \hat{V}^0 applied to response functions
      V_x_response = create_potential(world, x_response, xc, Rparams.print_level, "x");
      if(not Rparams.tda) V_y_response = create_potential(world, y_response, xc, Rparams.print_level, "y");

      // Load balance
      // Only balancing on x-components. Smart?
      if(world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0)) )
      {
         // Start a timer
         if(Rparams.print_level >= 1) start_timer(world); 
         if(world.rank() == 0) print(""); // Makes it more legible
 
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

         if(Rparams.print_level >= 1) end_timer(world, "Load balancing:");
      }

      // TDA approximation
      if(Rparams.tda)
      {
         // Constructing S
         S = expectation(world, x_response, x_response); 

         // Debugging output 
         if(Rparams.print_level >= 2 and world.rank() == 0) 
         {
            print("   Overlap matrix:");
            print(S);
         }
 
         // Constructing response matrix
         Tensor<double> A_x = create_response_matrix(world, x_fe, x_gamma, V_x_response, x_response, Gparams.orbitals, hamiltonian, Rparams.print_level, "x");

         // Augment S_x, A_x, x_gamma, x_response, V_x_response and x_gamma
         // if using a larger subspace and not iteration zero
         if(iteration < Rparams.larger_subspace and iteration > 0)
         {
            augment(world, S, A_x, x_gamma, x_response, V_x_response, x_fe, old_S, old_A, 
                    old_x_gamma, old_x_response, old_V_x_response, old_x_fe, Rparams.print_level);
         }

         // Solve Ax = Sxw
         // Just to be sure dimensions work out, clear omega
         omega.clear(); 
         diag_fock_matrix(world, A_x, x_response, V_x_response, x_gamma, x_fe, omega, S, 
                          FunctionDefaults<3>::get_thresh());

         // If larger subspace, need to "un-augment" everything 
         if(iteration < Rparams.larger_subspace)
         {  
            unaugment(world, m, iteration, omega, S, A_x, x_gamma, x_response, V_x_response, x_fe, old_S, 
                      old_A, old_x_gamma, old_x_response, old_V_x_response, old_x_fe, Rparams.print_level);          
         }
      }
      // Full TDHF
      else 
      {
         // Constructing S
         S = expectation(world, x_response, x_response) - expectation(world, y_response, y_response);

         // Debugging output 
         if(world.rank() == 0 and Rparams.print_level >= 2)
         {
            print("\n   Overlap Matrix:");
            print(S);
         }

         // Construct full response matrix 
         Tensor<double> A = create_full_response_matrix(world, x_gamma, V_x_response, B_x, x_fe, x_response, y_gamma, 
                                                        V_y_response, B_y, y_fe, y_response, Gparams.orbitals, hamiltonian, 
                                                        Rparams.small, FunctionDefaults<3>::get_thresh(), Rparams.print_level);

         // Larger subspace augmentation BROKEN!!!!!
         //if(iteration < Rparams.larger_subspace and iteration > 0)
         //{
         //   augment_full(world, S, A, 
         //                B_x, x_gamma, x_response, V_x_response, x_fe, 
         //                B_y, y_gamma, y_response, V_y_response, y_fe,
         //                old_S, old_A, 
         //                old_B_x, old_x_gamma, old_x_response, old_V_x_response, old_x_fe,
         //                old_B_y, old_y_gamma, old_y_response, old_V_y_response, old_y_fe,
         //                Rparams.print_level);
         //}

         // Diagonalize         
         // Just to be sure dimensions work out, clear omega
         omega.clear();
         Tensor<double> U = diag_full_response(world, S, A, 
                                               x_response, V_x_response, x_gamma, x_fe, B_x,
                                               y_response, V_y_response, y_gamma, y_fe, B_y,
                                               omega, FunctionDefaults<3>::get_thresh(), Rparams.print_level);

         // Larger subspace un-augmentation BROKEN!!!!
         //if(iteration < Rparams.larger_subspace)
         //{  
         //   unaugment_full(world, m, iteration, U, omega, S, A,
         //                  x_gamma, x_response, V_x_response, x_fe, B_x,
         //                  y_gamma, y_response, V_y_response, y_fe, B_y,
         //                  old_S, old_A,
         //                  old_x_gamma, old_x_response, old_V_x_response, old_x_fe, old_B_x,
         //                  old_y_gamma, old_y_response, old_V_y_response, old_y_fe, old_B_y,
         //                  Rparams.print_level);
         //}
      }

      // Basic output
      if(Rparams.print_level >= 1 and world.rank() == 0)
      {
         print("\n   Excitation Energies:");
         print(omega);       
      }

      // Calculate energy residual and update old_energy 
      energy_residuals = abs(omega - old_energy);
      old_energy = copy(omega);

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("   Energy residuals:");
         if(world.rank() == 0) print(energy_residuals);
      }

      // Analysis gets messed up if BSH is last thing applied
      // so exit early if last iteration
      if(iteration == Rparams.max_iter-1) 
      {
         end_timer(world," This iteration:");
         break;
      }

      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      x_shifts = create_shift(world, Gparams.energies, omega, Rparams.print_level, "x");
      if(not Rparams.tda) 
      {
         omega = -omega; // Negative here is so that these Greens functions are (eps - omega) 
         y_shifts = create_shift_target(world, Gparams.energies, omega, Gparams.energies[n-1], Rparams.print_level, "y");
         omega = -omega;
      }

      // Apply the shifts
      shifted_V_x_response = apply_shift(world, x_shifts, V_x_response, x_response);
      if(not Rparams.tda)
      {
         y_shifts = -y_shifts; 
         shifted_V_y_response = apply_shift(world, y_shifts, V_y_response, y_response);
      }

      // Construct RHS of equation
      ResponseFunction rhs_x = x_gamma + shifted_V_x_response;
      ResponseFunction rhs_y;
      if(not Rparams.tda)
      {
         // Add in coupling
         rhs_x = rhs_x + B_y;

         // And construct y
         rhs_y = y_gamma + shifted_V_y_response + B_x;
      }

      // Used to be localized orbital correction
      // but needs to happen in all cases
      ResponseFunction temp = scale_2d(world, x_response, ham_no_diag); 
      rhs_x = rhs_x - temp;

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms of off-diagonal hamiltonian correction for x components:");
         print_norms(world, temp);
      }

      if(not Rparams.tda)
      {
         temp = scale_2d(world, y_response, ham_no_diag);
         rhs_y = rhs_y - temp; 

         // Debugging output
         if(Rparams.print_level >= 2)
         {
            if(world.rank() == 0) print("   Norms of off-diagonal hamiltonian correction for y components:");
            print_norms(world, temp);
         }
      }      

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms of RHS of main equation:");
         if(world.rank() == 0) print("   x components:");
         print_norms(world, rhs_x);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y components:");
            print_norms(world, rhs_y);
         }
      }

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_operators = create_bsh_operators(world, x_shifts, Gparams.energies, omega, Rparams.small, FunctionDefaults<3>::get_thresh());
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_y_operators; 
      if(not Rparams.tda) 
      {
         omega = -omega;
         bsh_y_operators = create_bsh_operators(world, y_shifts, Gparams.energies, omega, Rparams.small, FunctionDefaults<3>::get_thresh());
         omega = -omega;
      }     

      // Save current into old
      old_x_response = x_response.copy();
      if(not Rparams.tda) old_y_response = y_response.copy();

      // Scale by -2.0 (coefficient in eq. 37 of reference paper)
      rhs_x = scale(rhs_x, -2.0);
      if(not Rparams.tda) rhs_y = scale(rhs_y, -2.0);     

      // Apply BSH and get updated response components 
      if(Rparams.print_level >= 1) start_timer(world);
      bsh_x_resp = apply(world, bsh_x_operators, rhs_x);
      if(not Rparams.tda) bsh_y_resp  = apply(world, bsh_y_operators, rhs_y);
      if(Rparams.print_level >= 1) end_timer(world, "Apply BSH:");

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms after application of BSH");
         if(world.rank() == 0) print("   x-components:");
         print_norms(world, bsh_x_resp);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y-components:");
            print_norms(world, bsh_y_resp);
         }
      }

      // Project out ground state
      for(int i = 0; i < m; i++) bsh_x_resp[i] = projector(bsh_x_resp[i]);
      if(not Rparams.tda) for(int i = 0; i < m; i++) bsh_y_resp[i] = projector(bsh_y_resp[i]);

      // Only update non-converged components
      for(int i = 0; i < m; i++)
      {
         if(not converged[i])
         {
            x_response[i] = bsh_x_resp[i];
            if(not Rparams.tda) y_response[i] = bsh_y_resp[i];
         }
      }

      // Get the difference between old and new
      x_differences = old_x_response - x_response;
      if(not Rparams.tda) y_differences = old_y_response - y_response;

      // Next calculate 2-norm of these vectors of differences
      // Remember: the entire vector is one state
      for(int i = 0; i < m; i++) x_norms(i) = norm2(world, x_differences[i]);
      if(not Rparams.tda) for(int i = 0; i < m; i++) y_norms(i) = norm2(world, y_differences[i]);

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) print("\n   2-norm of response function residuals:");
         if(world.rank() == 0) print("   x components:");
         if(world.rank() == 0) print(x_norms);

         if(not Rparams.tda)
         {
            if(world.rank() == 0) print("   y components:");
            if(world.rank() == 0) print(y_norms);
         }
      }

      // KAIN solver update 
      // Returns next set of components
      // If not kain, save the new components
      if(Rparams.kain)
      {
         if(not Rparams.tda)
         {
            // Add y functions to bottom of x functions
            // (for KAIN)
            for(int i = 0; i < m; i++)
            {
               x_response.push_back(y_response[i]);
               x_differences.push_back(y_differences[i]);
            }
         }

         // Do timers here (and not inside kain)
         start_timer(world);
         x_response = kain.update(x_response, x_differences, FunctionDefaults<3>::get_thresh(), 3.0);
         end_timer(world, " KAIN update:");

         // Add new functions back into y and
         // reduce x size back to original
         if(not Rparams.tda)
         {
            for(int i = 0; i < m; i++) y_response[i] = x_response[m + i];
            for(int i = 0; i < m; i++)
            { 
               x_response.pop_back();          
               x_differences.pop_back();          
            }
         }
      }

      // Apply mask
      for(int i = 0; i < m; i++) x_response[i] = mask * x_response[i];
      if(not Rparams.tda)
      { 
         for(int i = 0; i < m; i++) y_response[i] = mask * y_response[i];
      }

      // Only checking on X components even for full as Y are so small
      if(not relax)
      {
         for(int i = 0; i < m; i++)
         {
            if(iteration >= 1 && not converged[i] && fabs(x_norms[i]) < Rparams.dconv)
            {
               converged[i] = true;
               num_conv++;
               if(world.rank() == 0) print("   Response function", i, " has converged. Freezing it.");
            }
         }

         // Check if relaxing needs to start 
         if(num_conv == m)
         {
            relax_start = iteration;
            relax = true;
            if(world.rank() == 0) print("   All components converged. Unfreezing all states for final relaxation.");

            num_conv = 0;
            for(int i = 0; i < m; i++)
            {
               converged[i] = false;
            }
         }
      }
      else
      { 
         // Relaxing
         // Run at least 2 iterations
         if(iteration >= relax_start + 2)
         {
            // Check each state again
            for(int i = 0; i < m; i++)
            {
               if(not converged[i] && fabs(x_norms[i]) < Rparams.dconv)
               {
                  converged[i] = true;
                  num_conv++;                  
               }
            }
            if(num_conv == m) all_converged = true;
         }
      }
      
      // Update counter
      iteration += 1;

      // Done with the iteration.. truncate
      truncate(world, x_response);
      if(not Rparams.tda) truncate(world, y_response);

      // Save
      if(Rparams.save)
      {
         start_timer(world);
         save(world);
         end_timer(world, "Saving:");
      }

      // Basic output
      if(Rparams.print_level >= 1) end_timer(world, " This iteration:");

// TESTING
// get transition density
//if(world.rank() == 0) print("Making density.");
//std::vector<real_function_3d> densities = transition_density(world, Gparams.orbitals, x_response, y_response); 
//// Doing line plots along each axis
//if(world.rank() == 0) print("\n\nStarting plots");
//coord_3d lo,hi;
//char plotname[500];
//double Lp = std::min(Gparams.L, 24.0);
//if(world.rank() == 0) print("x:");
//// x axis 
//lo[0] = 0.0; lo[1] = 0.0; lo[2] = 0.0;
//hi[0] =  Lp; hi[1] = 0.0; hi[2] = 0.0;
////// plot ground state
////sprintf(plotname, "plot_ground_x.plt");
////plot_line(plotname, 5001, lo, hi, Gparams.orbitals[0]);
////
//// plot each x_k^p and the density
//for(int i = 0; i < m; i++)
//{
//   sprintf(plotname, "plot_orbital_%d_%d_x%d.plt", FunctionDefaults<3>::get_k(), i, iteration-1);
//   plot_line(plotname, 5001, lo, hi, x_response[i][0]);  
//}
//
//if(world.rank() == 0) print("y:");
//// y axis
//lo[0] = 0.0; lo[1] = 0.0; lo[2] = 0.0;
//hi[0] = 0.0; hi[1] =  Lp; hi[2] = 0.0;
//// plot ground state
//sprintf(plotname, "plot_ground1_y.plt");
//plot_line(plotname, 5001, lo, hi, Gparams.orbitals[0]);
//
//// plot each x_k^p and the density
//for(int i = 0; i < m; i++)
//{
//   sprintf(plotname, "plot_orbital_%d_%d_y%d.plt", FunctionDefaults<3>::get_k(), i, iteration-1);
//   plot_line(plotname, 5001, lo, hi, x_response[i][0]);
//}

//if(world.rank() == 0) print("z:");
// z axis
//lo[0] = 0.0; lo[1] = 0.0; lo[2] = -Lp;
//hi[0] = 0.0; hi[1] = 0.0; hi[2] =  Lp;
// plot ground state
//sprintf(plotname, "plot_ground1_z.plt");
//plot_line(plotname, 5001, lo, hi, Gparams.orbitals[0]);

// plot each x_k^p and the density
//for(int i = 0; i < m; i++)
//{
//   for(int j = 0; j < n; j++)
//   {
//      sprintf(plotname, "plot_orbital_%d_%d_%d_z%d.plt", FunctionDefaults<3>::get_k(), i, j, iteration-1);
//      plot_line(plotname, 20001, lo, hi, x_response[i][j]);
//  }
//}
//world.gop.fence();
//
// END TESTING
   }

   if(world.rank() == 0) print("\n");
   if(world.rank() == 0) print("   Finished TDHF Calculation ");
   if(world.rank() == 0) print("   ------------------------");
   if(world.rank() == 0) print("\n");

   // Did we converge?
   if(iteration == Rparams.max_iter && not all_converged)
   {
      if(world.rank() == 0) print("   Failed to converge. Reason:");
      if(world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
      if(world.rank() == 0) print("    Running analysis on current values.\n");
   }

   // Sort
   sort(world, omega, x_response);

   // Print final things 
   if(world.rank() == 0) 
   {
      print(" Final excitation energies:");
      print(omega);
      print(" Final energy residuals:");
      print(energy_residuals);
      print(" Final x-state response function residuals:");
      print(x_norms);
   
      if(not Rparams.tda)
      {
         if(world.rank() == 0) print(" Final y-state response function residuals:");
         if(world.rank() == 0) print(y_norms);
      }
   }

   // A little more detailed analysis
   analysis(world);

// TEST
// Doing line plots along z axis
//if(world.rank() == 0) print("\n\nStarting plots");
//coord_3d lo,hi;
//char plotname[500];
//// z axis 
//lo[0] = 0.0; lo[1] = 0.0; lo[2] = -Gparams.L;
//hi[0] = 0.0; hi[1] = 0.0; hi[2] =  Gparams.L;
//for(int i = 0; i < Rparams.states; i++) {
//  for(unsigned int j = 0; j < Gparams.num_orbitals; j++) {
//    sprintf(plotname, "plot_exX_%d_%d.plt", i, j);
//    plot_line(plotname, 500001, lo, hi, x_response[i][j]);  
//    sprintf(plotname, "plot_exY_%d_%d.plt", i, j);
//    plot_line(plotname, 500001, lo, hi, y_response[i][j]);  
//  }
//}
// END TEST


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

   // 'sort' these inner products within in each row
   Tensor<double> cpy = copy(x_norms);
   Tensor<int> x_order(m,n); Tensor<int> y_order(m,n);
   for(int i = 0; i < m; i++)
   {
      for(int j = 0; j < n; j++)
      {
         double x = cpy(i,_).max();
         int z = 0;
         while(x != cpy(i,z)) z++;
         cpy(i,z) = -100.0;
         x_order(i,j) = z;

         // Also sort y if full response
         if(not Rparams.tda)
         {
            y_order(i,j) = z;
         }
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
      
      // Normalization (negative?)
      dipoles(i,0) *= -sqrt(2.0);
      dipoles(i,1) *= -sqrt(2.0);
      dipoles(i,2) *= -sqrt(2.0);
   }

   // Calculate oscillator strength
   Tensor<double> oscillator(m);
   for(int i = 0; i < m; i++)
   {
      oscillator(i) = 2.0/3.0 * (dipoles(i,0)*dipoles(i,0) + dipoles(i,1)*dipoles(i,1) + dipoles(i,2)*dipoles(i,2)) * omega(i);
   }

   // Calculate transition quadrapole moments
   Tensor<double> quadrupoles(m,3,3);

   // Run over each excited state 
   for(int i = 0; i < m; i++)
   {
      // Add in contribution from each ground state
      for(int j = 0; j < n; j++)
      {
         quadrupoles(i,0,0) += inner(Gparams.orbitals[j], x * x * x_response[i][j]);
         quadrupoles(i,0,1) += inner(Gparams.orbitals[j], x * y * x_response[i][j]);
         quadrupoles(i,0,2) += inner(Gparams.orbitals[j], x * z * x_response[i][j]);
         quadrupoles(i,1,0) += inner(Gparams.orbitals[j], y * x * x_response[i][j]);
         quadrupoles(i,1,1) += inner(Gparams.orbitals[j], y * y * x_response[i][j]);
         quadrupoles(i,1,2) += inner(Gparams.orbitals[j], y * z * x_response[i][j]);
         quadrupoles(i,2,0) += inner(Gparams.orbitals[j], z * x * x_response[i][j]);
         quadrupoles(i,2,1) += inner(Gparams.orbitals[j], z * y * x_response[i][j]);
         quadrupoles(i,2,2) += inner(Gparams.orbitals[j], z * z * x_response[i][j]);

         if(not Rparams.tda)
         {
            quadrupoles(i,0,0) += inner(Gparams.orbitals[j], x * x * y_response[i][j]);
            quadrupoles(i,0,1) += inner(Gparams.orbitals[j], x * y * y_response[i][j]);
            quadrupoles(i,0,2) += inner(Gparams.orbitals[j], x * z * y_response[i][j]);
            quadrupoles(i,1,0) += inner(Gparams.orbitals[j], y * x * y_response[i][j]);
            quadrupoles(i,1,1) += inner(Gparams.orbitals[j], y * y * y_response[i][j]);
            quadrupoles(i,1,2) += inner(Gparams.orbitals[j], y * z * y_response[i][j]);
            quadrupoles(i,2,0) += inner(Gparams.orbitals[j], z * x * y_response[i][j]);
            quadrupoles(i,2,1) += inner(Gparams.orbitals[j], z * y * y_response[i][j]);
            quadrupoles(i,2,2) += inner(Gparams.orbitals[j], z * z * y_response[i][j]);
         }
      }
      // Normalization
      quadrupoles(i,0,0) *= sqrt(2.0);
      quadrupoles(i,0,1) *= sqrt(2.0);
      quadrupoles(i,0,2) *= sqrt(2.0);
      quadrupoles(i,1,0) *= sqrt(2.0);
      quadrupoles(i,1,1) *= sqrt(2.0);
      quadrupoles(i,1,2) *= sqrt(2.0);
      quadrupoles(i,2,0) *= sqrt(2.0);
      quadrupoles(i,2,1) *= sqrt(2.0);
      quadrupoles(i,2,2) *= sqrt(2.0);
   }

   // Now print?
   if(world.rank() == 0)
   {
      for(int i = 0; i < m; i++)
      {
         printf("   Response Function %d\t\t%7.8f a.u.", i, omega(i));
         print ("\n   --------------------------------------------");

         print("\n   Transition Dipole Moments");
         printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n", dipoles(i,0), dipoles(i,1), dipoles(i,2));

         printf("\n   Dipole Oscillator Strength: %7.8f\n", oscillator(i));

         print("\n   Transition Quadrupole Moments");
         printf("   %16s %16s %16s\n", "X", "Y", "Z");
         printf("   X %16.8f %16.8f %16.8f\n", quadrupoles(i,0,0), quadrupoles(i,0,1), quadrupoles(i,0,2));
         printf("   Y %16.8f %16.8f %16.8f\n", quadrupoles(i,1,0), quadrupoles(i,1,1), quadrupoles(i,1,2));
         printf("   Z %16.8f %16.8f %16.8f\n", quadrupoles(i,2,0), quadrupoles(i,2,1), quadrupoles(i,2,2));

         // Print contributions
         // Only print the top 5? 
         if(Rparams.tda)
         {
            print("\n   Dominant Contributions:");
            for(int j = 0; j < std::min(5,n); j++)
            {
               printf("   Occupied %d   %7.8f\n", x_order(i,j), x_norms(i,x_order(i,j)));
            }

            print("\n");
         }
         else
         {
            print("\n   Dominant Contributions:");
            print("                  x          y");
            for(int j = 0; j < std::min(5,n); j++)
            {
               printf("   Occupied %d   %7.8f %7.8f\n", x_order(i,j), x_norms(i,x_order(i,j)), y_norms(i,y_order(i,j)));
            }

            print("\n");

         }
      }
   }
}

// Simplified iterate scheme for guesses
void TDHF::iterate_guess(World & world,
                         ResponseFunction & guesses)
{
   // Variables needed to iterate
   int iteration = 0;                                         // Iteration counter
   QProjector<double, 3> projector(world, Gparams.orbitals);  // Projector to project out ground state
   int n = guesses[0].size();                                 // Number of ground state orbitals
   int m = guesses.size();;                                   // Number of excited components
   Tensor<double> shifts;                                     // Holds the shifted energy values
   ResponseFunction bsh_resp;      // Holds wave function corrections
   ResponseFunction gamma;         // Holds the perturbed two electron piece
   ResponseFunction fe;            // Holds the ground state-fock and energy scaled x response components 
   ResponseFunction V;             // Holds V^0 applied to response functions
   ResponseFunction shifted_V;     // Holds the shifted V^0 applied to response functions
   Tensor<double> S;               // Overlap matrix of response components for x states
   real_function_3d v_xc;          // For TDDFT

   // If DFT, initialize the XCOperator
   XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

   // Useful to have
   ResponseFunction zeros(world, m, n);

   // Now to iterate
   while( iteration < Rparams.guess_max_iter)
   {
      // Start a timer for this iteration
      start_timer(world); 

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) printf("\n   Guess Iteration %d at time %.1fs\n", iteration, wall_time());
         if(world.rank() == 0) print(" -------------------------------------");
      }

      // Load balance
      // Only balancing on x-components. Smart?
      if(world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0)) and iteration != 0)
      {
         // Start a timer
         if(Rparams.print_level >= 1) start_timer(world); 
         if(world.rank() == 0) print(""); // Makes it more legible
 
         LoadBalanceDeux<3> lb(world);
         for(int j = 0; j < n; j++)
         {
            for(int k = 0; k < Rparams.states; k++)
            {
               lb.add_tree(guesses[k][j], lbcost<double,3>(1.0,8.0),true);
               //lb.add_tree(V[k][j], lbcost<double,3>(1.0,8.0), true);
               //lb.add_tree(gamma[k][j], lbcost<double,3>(1.0,8.0), true);
            }
         }
         FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

         if(Rparams.print_level >= 1) end_timer(world, "Load balancing:");
      }

      // Project out ground state 
      for(int i = 0; i < m; i++) guesses[i] = projector(guesses[i]);

      // Truncate before doing expensive things
      truncate(world, guesses);

      // Normalize after projection
      if(Rparams.tda) normalize(world, guesses);

      // Create gamma
      gamma = create_gamma(world, guesses, zeros, Gparams.orbitals, Rparams.small, 
                           FunctionDefaults<3>::get_thresh(), Rparams.print_level, "x");

      // Create \hat{V}^0 applied to response functions
      V = create_potential(world, guesses, xc, Rparams.print_level, "x");

      // Constructing S
      S = expectation(world, x_response, x_response); 

      // Debugging output
      if(Rparams.print_level >= 2 and world.rank() == 0) 
      {
         print("\n   Overlap matrix:");
         print(S);
      }

      // Constructing response matrix
      Tensor<double> A = create_response_matrix(world, fe, gamma, V, guesses, Gparams.orbitals, hamiltonian, Rparams.print_level, "x");
      diag_fock_matrix(world, A, guesses, V, gamma, fe, omega, S, FunctionDefaults<3>::get_thresh());

      // Ensure right number of omegas
      if(omega.dim(0) != m)
      {
         if(world.rank() == 0) print("\n   Adding", m - omega.dim(0), "eigenvalue(s) (counters subspace size reduction in diagonalizatoin).");
         Tensor<double> temp(m);
         temp(Slice(0,omega.dim(0)-1)) = omega;
         for(int i = omega.dim(0); i < m; i++) temp[i] = 2.5 * i;
         omega = copy(temp);
      }

      // Basic output
      if(Rparams.print_level >= 1 and world.rank() == 0)
      {
         print("\n   Excitation Energies:");
         print(omega);       
      }

      // Only do BSH if not the last iteration
      if (iteration+1 < Rparams.guess_max_iter)
      {
         //  Calculates shifts needed for potential / energies
         //  If none needed, the zero tensor is returned
         shifts = create_shift(world, Gparams.energies, omega, Rparams.print_level, "x");

         // Apply the shifts
         shifted_V = apply_shift(world, shifts, V, guesses);

         // Construct RHS of equation
         ResponseFunction rhs = gamma + shifted_V;

         // Add in all off diagonal elements of ground state Fock matrix 
         ResponseFunction temp = scale_2d(world, guesses, ham_no_diag); 
         rhs = rhs - temp;

         // Construct BSH operators
         std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_operators = create_bsh_operators(world, shifts, Gparams.energies, omega, Rparams.small, FunctionDefaults<3>::get_thresh());

         // Scale by -2.0 (coefficient in eq. 37 of reference paper)
         rhs = scale(rhs, -2.0);

         // Apply BSH and get updated components 
         if(Rparams.print_level >= 1) start_timer(world);
         bsh_resp = apply(world, bsh_operators, rhs);
         if(Rparams.print_level >= 1) end_timer(world, "Apply BSH:");

         // Project out ground state
         for(int i = 0; i < m; i++) bsh_resp[i] = projector(bsh_resp[i]);

         // Save new components
         guesses = bsh_resp;

         // Apply mask
         for(int i = 0; i < m; i++) guesses[i] = mask * guesses[i];
      }     

      // Update counter
      iteration += 1;

      // Done with the iteration.. truncate
      truncate(world, guesses);

      // Basic output
      if(Rparams.print_level >= 1)
      {
         end_timer(world, " This iteration:");
      }
   }
}   // Done with iterate guess


// Create and diagonalize the CIS matrix for improved initial guess 
ResponseFunction TDHF::diagonalize_CIS_guess(World & world,
                                             std::vector<real_function_3d> & virtuals,
                                             Tensor<double> & omega,
                                             std::vector<real_function_3d> & orbitals,
                                             Tensor<double> & energies,
                                             double small,
                                             double thresh,
                                             int print_level)
{
   // Projecter for removing ground state
   QProjector<double, 3> Q(world, orbitals);

   // Diagonalize under ground state hamiltonian a few times
   // Create overlap
   compress(world, virtuals);
   Tensor<double> S = matrix_inner(world, virtuals, virtuals);

   // Create Fock matrix
   // -1 suppresses output
   Tensor<double> Fmat = create_ground_hamiltonian(world, virtuals, -1);

   // Diagonalize
   Tensor<double> U, evals, dummy(virtuals.size());
   U = get_fock_transformation(world, S, Fmat, evals, thresh); 

   // Transform and truncate functions
   virtuals = madness::transform(world, virtuals, U);
   truncate(world, virtuals, thresh, false);

   // filter out any ground state functions that crept in
   std::vector<real_function_3d> true_virtuals;
   for(unsigned int a = 0; a < virtuals.size(); a++)
   {
      if(evals(a) > 0.0) true_virtuals.push_back(virtuals[a]);
   }

   // Make sure we still have functions
   if(true_virtuals.empty()) MADNESS_EXCEPTION("Virtuals are empty: Too much overlap with occupied orbitals", 1);

   // Saving new components
   virtuals = Q(true_virtuals);

   // Debugging output
   if(print_level >= 2 and world.rank() == 0) print("   Remaining virtuals:", virtuals.size()); 

   // Now make the CIS matrix (copied from Jakob)
   if(world.rank() == 0) print("   Forming CIS matrix for improved initial guess.");

   // Start timer
   if(Rparams.print_level >= 1) start_timer(world);

   int I = -1; // combined index from i and a, start is -1 so that initial value is 0
   int J = -1; // combined index from j and b, start is -1 so that initial value is 0

   const int m = virtuals.size();
   const int n = orbitals.size();

   Tensor<double> MCIS(m*n, m*n);
   real_convolution_3d op = CoulombOperator(world, small, thresh);

   for(int i = 0; i < n; i++)
   {
      const real_function_3d brai = orbitals[i];
      const std::vector<real_function_3d> igv = apply(world, op, virtuals * brai);
      const std::vector<real_function_3d> igm = apply(world, op, orbitals * brai);

      for(int a = 0; a < m; a++)
      {
         I++;
         J = -1;
         for(int j = 0; j < n; j++)
         {
            const real_function_3d braj = orbitals[j];

            for(int b = 0; b < m; b++)
            {
               J++;
               double diag_element = 0.0;

               if(i == j and a == b) diag_element = Fmat(a,a) - energies(i);

               MCIS(I,J) = diag_element + 2.0 * inner(braj * virtuals[b], igv[a]) - inner(virtuals[a] * virtuals[b], igm[j]);
            }
         }
      }
   }

   // End timer
   if(print_level >= 1) end_timer(world, "Form CIS matrix:"); 

   // Debugging output
   if(print_level >= 2 and world.rank() == 0)
   {
      print("   CIS matrix:");
      print(MCIS);
   }

   // Diagonalize CIS matrix
   syev(MCIS, U, evals);

   // Always print this?
   if(world.rank() == 0)
   {
      print("   Initial CIS matrix eigenvalues:");
      print(evals); 
   }

   // Now construct the initial guess
   ResponseFunction f(world, U.dim(0), n);
   omega = Tensor<double>(U.dim(0));

   I = -1;
   for(int i = 0; i < n; i++)
   {
      for(int a = 0; a < m; a++)
      {
         I++;
         J = -1;
         if(evals(I) < 0.0)
         {
            if(world.rank() == 0) print("   Skipping negative root:", evals(I));
            continue;
         }
         for(int j = 0; j < n; j++)
         {
            for(int b = 0; b < m; b++)
            {
               J++;
               f[I][j] += U(J,I) * virtuals[b];
               omega(I) = evals(I);
            }
         }
      }
   }

   truncate(world, f);

   // Done. Whew.
   return f;
}

// Adds in random noise to a vector of vector of functions
ResponseFunction TDHF::add_randomness(World & world,
                                      ResponseFunction & f,
                                      double magnitude)

{
   // Copy input functions
   ResponseFunction f_copy = f.copy();

   // Lambda function to add in noise
   auto noise = [](const Key<3> & key, Tensor<double> & x) mutable
   {
      Tensor<double> y(x.size());
      y.fillrandom();
      //y.scale(magnitude);
      y.scale(1e3);
      x = x + y;
      //x(0,0,0) += y(0,0,0)-0.5;
   };

   // Go through each function in f_copy and add in random noise
   for(unsigned int i = 0; i < f_copy.size(); i++)
   {
      for(unsigned int j = 0; j < f_copy[0].size(); j++)
      {
         // Add in random noise using rng and a the defined lambda function
         f_copy[i][j].unaryop(noise);
      }
      
      // Apply mask to get boundary condition right
      f_copy[i] = mask * f_copy[i];
   }

   // Done
   return f_copy;
}

// Creates the ground state hamiltonian from given functions f
Tensor<double> TDHF::create_ground_hamiltonian(World & world,
                                               std::vector<real_function_3d> f,
                                               int print_level)
{
   // Basic output
   if(print_level >= 1) start_timer(world); 

   // Get sizes
   int m = f.size();

   // Debugging
   if(print_level > 2)
   {
      Tensor<double> S = matrix_inner(world, f, f);
      if(world.rank() == 0) print("   Ground state overlap:");
      if(world.rank() == 0) print(S);
   }

   // Calculate T
   // Make the derivative operators in each direction
   real_derivative_3d Dx(world, 0);
   real_derivative_3d Dy(world, 1);
   real_derivative_3d Dz(world, 2);

   // Apply derivatives once, and take inner products
   // according to this formula (faster / less noise):
   //  < f | \nabla^2 | f > = - < \nabla f | \nabla f >
   reconstruct(world, f);
   std::vector<real_function_3d> fx = apply(world, Dx, f);
   std::vector<real_function_3d> fy = apply(world, Dy, f);
   std::vector<real_function_3d> fz = apply(world, Dz, f);
   compress(world, fx, false);
   compress(world, fy, false);
   compress(world, fz, false);
   world.gop.fence();

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
   real_function_3d v_nuc = manager.vnuclear();
   v_nuc.truncate();

   // V_coul next
   // This does not include final multiplication of each orbital 
   // 2 is from integrating out spin
   real_function_3d v_coul = 2.0 * coulomb(world);

   // Clear old stored potentials
   stored_v_coul.clear();
   stored_v_nuc.clear();

   // If storing potentials, save them here
   if(Rparams.store_potential)
   {
      stored_v_nuc = copy(v_nuc);
      stored_v_coul = copy(v_coul);
   }

   // Sum coulomb (pre multiplied) and v_nuc
   // v_nuc comes out negative from potential manager, so add it
   real_function_3d v = v_coul + v_nuc;

   // Apply V to f functions
   std::vector<real_function_3d> vf = v * f;

   // Clear stored_potential
   stored_potential.clear();

   // ALWAYS DO THIS FOR THE STORED POTENTIAL!!
   // exchange last
   // 'small memory' algorithm from SCF.cc 
   real_convolution_3d op = CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());
   std::vector<real_function_3d> Kf = zero_functions_compressed<double,3>(world, m);
   for(int i=0; i<m; ++i)
   {
      std::vector<real_function_3d> psif = mul_sparse(world, f[i], f, FunctionDefaults<3>::get_thresh());
      truncate(world, psif);
      psif = apply(world, op, psif);
      truncate(world, psif);

      // Save the potential here if we are saving it
      if(Rparams.store_potential)
      {
         stored_potential.push_back(psif);
      }

      psif = mul_sparse(world, f[i], psif, FunctionDefaults<3>::get_thresh());
      gaxpy(world, 1.0, Kf, 1.0, psif);
   }

   // Only use the exchange above if HF:
   Tensor<double> V;
   real_function_3d v_xc;

   if(Gparams.xc == "hf")
   { 
      // Construct V
      V = matrix_inner(world, f, vf) - matrix_inner(world, f, Kf);
   }
   else // DFT
   {
      XCOperator xcop = create_xcoperator(world, f, Gparams.xc);

      real_function_3d v_xc = xcop.make_xc_potential();
      v = v + v_xc;
      std::vector<real_function_3d> vf = v * f;
      if((*xcop.xc).hf_exchange_coefficient() > 0.0)
      { 
         // XCOperator has member variable xc, which is an xcfunctional
         // which has the hf_exchange_coeff we need here
         gaxpy(world, 1.0, vf, -(*xcop.xc).hf_exchange_coefficient(), Kf);
      }
      V = matrix_inner(world, f, vf);
   }

   // Now create the hamiltonian
   hamiltonian = T + V;

   // Save a matrix that is
   // (T+V) - Lambda * eye 
   // Copy hamiltonian and zero the diagonal 
   ham_no_diag = copy(hamiltonian); 
   for(int i = 0; i < m; i++) ham_no_diag(i,i) = 0.0; 

   // Debug output
   if(print_level >= 2 and world.rank() == 0)
   {
      print("   Ground state hamiltonian:");
      print(hamiltonian);
   }

   // End timer
   if(print_level >=1) end_timer(world, "   Create grnd ham:");

   return hamiltonian;
}

// Creates the transition density
std::vector<real_function_3d> TDHF::transition_density(World& world,
                                                       std::vector<real_function_3d>& orbitals,
                                                       ResponseFunction& x,
                                                       ResponseFunction& y)
{
   // Get sizes
   int m = x.size();
   int n = x[0].size();

   // Return container 
   std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

   // Run over virtual...
   for(int i = 0; i < m; i++)
   {
      // Run over occupied...
      for(int j = 0; j < n; j++)
      {
         // y functions are zero if TDA is active
         densities[i] = densities[i] + orbitals[j] * x[i][j] + orbitals[j] * y[i][j];
      }
   }

   // Done!
   return densities;
}

template<std::size_t NDIM>
void TDHF::set_protocol(World & world, double thresh)
{
   int k;
   // Allow for imprecise conversion of threshold
   if(thresh >= 0.9e-2)
       k = 4;
   else if(thresh >= 0.9e-4)
       k = 6;
   else if(thresh >= 0.9e-6)
       k = 8;
   else if(thresh >= 0.9e-8)
       k = 10;
   else
       k = 12;
   
   // k defaults to make sense with thresh, override by providing k in input file
   if (Rparams.k > 0)
   {
       FunctionDefaults<NDIM>::set_k(Rparams.k);
   } 
   else FunctionDefaults<NDIM>::set_k(k);

   // MolDFT sets all these, so copying
   FunctionDefaults<NDIM>::set_thresh(thresh);
   FunctionDefaults<NDIM>::set_refine(true);
   FunctionDefaults<NDIM>::set_initial_level(2);
   FunctionDefaults<NDIM>::set_autorefine(false);
   FunctionDefaults<NDIM>::set_apply_randomize(false);
   FunctionDefaults<NDIM>::set_project_randomize(false);
   GaussianConvolution1DCache<double>::map.clear();

   // dconv defaults to thresh*100, overrirde by providing dconv in input file
   if (Rparams.dconv_set == false)
   {
      Rparams.dconv = thresh*100;
   }

   // Basic print
   if(world.rank() == 0)
   {
       print("\nSolving NDIM=",NDIM," with thresh", thresh, "    k",
             FunctionDefaults<NDIM>::get_k(), "  dconv", std::max(thresh, Rparams.dconv), "\n");
   }
}


void TDHF::check_k(World& world, 
                   double thresh,
                   int k)
{
   // Boolean to redo ground hamiltonian calculation if
   // ground state orbitals change
   bool redo = false;

   // Verify ground state orbitals have correct k
   if(FunctionDefaults<3>::get_k() != Gparams.orbitals[0].k())
   {
      // Re-read orbitals from the archive (assuming
      // the archive has orbitals stored at a higher
      // k value than what was previously computed 
      // with)
      Gparams.read(world, Rparams.archive);
      reconstruct(world, Gparams.orbitals);

      // Reset correct k (its set in Gparams.read)
      FunctionDefaults<3>::set_k(k);

      // Project each ground state to correct k
      for(unsigned int i = 0; i < Gparams.orbitals.size(); i++)
         Gparams.orbitals[i] = project(Gparams.orbitals[i], FunctionDefaults<3>::get_k(), thresh, false);
      world.gop.fence();

      // Clean up a bit
      truncate(world, Gparams.orbitals);

      // Ground state orbitals changed, clear old hamiltonian
      redo = true; 
   }

   // Recalculate ground state hamiltonian here
   if(redo or !hamiltonian.has_data())
   {
      hamiltonian = create_ground_hamiltonian(world, Gparams.orbitals, Rparams.print_level); 
   }

   // If we stored the potential, check that too
   if(Rparams.store_potential)
   {
      if(FunctionDefaults<3>::get_k() != stored_potential[0][0].k())
      { 
         // Project the potential into correct k
         for(unsigned int i = 0; i < stored_potential.size(); i++)
         {
            reconstruct(world, stored_potential[i]);
            for(unsigned int j = 0; j < stored_potential[0].size(); j++)
               stored_potential[i][j] = project(stored_potential[i][j], FunctionDefaults<3>::get_k(), thresh, false);
            world.gop.fence();
         }     
      }
      if(FunctionDefaults<3>::get_k() != stored_v_coul.k()) stored_v_coul = project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
      if(FunctionDefaults<3>::get_k() != stored_v_nuc.k()) stored_v_nuc = project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
   }

   // Verify response functions have correct k
   if(x_response.size() != 0)
   {
      if(FunctionDefaults<3>::get_k() != x_response[0][0].k())
      {
         // Project all x components into correct k
         for(unsigned int i = 0; i < x_response.size(); i++)
         {
            reconstruct(world, x_response[i]);
            for(unsigned int j = 0; j < x_response[0].size(); j++)
               x_response[i][j] = project(x_response[i][j], FunctionDefaults<3>::get_k(), thresh, false);
            world.gop.fence();
         }
         truncate(world, x_response); 

         // Do same for y components if applicable
         // (Always do this, as y will be zero 
         //  and still used in doing DFT and TDA)
         // Project all y components into correct k
         for(unsigned int i = 0; i < y_response.size(); i++)
         {
            reconstruct(world, y_response[i]);
            for(unsigned int j = 0; j < y_response[0].size(); j++)
               y_response[i][j] = project(y_response[i][j], FunctionDefaults<3>::get_k(), thresh, false);
            world.gop.fence();
         }
         truncate(world, y_response);
         
      } 
   }

   // Don't forget the mask function as well
   if(FunctionDefaults<3>::get_k() != mask.k())
   {      
      mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
   }

   // Make sure everything is done before leaving
   world.gop.fence();
}

// Creates random guess functions semi-intelligently(?)
ResponseFunction TDHF::create_random_guess(World & world, 
                                           int m, 
                                           int n,
                                           std::vector<real_function_3d> & grounds,
                                           Molecule & molecule)
{
   // Basic output 
   if(world.rank() == 0) print("   Using a random guess for initial response functions.\n");

   // Create empty container and add in randomness
   ResponseFunction f(world, m, n); 
   f = add_randomness(world, f, 1e3);

   // Create and apply a centered gaussian on each atom so that the randomness is localized around the atoms
   real_function_3d gaus = real_factory_3d(world);
   for(auto atom : molecule.get_atoms())
   {
      real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new GaussianGuess<3>(atom.get_coords(), 0.01, std::vector<int>{0,0,0})));
      gaus = gaus + x;
   }
   f = f * gaus;

   // Project out groundstate from guesses
   QProjector<double, 3> projector(world, grounds);
   for(unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);
  
   // Normalize
   normalize(world, f);

   return f;
}

// Creates random guess functions semi-intelligently(?)
std::vector<real_function_3d> TDHF::create_random_guess(World & world, 
                                                        int m, 
                                                        std::vector<real_function_3d> & grounds,
                                                        Molecule & molecule)
{
   // Basic output 
   if(world.rank() == 0) print("   Using a random guess for initial response functions.");

   // Create empty container and add in randomness
   std::vector<real_function_3d> f = zero_functions_compressed<double, 3>(world, m);

   // Create and apply a centered gaussian on each atom so that the randomness is localized around the atoms
   real_function_3d gaus = real_factory_3d(world);
   for(auto atom : molecule.get_atoms())
   {
      real_function_3d x = real_factory_3d(world).functor(real_functor_3d(new GaussianGuess<3>(atom.get_coords(), 0.01, std::vector<int>{0,0,0})));
      gaus = gaus + x;
   }

   // Lambda function to add in noise
   auto lambda = [](const Key<3> & key, Tensor<double> & x) mutable
   {
      Tensor<double> y(x.size());
      y.fillrandom();
      y.scale(1e3);
      x = x + y;
   };

   // Go through each function in f_copy and add in random noise
   for(unsigned int i = 0; i < f.size(); i++)
   {
      // Add in random noise using rng and a the defined lambda function
      f[i].unaryop(lambda);
  
      // Apply mask to get boundary condition right
      f[i] = mask * f[i] * gaus;
   }

   // Project out groundstate from guesses
   QProjector<double, 3> projector(world, grounds);
   for(unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);
  
   // Normalize
   for(unsigned int i = 0; i < f.size(); i++)
   {
      double norm = f[i].norm2();
      f[i].scale(1.0/norm);
   } 

   return f;
}

// Creates an initial guess function from nwchem output files
ResponseFunction TDHF::create_nwchem_guess(World & world, 
                                           int m)
{
   // Basic output 
   if(world.rank() == 0) print("   Creating an initial guess from NWChem file", Rparams.nwchem);

   // Create empty containers
   ResponseFunction f;

   // Create the nwchem reader
   slymer::NWChem_Interface nwchem(Rparams.nwchem, std::cout);
   
   // For parallel runs, silencing all but 1 slymer instance
   if(world.rank() != 0) {
      std::ostream dev_null(nullptr);
      nwchem.err = dev_null;
   }              

   // Read in basis set
   nwchem.read(slymer::Properties::Basis);

   // Read in the molecular orbital coefficients, energies,
   // and occupancies
   nwchem.read(slymer::Properties::MOs | slymer::Properties::Energies | slymer::Properties::Occupancies);
        
   // Create the nwchem orbitals as madness functions
   std::vector<real_function_3d> temp1;
   for(auto basis : slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set))
   {
       // Get the center of gaussian as its special point
       std::vector<coord_3d> centers;
       coord_3d r;
       r[0] = basis.get().center[0]; r[1] = basis.get().center[1]; r[2] = basis.get().center[2];
       centers.push_back(r);

       // Now make the function
       temp1.push_back(FunctionFactory<double,3>(world).functor(std::shared_ptr<FunctionFunctorInterface<double,3>>(new slymer::Gaussian_Functor(basis.get(), centers)))); 

       // Let user know something is going on
       if(temp1.size() % 10 == 0 and world.rank() == 0) print("Created", temp1.size(), "functions.");
   } 
   if(world.rank() == 0) print("Finished creating", temp1.size(), "functions.");

   // Normalize ao's
   madness::normalize(world, temp1);
   
   // Transform ao's now
   std::vector<real_function_3d> temp = madness::transform(world, temp1, nwchem.MOs, FunctionDefaults<3>::get_thresh(), true); 

   // Now save the unoccupied orbitals
   std::vector<real_function_3d> temp2;
   int num_virt = 0;
   for(unsigned int i = 0; i < temp1.size(); i++)
   {
       if(nwchem.occupancies[i] == 0)
       {
           temp2.push_back(temp[i]);
           num_virt++;
       }
   }

   // Create as many vectors of functions as we can from these nwchem virtual orbitals 
   // putting 1 virtual orbital from nwchem per vector 
   for(int i = 0; i < std::max(m, num_virt); i++)
   {
       // Create the vector to add the new function to
      std::vector<real_function_3d> v1 = zero_functions_compressed<double, 3>(world, Gparams.orbitals.size());

      // Put the "new" function into the vector
      v1[i % v1.size()] = temp2[i];

      // Add vector to return container
      f.push_back(v1);

      // See if we've made enough functions
      if(int(f.size()) >= m) break;
   }
   if(world.rank() == 0) print("Created", f.size(), "guess functions from provided NWChem data.");

   // If not enough functions have been made, start adding symmetry adapted functions
   int n = f.size(); 

   // If still not enough functions have been made, add in random guesses
   if(n < m)
   {
      // Tell user the bad news
      if(world.rank() == 0) print("\n   Only", n, "guess functions were provided by augmenting NWChem functions.\n   Augmenting with random functions.");

      // Create the random guess
      ResponseFunction rand = create_random_guess(world, m-n, Gparams.num_orbitals, Gparams.orbitals, Gparams.molecule);

      // Add to vector of functions
      for(unsigned int i = 0; i < rand.size(); i++) f.push_back(rand[i]);
   }

   // Project out groundstate from guesses
   QProjector<double, 3> projector(world, Gparams.orbitals);
   for(unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);
  
   // Truncate and normalize
   truncate(world, f);
   normalize(world, f);

   return f;
}

// Creates potentials using the ResponsePotential object
// Potentials are modified in place
void TDHF::create_all_potentials(World & world,
                                 ResponseFunction& x,
                                 ResponseFunction& x_gamma,
                                 ResponseFunction& x_V0,
                                 ResponsePotential& potentials,
                                 int print_level)
{
   // Intermediaries
   ResponseFunction gammaJ, gammaK, groundJ, groundK;
   
   // Calc. coulomb like terms
   if(print_level >= 1) start_timer(world);
   potentials.coulomb_terms(x, gammaK, groundJ);
   if(print_level >= 1) end_timer(world, "Coulomb terms:");
   
   // Calc. exchange like terms
   if(print_level >= 1) start_timer(world);
   potentials.exchange_terms(x, gammaJ, groundK);
   if(print_level >= 1) end_timer(world, "Exchange terms:");
   
   // Assemble pieces together
   x_gamma = gammaJ - gammaK;
   x_V0 = groundJ - groundK;
   
   // Debugging output
   if(print_level >= 2)
   {
      // Coulomb
      if(world.rank() == 0) printf("   Coulomb Deriv matrix:\n");
      Tensor<double> temp = expectation(world, x, gammaJ);
      if(world.rank() == 0) print(temp);

      // Exchange or VXC
      if(Rparams.xc == "hf" and world.rank() == 0) printf("   Exchange Deriv matrix:\n");
      if(Rparams.xc != "hf" and world.rank() == 0) printf("   Negative of XC Deriv matrix:\n");
      temp = expectation(world, x, gammaK);      
      if(world.rank() == 0) print(temp);

      // Total Gamma
      if(world.rank() == 0) printf("   Gamma matrix:\n");
      temp = expectation(world, x, x_gamma);
      if(world.rank() == 0) print(temp);
   
      // Coulomb (ground) 
      if(world.rank() == 0) printf("   Coulomb + Nuclear potential matrix:\n");
      temp = expectation(world, x, groundJ);
      if(world.rank() == 0) print(temp);

      // Exchange or VXC (ground)
      if(Rparams.xc == "hf" and world.rank() == 0) printf("   Exchange potential matrix:\n");
      if(Rparams.xc != "hf" and world.rank() == 0) printf("   XC potential matrix:\n");
      temp = expectation(world, x, groundK);
      if(world.rank() == 0) print(temp);
      if(world.rank() == 0) printf("   Total Potential Energy matrix:\n");
      temp = expectation(world, x, x_V0);
      if(world.rank() == 0) print(temp);
   }
}

// Main function, makes sure everything happens in correcct order
// Solves for response components
void TDHF::solve(World & world)
{
   // Get start time
   start_timer(world); 

   // Plotting input orbitals
   if(Rparams.plot_initial)
   {
      if(world.rank() == 0) print("\n   Plotting ground state densities.\n");
      if(Rparams.plot_L > 0.0) do_vtk_plots(world, Rparams.plot_pts, Rparams.plot_L, 0, Gparams.num_orbitals, Gparams.molecule, square(world, Gparams.orbitals), "ground");
      else do_vtk_plots(world, Rparams.plot_pts, Gparams.L/2.0, 0, Gparams.num_orbitals, Gparams.molecule, square(world, Gparams.orbitals), "ground");
   }

   // Warm and fuzzy 
   if(world.rank() == 0)
   {
      print("\n\n     Response Calculation");
      print("   ------------------------");
   }

   // Ready to iterate! 
   for(unsigned int proto = 0; proto < Rparams.protocol_data.size(); proto++)
   {
      // Set defaults inside here
      set_protocol<3>(world, Rparams.protocol_data[proto]);

      // Do something to ensure all functions have same k value
      check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());

      // Create the active subspace (select which ground state orbitals to calculate excitations from)
      //if(Rparams.e_window) select_active_subspace(world);

      if(proto == 0)
      {
         if(Rparams.restart)
         {
            if(world.rank() == 0) print("   Restarting from file:", Rparams.restart_file);
            load(world, Rparams.restart_file);
            check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());
         }
         else
         {
            // Create trial functions by...
            // (Always creating (at least) twice the amount requested for initial diagonalization)
            if(world.rank() == 0)  print("\n   Creating trial functions.\n");
            if(Rparams.random)
            {
               // Random guess
               x_response = create_random_guess(world, 2*Rparams.states, Gparams.num_orbitals, Gparams.orbitals, Gparams.molecule);
            }
            else if(Rparams.nwchem != "")
            {
               // Virtual orbitals from NWChem
               x_response = create_nwchem_guess(world, 2*Rparams.states);
            }
            else
            {
               // Use a symmetry adapted operator on ground state functions
               x_response = create_trial_functions(world, 2*Rparams.states, Gparams.orbitals, Rparams.print_level);
            }

            // Load balance
            // Only balancing on x-components. Smart?
            if(world.size() > 1)
            {
               // Start a timer
               if(Rparams.print_level >= 1) start_timer(world); 
               if(world.rank() == 0) print(""); // Makes it more legible
 
               LoadBalanceDeux<3> lb(world);
               for(int j = 0; j < Rparams.states; j++)
               {
                  for(unsigned int k = 0; k < Gparams.num_orbitals; k++)
                  {
                     lb.add_tree(x_response[j][k], lbcost<double,3>(1.0,8.0),true);
                  }
               }
               for(unsigned int j = 0; j < Gparams.num_orbitals; j++)
               {
                  lb.add_tree(Gparams.orbitals[j], lbcost<double,3>(1.0,8.0),true);
               }
               FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

               if(Rparams.print_level >= 1) end_timer(world, "Load balancing:");
            }
            
            // Project out groundstate from guesses
            QProjector<double, 3> projector(world, Gparams.orbitals);
            for(unsigned int i = 0; i < x_response.size(); i++) x_response[i] = projector(x_response[i]);

            // Ensure orthogonal guesses
            for(int i = 0; i < 2; i++)
            {
               start_timer(world);
               // Orthog
               x_response = gram_schmidt(world, x_response);
               end_timer(world, "orthog");

               start_timer(world); 
               // Normalize
               normalize(world, x_response);
               end_timer(world, "normalize");
            } 

            // Diagonalize guess
            if(world.rank() == 0) print("\n   Iterating trial functions for an improved initial guess.\n");
            iterate_guess(world, x_response);

            // Sort
            sort(world, omega, x_response);

            // Basic output
            if(Rparams.print_level >= 1 and world.rank() == 0)
            {
               print("\n   Final initial guess excitation energies:");
               print(omega);
            }

            // Select lowest energy functions from guess   
            x_response = select_functions(world, x_response, omega, Rparams.states, Rparams.print_level);

            // Initial guess for y are zero functions
            y_response = ResponseFunction(world, Rparams.states, Gparams.num_orbitals);
         }
      }

      // Now actually ready to iterate...
      iterate(world);
   }   

   // Plot the response function if desired
   if(Rparams.plot)
   {
      // Need to get densities first
      std::vector<real_function_3d> densities = transition_density(world, Gparams.orbitals, x_response, y_response); 

      // For the instance where we don't plot all the orbitals
      std::vector<real_function_3d> plot_densities;
      for(int i : Rparams.plot_data)
      {
         plot_densities.push_back(densities[i]);
      }

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
   end_timer(world, "total:");
}

// Iterates the response functions until converged or out of iterations
void TDHF::iterate_polarizability(World & world,
                                  ResponseFunction &dipoles)
{
   // Variables needed to iterate
   int iteration = 0;                                         // Iteration counter
   QProjector<double, 3> projector(world, Gparams.orbitals);  // Projector to project out ground state
   int n = Gparams.num_orbitals;                              // Number of ground state orbitals
   int m = Rparams.states;                                    // Number of excited states
   Tensor<double> x_norms(m);                                 // Holds the norms of x function residuals (for convergence)
   Tensor<double> y_norms(m);                                 // Holds the norms of y function residuals (for convergence)
   Tensor<double> x_shifts(m);                                // Holds the shifted energy values
   Tensor<double> y_shifts(m);                                // Holds the shifted energy values
   ResponseFunction bsh_x_resp;                               // Holds wave function corrections
   ResponseFunction bsh_y_resp;                               // Holds wave function corrections
   ResponseFunction x_differences;                            // Holds wave function corrections
   ResponseFunction y_differences;                            // Holds wave function corrections
   ResponseFunction x_gamma;                                  // Holds the perturbed two electron piece
   ResponseFunction y_gamma;                                  // Holds the perturbed two electron piece
   ResponseFunction x_fe;                                     // Holds the ground state-fock and energy scaled x response components 
   ResponseFunction y_fe;                                     // Holds the ground state-fock and energy scaled y response components 
   ResponseFunction V_x_response;                             // Holds V^0 applied to response functions
   ResponseFunction V_y_response;                             // Holds V^0 applied to response functions
   ResponseFunction B_x;                                      // Holds the off diagonal perturbed piece of y equation
   ResponseFunction B_y;                                      // Holds the off diagonal perturbed piece of x equation
   ResponseFunction shifted_V_x_response;                     // Holds the shifted V^0 applied to response functions
   ResponseFunction shifted_V_y_response;                     // Holds the shifted V^0 applied to response functions
   ResponseFunction old_x_response;                           // Holds the old x_response vector of vectors
   ResponseFunction old_y_response;                           // Holds the old y_response vector of vectors
   real_function_3d v_xc;                                     // For TDDFT
   bool converged = false;                                    // Converged flag


   // If DFT, initialize the XCOperator
   XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

   // The KAIN solver
   XNonlinearSolver<ResponseFunction, double, TDHF_allocator> kain(TDHF_allocator(world, (Rparams.omega != 0.0) ? 2*m : m, n), false); 

   // Setting max sub size for KAIN solver
   if(Rparams.kain) kain.set_maxsub(Rparams.maxsub);

   // Set omega (its constant here,
   // and has only 1 entry for each axis)
   omega = Tensor<double>(3);
   omega = Rparams.omega; 

   // Verify if any shift is needed (NEEDS CHECKING)
   if((Gparams.energies[n-1] + Rparams.omega) > 0.0)
   {
      // Calculate minimum shift needed such that \eps + \omega + shift < 0 
      // for all \eps, \omega
      x_shifts = create_shift(world, Gparams.energies, omega, Rparams.print_level, "x");
      y_shifts = Gparams.energies[n-1] + Rparams.omega + 0.05;
   }

   // Construct BSH operators
   std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_operators = create_bsh_operators(world, x_shifts, Gparams.energies, omega, Rparams.small, FunctionDefaults<3>::get_thresh());
   std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_y_operators; 

   // Negate omega to make this next set of BSH operators \eps - omega
   if(Rparams.omega != 0.0)
   {
      omega = -omega;
      bsh_y_operators = create_bsh_operators(world, y_shifts, Gparams.energies, omega, Rparams.small, FunctionDefaults<3>::get_thresh());
   }

   // Now to iterate
   while(iteration < Rparams.max_iter and !converged)
   {
      // Start a timer for this iteration
      start_timer(world); 

      // Basic output
      if(Rparams.print_level >= 1)
      {
         if(world.rank() == 0) printf("\n   Iteration %d at time %.1fs\n", iteration, wall_time());
         if(world.rank() == 0) print(" -------------------------------");
      }

      // If omega = 0.0, x = y
      if(Rparams.omega == 0.0) y_response = x_response.copy();

      // Save current to old
      old_x_response = x_response.copy();
      if(Rparams.omega != 0.0) old_y_response = y_response.copy(); 
//      world.gop.fence(); // Norm calc. below sometimes hangs without this (?)

      // Get norms 
      for(int i = 0; i < m; i++) x_norms[i] = sqrt(inner(x_response[i], x_response[i]) - inner(y_response[i], y_response[i])); 

      // Scale x and y 
      Tensor<double> rec_norms(m);
      for(int i = 0; i < m; i++) rec_norms(i) = 1.0/std::max(1.0, x_norms(i));
      x_response.scale(rec_norms); y_response.scale(rec_norms);

      // Create gamma
      x_gamma = create_gamma(world, x_response, y_response, Gparams.orbitals, Rparams.small, 
                             FunctionDefaults<3>::get_thresh(), Rparams.print_level, "x");
      if(Rparams.omega != 0.0) y_gamma = create_gamma(world, y_response, x_response, Gparams.orbitals, Rparams.small, 
                                                      FunctionDefaults<3>::get_thresh(), Rparams.print_level, "y");

      // Create \hat{V}^0 applied to response functions
      V_x_response = create_potential(world, x_response, xc, Rparams.print_level, "x");
      if(Rparams.omega != 0.0) V_y_response = create_potential(world, y_response, xc, Rparams.print_level, "y");

      // Apply shift
      V_x_response = apply_shift(world, x_shifts, V_x_response, x_response);
      if(Rparams.omega != 0.0) V_y_response = apply_shift(world, y_shifts, V_y_response, y_response);

      // Create \epsilon applied to response functions
      x_fe = scale_2d(world, x_response, ham_no_diag); 
      if(Rparams.omega != 0.0) y_fe = scale_2d(world, y_response, ham_no_diag); 
      if(Rparams.print_level >= 2) 
      {
         Tensor<double> t = expectation(world, x_response, x_fe);
         if(world.rank() == 0)
         {
            print("   Energy scaled response orbitals for x components:");
            print(t);
         }

         if(Rparams.omega != 0.0)
         {
            t = expectation(world, y_response, y_fe);
            if(world.rank() == 0)
            {
               print("   Energy scaled response orbitals for y components:");
               print(t);
            }
         }
      }

      // Load balance
      // Only balancing on x-components. Smart?
      // Only balance on first two iterations or every 5th iteration
      if(world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0)) )
      {
         // Start a timer
         if(Rparams.print_level >= 1) start_timer(world); 
         if(world.rank() == 0) print(""); // Makes it more legible
 
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

         if(Rparams.print_level >= 1) end_timer(world, "Load balancing:");
      }

      // Calculate coupling terms
      B_y = create_B(world, y_response, x_response, Gparams.orbitals, Rparams.small, FunctionDefaults<3>::get_thresh(), Rparams.print_level);
      if(Rparams.omega != 0.0) B_x = create_B(world, x_response, y_response, Gparams.orbitals, Rparams.small, FunctionDefaults<3>::get_thresh(), Rparams.print_level);

      // Scale dipoles by same value
      ResponseFunction dip_copy(dipoles);
      dip_copy.scale(rec_norms);

      // Construct RHS of equation
      ResponseFunction rhs_x, rhs_y;
      rhs_x = V_x_response - x_fe + dip_copy + x_gamma + B_y;
      if(Rparams.omega != 0.0) rhs_y = V_y_response - y_fe + dip_copy + y_gamma + B_x; 

      // Project out ground state
      for(int i = 0; i < m; i++) rhs_x[i] = projector(rhs_x[i]);
      if(Rparams.omega != 0.0) for(int i = 0; i < m; i++) rhs_y[i] = projector(rhs_y[i]);

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms of RHS of main equation x components:");
         print_norms(world, rhs_x);

         if(Rparams.omega != 0.0)
         {
            if(world.rank() == 0) print("   Norms of RHS of main equation y components:");
            print_norms(world, rhs_y);
         }
      }

      // Apply BSH and get updated response components 
      if(Rparams.print_level >= 1) start_timer(world);
      bsh_x_resp = apply(world, bsh_x_operators, rhs_x);
      if(Rparams.omega != 0.0) bsh_y_resp  = apply(world, bsh_y_operators, rhs_y);
      if(Rparams.print_level >= 1) end_timer(world, "Apply BSH:");

      // Scale by -2.0 (coefficient in eq. 37 of reference paper)
      for(int i = 0; i < m; i++) bsh_x_resp[i] = bsh_x_resp[i] * (std::max(1.0, x_norms[i]) * -2.0); 
      if(Rparams.omega != 0.0) for(int i = 0; i < m; i++) bsh_y_resp[i] = bsh_y_resp[i] * (std::max(1.0, x_norms[i]) * -2.0); 

      // Debugging output
      if(Rparams.print_level >= 2)
      {
         if(world.rank() == 0) print("   Norms after application of BSH to x components:");
         print_norms(world, bsh_x_resp);

         if(Rparams.omega != 0.0)
         {
            if(world.rank() == 0) print("   Norms after application of BSH to y components:");
            print_norms(world, bsh_y_resp);
         }
      }

      // Update orbitals
      x_response = bsh_x_resp;
      if(Rparams.omega != 0.0) y_response = bsh_y_resp;

      // Get the difference between old and new
      x_differences = old_x_response - x_response;
      if(Rparams.omega != 0.0) y_differences = old_y_response - y_response;

      // Next calculate 2-norm of these vectors of differences
      // Remember: the entire vector is one state
      for(int i = 0; i < m; i++) x_norms(i) = norm2(world, x_differences[i]);
      if(Rparams.omega != 0.0) for(int i = 0; i < m; i++) y_norms(i) = norm2(world, y_differences[i]);

      // Basic output
      if(Rparams.print_level >= 1 and world.rank() == 0)
      {
         print("\n   2-norm of response function residuals of x components:");
         print(x_norms);
         
         if(Rparams.omega != 0.0)
         {
            print("   2-norm of response function residuals of y components:");
            print(y_norms);
         }
      }

      // Check convergence
      if(std::max(x_norms.absmax(), y_norms.absmax()) < Rparams.dconv and iteration > 0) 
      {
          if(Rparams.print_level >= 1) end_timer(world, "This iteration:");
          if(world.rank() == 0) print("\n   Converged!");
          converged = true;
          break;
      }

      // KAIN solver update 
      // Returns next set of components
      // If not kain, save the new components
      if(Rparams.kain)
      {
         if(Rparams.omega != 0.0)
         {
            // Add y functions to bottom of x functions
            // (for KAIN)
            for(int i = 0; i < m; i++)
            {
               x_response.push_back(y_response[i]);
               x_differences.push_back(y_differences[i]);
            }
         }

         start_timer(world);
         x_response = kain.update(x_response, x_differences, FunctionDefaults<3>::get_thresh(), 3.0);
         end_timer(world, " KAIN update:");

         if(Rparams.omega != 0.0)
         {
            // Add new functions back into y and
            // reduce x size back to original
            for(int i = 0; i < m; i++) y_response[i] = x_response[m + i];
            for(int i = 0; i < m; i++)
            { 
               x_response.pop_back();          
               x_differences.pop_back();          
            }
         }
      }

      // Apply mask
      for(int i = 0; i < m; i++) x_response[i] = mask * x_response[i];
      if(Rparams.omega != 0.0) for(int i = 0; i < m; i++) y_response[i] = mask * y_response[i];
 
      // Update counter
      iteration += 1;

      // Done with the iteration.. truncate
      truncate(world, x_response);
      if(Rparams.omega != 0.0) truncate(world, y_response);

      // Save
      if(Rparams.save)
      { 
         start_timer(world);
         save(world);
         end_timer(world, "Save:");
      }
      // Basic output
      if(Rparams.print_level >= 1) end_timer(world, " This iteration:");

   }
}  // Done with iterate_polarizability 

// Calculates polarizability according to
// alpha_ij(\omega) = -sum_{ directions } < x_j | r_i | 0 > + < 0 | r_i | y_j >
void TDHF::polarizability(World& world,
                          Tensor<double> polar)
{
   // Get transition density
   //std::vector<real_function_3d> rhos = transition_density(world, Gparams.orbitals, x_response, y_response); 
   std::vector<real_function_3d> rhos;
   if(Rparams.omega == 0) rhos = transition_density(world, Gparams.orbitals, x_response, x_response); 
   else rhos = transition_density(world, Gparams.orbitals, x_response, y_response); 

   // For each r_axis
   for(int axis = 0; axis < 3; axis++)
   {
      real_function_3d drho = rhos[axis];

      // Run over axis and calc.
      // the polarizability
      for(int i = 0; i < 3; i++)
      {
         // Create dipole operator in the 'i' direction
         std::vector<int> f(3,0);
         f[i] = 1;
         real_function_3d dip = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
         polar(axis, i) = -2.0 * dip.inner(drho);
      }
   }
}

// Main function, makes sure everything happens in correct order
// Solves for polarizability
void TDHF::solve_polarizability(World & world)
{
   // Get start time
   start_timer(world); 

   // Warm and fuzzy 
   if(world.rank() == 0)
   {
      print("\n\n    Response Calculation");
      print("   ------------------------");
   }

   // Create the polarizability tensor
   Tensor<double> polar_tensor(3,3);

   // Keep a copy of dipoles * MO (needed explicitly in eq.)
   ResponseFunction dipoles;

   // For each protocol
   for(unsigned int proto = 0; proto < Rparams.protocol_data.size(); proto++)
   {
      // Set defaults inside here
      set_protocol<3>(world, Rparams.protocol_data[proto]);

      // Do something to ensure all functions have same k value
      check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());

      // Create guesses if no response functions
      // If restarting, load here
      if(proto == 0)
      {
         if(Rparams.restart)
         {
            if(world.rank() == 0) print("   Initial guess from file:", Rparams.restart_file);
            load(world, Rparams.restart_file);
            check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());
         }
         else // Dipole guesses
         { 
            x_response = dipole_guess(world, Gparams.orbitals); 
            y_response = x_response.copy(); 
         }
      }

      // Set the dipoles (ground orbitals are probably
      // more accurate now, so recalc the dipoles)
      dipoles = dipole_guess(world, Gparams.orbitals); 

      // Now actually ready to iterate...
      iterate_polarizability(world, dipoles);
   }

   // Have response function, now calculate polarizability for this axis
   polarizability(world, polar_tensor); 

   // Final polarizability analysis
   //diagonalize
   Tensor<double> V, epolar;
   syev(polar_tensor, V, epolar);
   double Dpolar_average = 0.0;
   double Dpolar_iso = 0.0;
   for(unsigned int i=0; i<3; ++i) Dpolar_average = Dpolar_average + epolar[i];
   Dpolar_average = Dpolar_average /3.0;
   Dpolar_iso= sqrt(.5)*sqrt( std::pow(polar_tensor(0,0) -  polar_tensor(1,1),2) +
                     std::pow(polar_tensor(1,1) -  polar_tensor(2,2),2) +
                     std::pow(polar_tensor(2,2) -  polar_tensor(0,0),2));

   if (world.rank() == 0)
   {
       print("\nTotal Dynamic Polarizability Tensor");
       printf("\nFrequency  = %.6f a.u.\n\n", Rparams.omega);
       //printf("\nWavelength = %.6f a.u.\n\n", Rparams.omega * ???); 
       print(polar_tensor);
       printf("\tEigenvalues = ");
       printf("\t %.6f \t %.6f \t %.6f \n", epolar[0], epolar[1], epolar[2]);
       printf("\tIsotropic   = \t %.6f \n", Dpolar_average);
       printf("\tAnisotropic = \t %.6f \n", Dpolar_iso);
       printf("\n");
   }

   // Print total time
   // Precision is set to 10 coming in, drop it to 2
   std::cout.precision(2);
   std::cout << std::fixed;

   // Get start time
   end_timer(world, "total:");
} 
// End solve_polar


// Exactam eam
// Deuces
