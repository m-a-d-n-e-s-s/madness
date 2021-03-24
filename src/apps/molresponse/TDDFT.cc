/*
 * Copyright 2021 Adrian Hurtado
 *
 *

 *
 *   Written by: bsundahl and molresponse
 *   Date: A long time ago... and today
 *
 */

#include "TDDFT.h"

#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "NWChem.h"  // For nwchem interface
#include "Plot_VTK.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

// KAIN allocator for vectorfunctions
struct TDHF_allocator {
  // Member variables
  World& world;
  const size_t num_vir;
  const size_t num_occ;  // Constructor
  TDHF_allocator(World& world, const int num_vir, const int num_occ)
      : world(world), num_vir(num_vir), num_occ(num_occ) {}

  // Overloading () operator
  response_space operator()() {
    response_space f(world, num_vir, num_occ);

    return f;
  }

  // Copy constructor
  TDHF_allocator operator=(const TDHF_allocator& other) {
    TDHF_allocator tmp(world, other.num_occ, other.num_vir);
    return tmp;
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

  x = (x * x * (3. - 2. * x));
  return x;
}

static double mask3(const coord_3d& ruser) {
  coord_3d rsim;
  user_to_sim(ruser, rsim);
  double x = rsim[0], y = rsim[1], z = rsim[2];
  double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
  double rlo = 1.0 / lo;

  if (x < lo)
    result *= mask1(x * rlo);
  else if (x > hi)
    result *= mask1((1.0 - x) * rlo);
  if (y < lo)
    result *= mask1(y * rlo);
  else if (y > hi)
    result *= mask1((1.0 - y) * rlo);
  if (z < lo)
    result *= mask1(z * rlo);
  else if (z > hi)
    result *= mask1((1.0 - z) * rlo);

  return result;
}

// Collective constructor
TDDFT::TDDFT(World& world, const char* filename)
    : TDDFT(world,
            (world.rank() == 0 ? std::make_shared<std::ifstream>(filename)
                               : nullptr)) {}

// Constructor that actually does stuff
TDDFT::TDDFT(World& world, std::shared_ptr<std::istream> input) {
  // Start the timer
  molresponse::start_timer(world);

  // Try and open input file
  if (world.rank() == 0) {
    if (input->fail())
      MADNESS_EXCEPTION("Response failed to open input stream", 0);
    // Welcome user (future ASCII art of Robert goes here)
    print("\n   Preparing to solve the TDHF equations.\n");
    // Read input files
    Rparams.read(*input);
    // Print out what was read in
  }
  Gparams.read(world, Rparams.archive);
  if (world.rank() == 0) {
    Gparams.print_params();
    print_molecule(world, Gparams);
  }
  // if a proerty calculation set the number of states
  if (Rparams.property) {
    Rparams.SetNumberOfStates(Gparams.molecule);
  }

  // print params
  if (world.rank() == 0) {
    Rparams.print_params();
  }
  // Broadcast to all other nodes
  world.gop.broadcast_serializable(Rparams, 0);

  // Read in archive
  // Create the projector Qhat to be used in any calculation

  // Set some function defaults
  FunctionDefaults<3>::set_cubic_cell(-Gparams.L, Gparams.L);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // Initialize response state xcfuntion object
  xcf.initialize(Rparams.xc, false, world, true);

  // Create the masking function
  mask = real_function_3d(
      real_factory_3d(world).f(mask3).initial_level(4).norefine());

  if (world.size() > 1) {
    // Start a timer
    if (Rparams.print_level >= 1) molresponse::start_timer(world);
    if (world.rank() == 0) print("");  // Makes it more legible

    LoadBalanceDeux<3> lb(world);
    for (unsigned int j = 0; j < Gparams.num_orbitals; j++) {
      lb.add_tree(Gparams.orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
    }
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

    if (Rparams.print_level >= 1)
      molresponse::end_timer(world, "Load balancing:");
  }
}

// Constructor that actually does stuff
TDDFT::TDDFT(World& world,
             ResponseParameters rparams,
             GroundParameters gparams) {
  // Start the timer
  this->Rparams = rparams;
  this->Gparams = gparams;

  if (rparams.response_type.compare("excited_state") == 0) {
    this->omega = Tensor<double>(rparams.states);
  } else {
    this->omega = Tensor<double>(1);
  }

  molresponse::start_timer(world);

  // Broadcast to all other nodes
  world.gop.broadcast_serializable(Rparams, 0);

  // Read in archive
  // Create the projector Qhat to be used in any calculation

  // Set some function defaults
  FunctionDefaults<3>::set_cubic_cell(-Gparams.L, Gparams.L);
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // Initialize response state xcfuntion object
  xcf.initialize(Rparams.xc, false, world, true);

  // Create the masking function
  mask = real_function_3d(
      real_factory_3d(world).f(mask3).initial_level(4).norefine());

  if (world.size() > 1) {
    // Start a timer
    if (Rparams.print_level >= 1) molresponse::start_timer(world);
    if (world.rank() == 0) print("");  // Makes it more legible

    LoadBalanceDeux<3> lb(world);
    for (unsigned int j = 0; j < Gparams.num_orbitals; j++) {
      lb.add_tree(Gparams.orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
    }
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

    if (Rparams.print_level >= 1)
      molresponse::end_timer(world, "Load balancing:");
  }
}
response_space& TDDFT::GetResponseFunctions(std::string xy) {
  if (xy == "x") {
    return x_response;
  } else if (xy == "y") {
    return y_response;
  } else {
    MADNESS_EXCEPTION("not a valid response state", 0);
  }
}
response_space& TDDFT::GetPVector() { return P; }
response_space& TDDFT::GetQVector() { return Q; }
// Get response parameters
ResponseParameters TDDFT::GetResponseParameters() { return Rparams; }
GroundParameters TDDFT::GetGroundParameters() { return Gparams; }
Property TDDFT::GetPropertyObject() { return p; }
// Get Frequencies Omega
Tensor<double> TDDFT::GetFrequencyOmega() {
  print("Frequencies : ", omega);
  return omega;
}
// Save the current response calculation
void TDDFT::save(World& world, std::string name) {
  // Archive to write everything to
  archive::ParallelOutputArchive ar(world, name.c_str(), 1);
  // Just going to enforce 1 io server

  // Saving, in this order;
  //  string           ground-state archive name (garch_name)
  //  bool             TDA flag
  // size_t                number of ground state orbitals (n)
  // size_t                number of excited state orbitals (m)
  //  Tensor<double>   energies of m x-components
  //  for i from 0 to m-1
  //     for j from 0 to n-1
  //        Function<double,3> x_response[i][j]
  //  (If TDA flag == True)
  //  (Tensor<double>  energies of m y-components    )
  //  (for i from 0 to m-1                       )
  //  (   for j from 0 to n-1                    )
  //  (      Function<double,3> y_response[i][j] )
  ar& Gparams.inFile;
  ar& Rparams.tda;
  ar& Gparams.num_orbitals;
  ar& Rparams.states;
  ar& omega;

  for (size_t i = 0; i < Rparams.states; i++)
    for (size_t j = 0; j < Gparams.num_orbitals; j++) ar& x_response[i][j];
  if (not Rparams.tda) {
    for (size_t i = 0; i < Rparams.states; i++)
      for (size_t j = 0; j < Gparams.num_orbitals; j++) ar& y_response[i][j];
  }
}

// Load a response calculation
void TDDFT::load(World& world, std::string name) {
  // The archive to read from
  archive::ParallelInputArchive ar(world, name.c_str());

  // Reading in, in this order;
  //  string           ground-state archive name (garch_name)
  //  bool             TDA flag
  // size_t                number of ground state orbitals (n)
  // size_t                number of excited state orbitals (m)
  //  Tensor<double>   energies of m x-components
  //  for i from 0 to m-1
  //     for j from 0 to n-1
  //        Function<double,3> x_response[i][j]
  //  (If TDA flag == True)
  //  (Tensor<double>  energies of m y-components    )
  //  (for i from 0 to m-1                       )
  //  (   for j from 0 to n-1                    )
  //  (      Function<double,3> y_response[i][j] )

  ar& Rparams.archive;
  ar& Rparams.tda;
  ar& Gparams.num_orbitals;
  ar& Rparams.states;
  ar& omega;

  x_response = response_space(world, Rparams.states, Gparams.num_orbitals);

  for (size_t i = 0; i < Rparams.states; i++)
    for (size_t j = 0; j < Gparams.num_orbitals; j++) ar& x_response[i][j];
  world.gop.fence();

  y_response = response_space(world, Rparams.states, Gparams.num_orbitals);
  if (not Rparams.tda) {
    for (size_t i = 0; i < Rparams.states; i++)
      for (size_t j = 0; j < Gparams.num_orbitals; j++) ar& y_response[i][j];
    world.gop.fence();
  }
}

// (Each state's norm should be 1, not the
// individual functions norms)
void TDDFT::normalize(World& world, response_space& f) {
  // Run over rows
  for (unsigned int i = 0; i < f.size(); i++) {
    // Get the normalization constant
    // (Sum included inside inner)
    double norm = inner(f[i], f[i]);
    norm = sqrt(norm);
    // Doing this to deal with zero functions.
    // Maybe not smrt.
    if (norm == 0) continue;

    // And scale
    f[i] = f[i] * (1.0 / norm);
  }
}

// (Each state's norm should be 1, not the
// individual functions norms)
//  non-standard normalization for eigen value problem
void TDDFT::normalize(World& world, response_space& f, response_space& g) {
  // Run over rows
  for (size_t i = 0; i < f.size(); i++) {
    // Get the normalization constant
    // (Sum included inside inner)
    double normf = inner(f[i], f[i]);
    double normg = inner(g[i], g[i]);
    double norm = sqrt(normf - normg);

    // Doing this to deal with zero functions.
    // Maybe not smrt.
    if (norm == 0) continue;
    // And scale
    scale(world, f[i], (1.0 / norm));
    scale(world, g[i], (1.0 / norm));
  }
}

// Prints norms of the given vector of vector of functions
void TDDFT::print_norms(World& world, response_space f) {
  // Container
  Tensor<double> norms(f.size(), f[0].size());

  // Calc the norms
  for (unsigned int i = 0; i < f.size(); i++) {
    for (unsigned int j = 0; j < f[0].size(); j++) {
      norms(i, j) = f[i][j].norm2();
    }
  }

  // Print em in a smart way
  if (world.rank() == 0) print(norms);
}

// Small function to print geometry of a molecule nicely

// Radial function
static double kronecker(size_t l, size_t n) {
  if (l == n) return 1.0;
  return 0.0;
}
// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
std::map<std::vector<int>, real_function_3d> TDDFT::solid_harmonics(
    World& world,
    int n) {
  // Container to return
  std::map<std::vector<int>, real_function_3d> result;

  // Create the basic x, y, z, constant and zero
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));
  real_function_3d c = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 0})));
  real_function_3d zero = real_factory_3d(world);

  // Add in first few, since they're simple
  // Assuming n >= 1
  result[std::vector<int>{0, 0}] = copy(c);
  result[std::vector<int>{0, -1}] = zero;
  result[std::vector<int>{0, 1}] = zero;
  result[std::vector<int>{-1, 0}] = zero;

  // Generate the solid harmonics recursively from here
  for (int l = 0; l < n; l++) {
    // Calculate ends of this row first
    result[std::vector<int>{l + 1, l + 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (x * result[std::vector<int>{l, l}] -
         (1 - kronecker(l, 0) * y * result[std::vector<int>{l, -l}]));
    result[std::vector<int>{l + 1, -l - 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (y * result[std::vector<int>{l, l}] +
         (1 - kronecker(l, 0) * x * result[std::vector<int>{l, -l}]));

    // Formula below calls for some functions that don't exist.
    // Need zeroes where that would occur
    result[std::vector<int>{l + 1, l + 2}] = zero;
    result[std::vector<int>{l + 1, -l - 2}] = zero;

    // Run over quantum number m
    for (int m = -l; m < l + 1; m++) {
      // Calculate remaining terms
      result[std::vector<int>{l + 1, m}] =
          1.0 / std::sqrt((l + m + 1) * (l - m + 1)) *
          ((2 * l + 1) * z * result[std::vector<int>{l, m}] -
           sqrt((l + m) * (l - m)) * (x * x + y * y + z * z) *
               result[std::vector<int>{l - 1, m}]);
    }
  }

  // Get rid of any zero functions we added
  for (auto it = result.begin(); it != result.end();) {
    if (it->second.norm2() == 0)
      it = result.erase(it);
    else
      ++it;
  }

  // Also get rid of the constant
  result.erase(std::vector<int>{0, 0});

  // Done
  return result;
}

// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
std::map<std::vector<int>, real_function_3d> TDDFT::simple_spherical_harmonics(
    World& world,
    int n) {
  // Container to return
  std::map<std::vector<int>, real_function_3d> result;

  // Create the basic x, y, z, constant and zero
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));
  real_function_3d c = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 0})));
  real_function_3d zero = real_factory_3d(world);

  real_function_3d rfunc = (x * x + y * y + z * z);
  double r = rfunc.norm2();

  std::vector<real_function_3d> funcs;
  funcs[0] = c;

  funcs[1] = y.scale(1 / r);
  funcs[2] = z.scale(1 / r);
  funcs[3] = x.scale(1 / r);

  funcs[4] = x * y;
  funcs[5] = y * z;
  funcs[6] = 2 * z * z - x * x - y * y;
  funcs[7] = z * x;
  funcs[8] = x * x - y * y;

  funcs[4] = funcs[4].scale(1 / (r * r));
  funcs[5] = funcs[5].scale(1 / (r * r));
  funcs[6] = funcs[6].scale(1 / (r * r));
  funcs[7] = funcs[7].scale(1 / (r * r));
  funcs[8] = funcs[8].scale(1 / (r * r));

  size_t num = 0;
  for (int l = 0; l < n; l++) {
    for (int m = -l; m <= l; m++) {
      result[vector<int>{l, m}] = funcs[num];
      num += 1;
    }
  }

  // Done
  return result;
}

// Returns initial guess functions as
// ground MO * solid harmonics
response_space TDDFT::create_trial_functions(
    World& world,
    size_t k,
    std::vector<real_function_3d>& orbitals,
    size_t print_level) {
  // Get size
  // /
  size_t n = orbitals.size();

  // Create solid harmonics such that num. solids * num. orbitals > k.
  // The total number of solid harmonics that exist up to level n is
  // (n+1)^2 (because we count from zero)
  // Always do at least 8 (through the d orbital angular momentum functions,
  // minus )
  std::map<std::vector<int>, real_function_3d> solids =
      solid_harmonics(world, std::max(2.0, ceil(sqrt(k / n) - 1)));

  // Useful info.
  if (world.rank() == 0)
    print("   Created", solids.size(), "solid harmonics.\n");

  // Container to return
  response_space trials;

  // Counter for number of trials created
  size_t count = 0;

  // Multiply each solid harmonic onto a ground state orbital
  for (size_t i = 0; i < n; i++) {
    // For each solid harmonic
    for (auto key : solids) {
      // Temp zero functions
      std::vector<real_function_3d> temp =
          zero_functions_compressed<double, 3>(world, n);

      // Create one non-zero function and add to trials
      temp[count % n] = key.second * orbitals[n - count % n - 1];
      trials.push_back(temp);
      count++;
    }

    // Stop when we first get beyond k components
    if (count >= k) break;
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) print("   Norms of guess functions:");
    print_norms(world, trials);
  }

  // Truncate
  madness::truncate(
      world, trials, madness::FunctionDefaults<3>::get_thresh(), true);

  // Done
  return trials;
}

// Returns initial guess functions as
// ground MO * <x,y,z>
response_space TDDFT::create_trial_functions2(
    World& world,
    std::vector<real_function_3d>& orbitals,
    size_t print_level) {
  // Get size
  size_t n = orbitals.size();
  size_t directions = 3;
  // (n+1)^2 (because we count from zero)
  // adsf
  //
  std::vector<real_function_3d> xyz = createDipoleFunctionMap(world);
  // create 3 x n orbital functions
  std::vector<vector<real_function_3d>> functions;

  print("Debug Norms", xyz[0].norm2());
  print("Debug Norms", xyz[1].norm2());
  print("Debug Norms", xyz[2].norm2());

  for (size_t d = 0; d < directions; d++) {
    vector<real_function_3d> temp;
    for (size_t i = 0; i < n; i++) {
      // create x functions then y.. then z ..
      temp.push_back(orbitals[i] * xyz[d]);
      print("Debug Norms of temp ", i, "=", temp[i].norm2());
      print("Debug Norms of orbitals ", i, "=", orbitals[i].norm2());
    }
    // all the x then the y then the z
    functions.push_back(temp);
  }
  print("number of orbitals: ", orbitals.size());

  // Container to return
  size_t count = 0;

  response_space trials(world, 3 * n * n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t d = 0; d < directions; d++) {
      for (size_t o = 0; o < n; o++) {
        //        trials[i + j + o][o] = functions[i][j];
        trials[count][o] = copy(functions.at(d).at(o));
        count++;
      }
    }
  }
  // The above generates response function as follows
  // all functions start off as zeros
  // 1  [x1 0 0 ]
  // 2  [0 x1 0 ]
  // 3  [0 0 x1 ]
  // 4  [y1 0 0 ]
  // 5  [0 y1 0 ]
  // 6  [0 0 y1 ]
  // 7  [z1 0 0 ]
  // 8  [0 z1 0 ]
  // 9  [0 0 z1 ]
  // 10 [x2 0 0 ]
  // 11 [0 x2 0 ]
  // 12 [0 0 x2 ]
  // 13 [y2 0 0 ]
  // 14 [0 y2 0 ]
  // 15 [0 0 y2 ]
  // 16 [z2 0 0 ]
  // 17 [0 z2 0 ]
  // 18 [0 0 z2 ]
  // 19 [x3 0 0 ]
  // 20 [0 x3 0 ]
  // 21 [0 0 x3 ]
  // 22 [y3 0 0 ]
  // 23 [0 y3 0 ]
  // 24 [0 0 y3 ]
  // 25 [z3 0 0 ]
  // 26 [0 z3 0 ]
  // 27 [0 0 z3 ]
  // for each orbital for each direction
  // Counter for number of trials created
  // Multiply each solid harmonic onto a ground state orbital
  //  for each orbital we

  // For each solid harmonic
  // Temp zero functions

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) print("   Norms of guess functions:");
    print_norms(world, trials);
  }

  // Truncate
  madness::truncate(world, trials);

  // Done
  return trials;
}

// Returns a list of solid harmonics such that:
// solid_harm.size() * num_ground_orbs > 2 * num. resp. components
std::vector<real_function_3d> TDDFT::createDipoleFunctionMap(World& world) {
  // Container to return

  // Create the basic x, y, z, constant and zero
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));

  //  real_function_3d rfunc = (x * x + y * y + z * z);
  // double r = rfunc.norm2();

  std::vector<real_function_3d> funcs;
  funcs.push_back(x);
  funcs.push_back(y);
  funcs.push_back(z);
  // Done
  return funcs;
}

typedef Tensor<double> tensorT;
typedef Function<double, 3> functionT;
typedef std::shared_ptr<FunctionFunctorInterface<double, 3>> functorT;
typedef FunctionFactory<double, 3> factoryT;

response_space TDDFT::PropertyRHS(World& world, Property& p) const {
  if (Rparams.print_level >= 1) {
    molresponse::start_timer(world);
  }

  print("Creating RHS for ", p.property, "operator");
  response_space rhs(world, p.num_operators, Gparams.num_orbitals);

  reconstruct(world, Gparams.orbitals);
  QProjector<double, 3> Qhat(world, Gparams.orbitals);
  // Set the dipoles (ground orbitals are probably
  // more accurate now, so recalc the dipoles)
  // why is it called dipole guess.
  // This is just orbitals times dipole operator
  std::vector<real_function_3d> orbitals = Gparams.orbitals;

  print("num operators ", p.num_operators);
  for (size_t i = 0; i < p.num_operators; i++) {
    // question here....MolecularDerivativeFunctor takes derivative with
    // respect to axis atom and axis
    // here we save
    // need to project

    rhs[i] =
        mul(world, p.operator_vector.at(i), Gparams.orbitals, Rparams.small);

    truncate(world, rhs[i]);
    // rhs[i].truncate_vec();

    // project rhs vectors for state
    rhs[i] = Qhat(rhs[i]);
    // truncate(world, rhs[i], true);
    for (size_t j = 0; j < orbitals.size(); j++) {
      print("RHS norm for after orbital ",
            j,
            "Response state  ",
            i,
            ": ",
            rhs[i][j].norm2());
    }

    world.gop.fence();
    // core projector contribution
  }

  // if (world.rank() ==dipole 0) print("derivatives:\n", r, ru, rc, ra);
  molresponse::end_timer(world, "rhs vectors");
  return rhs;
}
// Returns the derivative of the coulomb operator, applied to ground state
// orbitals Returns the Electron Interaction Gamma Response Function for all
// states. This function assumes all orbitals and response functions are real f
// and g are the response functions phi are the ground state orbitals small and
// thresh are accuracy parameters for the creating of the coulomb operator
response_space TDDFT::CreateCoulombDerivativeRF(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();     // number of resposne states or frequencies
  size_t n = f[0].size();  // number of ground states  x[m][n]
  // Zero function, to be returned
  response_space deriv_J(world, m, n);  // J_p--Jderivative
  // Need the coulomb operator
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  // Temperary storage
  real_function_3d f_density = real_function_3d(world);
  // Need to run over each state
  for (size_t k = 0; k < m; k++) {
    // transition_density = dot(world, f[k] + g[k], orbitals); //sum the vector
    // of functions
    // This works because we assume x,y,phi_i all to be real
    // Apply coulomb operator
    f_density = apply(op, dot(world, f[k], phi));
    // transition_density = apply(op, rho);
    for (size_t p = 0; p < n; p++) {
      // Multiply by ground state orbital p
      // and save the result
      deriv_J[k][p] = f_density * phi[p];
    }
  }
  return deriv_J;
}
// Returns the derivative of the conjugate couloumb derivative operator, applied
// to to the groundstate orbitals.  (TODO: set up imaginary functions)
response_space TDDFT::CreateCoulombDerivativeRFDagger(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();     // number of resposne states or frequencies
  size_t n = f[0].size();  // number of ground states  x[m][n]
  // Zero function, to be returned
  response_space deriv_J_dagger(world, m, n);  // J_p--Jderivative
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  real_function_3d f_density = real_function_3d(world);
  for (size_t k = 0; k < m; k++) {  // for each of the m response states
    // dot vector of response functions with orbitals phi
    f_density = apply(op, dot(world, phi, f[k]));
    // f_density = apply(op,dot(world,dagger(phi),f[k])));
    // TODO write or find a dagger function
    //
    // apply to each orbital to make up jdaggerKP
    for (size_t p = 0; p < n; p++) {
      deriv_J_dagger[k][p] = f_density * phi[p];
    }
  }
  return deriv_J_dagger;
}

// Does what it sounds like it does
response_space TDDFT::CreateExchangeDerivativeRF(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Zero function, to be returned
  response_space deriv_k(world, m, n);

  // Need the coulomb operator
  real_convolution_3d op = CoulombOperator(world, small, thresh);

  // Potential is not stored by default
  // Need to run over occupied orbitals
  // Need to run over all virtual orbitals originating from orbital p
  // Need to sum over occupied orbitals
  if (Rparams.store_potential) {
    for (size_t p = 0; p < n; p++) {
      for (size_t k = 0; k < m; k++) {
        for (size_t i = 0; i < n; i++) {
          deriv_k[k][p] += stored_potential[i][p] * f[k][i];
        }
      }
    }
  } else {                            // But the storage can be turned off...{
    for (size_t p = 0; p < n; p++) {  //
      for (size_t k = 0; k < m; k++) {
        for (size_t i = 0; i < n; i++) {
          // and add to total
          real_function_3d rho = phi[i] * phi[p];
          // Apply coulomb operator
          rho = apply(op, rho);
          // Multiply by response function (k,i)
          // and add to total
          deriv_k[k][p] += rho * f[k][i];
        }
      }
    }
  }
  return deriv_k;
}

// Does what it sounds like it does
response_space TDDFT::CreateExchangeDerivativeRFDagger(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();
  // Zero function, to be returned
  response_space deriv_k_dagger(world, m, n);
  // Need the coulomb operator
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  // Need to run over occupied orbitals
  for (size_t p = 0; p < n; p++) {
    // Need to run over all virtual orbitals originating from orbital p
    for (size_t k = 0; k < m; k++) {
      // Need to sum over occupied orbitals
      for (size_t i = 0; i < n; i++) {
        // Get density (ground state orbitals)
        real_function_3d rho = f[k][i] * phi[p];
        // real_function_3d rho = dagger(f[k][i]) * phi[p];TODO:DAGGER()
        // Apply coulomb operator
        rho = apply(op, rho);
        // and add to total
        deriv_k_dagger[k][p] += rho * phi[i];
      }
    }
  }
  return deriv_k_dagger;
}

response_space TDDFT::CreateXCDerivativeRF(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Initialize response function
  response_space deriv_XC(world, m, n);
  // Get WF for each resonse function
  std::vector<real_function_3d> WxconF = GetWxcOnFDensities(world, phi, f);
  // apply xc kernel to ground staate orbitals
  for (size_t i = 0; i < m; i++) {
    deriv_XC[i] = mul_sparse(world, WxconF[i], phi, thresh, false);
  }
  world.gop.fence();
  return deriv_XC;
}

response_space TDDFT::CreateXCDerivativeRFDagger(
    World& world,
    const response_space& f,
    const std::vector<real_function_3d>& phi,
    double small,
    double thresh) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Initialize response function
  response_space deriv_XC(world, m, n);
  // Get WF for each resonse function
  std::vector<real_function_3d> conjWxconF =
      GetConjugateWxcOnFDensities(world, phi, f);
  // apply xc kernel to ground staate orbitals
  for (size_t i = 0; i < m; i++) {
    deriv_XC[i] = mul_sparse(world, conjWxconF[i], phi, thresh, false);
  }
  world.gop.fence();
  return deriv_XC;
}

// Creates diagonal (letter A) portions of response matrix
response_space TDDFT::createAf(World& world,
                               response_space& Vf,
                               response_space& F0_f,
                               response_space& Epsilonf,
                               response_space& Hf,
                               response_space& f,
                               std::vector<real_function_3d>& orbitals,
                               size_t print_level,
                               std::string xy) {
  size_t m = f.size();
  size_t n = f[0].size();

  response_space Af(world, m, n);
  // Create the ground-state fock operator on response components
  // Create F0 on x or y response  (Fx-xF)
  // Af = F0*xp -Fpp*xp +Hf,p -\sum_{neq p} xi*Fip

  // ResponseFunction F0_f = CreateFock(world, Vf, f, print_level, xy); // Fx
  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0)
      printf("   Ground Fock matrix for %s components:\n", xy.c_str());
    Tensor<double> temp2 = expectation(world, f, F0_f);
    if (world.rank() == 0) print(temp2);
    printf("   Testing Expectation Ground Fock matrix for %s components:\n",
           xy.c_str());
    //if (print_level >= 1) molresponse::start_timer(world);
   // Tensor<double> temp3 = expectation2(world, f, F0_f);
    //if (print_level >= 1) molresponse::end_timer(world, " Expectation 2");
    //if (world.rank() == 0) print(temp3);
  }
  //
  Af = Hf + F0_f;
  Af = Af - Epsilonf;
  // Need to project
  // It actually should not be necessary here if we project G and H
  // lets not project and have G and H projected
  /*
  QProjector<double, 3> projector(world, orbitals);
  for (size_t i = 0; i < m; i++)
    Af[i] = projector(Af[i]);

    */
  return Af;

  // And return the sum
}

// Creates the off diagonal (letter B) portions of response matrix
// Simply projects out ground state from Gf response functions
response_space TDDFT::createBf(World& world,
                               response_space& Gf,
                               std::vector<real_function_3d>& orbitals,
                               size_t print_level) {
  // Start a timer
  // if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  // size_t m = Gf.size();
  // size_t n = Gf[0].size();
  // ResponseFunction Bf(world, m, n);
  // Project out the ground state
  // QProjector<double, 3> projector(world, orbitals);
  // for (size_t i = 0; i < m; i++) Bf[i] = projector(Gf[i]);

  // End timer
  // if (print_level >= 1) molresponse::end_timer(world, "   Creating Bf:");

  // Done
  return Gf;
}
// Calculates ground state coulomb potential
real_function_3d TDDFT::Coulomb(World& world) {
  // Coulomb operator
  real_convolution_3d op =
      CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());

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
response_space TDDFT::exchange(World& world, response_space& f) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Coulomb operator
  real_convolution_3d op =
      CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());

  // Container for results and others
  response_space result(world, m, n);
  real_function_3d psif = real_function_3d(world);

  // Modified 'small memory' algorithm from SCF.cc
  f.reconstruct_rf();

  // Run over each excited state
  for (size_t k = 0; k < m; k++) {
    // And run over each occupied state
    for (size_t j = 0; j < n; j++) {
      // Get a vector of transition densities
      auto phix = mul_sparse(
          world, Gparams.orbitals[j], f[k], FunctionDefaults<3>::get_thresh());

      // Clean up
      truncate(world, phix);

      // Apply operator to each member of vector
      phix = apply(world, op, phix);

      // Clean up
      truncate(world, phix);

      // Final multiplication of each member of vector by a single function
      phix = mul_sparse(
          world, Gparams.orbitals[j], phix, FunctionDefaults<3>::get_thresh());

      // Add the vector to result
      gaxpy(world, 1.0, result[k], 1.0, phix);
    }
  }

  // Truncate
  madness::truncate(world, result);

  // Done!
  return result;
}

// Returns the ground state potential applied to functions f
response_space TDDFT::CreatePotential(World& world,
                                      response_space& f,
                                      XCOperator xc,
                                      size_t print_level,
                                      std::string xy) {
  // Start a timer
  if (print_level >= 3) molresponse::start_timer(world);

  // Return container
  response_space V_x_resp(world, f.size(), f[0].size());

  // Computing \hat{V}^0 = v_nuc + v_coul + v_exch
  // v_nuc first
  real_function_3d v_nuc, v_coul;
  if (not Rparams.store_potential) {
    // "a" is the core type
    PotentialManager manager(Gparams.molecule, "a");
    manager.make_nuclear_potential(world);
    // v_nuc = manager.vnuclear().truncate();
    v_nuc = manager.vnuclear();
    v_nuc.truncate();

    // v_coul next
    // This does not include final multiplication of each orbital
    // 2.0 scale is from spin integration
    v_coul = Coulomb(world);
    v_coul.scale(2.0);
  } else {  // Already pre-computed
    v_nuc = stored_v_nuc;
    v_coul = stored_v_coul;
  }

  // Intermediaries

  response_space v_exch(world, f.size(), f[0].size());
  real_function_3d v_xc = real_factory_3d(world).fence(true);

  // If including any exact HF exchange
  if (xcf.hf_exchange_coefficient()) {
    // V_exch last
    // Multiplication by f functions is included in construction
    v_exch = exchange(world, f);
  }
  if (xcf.hf_exchange_coefficient() != 1.0) {
    // Calculate DFT potential
    v_xc = xc.make_xc_potential();
  }

  // Assemble all the pieces for V_x
  V_x_resp = (f * (v_coul + v_nuc + v_xc));
  V_x_resp = V_x_resp - (v_exch * xcf.hf_exchange_coefficient());

  // Debugging output
  if (print_level >= 2) {
    // Print potential energy matrices
    if (world.rank() == 0)
      printf("   Nuclear potential matrix for %s components:\n", xy.c_str());
    response_space temp1 = f * v_nuc;
    Tensor<double> temp = expectation(world, f, temp1);
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0)
      printf("   Coulomb potential matrix for %s components:\n", xy.c_str());
    response_space temp2 = f * v_coul;
    temp = expectation(world, f, temp2);
    if (world.rank() == 0) print(temp);
    if (xcf.hf_exchange_coefficient()) {
      if (world.rank() == 0)
        printf("   Exchange potential matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, v_exch);
    } else if (xcf.hf_exchange_coefficient() != 1.0) {
      if (world.rank() == 0)
        printf("   XC potential matrix for %s components:\n", xy.c_str());
      v_exch = f * v_xc;
      temp = expectation(world, f, v_exch);
    }
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0)
      printf("   Total Potential Energy matrix for %s components:\n",
             xy.c_str());
    temp = expectation(world, f, V_x_resp);
    if (world.rank() == 0) print(temp);
  }
  V_x_resp.truncate_rf();

  // Basic output
  if (print_level >= 3) molresponse::end_timer(world, "Creating V0 * x:");

  // Done
  return V_x_resp;
}

void TDDFT::computeElectronResponse(World& world,
                                    ElectronResponseFunctions& I,
                                    response_space& x,
                                    response_space& y,
                                    std::vector<real_function_3d>& orbitals,
                                    XCOperator xc,
                                    Tensor<double>& hamiltonian,
                                    Tensor<double>& ham_no_diag,
                                    double small,
                                    double thresh,
                                    size_t print_level,
                                    std::string xy) {
  // Start a timer
  if (print_level >= 1) molresponse::start_timer(world);

  I.Vx = CreatePotential(world, x, xc, print_level, "x");
  I.F0_x = CreateFock(world, I.Vx, x, print_level, "x");
  // epsilon with diag for FullR matrix
  I.EpsilonX = x * hamiltonian;  // scale_2d(world, x, hamiltonian);
  I.EpsilonXNoDiag =
      x * ham_no_diag;  // scale_2d(world, x, ham_no_diag);  // for rhs
  // compute Electron Interaction Terms for this Iteration
  I.Hx = ComputeHf(world, x, orbitals, small, thresh, print_level, "x");
  // print(Hx);
  // else Compute everything
  if (not Rparams.tda) {  // not sure why this is the condition
    I.Gy = ComputeGf(world, y, orbitals, small, thresh, print_level, "y");

    I.Vy = CreatePotential(world, y, xc, print_level, "y");
    I.F0_y = CreateFock(world, I.Vy, y, print_level, "y");
    I.EpsilonY = y * hamiltonian;        // scale_2d(world, y, hamiltonian);
    I.EpsilonYNoDiag = y * ham_no_diag;  // scale_2d(world, y, ham_no_diag);
    I.Hy = ComputeHf(world, y, orbitals, small, thresh, print_level, "y");
    I.Gx = ComputeGf(world, x, orbitals, small, thresh, print_level, "x");
  }
  molresponse::end_timer(world, "Creating Electron Responses");

  // Create \hat{V}^0 applied to response functions
}
// Returns a Tensor of inner products, where
// result(i,j) = inner(a[i],b[j]).sum()
Tensor<double> TDDFT::expectation(World& world,
                                  const response_space& A,
                                  const response_space& B) {
  // Get sizes
  size_t dim_1 = A.size();
  size_t dim_2 = A[0].size();
  // Need to take transpose of each input ResponseFunction
  response_space A_t(world, dim_2, dim_1);
  response_space B_t(world, dim_2, dim_1);
  for (size_t i = 0; i < dim_1; i++) {
    for (size_t j = 0; j < dim_2; j++) {
      A_t[j][i] = A[i][j];
      B_t[j][i] = B[i][j];
    }
  }
  // Container for result
  Tensor<double> result(dim_1, dim_1);
  /**
   * @brief
   * [x1 x2 x3]T[x1 x2 x3]
   *
   */
  // Run over dimension two
  // each vector in orbital has dim_1 response functoins associated
  for (size_t p = 0; p < dim_2; p++) {
    result += matrix_inner(world, A_t[p], B_t[p]);
  }

  // Done
  return result;
}

Tensor<double> TDDFT::expectation2(World& world,
                                   const response_space& a,
                                   const response_space& b) {
  MADNESS_ASSERT(a.size() > 0);
  MADNESS_ASSERT(a.size() == b.size());
  MADNESS_ASSERT(a[0].size() > 0);
  MADNESS_ASSERT(a[0].size() == b[0].size());
  // Get sizes
  size_t dim_1 = a.size();

  Tensor<double> result(a.size(), a.size());
  // hold the results
  response_space c(world, dim_1, dim_1);

  for (size_t i(0); i < dim_1; i++) {
    for (size_t j(0); j < dim_1; j++) {
      c[i][j] = dot(world, a[i], b[j]);
      result(i, j) = c[i][j].trace();
    }
  }
  /**
   * @brief
   * [x1 x2 x3]T[x1 x2 x3]
   *
   */
  // Run over dimension two
  return result;
}
void TDDFT::PrintRFExpectation(World& world,
                               response_space f,
                               response_space g,
                               std::string fname,
                               std::string gname) {
  Tensor<double> t = expectation(world, f, g);
  if (world.rank() == 0) {
    print(" Expectation between ", fname, " and ", gname);
    print(t);
  }
}
void TDDFT::PrintResponseVectorNorms(World& world,
                                     response_space f,
                                     std::string fname) {
  if (world.rank() == 0) {
    print(" Norms of ResVector: ", fname);
    print(f.norm2());
  }
}
// Returns the ground state fock operator applied to functions f
response_space TDDFT::CreateFock(World& world,
                                 response_space& Vf,
                                 response_space& f,
                                 size_t print_level,
                                 std::string xy) {
  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0)
      printf("   Creating perturbed fock matrix for %s components\n",
             xy.c_str());
  }
  // Container to return
  response_space fock;  // Fock = (T + V) * orbitals
                        // Already have V (input parameter)
  // Create T
  // Make the derivative operators in each direction
  real_derivative_3d Dx(world, 0);
  real_derivative_3d Dy(world, 1);
  real_derivative_3d Dz(world, 2);
  // Apply derivatives to orbitals
  f.reconstruct_rf();
  response_space dvx = apply(world, Dx, f);
  response_space dvy = apply(world, Dy, f);
  response_space dvz = apply(world, Dz, f);
  // Apply again for 2nd derivatives
  response_space dvx2 = apply(world, Dx, dvx);
  response_space dvy2 = apply(world, Dy, dvy);
  response_space dvz2 = apply(world, Dz, dvz);
  // Add together derivatives
  fock = (dvx2 + dvy2 + dvz2) * (-0.5);
  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0)
      printf("   Kinetic energy matrix for %s components:\n", xy.c_str());
    Tensor<double> temp = expectation(world, f, fock);
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0)
      printf("   Potential energy matrix for %s components:\n", xy.c_str());
    temp = expectation(world, f, Vf);
    if (world.rank() == 0) print(temp);
  }
  // Add in potential
  fock = fock + Vf;
  fock.truncate_rf();
  // Done
  return fock;
}

Zfunctions TDDFT::ComputeZFunctions(World& world,
                                    std::vector<real_function_3d> rho_omega,
                                    response_space orbital_products,
                                    response_space& x,
                                    response_space& y,
                                    XCOperator xc,
                                    double x_shifts,
                                    double y_shifts,
                                    const GroundParameters& Gparams,
                                    const ResponseParameters& Rparams,
                                    Tensor<double> ham_no_diagonal) {
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();

  Zfunctions Z;
  // x functions
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  Z.v0_x = CreatePotential(world, x, xc, Rparams.print_level, "x");

  // + \Delta xp
  Z.v0_x += x * x_shifts;  // scale(Z.v0_x, x_shifts);

  Z.x_f_no_diag = x * ham_no_diag;  // scale_2d(world, x, ham_no_diagonal);

  GammaResponseFunctions gamma = ComputeGammaFunctions(
      world,  x, y, xc, Gparams, Rparams,Rparams.omega != 0.0);
  // y functions
  if (Rparams.omega != 0.0) {
    Z.v0_y = CreatePotential(world, y, xc, Rparams.print_level, "y");
    // no need to apply shift in y case
    // Z.v0_y = apply_shift(world, y_shifts, Z.v0_y, y);

    Z.y_f_no_diag = y * ham_no_diag;
    // scale_2d(world, y, ham_no_diagonal);
  }

  Z.Z_x = Z.v0_x - Z.x_f_no_diag + gamma.gamma;
  if (Rparams.omega != 0.0) {
    Z.Z_y = Z.v0_y - Z.y_f_no_diag + gamma.gamma_conjugate;
  }

  if (Rparams.print_level >= 2) {
    Tensor<double> t = expectation(world, x, Z.v0_x);
    if (world.rank() == 0) {
      print("   Energy scaled response orbitals for x components:");
      print(t);
    }

    if (Rparams.omega != 0.0) {
      t = expectation(world, y, Z.v0_y);
      if (world.rank() == 0) {
        print("   Energy scaled response orbitals for y components:");
        print(t);
      }
    }
  }

  world.gop.fence();
  return Z;
}
X_space TDDFT::ComputeResponseResidual(World& world,
                                       std::vector<real_function_3d> rho_omega,
                                       response_space orbital_products,
                                       response_space& x,
                                       response_space& y,
                                       response_space rhs_x,
                                       response_space rhs_y,
                                       XCOperator xc,
                                       const GroundParameters& Gparams,
                                       const ResponseParameters& Rparams,
                                       Tensor<double> hamiltonian,
                                       double omega,
                                       size_t iteration) {
  // compute
  size_t m = Rparams.states;
  size_t n = Gparams.num_orbitals;
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();

  real_convolution_3d op = CoulombOperator(world, small, thresh);
  response_space v0x(world, m, n);
  response_space v0y(world, m, n);
  response_space F0x(world, m, n);
  response_space F0y(world, m, n);
  // x scaled by hamiltonian
  response_space xham(world, m, n);
  response_space yham(world, m, n);
  // 2 electron terms
  response_space Hx(world, m, n);
  response_space Gy(world, m, n);
  response_space Hy(world, m, n);
  response_space Gx(world, m, n);

  response_space omegaX(world, m, n);
  response_space omegaY(world, m, n);

  ResidualResponseVectors residual(world, m, n);
  /* We first compute the 1 electron potentials
   * We compute the the pieces that depend on x response functions
   */
  // x functions
  // V0 applied to x response function
  v0x = CreatePotential(world, x, xc, Rparams.print_level, "x");
  F0x = CreateFock(world, v0x, x, Rparams.print_level, "x");
  F0x.truncate_rf();

  // x response scaled by off diagonal ham
  xham = x * hamiltonian;  // scale_2d(world, x, hamiltonian);
  omegaX = x_response * omega;
  if (Rparams.print_level == 3) {
    print("norms of x scaled by ham no diag");
    print(xham.norm2());
  }
  // If not static we compute the y components
  if (Rparams.omega != 0.0) {
    v0y = CreatePotential(world, y, xc, Rparams.print_level, "y");
    F0y = CreateFock(world, v0y, y, Rparams.print_level, "y");
    F0y.truncate_rf();
    yham = y * hamiltonian;  // scale_2d(world, y, hamiltonian);
    omegaY = y_response * omega;
  }
  // Some printing for debugging
  if (Rparams.print_level >= 2) {
    { PrintRFExpectation(world, x, v0x, "x", "V0X"); }

    if (Rparams.omega != 0.0) {
      PrintRFExpectation(world, y, v0y, "y", "V0Y");
    }
  }
  // Last we compute the 2-electron peices
  //
  //
  if (Rparams.old_two_electron) {
    Hx = ComputeHf(
        world, x, Gparams.orbitals, small, thresh, Rparams.print_level, "x");

    Gy = ComputeGf(
        world, y, Gparams.orbitals, small, thresh, Rparams.print_level, "y");
    if (Rparams.omega != 0.0) {
      Hy = ComputeHf(
          world, y, Gparams.orbitals, small, thresh, Rparams.print_level, "x");

      Gx = ComputeGf(
          world, x, Gparams.orbitals, small, thresh, Rparams.print_level, "y");
    }

    // Z.Z_x = (Z.v0_x - Z.x_f_no_diag + rhs_x) * -2;
    residual.x = (F0x - xham + Hx + Gy + rhs_x - omegaX);
    if (Rparams.omega != 0.0) {
      residual.y = (F0y - yham + Hy + Gx + rhs_y + omegaY);
      // Z.Z_y = (Z.v0_y - Z.y_f_no_diag + rhs_y) * -2;
    }
  } else {
    GammaResponseFunctions gamma = ComputeGammaFunctions(
        world,  x, y, xc, Gparams, Rparams,Rparams.omega != 0.0);
    // We can use the old algorithm here for testings
    // we then assemble the right hand side vectors
    residual.x = (F0x - omegaX - xham + gamma.gamma + rhs_x);
    // Z.Z_x = Z.v0_x - Z.x_f_no_diag + rhs_x;
    if (Rparams.omega != 0.0) {
      residual.y = (F0y + omegaY - yham + gamma.gamma_conjugate);
      // Z.Z_y = Z.v0_y - Z.y_f_no_diag + rhs_y;
    }
  }
  residual.x.truncate_rf();
  residual.y.truncate_rf();
  return X_space(residual.x, residual.y);
}

void TDDFT::xy_from_XVector(response_space& x,
                            response_space& y,
                            std::vector<X_vector>& Xvectors) {
  MADNESS_ASSERT(x.size() == Xvectors.size());
  MADNESS_ASSERT(y.size() == Xvectors.size());
  MADNESS_ASSERT(x[0].size() == size_orbitals(Xvectors[0]));
  MADNESS_ASSERT(y[0].size() == size_orbitals(Xvectors[0]));

  for (size_t b = 0; b < x.size(); b++) {
    std::cout << "moving state " << b << std::endl;
    for (auto xs : x[b]) {
      std::cout << "norm xs before move" << xs.norm2() << std::endl;
    }
    for (auto xs : Xvectors[b].X[0]) {
      std::cout << "norm Xvector before move" << xs.norm2() << std::endl;
    }
    x[b] = Xvectors[b].X[0];
    y[b] = Xvectors[b].Y[0];
    for (auto xs : x[b]) {
      std::cout << "norm xs after move" << xs.norm2() << std::endl;
    }
  }
}

void TDDFT::vector_stats(const std::vector<double>& v,
                         double& rms,
                         double& maxabsval) const {
  rms = 0.0;
  maxabsval = v[0];
  for (size_t i = 0; i < v.size(); ++i) {
    rms += v[i] * v[i];
    maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
  }
  rms = sqrt(rms / v.size());
}
double TDDFT::do_step_restriction(World& world,
                                  const vecfuncT& x,
                                  vecfuncT& x_new,
                                  std::string spin) const {
  std::vector<double> anorm = norm2s(world, sub(world, x, x_new));
  size_t nres = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    print("anorm ", i, " : ", anorm[i]);
    if (anorm[i] > Rparams.maxrotn) {
      double s = Rparams.maxrotn / anorm[i];
      ++nres;
      if (world.rank() == 0) {
        if (nres == 1 and (Rparams.print_level > 1))
          printf("  restricting step for %s orbitals:", spin.c_str());
        printf(" %d", i);
      }
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
    }
  }
  if (nres > 0 && world.rank() == 0 and (Rparams.print_level > 1)) printf("\n");

  world.gop.fence();
  double rms, maxval;
  vector_stats(anorm, rms, maxval);
  if (world.rank() == 0 and (Rparams.print_level > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}
// Construct the Hamiltonian
Tensor<double> TDDFT::CreateResponseMatrix(
    World& world,
    response_space& x,
    ElectronResponseFunctions& I,
    std::vector<real_function_3d>& ground_orbitals,
    size_t print_level,
    std::string xy) {
  // Start a timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Construct intermediary
  // Sets fe to be (\hat{fock} - eps)*f
  response_space Ax = createAf(world,
                               I.Vx,
                               I.F0_x,
                               I.EpsilonX,
                               I.Hx,
                               x,
                               ground_orbitals,
                               print_level,
                               xy);

  // Make the matrix
  Tensor<double> xAx = response_space_inner(x, Ax);

  // Enforce symmetry....not sure if this is a good idea
  //
  for (int64_t i = 0; i < xAx.dim(0); i++) {
    for (int64_t j = i + 1; j < xAx.dim(1); j++) {
      //      print(i, j);
      //      print(xAx(i, j));
      //     print(xAx(j, i));
      xAx(j, i) = xAx(i, j);// TODO Complex conjugate
    }
  }

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    printf("\n   Response matrix for %s components:\n", xy.c_str());
    print(xAx);
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "Create resp. matrix:");

  // Done
  return xAx;
}

// Constructs the matrix, so really it does
// [ X Y ] [ A  B ] [ X ]
//         [ B  A ] [ Y ]
Tensor<double> TDDFT::CreateFullResponseMatrix(
    World& world,
    response_space& x,  // x response functions
    response_space& y,  // y response functions
    ElectronResponseFunctions& I,
    std::vector<real_function_3d>& ground_orbitals,  // ground state orbitals
    double small,
    double thresh,
    size_t print_level) {
  // Start timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Create the A pieces
  // (Sets fe_x and fe_y to be (\hat{F}-eps) * resp. funcs.
  response_space A_x = createAf(world,
                                I.Vx,
                                I.F0_x,
                                I.EpsilonX,
                                I.Hx,
                                x,
                                ground_orbitals,
                                print_level,
                                "x");
  response_space A_y = createAf(world,
                                I.Vy,
                                I.F0_y,
                                I.EpsilonY,
                                I.Hy,
                                y,
                                ground_orbitals,
                                print_level,
                                "y");

  // Create the B pieces
  response_space B_x = createBf(world, I.Gx, ground_orbitals, print_level);
  response_space B_y = createBf(world, I.Gy, ground_orbitals, print_level);

  // Debugging output
  if (print_level >= 2) {
    Tensor<double> t = expectation(world, x, A_x);
    if (world.rank() == 0) {
      print("\n   A_x:");
      print(t);
    }

    t = expectation(world, x, B_y);
    if (world.rank() == 0) {
      print("\n   B_y:");
      print(t);
    }

    t = expectation(world, y, A_y);
    if (world.rank() == 0) {
      print("\n   A_y:");
      print(t);
    }

    t = expectation(world, y, B_x);
    if (world.rank() == 0) {
      print("\n   B_x:");
      print(t);
    }
  }

  // Add corresponding pieces together
  A_x = A_x + B_y;
  A_y = B_x + A_y;  // Needs adjoint if complex

  // Construct matrix to be returned
  Tensor<double> response_matrix =
      expectation(world, x, A_x) + expectation(world, y, A_y);

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("\n   Full Coupled Response Matrix:");
    print(response_matrix);
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "Create resp. matrix:");

  // Done
  return response_matrix;
}

// Returns the shift needed to make sure that
// -2.0 * (ground_state_energy + excited_state_energy)
// is negative. Please note: The same shift needs to
// be applied to the potential.
Tensor<double> TDDFT::create_shift(World& world,
                                   Tensor<double>& ground,
                                   Tensor<double>& omega,
                                   size_t print_level,
                                   std::string xy) {
  // Start a timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = omega.size();
  size_t n = ground.size();

  // Container to hold shift
  Tensor<double> result(m, n);

  // Run over excited components
  for (size_t k = 0; k < m; k++) {
    // Run over ground components
    for (size_t p = 0; p < n; p++) {
      if (ground(p) + omega(k) > 0) {
        // Calculate the shift needed to get energy to -0.05,
        // which was arbitrary (same as moldft)
        result(k, p) = -(ground(p) + omega(k) + 0.05);

        // Basic output
        if (print_level >= 2) {
          if (world.rank() == 0)
            printf(
                "   Shift needed for transition from ground orbital %d to "
                "response %s state %d\n",
                static_cast<int>(p),
                xy.c_str(),
                static_cast<int>(k));
          if (world.rank() == 0) print("   Ground energy =", ground(p));
          if (world.rank() == 0) print("   Excited energy =", omega(k));
          if (world.rank() == 0) print("   Shifting by", result(k, p));
          if (world.rank() == 0) print("");
        }
      }
    }
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "Create shift:");

  // Done
  return result;
}

// Returns the shift needed to make sure that
// (ground_state_energy + excited_state_energy + shift) = target
// Please note: The same shift needs to be applied to the potential.
Tensor<double> TDDFT::create_shift_target(World& world,
                                          Tensor<double>& ground,
                                          Tensor<double>& omega,
                                          double target,
                                          size_t print_level,
                                          std::string xy) {
  // Start a timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = omega.size();
  size_t n = ground.size();

  // Container to hold shift
  Tensor<double> result(m, n);

  // Run over excited components
  for (size_t k = 0; k < m; k++) {
    // Run over ground components
    for (size_t p = 0; p < n; p++) {
      // Calculate the shift needed to get energy to target
      result(k, p) = -(ground(p) + omega(k) - target);

      // Basic output
      if (print_level >= 2) {
        if (world.rank() == 0)
          printf(
              "   Shift needed for transition from ground orbital %d to "
              "response %s state %d\n",
              static_cast<int>(p),
              xy.c_str(),
              static_cast<int>(k));
        if (world.rank() == 0) print("   Ground energy =", ground(p));
        if (world.rank() == 0) print("   Excited energy =", omega(k));
        if (world.rank() == 0) print("   Shifting by", result(k, p));
        if (world.rank() == 0) print("");
      }
    }
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "Create shift:");

  // Done
  return result;
}

// Returns the given shift applied to the given potential
response_space TDDFT::apply_shift(World& world,
                                  Tensor<double>& shifts,
                                  response_space& V,
                                  response_space& f) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Sizes inferred from V
  size_t n = V[0].size();
  size_t m = V.size();

  // Container to return
  response_space shifted_V(world, m, n);

  // Run over occupied
  for (size_t k = 0; k < m; k++) {
    // Run over virtual
    for (size_t p = 0; p < n; p++) {
      shifted_V[k][p] = V[k][p] + shifts(k, p) * f[k][p];
    }
  }

  shifted_V.truncate_rf();

  // End timer
  if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply shift:");

  // Done
  return shifted_V;
}

// Returns the given shift applied to the given potential
response_space TDDFT::apply_shift(World& world,
                                  double& shift,
                                  response_space& V,
                                  response_space& f) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Sizes inferred from V
  size_t n = V[0].size();
  size_t m = V.size();

  // Container to return
  response_space shifted_V(world, m, n);

  // Run over occupied
  for (size_t k = 0; k < m; k++) {
    // Run over virtual
    for (size_t p = 0; p < n; p++) {
      shifted_V[k][p] = V[k][p] + shift * f[k][p];
    }
  }

  shifted_V.truncate_rf();

  // End timer
  if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply shift:");

  // Done
  return shifted_V;
}

// Function to make a vector of BSH operators using ground and excited
// state energies
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
TDDFT::create_bsh_operators(World& world,
                            Tensor<double>& shift,
                            Tensor<double>& ground,
                            Tensor<double>& omega,
                            double small,
                            double thresh) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Sizes inferred from ground and omega
  size_t n = ground.size();
  size_t m = omega.size();

  // Make the vector
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> operators;

  // Make a BSH operator for each response function
  // Run over excited components
  for (size_t k = 0; k < m; k++) {
    // Container for intermediary
    std::vector<std::shared_ptr<real_convolution_3d>> temp(n);

    // Run over occupied components
    for (size_t p = 0; p < n; p++) {
      temp[p] = std::shared_ptr<SeparatedConvolution<double, 3>>(
          BSHOperatorPtr3D(world,
                           sqrt(-2.0 * (ground(p) + omega(k) + shift(k, p))),
                           small,
                           thresh));
    }

    // Add intermediary to return container
    operators.push_back(temp);
  }

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Creating BSH ops:");

  // Done
  return operators;
}

// Function to make a vector of BSH operators using ground and excited
// state energies
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
TDDFT::CreateBSHOperatorPropertyVector(World& world,
                                       Tensor<double>& shift,
                                       Tensor<double>& ground,
                                       Tensor<double>& omega,
                                       double small,
                                       double thresh) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Sizes inferred from ground and omega
  size_t n = ground.size();  // number of orbitals
  size_t num_states = Rparams.states;
  size_t num_freq = omega.size();  // number of frequency states
  // print("num of freq", num_freq);

  // Make the vector
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> operators;

  // Make a BSH operator for each response function
  // Run over excited components
  // print("num of states bsh step", num_states);
  for (size_t k = 0; k < num_freq; k++) {
    // Container for intermediary
    std::vector<std::shared_ptr<real_convolution_3d>> temp(n);
    for (size_t state = 0; state < num_states; state++) {
      // Run over occupied components
      for (size_t p = 0; p < n; p++) {
        temp[p] = std::shared_ptr<SeparatedConvolution<double, 3>>(
            BSHOperatorPtr3D(world,
                             sqrt(-2.0 * (ground(p) + omega(k) + shift(k, p))),
                             small,
                             thresh));
      }
      operators.push_back(temp);
    }

    // Add intermediary to return container
  }

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Creating BSH ops:");

  // Done
  return operators;
}
// creating a shift in a property calculation requires only one double for the
// shift
std::vector<std::shared_ptr<real_convolution_3d>>
TDDFT::CreateBSHOperatorPropertyVector(World& world,
                                       double& shift,
                                       Tensor<double>& ground,
                                       double& omega,
                                       double lo,
                                       double eps) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Sizes inferred from ground and omega
  size_t num_ground_states = ground.size();  // number of orbitals
  // print("num of freq", num_freq);

  // Make the vector
  std::vector<std::shared_ptr<real_convolution_3d>> ghat_operators(
      num_ground_states);

  // Make a BSH operator for each response function
  // Run over excited components
  // print("num of states bsh step", num_states);
  // Container for intermediary
  // Run over occupied components
  for (size_t p = 0; p < num_ground_states; p++) {
    double mu = sqrt(-2.0 * (ground(p) + omega + shift));
    ghat_operators[p] = std::shared_ptr<SeparatedConvolution<double, 3>>(
        BSHOperatorPtr3D(world, mu, lo, eps));
  }
  // Add intermediary to return container

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Creating BSH ops:");

  // Done
  return ghat_operators;
}
// Returns the second order update to the energies of the excited components
// Not currently used.
Tensor<double> TDDFT::calculate_energy_update(World& world,
                                              response_space& rhs,
                                              response_space& f_residuals,
                                              response_space& new_f,
                                              size_t print_level,
                                              std::string xy) {
  /*
   *  The correction is:
   *      \delta \omega^{(k)} = - \frac{ \sum_p\left< \hat{V}^0 x_p^{(k)}(r) +
   * (1 - \hat{\rho}^0) \Gamma_p^{(k)}(r)\right| \left. x_p^{(k)} -
   * \~{x}_p^{(k)} \right> } { \sum_p \left| \left| \~{x}_p^{(k)} \right|
   * \right|^2 }
   */

  // Basic output
  if (print_level >= 1) {
    if (world.rank() == 0)
      printf("   Calculating energy residuals for %s components\n", xy.c_str());
  }

  // Size inferred
  size_t m = rhs.size();

  // Container for updates
  Tensor<double> updates(m);

  // Need to run over all functions in rhs and calculate inner products.
  // rhs contains the bra in the braket notation above, and f_residuals
  // is the ket.

  // Run over excited components
  for (size_t k = 0; k < m; k++) {
    // vmra.h function, line 627
    // Sum is included inside function call
    updates(k) = inner(f_residuals[k], rhs[k]);

    // Normalize update function
    // The -1.0 is the leading coefficient in the update formula
    // the 1/2 is to undo the scaling of V
    updates(k) = -1.0 / 2.0 * updates(k) / inner(new_f[k], new_f[k]);
  }

  if (print_level >= 1) {
    // Print energy deltas
    if (world.rank() == 0)
      printf("   Energy residuals for %s components:\n", xy.c_str());
    if (world.rank() == 0) print("Er: ", abs(updates));
  }

  // Done?
  return updates;
}

// Specialized for response calculations that returns orthonormalized
// functions
response_space TDDFT::gram_schmidt(World& world, response_space& f) {
  // Sizes inferred
  size_t m = f.size();

  // Return container
  response_space result = f.copy();

  // Orthogonalize
  for (size_t j = 0; j < m; j++) {
    // Need to normalize the row
    double norm = norm2(world, result[j]);

    // Now scale each entry
    scale(world, result[j], 1.0 / norm);

    // Project out from the rest of the vectors
    for (size_t k = j + 1; k < m; k++) {
      // Temp function to hold the sum
      // of inner products
      // vmra.h function, line 627
      double temp = inner(result[j], result[k]);

      // Now subtract
      gaxpy(world, 1.0, result[k], -temp, result[j]);
    }
  }

  result.truncate_rf();

  // Done
  return result;
}

// Specialized for response calculations that returns orthonormalized
// functions
void TDDFT::gram_schmidt(World& world, response_space& f, response_space& g) {
  // Sizes inferred
  size_t m = f.size();

  // Orthogonalize
  for (size_t j = 0; j < m; j++) {
    // Need to normalize the row
    double norm = inner(f[j], f[j]) - inner(g[j], g[j]);

    // Now scale each entry
    scale(world, f[j], 1.0 / sqrt(norm));
    scale(world, g[j], 1.0 / sqrt(norm));

    // Project out from the rest of the vectors
    for (size_t k = j + 1; k < m; k++) {
      // Temp function to hold the sum
      // of inner products
      // vmra.h function, line 627
      double temp = inner(f[j], f[k]) - inner(g[j], g[k]);

      // Now subtract
      gaxpy(world, 1.0, f[k], -temp, f[j]);
      gaxpy(world, 1.0, g[k], -temp, g[j]);
    }
  }

  f.truncate_rf();
  g.truncate_rf();
}

// Returns the max norm of the given vector of functions
double TDDFT::calculate_max_residual(World& world, response_space& f) {
  // Container for max
  double max = 0.0;

  // Run over all functions in f
  for (unsigned int i = 0; i < f.size(); i++) {
    double temp = 0.0;

    for (unsigned int j = 0; j < f[0].size(); j++) {
      temp += pow(f[i][j].norm2(), 2);
    }

    temp = sqrt(temp);

    if (temp > max) max = temp;
  }

  // Done
  return max;
}

// Selects the 'active' orbitals from ground state orbitals to be used in the
// calculation (based on energy distance from the HOMO). Function needs
// knowledge of Gparams.orbitals and Gparams.energies. Function sets
// act_orbitals and num_act_orbitals.
void TDDFT::select_active_subspace(World& world) {
  // Default output
  if (Rparams.print_level >= 0) {
    // Set print output to something reasonable
    std::cout.precision(2);
    std::cout << std::fixed;

    if (world.rank() == 0)
      print(
          "   Selecting ground state subspace to excite from for "
          "components.");
    if (world.rank() == 0)
      print("   This is all orbitals between",
            Rparams.range_low,
            "and",
            Rparams.range_high,
            "\n");

    // Reset precision
    std::cout.precision(10);
    std::cout << std::scientific;
  }

  // Determine active orbitals based on energy differences
  // from HOMO
  for (unsigned int i = 0; i < Gparams.num_orbitals; i++) {
    if (Rparams.range_low < Gparams.energies(i) and
        Gparams.energies(i) < Rparams.range_high) {
      // This orbital should be active, so add to list
      active.push_back(i);
    }
  }

  // Make sure we have at least one ground state orbital to excite from
  MADNESS_ASSERT(active.size() > 0);

  // Now that we know size, allocate act_ground_energies
  act_ground_energies = Tensor<double>(active.size());

  // Now to pull the functions and energies and add to act_orbitals and
  // act_ground_energies
  for (unsigned int i = 0; i < active.size(); i++) {
    act_orbitals.push_back(Gparams.orbitals[active[i]]);
    act_ground_energies(i) =
        Gparams.energies(active[i]);  // Put energies on diagonal
  }

  // Also set the active size
  act_num_orbitals = act_orbitals.size();

  print("Found", act_num_orbitals, "active orbitals.");
}

// Selects from a list of functions and energies the k functions with the
// lowest energy
response_space TDDFT::select_functions(World& world,
                                       response_space& f,
                                       Tensor<double>& energies,
                                       size_t k,
                                       size_t print_level) {
  // Container for result
  response_space answer;

  // Debugging output
  if (print_level >= 1) {
    if (world.rank() == 0)
      print("\n   Selecting the", k, "lowest excitation energy components.\n");
  }

  // Get rid of extra functions and save
  // the first k
  while (f.size() > k) f.pop_back();
  answer = f;
  answer.truncate_rf();

  // Get rid of extra energies and save
  // the first k
  energies = energies(Slice(0, k - 1));

  // Basic output
  if (print_level >= 1) {
    if (world.rank() == 0)
      print("   The selected components have excitation energies:");
    if (world.rank() == 0) print(energies);
  }

  // Done
  return answer;
}

// Calculate the exponentiation of a matrix through first order (I think)
Tensor<double> TDDFT::matrix_exponential(const Tensor<double>& A) {
  const double tol = 1e-13;
  MADNESS_ASSERT(A.dim((0) == A.dim(1)));

  // Scale A by a power of 2 until it is "small"
  double anorm = A.normf();
  size_t n = 0;
  double scale = 1.0;
  while (anorm * scale > 0.1) {
    ++n;
    scale *= 0.5;
  }
  Tensor<double> B = scale * A;  // B = A*2^-n

  // Compute exp(B) using Taylor series
  Tensor<double> expB = Tensor<double>(2, B.dims());
  for (int64_t i = 0; i < expB.dim(0); ++i) expB(i, i) = 1.0;

  size_t k = 1;
  Tensor<double> term = B;
  while (term.normf() > tol) {
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
/// @param[in]  thresh_degenerate       threshold for orbitals being
/// degenerate
/// @return             the unitary matrix U: U^T F U = evals
Tensor<double> TDDFT::get_fock_transformation(World& world,
                                              Tensor<double>& overlap,
                                              Tensor<double>& fock,
                                              Tensor<double>& evals,
                                              const double thresh_degenerate) {
  // Run an SVD on the overlap matrix and ignore values
  // less than thresh_degenerate
  Tensor<double> r_vecs;
  Tensor<double> s_vals;
  Tensor<double> l_vecs;
  Tensor<double> overlap_copy = copy(overlap);
  svd(overlap_copy, l_vecs, s_vals, r_vecs);

  // Debugging output
  if (Rparams.print_level >= 2 and world.rank() == 0) {
    print("\n   Singular values of overlap matrix:");
    print(s_vals);
    print("   Left singular vectors of overlap matrix:");
    print(l_vecs);
  }

  // Check how many singular values are less than 10*thresh_degen
  size_t num_sv = 0;
  for (int64_t i = 0; i < s_vals.dim(0); i++) {
    if (s_vals(i) < 10 * thresh_degenerate) {
      if (world.rank() == 0 and num_sv == 0) print("");
      if (world.rank() == 0)
        printf(
            "   Detected singular value (%.8f) below threshold (%.8f). "
            "Reducing subspace size.\n",
            s_vals(i),
            10 * thresh_degenerate);
      num_sv++;
    }
    if (world.rank() == 0 and i == s_vals.dim(0) - 1 and num_sv > 0) print("");
  }

  // Going to use these a lot here, so just calculate them
  size_t size_l = s_vals.dim(0);
  size_t size_s = size_l - num_sv;
  Tensor<double> l_vecs_s(size_l, num_sv);

  // Transform into this smaller space if necessary
  if (num_sv > 0) {
    // Cut out the singular values that are small
    // (singular values come out in descending order)
    overlap = Tensor<double>(size_s, size_s);
    for (size_t i = 0; i < size_s; i++) overlap(i, i) = s_vals(i);

    // Copy the active vectors to a smaller container
    l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

    // Debugging output
    if (Rparams.print_level >= 2 and world.rank() == 0) {
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

  int64_t nmo = fock.dim(0);  // NOLINT

  bool switched = true;
  while (switched) {
    switched = false;
    for (int64_t i = 0; i < nmo; i++) {
      for (int64_t j = i + 1; j < nmo; j++) {
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
  for (int64_t i = 0; i < nmo; ++i)  // NOLINT
    if (U(i, i) < 0.0) U(_, i).scale(-1.0);

  // Rotations between effectively degenerate components confound
  // the non-linear equation solver ... undo these rotations
  int64_t ilo = 0;  // first element of cluster NOLINT
  while (ilo < nmo - 1) {
    int64_t ihi = ilo;  // NOLINT
    while (fabs(evals[ilo] - evals[ihi + 1]) <
           thresh_degenerate * 100.0 * std::max(fabs(evals[ilo]), 1.0)) {
      ++ihi;
      if (ihi == nmo - 1) break;
    }
    int64_t nclus = ihi - ilo + 1;  // NOLINT
    if (nclus > 1) {
      Tensor<double> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));

      // Polar Decomposition
      Tensor<double> VH(nclus, nclus);
      Tensor<double> W(nclus, nclus);
      Tensor<double> sigma(nclus);

      svd(q, W, sigma, VH);
      q = transpose(inner(W, VH));  // Should be conj. tranpose if complex
      U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
    }
    ilo = ihi + 1;
  }

  fock = 0;
  for (unsigned int i = 0; i < nmo; ++i) fock(i, i) = evals(i);

  // If we transformed into the smaller subspace, time to transform back
  if (num_sv > 0) {
    // Temp. storage
    Tensor<double> temp_U(size_l, size_l);
    Tensor<double> U2(size_l, size_l);

    // Copy U back to larger size
    temp_U(Slice(0, size_s - 1), Slice(0, size_s - 1)) = copy(U);
    for (size_t i = size_s; i < size_l; i++) temp_U(i, i) = 1.0;

    // Transform back
    mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());

    U = copy(U2);
  }

  return U;
}

// Sorts the given Tensor of energies
Tensor<int> TDDFT::sort_eigenvalues(World& world,
                                    Tensor<double>& vals,
                                    Tensor<double>& vecs) {
  // Get relevant sizes
  size_t k = vals.size();

  // Tensor to hold selection order
  Tensor<int> selected(k);

  // Copy everything...
  std::vector<double> vals_copy;
  for (size_t i = 0; i < k; i++) vals_copy.push_back(vals[i]);
  Tensor<double> vals_copy2 = copy(vals);
  Tensor<double> vecs_copy = copy(vecs);

  // Now sort vals_copy
  std::sort(vals_copy.begin(), vals_copy.end());

  // Now sort the rest of the things, using the sorted energy list
  // to find the correct indices
  for (size_t i = 0; i < k; i++) {
    // Find matching index in sorted vals_copy
    size_t j = 0;
    while (fabs(vals_copy[i] - vals_copy2[j]) > 1e-8 && j < k) j++;

    // Add in to list which one we're taking
    selected(i) = j;

    // Put corresponding things in the correct place
    vals(i) = vals_copy[i];
    vecs(_, i) = vecs_copy(_, j);

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
Tensor<double> TDDFT::diagonalizeFockMatrix(World& world,
                                            Tensor<double>& fock,
                                            response_space& psi,
                                            ElectronResponseFunctions& I,
                                            Tensor<double>& evals,
                                            Tensor<double>& overlap,
                                            const double thresh) {
  // compute the unitary transformation matrix U that diagonalizes
  // the fock matrix
  Tensor<double> U =
      get_fock_transformation(world, overlap, fock, evals, thresh);

  // Sort into ascending order
  Tensor<int> selected = sort_eigenvalues(world, evals, U);

  // Debugging output
  if (Rparams.print_level >= 2 and world.rank() == 0) {
    print("   U:");
    print(U);
  }

  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // transform the orbitals and the potential
  // Truncate happens inside here
  I.Vx = transform(world, I.Vx, U);
  I.Hx = transform(world, I.Hx, U);
  I.EpsilonX = transform(world, I.EpsilonX, U);
  I.EpsilonXNoDiag = transform(world, I.EpsilonXNoDiag, U);
  psi = transform(world, psi, U);

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Transform orbs.:");

  // Normalize x
  normalize(world, psi);

  // Debugging output
  if (Rparams.print_level >= 2 and world.rank() == 0) {
    print("   Eigenvector coefficients from diagonalization:");
    print(U);
  }

  return U;
}

// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
response_space TDDFT::transform(World& world,
                                response_space& f,
                                Tensor<double>& U) {
  // Return container
  response_space result;

  // Go element by element
  for (unsigned int i = 0; i < f.size(); i++) {
    // Temp for the result of one row
    std::vector<real_function_3d> temp =
        zero_functions_compressed<double, 3>(world, f[0].size());

    for (unsigned int j = 0; j < f.size(); j++) {
      gaxpy(world, 1.0, temp, U(j, i), f[j]);
    }

    // Add to temp to result
    result.push_back(temp);
  }

  result.truncate_rf();

  // Done
  return result;
}

// If using a larger subspace to diagonalize in, this will put everything in
// the right spot
/**
 * @brief Diagonolize in larger subspace x
 * 
 * @param world 
 * @param S <x|x> 
 * @param A <x|Ax>
 * @param Current Ax
 * @param Last Axold
 * @param x_response 
 * @param old_S 
 * @param old_A 
 * @param old_x_response 
 * @param print_level 
 */
void TDDFT::augment(World& world,
                    Tensor<double>& S,
                    Tensor<double>& A,
                    ElectronResponseFunctions& Current,
                    ElectronResponseFunctions& Last,
                    response_space& x_response,
                    // Contains fock and energy scaled orbitals
                    Tensor<double>& old_S,
                    Tensor<double>& old_A,
                    response_space& old_x_response,
                    size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = x_response.size();
  // Create work space, will overwrite S and A in the end
  Tensor<double> temp_S(2 * m, 2 * m);
  Tensor<double> temp_A(2 * m, 2 * m);
  /**
   * @brief Need to create off diagonal blocks of A
   *  A=
   *  [xAx      xAx_old ]
   *  [xoldAx  xoldAxold]
   * 
   */
  // Need to create off diagonal blocks of A, so
  response_space TempCurrentAf = Current.Hx + Current.F0_x - Current.EpsilonX;// Ax
  response_space TempLastAf = Last.Hx + Last.F0_x - Last.EpsilonX;//Ax_old
  // Calculate correct inner products of upper off diagonal
  Tensor<double> off = expectation(world, x_response, TempLastAf);//xAxold
  temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);//top right
  // Now for lower off diagonal
  off = expectation(world, old_x_response, TempCurrentAf);//x_oldAx
  temp_A(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);//bottom left
  temp_A(Slice(0, m - 1), Slice(0, m - 1)) = copy(A);// xAx top left
  temp_A(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_A);// xoldAxold bottom right
  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("   Before symmeterizing A:");
    print(temp_A);
  }
  // Save temp_A as A_x
  // Need to symmeterize A as well (?)
  A = 0.5 * (temp_A + transpose(temp_A));
  /**
   * @brief Creating S
   * S= [<x|x>    <x|xold>   ]
   *    [<xold|x> <xold|xold>]
   */
  // Now create upper off diagonal block of S
  off = expectation(world, x_response, old_x_response);
  // Use slicing to put in correct spot
  temp_S(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);// top right <x|xold>
  // Now the lower off diagonal block
  // (Go ahead and cheat and use the transpose...)
  off = transpose(off);// just transpose <xold|x>
  // Use slicing to put in correct spot
  temp_S(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);// bottom right <xold|x>
  // Put together the rest of S
  temp_S(Slice(0, m - 1), Slice(0, m - 1)) = copy(S);//top left <x|x>
  temp_S(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_S);//<xold|xold>
  // Save temp_S as S_x
  S = copy(temp_S);
  // Add in old vectors to current vectors for the appropriate ones
  for (size_t i = 0; i < m; i++) {
    Current.Hx.push_back(Last.Hx[i]);
    x_response.push_back(old_x_response[i]);
    Current.EpsilonX.push_back(Last.EpsilonX[i]);
    Current.F0_x.push_back(Last.F0_x[i]);
    Current.Vx.push_back(Last.Vx[i]);
  }

  // End the timer
  if (print_level >= 1) molresponse::end_timer(world, "Aug. resp. matrix:");

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("\n   Augmented response matrix:");
    print(A);
  }

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("   Augmented overlap matrix:");
    print(S);
  }

  // SUPER debugging
  if (print_level >= 3) {
    if (world.rank() == 0)
      print("   Calculating condition number of aug. response matrix");
  }
}

// If using a larger subspace to diagonalize in, this will put everything in
// the right spot
void TDDFT::augment_full(World& world,
                         Tensor<double>& S,
                         Tensor<double>& A,
                         response_space& B_x,
                         response_space& x_gamma,
                         response_space& x_response,
                         response_space& V_x_response,
                         response_space& x_fe,  // Contains V_x_response
                         response_space& B_y,
                         response_space& y_gamma,
                         response_space& y_response,
                         response_space& V_y_response,
                         response_space& y_fe,  // Contains V_y_response
                         Tensor<double>& old_S,
                         Tensor<double>& old_A,
                         response_space& old_B_x,
                         response_space& old_x_gamma,
                         response_space& old_x_response,
                         response_space& old_V_x_response,
                         response_space& old_x_fe,
                         response_space& old_B_y,
                         response_space& old_y_gamma,
                         response_space& old_y_response,
                         response_space& old_V_y_response,
                         response_space& old_y_fe,
                         size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = x_gamma.size();

  // Create work space, will overwrite S and A in the end
  Tensor<double> temp_S(2 * m, 2 * m);
  Tensor<double> temp_A(2 * m, 2 * m);

  // Need to create off diagonal blocks of A, so
  // create temps that are the sums of current and
  // old components respectively
  response_space temp_cur_x = x_gamma + x_fe + B_y;
  response_space temp_cur_y = y_gamma + y_fe + B_x;
  response_space temp_old_x = old_x_gamma + old_x_fe + old_B_y;
  response_space temp_old_y = old_y_gamma + old_y_fe + old_B_x;

  // Calculate correct inner products of upper off diagonal
  Tensor<double> off = expectation(world, x_response, temp_old_x) +
                       expectation(world, y_response, temp_old_y);

  // Use slicing to put in correct spot
  temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);
  // temp_A(Slice(m, 2*m-1), Slice(0, m-1)) = copy(off);

  // Now for lower off diagonal
  off = expectation(world, old_x_response, temp_cur_x) +
        expectation(world, old_y_response, temp_cur_y);

  // Use slicing to put in correct spot
  temp_A(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);
  // temp_A(Slice(0, m-1), Slice(m, 2*m-1)) = copy(off);

  // Put together the rest of A
  temp_A(Slice(0, m - 1), Slice(0, m - 1)) = copy(A);
  temp_A(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_A);

  // Save temp_A into A
  A = copy(temp_A);

  // Now create upper off diagonal block of S
  off = expectation(world, x_response, old_x_response) -
        expectation(world, y_response, old_y_response);

  // Use slicing to put in correct spot
  temp_S(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);

  // Now the lower off diagonal block
  // (Cheating by using the transpose...)
  off = transpose(off);

  // Use slicing to put in correct spot
  temp_S(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);

  // Put together the rest of S
  temp_S(Slice(0, m - 1), Slice(0, m - 1)) = copy(S);
  temp_S(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_S);

  // Save temp_S as S
  S = copy(temp_S);

  // Finally, add in old vectors to current vectors for the appropriate ones
  for (size_t i = 0; i < m; i++) {
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
  if (print_level >= 1) molresponse::end_timer(world, "Aug. resp. matrix:");

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("\n   Augmented response matrix:");
    print(A);
  }

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("   Augmented overlap matrix:");
    print(S);
  }
}

// If using a larger subspace to diagonalize in, after diagonalization this
// will put everything in the right spot
void TDDFT::unaugment(World& world,
                      size_t num_states,
                      size_t iter,
                      Tensor<double>& omega,
                      Tensor<double>& S_x,
                      Tensor<double>& A_x,
                      ElectronResponseFunctions& Current,
                      ElectronResponseFunctions& Last,
                      response_space& x_response,
                      Tensor<double>& old_S,
                      Tensor<double>& old_A,
                      response_space& old_x_response,
                      size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Note: the eigenvalues and vectors were sorted after diagonalization
  // and hence all the functions are sorted in ascending order of energy

  // Quick copy of m lowest eigenvalues
  omega = omega(Slice(0, num_states - 1));

  // Pop off the "m" vectors off the back end of appropriate vectors
  // (only after first iteration)
  if (iter > 0) {
    for (size_t i = 0; i < num_states; i++) {
      //  x_fe.pop_back();
      // V_x_response.pop_back();
      // x_gamma.pop_back();
      // x_response.pop_back();
      Current.Hx.pop_back();
      x_response.pop_back();
      Current.EpsilonX.pop_back();
      Current.F0_x.pop_back();
      Current.Vx.pop_back();
    }
  }
  Last.Hx = Current.Hx.copy();
  Last.EpsilonX = Current.EpsilonX.copy();
  Last.F0_x = Current.F0_x.copy();
  Last.Vx = Current.Vx.copy();

  // Save the "current" into "old"
  // old_x_fe = x_fe.copy();
  // old_x_gamma = x_gamma.copy();
  // old_V_x_response = V_x_response.copy();

  // Project out ground state
  // QProjector<double, 3> projector(world, Gparams.orbitals);
  // for(size_t i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
  old_x_response = x_response.copy();

  // Copy S into old_S
  // S is identity
  old_S = expectation(world, x_response, x_response);

  // Copy A into old_A
  // A is (nearly?) diagonal
  old_A = Tensor<double>(num_states, num_states);
  for (size_t i = 0; i < num_states; i++) old_A(i, i) = omega(i);

  // End the timer
  if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}

// If using a larger subspace to diagonalize in, after diagonalization this
// will put everything in the right spot
void TDDFT::unaugment_full(World& world,
                           size_t m,
                           size_t iter,
                           Tensor<double>& U,
                           Tensor<double>& omega,
                           Tensor<double>& S,
                           Tensor<double>& A,
                           response_space& x_gamma,
                           response_space& x_response,
                           response_space& V_x_response,
                           response_space& x_fe,
                           response_space& B_x,
                           response_space& y_gamma,
                           response_space& y_response,
                           response_space& V_y_response,
                           response_space& y_fe,
                           response_space& B_y,
                           Tensor<double>& old_S,
                           Tensor<double>& old_A,
                           response_space& old_x_gamma,
                           response_space& old_x_response,
                           response_space& old_V_x_response,
                           response_space& old_x_fe,
                           response_space& old_B_x,
                           response_space& old_y_gamma,
                           response_space& old_y_response,
                           response_space& old_V_y_response,
                           response_space& old_y_fe,
                           response_space& old_B_y,
                           size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Note: the eigenvalues and vectors were sorted after diagonalization
  // and hence all the functions are sorted in ascending order of energy

  // Quick copy of m lowest eigenvalues
  omega = omega(Slice(0, m - 1));

  // Pop off the "m" vectors off the back end of appropriate vectors
  // (only after first iteration)
  if (iter > 0) {
    for (size_t i = 0; i < m; i++) {
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
  // QProjector<double, 3> projector(world, Gparams.orbitals);
  // for(size_t i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
  // for(size_t i = 0; i < m; i++) y_response[i] = projector(y_response[i]);
  old_x_response = x_response.copy();
  old_y_response = y_response.copy();

  // Now to pull out correct values from S and A (both are size 2*m by 2*m,
  // and only want m by m values)
  old_S = expectation(world, x_response, x_response) -
          expectation(world, y_response, y_response);

  // And construct old_A
  response_space t1 = old_x_fe + old_x_gamma + old_B_y;
  response_space t2 = old_y_fe + old_y_gamma + old_B_x;
  old_A = expectation(world, old_x_response, t1) +
          expectation(world, old_y_response, t2);

  // End the timer
  if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}

// Diagonalize the full response matrix, taking care of degenerate components
// Why diagonalization and then transform the x_fe vectors

Tensor<double> TDDFT::diagonalizeFullResponseMatrix(
    World& world,
    Tensor<double>& S,
    Tensor<double>& A,
    response_space& x,
    response_space& y,
    ElectronResponseFunctions& I,
    Tensor<double>& omega,
    const double thresh,
    size_t print_level) {
  // compute the unitary transformation matrix U that diagonalizes
  // the response matrix
  Tensor<double> vecs =
      GetFullResponseTransformation(world, S, A, omega, thresh);

  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Transform the vectors of functions
  // Truncate happens in here
  I.Hx = transform(world, I.Hx, vecs);
  I.Hy = transform(world, I.Hy, vecs);
  // mixing terms
  I.Gx = transform(world, I.Gx, vecs);
  I.Gy = transform(world, I.Gy, vecs);
  // Potentials
  I.Vx = transform(world, I.Vx, vecs);
  I.Vy = transform(world, I.Vy, vecs);
  // functions
  x = transform(world, x, vecs);
  y = transform(world, y, vecs);
  //
  I.F0_x = transform(world, I.F0_x, vecs);
  I.F0_y = transform(world, I.F0_y, vecs);
  // epsilonX
  I.EpsilonX = transform(world, I.EpsilonX, vecs);
  I.EpsilonY = transform(world, I.EpsilonY, vecs);
  I.EpsilonXNoDiag = transform(world, I.EpsilonXNoDiag, vecs);
  I.EpsilonYNoDiag = transform(world, I.EpsilonYNoDiag, vecs);
  // we do transform here
  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Transform orbs.:");

  // Normalize x and y
  normalize(world, x, y);

  // Debugging output
  if (world.rank() == 0 and print_level >= 2) {
    print("   Eigenvector coefficients from diagonalization:");
    print(vecs);
  }

  // Return the selected functions
  return vecs;
}

// Similar to what robert did above in "get_fock_transformation"
Tensor<double> TDDFT::GetFullResponseTransformation(
    World& world,
    Tensor<double>& S,
    Tensor<double>& A,
    Tensor<double>& evals,
    const double thresh_degenerate) {
  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  // Get size
  size_t m = S.dim(0);

  // Run an SVD on the overlap matrix and ignore values
  // less than thresh_degenerate
  Tensor<double> r_vecs, s_vals, l_vecs;
  Tensor<double> S_copy = copy(S);
  /**
   * @brief SVD on overlap matrix S
   * S=UsVT
   * S, U, s , VT
   * 
   */
  svd(S_copy, l_vecs, s_vals, r_vecs);

  // Debugging output
  if (Rparams.print_level >= 2 and world.rank() == 0) {
    print("\n   Singular values of overlap matrix:");
    print(s_vals);
    print("   Left singular vectors of overlap matrix:");
    print(l_vecs);
  }

  // Check how many singular values are less than 10*thresh_degen
  size_t num_sv = 0;
  for (int64_t i = 0; i < s_vals.dim(0); i++) {
    if (s_vals(i) < 10 * thresh_degenerate) {
      if (world.rank() == 0 and num_sv == 0) print("");
      if (world.rank() == 0)
        printf(
            "   Detected singular value (%.8f) below threshold (%.8f). "
            "Reducing subspace size.\n",
            s_vals(i),
            10 * thresh_degenerate);
      num_sv++;
    }
    if (world.rank() == 0 and i == s_vals.dim(0) - 1 and num_sv > 0) print("");
  }

  // Going to use these a lot here, so just calculate them
  size_t size_l = s_vals.dim(0);    // number of singular values
  size_t size_s = size_l - num_sv;  // smaller subspace size
  /**
   * @brief l_vecs_s(m,1)
   * 
   * @return Tensor<double> 
   */
  Tensor<double> l_vecs_s(
      size_l,
      num_sv);                     // number of sv by number smaller than thress
  Tensor<double> copyA = copy(A);  // we copy xAx

  // Transform into this smaller space if necessary
  if (num_sv > 0) {
    // Cut out the singular values that are small
    // (singular values come out in descending order)

    // S(m-sl,m-sl)
    S = Tensor<double>(size_s, size_s);  // create size of new size
    for (size_t i = 0; i < size_s; i++) S(i, i) = s_vals(i);
    // Copy the active vectors to a smaller container
    // left vectors [m,m-sl]
    l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

    // Debugging output
    if (Rparams.print_level >= 2 and world.rank() == 0) {
      print("   Reduced size left singular vectors of overlap matrix:");
      print(l_vecs_s);
    }

    // Transform
    // Work(m,m-sl)
    Tensor<double> work(size_l, size_s);
        /*
          c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)
          
          where it is assumed that the last index in each array is has unit
          stride and the dimensions are as provided.
          
          4-way unrolled k loop ... empirically fastest on PIII
          compared to 2/3 way unrolling (though not by much).
        */
   // dimi,dimj,dimk,c,a,b 
    mxm(size_l, size_s, size_l, work.ptr(), A.ptr(), l_vecs_s.ptr());
    // A*left
    copyA = Tensor<double>(size_s, size_s);
    Tensor<double> l_vecs_t = transpose(l_vecs);
    // s s l, copyA=lvect_t*A*left
    mxm(size_s, size_s, size_l, copyA.ptr(), l_vecs_t.ptr(), work.ptr());

    // Debugging output
    if (Rparams.print_level >= 2 and world.rank() == 0) {
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
  if (world.rank() == 0 and Rparams.print_level >= 2)
    print("\n   Max imaginary component of eigenvalues:", max_imag, "\n");
  MADNESS_ASSERT(max_imag <= 1e-5);  // MUST BE REAL!
  evals = real(omega);

  // Easier to just resize here
  m = evals.dim(0);

  bool switched = true;
  while (switched) {
    switched = false;
    for (size_t i = 0; i < m; i++) {
      for (size_t j = i + 1; j < m; j++) {
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
  for (size_t i = 0; i < m; ++i)
    if (U(i, i) < 0.0) U(_, i).scale(-1.0);

  // Rotations between effectively degenerate components confound
  // the non-linear equation solver ... undo these rotations
  size_t ilo = 0;  // first element of cluster
  while (ilo < m - 1) {
    size_t ihi = ilo;
    while (fabs(evals[ilo] - evals[ihi + 1]) <
           thresh_degenerate * 10.0 * std::max(fabs(evals[ilo]), 1.0)) {
      ++ihi;
      if (ihi == m - 1) break;
    }
    int64_t nclus = ihi - ilo + 1;
    if (nclus > 1) {
      Tensor<double> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));

      // Polar Decomposition
      Tensor<double> VH(nclus, nclus);
      Tensor<double> W(nclus, nclus);
      Tensor<double> sigma(nclus);

      svd(q, W, sigma, VH);
      q = transpose(inner(W, VH));  // Should be conj. tranpose if complex
      U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
    }
    ilo = ihi + 1;
  }

  // If we transformed into the smaller subspace, time to transform back
  if (num_sv > 0) {
    // Temp. storage
    Tensor<double> temp_U(size_l, size_l);
    Tensor<double> U2(size_l, size_l);

    // Copy U back to larger size
    temp_U(Slice(0, size_s - 1), Slice(0, size_s - 1)) = copy(U);
    for (size_t i = size_s; i < size_l; i++) temp_U(i, i) = 1.0;

    // Transform U back
    mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());
    U = copy(U2);
  }

  // Sort into ascending order
  Tensor<int> selected = sort_eigenvalues(world, evals, U);

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "Diag. resp. mat.");

  return U;
}

// Sorts the given tensor of eigenvalues and
// response functions
void TDDFT::sort(World& world, Tensor<double>& vals, response_space& f) {
  // Get relevant sizes
  size_t k = vals.size();

  // Copy everything...
  response_space f_copy(f);
  Tensor<double> vals_copy = copy(vals);
  Tensor<double> vals_copy2 = copy(vals);

  // Now sort vals_copy
  std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

  // Now sort the rest of the things, using the sorted energy list
  // to find the correct indices
  for (size_t i = 0; i < k; i++) {
    // Find matching index in sorted vals_copy
    size_t j = 0;
    while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

    // Put corresponding function, difference function, value residual and
    // value in the correct place
    f[i] = f_copy[j];
    vals(i) = vals_copy(i);

    // Change the value of vals_copy2[j] to help deal with duplicates?
    vals_copy2(j) = 10000.0;
  }
}
/**
 * @brief Diagonalize AX=SX omega
 * 
 * @param world 
 * @param S 
 * @param old_S 
 * @param old_A 
 * @param x_response 
 * @param old_x_response 
 * @param ElectronResponses 
 * @param OldElectronResponses 
 * @param omega 
 * @param iteration 
 * @param m 
 */
void TDDFT::deflateTDA(World& world,
                       Tensor<double>& S,
                       Tensor<double> old_S,
                       Tensor<double> old_A,
                       response_space& x_response,
                       response_space& old_x_response,
                       ElectronResponseFunctions& ElectronResponses,
                       ElectronResponseFunctions& OldElectronResponses,
                       Tensor<double>& omega,
                       size_t& iteration,
                       size_t& m) {

  S = response_space_inner(x_response, x_response);

  // Debugging output
  if (Rparams.print_level >= 2 and world.rank() == 0) {
    print("   Overlap matrix:");
    print(S);
  }
  // Constructing response matrix
  Tensor<double> XAX = CreateResponseMatrix(world,
                                            x_response,
                                            ElectronResponses,
                                            Gparams.orbitals,
                                            Rparams.print_level,
                                            "x");

  // Augment S_x, A_x, x_gamma, x_response, V_x_response and x_gamma
  // if using a larger subspace and not iteration zero (TODO ---Gotta
  // look at this and make sure it uses my new functions molresponse )
  // by default Rparams.larger_subspace = 0 therefore never uses this
  if (iteration < Rparams.larger_subspace and iteration > 0) {
    print("Using augmented subspace");
    augment(world,
            S,
            XAX,
            ElectronResponses,
            OldElectronResponses,
            x_response,
            old_S,
            old_A,
            old_x_response,
            Rparams.print_level);
  }

  // Solve Ax = Sxw
  // Just to be sure dimensions work out, clear omega
  omega.clear();
  diagonalizeFockMatrix(world,
                        XAX,
                        x_response,
                        ElectronResponses,
                        omega,
                        S,
                        FunctionDefaults<3>::get_thresh());

  // If larger subspace, need to "un-augment" everything
  if (iteration < Rparams.larger_subspace) {
    print("Unaugmenting subspace");
    unaugment(world,
              m,
              iteration,
              omega,
              S,
              XAX,
              ElectronResponses,
              OldElectronResponses,
              x_response,
              old_S,
              old_A,
              old_x_response,
              Rparams.print_level);
  }
}

void TDDFT::deflateFull(World& world,
                        Tensor<double>& S,
                        Tensor<double> old_S,
                        Tensor<double> old_A,
                        response_space& x_response,
                        response_space& y,
                        response_space& old_x_response,
                        response_space& old_y_response,
                        ElectronResponseFunctions& ElectronResponses,
                        ElectronResponseFunctions& OldElectronResponses,
                        Tensor<double>& omega,
                        size_t& iteration,
                        size_t& m) {
  S = expectation(world, x_response, x_response) -
      expectation(world, y_response, y_response);

  // Debugging output
  if (world.rank() == 0 and Rparams.print_level >= 2) {
    print("\n   Overlap Matrix:");
    print(S);
  }

  // Construct full response matrix
  Tensor<double> A = CreateFullResponseMatrix(world,
                                              x_response,
                                              y_response,
                                              ElectronResponses,
                                              Gparams.orbitals,
                                              Rparams.small,
                                              FunctionDefaults<3>::get_thresh(),
                                              Rparams.print_level);

  // Larger subspace augmentation BROKEN!!!!!
  // if(iteration < Rparams.larger_subspace and iteration > 0)
  //{
  //   augment_full(world, S, A,
  //                B_x, x_gamma, x_response, V_x_response, x_fe,
  //                B_y, y_gamma, y_response, V_y_response, y_fe,
  //                old_S, old_A,
  //                old_B_x, old_x_gamma, old_x_response,
  //                old_V_x_response, old_x_fe, old_B_y, old_y_gamma,
  //                old_y_response, old_V_y_response, old_y_fe,
  //                Rparams.print_level);
  //}

  // Diagonalize
  // Just to be sure dimensions work out, clear omega
  omega.clear();

  Tensor<double> U =
      diagonalizeFullResponseMatrix(world,
                                    S,
                                    A,
                                    x_response,
                                    y_response,
                                    ElectronResponses,
                                    omega,
                                    FunctionDefaults<3>::get_thresh(),
                                    Rparams.print_level);
  // Larger subspace un-augmentation BROKEN!!!!
  // if(iteration < Rparams.larger_subspace)
  //{
  //   unaugment_full(world, m, iteration, U, omega, S, A,
  //                  x_gamma, x_response, V_x_response, x_fe, B_x,
  //                  y_gamma, y_response, V_y_response, y_fe, B_y,
  //                  old_S, old_A,
  //                  old_x_gamma, old_x_response, old_V_x_response,
  //                  old_x_fe, old_B_x, old_y_gamma, old_y_response,
  //                  old_V_y_response, old_y_fe, old_B_y,/
  //                  Rparams.print_level);
  //}
}
// const double thresh, int print_level) {
// Creates the XCOperator object and initializes it with correct
// parameters
XCOperator TDDFT::create_xcoperator(World& world,
                                    std::vector<real_function_3d> orbitals,
                                    std::string xc) {
  // First calculate the ground state density
  std::vector<real_function_3d> vsq =
      square(world, Gparams.orbitals);  // we square each orbital
  compress(world, vsq);  // compress into multiwavelet representation
  real_function_3d rho = real_factory_3d(world);  // create function rho
  rho.compress();
  for (unsigned int i = 0; i < vsq.size(); ++i) {
    rho.gaxpy(1.0, vsq[i], 1.0, false);
  }
  world.gop.fence();

  // And create the object using Gparams.xc
  XCOperator xcop(world,
                  xc,
                  false,
                  rho,
                  rho);  // world,which xc, spin_polarized? ,spinup, spindown

  return xcop;
}

// Uses an XCOperator to construct v_xc for the ground state density
// Returns d^2/d rho^2 E_xc[rho]
std::vector<real_function_3d> TDDFT::create_fxc(
    World& world,
    std::vector<real_function_3d>& orbitals,
    response_space& f,
    response_space& g) {
  // Create the xcop
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  // Next need the perturbed density
  std::vector<real_function_3d> drho =
      transition_density(world, orbitals, f, g);

  // Return container
  std::vector<real_function_3d> vxc;

  // Finally create the functions we want, one per response state
  // (xc_args_prep_response happens inside this call)
  for (unsigned int i = 0; i < f.size(); i++) {
    vxc.push_back(xc.apply_xc_kernel(drho[i]));
  }
  // for each density apply xckernel

  return vxc;
}

// Uses an XCOperator to construct v_xc for the ground state density
// Returns d^2/d rho^2 E_xc[rho]
std::vector<real_function_3d> TDDFT::GetWxcOnFDensities(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Create the xcop
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);
  // Next need the perturbed density
  std::vector<real_function_3d> drhoM =
      GetTransitionDensities(world, orbitals, f);
  // Return container
  std::vector<real_function_3d> Wfxc;
  // Finally create the functions we want, one per response state
  // (xc_args_prep_response happens inside this call)
  for (unsigned int i = 0; i < f.size(); i++) {
    Wfxc.push_back(xc.apply_xc_kernel(drhoM[i]));
  }
  // for each density apply xckernel

  return Wfxc;
}

std::vector<real_function_3d> TDDFT::GetConjugateWxcOnFDensities(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Create the xcop
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  // Next need the perturbed density
  std::vector<real_function_3d> drhoM =
      GetConjugateTransitionDensities(world, orbitals, f);
  // Return container
  std::vector<real_function_3d> conjugateWfvxc;

  // Finally create the functions we want, one per response state
  // (xc_args_prep_response happens inside this call)
  for (unsigned int i = 0; i < f.size(); i++) {
    conjugateWfvxc.push_back(xc.apply_xc_kernel(drhoM[i]));
  }
  // for each density apply xckernel

  return conjugateWfvxc;
}

std::vector<real_function_3d> TDDFT::CreateXCDerivative(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Create the xcop
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);
  size_t m = f.size();     // get the number of response functions
  size_t n = f[0].size();  // get the number of orbitals function

  std::vector<real_function_3d> drho = zero_functions<double, 3>(world, m);
  // Run over virtual...
  for (size_t i = 0; i < m; i++) {
    // Run over occupied...
    for (size_t j = 0; j < n; j++) {
      // y functions are zero if TDA is active
      drho[i] = drho[i] + orbitals[j] * f[i][j];  //+ orbitals[j] * y[i][j];
    }
  }
  // Return container
  std::vector<real_function_3d> vxc;

  // Finally create the functions we want, one per response state
  // (xc_args_prep_response happens inside this call)
  for (unsigned int i = 0; i < f.size(); i++) {
    vxc.push_back(xc.apply_xc_kernel(drho[i]));
  }
  // for each density apply xckernel

  return vxc;
}

// Iterates the response functions until converged or out of iterations
//

// Simplified iterate scheme for guesses
void TDDFT::IterateGuess(World& world, response_space& guesses) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
  QProjector<double, 3> projector(
      world, Gparams.orbitals);  // Projector to project out ground state
  size_t n = guesses[0].size();  // Number of ground state orbitals
  size_t m = guesses.size();     // number initial guess orbitals
  Tensor<double> shifts;         // Holds the shifted energy values
  ElectronResponseFunctions I;
  response_space bsh_resp;  // Holds wave function corrections
  response_space gamma;     // Holds the perturbed two electron piece
  response_space rhs;
  response_space fe;  // Holds the ground state-fock and energy scaled x
  // response components
  response_space V;  // Holds V^0 applied to response functions
  response_space
      shifted_V;          // Holds the shifted V^0 applied to response functions
  Tensor<double> S;       // Overlap matrix of response components for x states
  real_function_3d v_xc;  // For TDDFT

  // If DFT, initialize the XCOperator
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  // Useful to have
  response_space zeros(world, m, n);

  // Now to iterate
  while (iteration < Rparams.guess_max_iter) {
    // Start a timer for this iteration
    molresponse::start_timer(world);
    // (Eventually take symmtry into consideration
    //  We need a number of points that takes into account the number
    //  of degenerate states)
    // reduce the space progressively
    //  Starting number of response functions
    //
    size_t N0 = guesses.size();
    //  number of response functions after p-1 iterations
    //  I believe we don't iterate the last time
    size_t Np = 2 * Rparams.states;
    // this controls the final number of points to choose from
    // number of iterations
    size_t p = Rparams.guess_max_iter - 1;
    // Number of points after i iterations
    size_t Ni = N0;

    // No*exp(log(Np/N0)/p*t) to exponential decay
    // the number of states down to 2*Rparams.states
    if (iteration > 1 && Rparams.guess_xyz && Ni > Np) {
      Ni = std::ceil(N0 * std::exp(std::log(static_cast<double>(Np) /
                                            static_cast<double>(N0)) /
                                   static_cast<double>(p) * iteration));
      sort(world, omega, guesses);
      print(omega);
      for (size_t i = 0; i < Ni; i++) {
        Ni = 0;
        if (omega[i] < 1) {
          Ni++;
        }
      }
      // this function selects k functions
      guesses =
          select_functions(world, guesses, omega, Ni, Rparams.print_level);
    }

    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0)
        printf("\n   Guess Iteration %d at time %.1fs\n",
               static_cast<int>(iteration),
               wall_time());
      if (world.rank() == 0) print(" -------------------------------------");
    }

    // Load balance
    // Only balancing on x-components. Smart?
    if (world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0)) and

        iteration != 0) {
      // Start a timer
      if (Rparams.print_level >= 1) molresponse::start_timer(world);
      if (world.rank() == 0) print("");  // Makes it more legible

      LoadBalanceDeux<3> lb(world);
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < Rparams.states; k++) {
          lb.add_tree(guesses[k][j], lbcost<double, 3>(1.0, 8.0), true);
          // lb.add_tree(V[k][j], lbcost<double,3>(1.0,8.0), true);
          // lb.add_tree(gamma[k][j], lbcost<double,3>(1.0,8.0), true);
        }
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, "Load balancing:");
    }

    // Project out ground state
    for (size_t i = 0; i < Ni; i++) guesses[i] = projector(guesses[i]);

    // Truncate before doing expensive things
    guesses.truncate_rf();

    // Normalize after projection
    if (Rparams.tda) normalize(world, guesses);
    // (TODO why not normalize if not tda)

    // Create gamma
    //    gamma = CreateGamma(world, guesses, zeros, Gparams.orbitals,
    //    Rparams.small,
    //                        FunctionDefaults<3>::get_thresh(),
    //                        Rparams.print_level,
    //                       "x");
    computeElectronResponse(world,
                            I,
                            guesses,
                            zeros,
                            Gparams.orbitals,
                            xc,
                            hamiltonian,
                            ham_no_diag,
                            Rparams.small,
                            FunctionDefaults<3>::get_thresh(),
                            Rparams.print_level,
                            "x");
    // Constructing S
    S = expectation(world, x_response, x_response);
    // Debugging output
    if (Rparams.print_level >= 2 and world.rank() == 0) {
      print("\n   Overlap matrix:");
      print(S);
    }
    Tensor<double> A = CreateResponseMatrix(
        world, guesses, I, Gparams.orbitals, Rparams.print_level, "x");

    diagonalizeFockMatrix(
        world, A, guesses, I, omega, S, FunctionDefaults<3>::get_thresh());

    // Ensure right number of omegas
    if (size_t(omega.dim(0)) != Ni) {
      if (world.rank() == 0)
        print("\n   Adding",
              Ni - omega.dim(0),
              "eigenvalue(s) (counters subspace size "
              "reduction in "
              "diagonalizatoin).");
      Tensor<double> temp(Ni);
      temp(Slice(0, omega.dim(0) - 1)) = omega;
      for (size_t i = omega.dim(0); i < Ni; i++) temp[i] = 2.5 * i;
      omega = copy(temp);
    }

    // Basic output
    if (Rparams.print_level >= 1 and world.rank() == 0) {
      print("\n   Excitation Energies:");
      print("gi=", iteration, " roots: ", omega);
    }

    // Only do BSH if not the last iteration
    if (iteration + 1 < Rparams.guess_max_iter) {
      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      shifts = create_shift(
          world, Gparams.energies, omega, Rparams.print_level, "x");

      // Apply the shifts
      shifted_V = apply_shift(world, shifts, I.Vx, guesses);

      // Construct RHS of equation
      //      ResponseFunction rhs = gamma + shifted_V;
      rhs = I.Hx;
      for (size_t i = 0; i < Ni; i++) {
        rhs[i] = projector(rhs[i]);
      }

      rhs = rhs - I.EpsilonXNoDiag + shifted_V;

      // Add in all off diagonal elements of ground state Fock
      // matrix
      //           rhs = rhs - g_fe2;

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
          bsh_operators =
              create_bsh_operators(world,
                                   shifts,
                                   Gparams.energies,
                                   omega,
                                   Rparams.small,
                                   FunctionDefaults<3>::get_thresh());

      // Scale by -2.0 (coefficient in eq. 37 of reference
      // paper)

      // Apply BSH and get updated components
      if (Rparams.print_level >= 1) molresponse::start_timer(world);
      bsh_resp = apply(world, bsh_operators, rhs);
      if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply BSH:");

      // Project out ground state
      for (size_t i = 0; i < Ni; i++) bsh_resp[i] = projector(bsh_resp[i]);

      // Save new components
      guesses = bsh_resp;
      guesses = guesses * -2.0;  // scale(guesses, -2.0);
      // Apply mask
      for (size_t i = 0; i < Ni; i++) guesses[i] = mask * guesses[i];
    }

    // Update counter
    iteration += 1;

    // Done with the iteration.. truncate
    guesses.truncate_rf();

    // Basic output
    if (Rparams.print_level >= 1) {  //
      molresponse::end_timer(world, " This iteration:");
    }
  }
}  // Done with iterate guess

// Create and diagonalize the CIS matrix for improved initial guess
response_space TDDFT::diagonalize_CIS_guess(
    World& world,
    std::vector<real_function_3d>& virtuals,
    Tensor<double>& omega,
    std::vector<real_function_3d>& orbitals,
    Tensor<double>& energies,
    double small,
    double thresh,
    size_t print_level) {
  // Projecter for removing ground state
  QProjector<double, 3> Q(world, orbitals);

  // Diagonalize under ground state hamiltonian a few times
  // Create overlap
  compress(world, virtuals);
  Tensor<double> S = matrix_inner(world, virtuals, virtuals);

  // Create Fock matrix
  // -1 suppresses output
  Tensor<double> Fmat = CreateGroundHamiltonian(world, virtuals, -1);

  // Diagonalize
  Tensor<double> U, evals, dummy(virtuals.size());
  U = get_fock_transformation(world, S, Fmat, evals, thresh);

  // Transform and truncate functions
  virtuals = madness::transform(world, virtuals, U);
  truncate(world, virtuals, thresh, false);

  // filter out any ground state functions that crept in
  std::vector<real_function_3d> true_virtuals;
  for (unsigned int a = 0; a < virtuals.size(); a++) {
    if (evals(a) > 0.0) true_virtuals.push_back(virtuals[a]);
  }

  // Make sure we still have functions
  if (true_virtuals.empty())
    MADNESS_EXCEPTION(
        "Virtuals are empty: Too much overlap with occupied orbitals", 1);

  // Saving new components
  virtuals = Q(true_virtuals);

  // Debugging output
  if (print_level >= 2 and world.rank() == 0)
    print("   Remaining virtuals:", virtuals.size());

  // Now make the CIS matrix (copied from Jakob)
  if (world.rank() == 0)
    print("   Forming CIS matrix for improved initial guess.");

  // Start timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  size_t I = -1;  // combined index from i and a, start is -1 so that initial
                  // value
  // is 0
  size_t J = -1;  // combined index from j and b, start is -1 so that initial
                  // value
  // is 0

  const size_t m = virtuals.size();
  const size_t n = orbitals.size();

  Tensor<double> MCIS(m * n, m * n);
  real_convolution_3d op = CoulombOperator(world, small, thresh);

  for (size_t i = 0; i < n; i++) {
    const real_function_3d brai = orbitals[i];
    const std::vector<real_function_3d> igv = apply(world, op, virtuals * brai);
    const std::vector<real_function_3d> igm = apply(world, op, orbitals * brai);

    for (size_t a = 0; a < m; a++) {
      I++;
      J = -1;
      for (size_t j = 0; j < n; j++) {
        const real_function_3d braj = orbitals[j];

        for (size_t b = 0; b < m; b++) {
          J++;
          double diag_element = 0.0;

          if (i == j and a == b) diag_element = Fmat(a, a) - energies(i);

          MCIS(I, J) = diag_element + 2.0 * inner(braj * virtuals[b], igv[a]) -
                       inner(virtuals[a] * virtuals[b], igm[j]);
        }
      }
    }
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "Form CIS matrix:");

  // Debugging output
  if (print_level >= 2 and world.rank() == 0) {
    print("   CIS matrix:");
    print(MCIS);
  }

  // Diagonalize CIS matrix
  syev(MCIS, U, evals);

  // Always print this?
  if (world.rank() == 0) {
    print("   Initial CIS matrix eigenvalues:");
    print(evals);
  }

  // Now construct the initial guess
  response_space f(world, U.dim(0), n);
  omega = Tensor<double>(U.dim(0));

  I = -1;
  for (size_t i = 0; i < n; i++) {
    for (size_t a = 0; a < m; a++) {
      I++;
      J = -1;
      if (evals(I) < 0.0) {
        if (world.rank() == 0) print("   Skipping negative root:", evals(I));
        continue;
      }
      for (size_t j = 0; j < n; j++) {
        for (size_t b = 0; b < m; b++) {
          J++;
          f[I][j] += U(J, I) * virtuals[b];
          omega(I) = evals(I);
        }
      }
    }
  }

  f.truncate_rf();

  // Done. Whew.
  return f;
}

// Adds in random noise to a vector of vector of functions
response_space TDDFT::add_randomness(World& world,
                                     response_space& f,
                                     double magnitude) {
  // Copy input functions
  response_space f_copy = f.copy();

  // Lambda function to add in noise
  auto noise = [](const Key<3>& key, Tensor<double>& x) mutable {
    Tensor<double> y(x.size());
    y.fillrandom();
    // y.scale(magnitude);
    y.scale(1e3);
    x = x + y;
    // x(0,0,0) += y(0,0,0)-0.5;
  };
  // TODO
  // Go through each function in f_copy and add in random noise
  for (unsigned int i = 0; i < f_copy.size(); i++) {
    for (unsigned int j = 0; j < f_copy[0].size(); j++) {
      // Add in random noise using rng and a the defined lambda function
      f_copy[i][j].unaryop(noise);  // not sure how this workes
      // i see we pass lambda function noise
    }

    // Apply mask to get boundary condition right
    f_copy[i] = mask * f_copy[i];  // apply mask
  }

  // Done
  return f_copy;
}

// Creates the ground state hamiltonian from given functions f
Tensor<double> TDDFT::CreateGroundHamiltonian(World& world,
                                              std::vector<real_function_3d> f,
                                              size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);
  // Get sizes
  size_t m = f.size();
  // Debugging
  if (print_level > 2) {
    Tensor<double> S = matrix_inner(world, f, f);
    if (world.rank() == 0) print("   Ground state overlap:");
    if (world.rank() == 0) print(S);
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
  Tensor<double> T =
      1.0 / 2.0 *
      (matrix_inner(world, fx, fx) + matrix_inner(world, fy, fy) +
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
  real_function_3d v_coul = 2.0 * Coulomb(world);

  // Clear old stored potentials
  stored_v_coul.clear();
  stored_v_nuc.clear();

  // If storing potentials, save them here
  if (Rparams.store_potential) {
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
  real_convolution_3d op =
      CoulombOperator(world, Rparams.small, FunctionDefaults<3>::get_thresh());
  std::vector<real_function_3d> Kf =
      zero_functions_compressed<double, 3>(world, m);
  for (size_t i = 0; i < m; ++i) {
    std::vector<real_function_3d> psif =
        mul_sparse(world, f[i], f, FunctionDefaults<3>::get_thresh());
    truncate(world, psif);
    psif = apply(world, op, psif);
    truncate(world, psif);

    // Save the potential here if we are saving it
    if (Rparams.store_potential) {
      stored_potential.push_back(psif);
    }

    psif = mul_sparse(world, f[i], psif, FunctionDefaults<3>::get_thresh());
    gaxpy(world, 1.0, Kf, 1.0, psif);
  }

  // Only use the exchange above if HF:
  Tensor<double> V;
  real_function_3d v_xc;

  if (Gparams.xc == "hf") {
    // Construct V
    V = matrix_inner(world, f, vf) - matrix_inner(world, f, Kf);
  } else {  // DFT

    XCOperator xcop = create_xcoperator(world, f, Gparams.xc);

    real_function_3d v_xc = xcop.make_xc_potential();
    v = v + v_xc;
    std::vector<real_function_3d> vf = v * f;
    if ((*xcop.xc).hf_exchange_coefficient() > 0.0) {
      // XCOperator has member variable xc, which is an xcfunctional
      // which has the hf_exchange_coeff we need here
      gaxpy(world, 1.0, vf, -(*xcop.xc).hf_exchange_coefficient(), Kf);
    }
    V = matrix_inner(world, f, vf);
  }

  // Now create the hamiltonian
  hamiltonian = T + V;

  for (int64_t i = 0; i < hamiltonian.dim(0); i++) {
    for (int64_t j = i + 1; j < hamiltonian.dim(1); j++) {
      //      print(i, j);
      //      print(xAx(i, j));
      //     print(xAx(j, i));
      hamiltonian(j, i) = hamiltonian(i, j);
    }
  }
  double traceOfHamiltonian(0);
  for (int64_t i = 0; i < hamiltonian.dim(0); i++) {
    traceOfHamiltonian += hamiltonian(i, i);
  }
  print("Trace of Hamiltonian");
  print(traceOfHamiltonian);
  // Save a matrix that is
  // (T+V) - Lambda * eye
  // Copy hamiltonian and zero the diagonal
  ham_no_diag = copy(hamiltonian);
  for (size_t i = 0; i < m; i++) ham_no_diag(i, i) = 0.0;

  // Debug output
  if (print_level >= 2 and world.rank() == 0) {
    print("   Ground state hamiltonian:");
    print(hamiltonian);
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "   Create grnd ham:");

  return hamiltonian;
}

// Creates the transition densities
std::vector<real_function_3d> TDDFT::transition_density(
    World& world,
    std::vector<real_function_3d> const& orbitals,
    response_space& x,
    response_space& y) {
  // Get sizes
  size_t m = x.size();

  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);
  // (REMARK)  implementation with imaginary wavefunctions
  // we need to create functions for conjugates
  // Run over virtual...
  /*
  x.reconstruct_rf();
  y.reconstruct_rf();
  reconstruct(world, Gparams.orbitals);
  */
  // print("Thresh before truncate in t-density");
  // print(FunctionDefaults<3>::get_thresh());
  x.truncate_rf();
  y.truncate_rf();
  truncate(world, Gparams.orbitals);
  for (size_t b = 0; b < m; b++) {
    // Run over occupied...
    // y functions are zero if TDA is active
    densities[b] =
        dot(world, x[b], Gparams.orbitals) + dot(world, y[b], Gparams.orbitals);
  }

  truncate(world, densities);
  world.gop.fence();
  // Done!
  return densities;
}
// Creates the transition density
std::vector<real_function_3d> TDDFT::GetTransitionDensities(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();
  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

  // Run over virtual...
  for (size_t i = 0; i < m; i++) {
    // Run over occupied...
    for (size_t j = 0; j < n; j++) {
      // y functions are zero if TDA is active
      densities[i] = densities[i] + orbitals[j] * f[i][j];
      // densities[i] =densities[i] + orbitals[j] * dagger(f[i][j]);TODO:
      // DAGGER
    }
  }

  // Done!
  return densities;
}

std::vector<real_function_3d> TDDFT::GetConjugateTransitionDensities(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();
  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);
  // Run over virtual...
  for (size_t i = 0; i < m; i++) {
    // Run over occupied...
    for (size_t j = 0; j < n; j++) {
      // y functions are zero if TDA is active
      densities[i] = densities[i] + orbitals[j] * f[i][j];
      // densities[i] + orbitals[j] * dagger(f[i][j]) TODO: DAGGER;
    }
  }
  // Done!
  return densities;
}

template <std::size_t NDIM>
void TDDFT::set_protocol(World& world, double thresh) {
  size_t k;
  // Allow for imprecise conversion of threshold
  if (thresh >= 0.9e-2)
    k = 4;
  else if (thresh >= 0.9e-4)
    k = 6;
  else if (thresh >= 0.9e-6)
    k = 8;
  else if (thresh >= 0.9e-8)
    k = 10;
  else
    k = 12;

  // k defaults to make sense with thresh, override by providing k in
  // input file
  if (Rparams.k > 0) {
    FunctionDefaults<NDIM>::set_k(Rparams.k);
  } else {
    FunctionDefaults<NDIM>::set_k(k);
  }

  // MolDFT sets all these, so copying
  FunctionDefaults<NDIM>::set_thresh(thresh);
  FunctionDefaults<NDIM>::set_refine(true);
  FunctionDefaults<NDIM>::set_initial_level(2);
  FunctionDefaults<NDIM>::set_autorefine(false);
  FunctionDefaults<NDIM>::set_apply_randomize(false);
  FunctionDefaults<NDIM>::set_project_randomize(false);
  // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
  // is this? Do i need it?)

  // dconv defaults to thresh*100, overrirde by providing dconv in input
  // file
  if (Rparams.dconv_set == false) {
    Rparams.dconv = thresh * 100;
  }

  // Basic print
  if (world.rank() == 0) {
    print("\nSolving NDIM=",
          NDIM,
          " with thresh",
          thresh,
          "    k",
          FunctionDefaults<NDIM>::get_k(),
          "  dconv",
          std::max(thresh, Rparams.dconv),
          "\n");
  }
}

void TDDFT::check_k(World& world, double thresh, size_t k) {
  // Boolean to redo ground hamiltonian calculation if
  // ground state orbitals change
  bool redo = false;

  // Verify ground state orbitals have correct k
  if (FunctionDefaults<3>::get_k() != Gparams.orbitals[0].k()) {
    // Re-read orbitals from the archive (assuming
    // the archive has orbitals stored at a higher
    // k value than what was previously computed
    // with)
    Gparams.read(world, Rparams.archive);
    reconstruct(world, Gparams.orbitals);

    // Reset correct k (its set in Gparams.read)
    FunctionDefaults<3>::set_k(k);

    // Project each ground state to correct k
    for (unsigned int i = 0; i < Gparams.orbitals.size(); i++)
      Gparams.orbitals[i] = project(
          Gparams.orbitals[i], FunctionDefaults<3>::get_k(), thresh, false);
    world.gop.fence();

    // Clean up a bit
    truncate(world, Gparams.orbitals);

    // Ground state orbitals changed, clear old hamiltonian
    redo = true;
  }

  // Recalculate ground state hamiltonian here
  if (redo or !hamiltonian.has_data()) {
    hamiltonian =
        CreateGroundHamiltonian(world, Gparams.orbitals, Rparams.print_level);
  }

  // If we stored the potential, check that too
  if (Rparams.store_potential) {
    if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k()) {
      // Project the potential into correct k
      for (unsigned int i = 0; i < stored_potential.size(); i++) {
        reconstruct(world, stored_potential[i]);
        for (unsigned int j = 0; j < stored_potential[0].size(); j++)
          stored_potential[i][j] = project(stored_potential[i][j],
                                           FunctionDefaults<3>::get_k(),
                                           thresh,
                                           false);
        world.gop.fence();
      }
    }
    if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
      stored_v_coul =
          project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
    if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
      stored_v_nuc =
          project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
  }

  // Verify response functions have correct k
  if (x_response.size() != 0) {
    if (FunctionDefaults<3>::get_k() != x_response[0][0].k()) {
      // Project all x components into correct k
      for (unsigned int i = 0; i < x_response.size(); i++) {
        reconstruct(world, x_response[i]);
        for (unsigned int j = 0; j < x_response[0].size(); j++)
          x_response[i][j] = project(
              x_response[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      x_response.truncate_rf();

      // Do same for y components if applicable
      // (Always do this, as y will be zero
      //  and still used in doing DFT and TDA)
      // Project all y components into correct k
      for (unsigned int i = 0; i < y_response.size(); i++) {
        reconstruct(world, y_response[i]);
        for (unsigned int j = 0; j < y_response[0].size(); j++)
          y_response[i][j] = project(
              y_response[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      y_response.truncate_rf();
    }
  }

  // Don't forget the mask function as well
  if (FunctionDefaults<3>::get_k() != mask.k()) {
    mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
  }
  // Don't forget right hand side

  // Verify response functions have correct k
  if (P.size() != 0) {
    if (FunctionDefaults<3>::get_k() != P[0][0].k()) {
      // Project all x components into correct k
      for (unsigned int i = 0; i < x_response.size(); i++) {
        reconstruct(world, P[i]);
        for (unsigned int j = 0; j < P[0].size(); j++)
          P[i][j] =
              project(P[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      x_response.truncate_rf();

      // Do same for y components if applicable
      // (Always do this, as y will be zero
      //  and still used in doing DFT and TDA)
      // Project all y components into correct k
      for (unsigned int i = 0; i < Q.size(); i++) {
        reconstruct(world, Q[i]);
        for (unsigned int j = 0; j < Q[0].size(); j++)
          Q[i][j] =
              project(Q[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      y_response.truncate_rf();
    }
  }
  // Make sure everything is done before leaving
  world.gop.fence();
}

// Creates random guess functions semi-intelligently(?)
response_space TDDFT::create_random_guess(
    World& world,
    size_t m,                                // m response states
    size_t n,                                // n ground states
    std::vector<real_function_3d>& grounds,  // guess should have size n
    Molecule& molecule) {
  // Basic output
  if (world.rank() == 0)
    print("   Using a random guess for initial response functions.\n");

  // Create empty container and add in randomness
  response_space f(world, m, n);
  f = add_randomness(world, f, 1e3);  // noise all over the world

  // Create and apply a centered gaussian on each atom so that the
  // randomness is localized around the atoms
  real_function_3d gaus = real_factory_3d(world);
  for (auto atom : molecule.get_atoms()) {
    real_function_3d x =
        real_factory_3d(world).functor(real_functor_3d(new GaussianGuess<3>(
            atom.get_coords(), 0.01, std::vector<int>{0, 0, 0})));
    gaus = gaus + x;
  }
  f = f * gaus;

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, grounds);
  for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

  // Normalize
  normalize(world, f);

  return f;
}

// Creates random guess functions semi-intelligently(?)
std::vector<real_function_3d> TDDFT::create_random_guess(
    World& world,
    size_t m,
    std::vector<real_function_3d>& grounds,
    Molecule& molecule) {
  // Basic output
  if (world.rank() == 0)
    print("   Using a random guess for initial response functions.");

  // Create empty container and add in randomness
  std::vector<real_function_3d> f =
      zero_functions_compressed<double, 3>(world, m);
  // create vector of m functions:w

  // Create and apply a centered gaussian on each atom so that the
  // randomness is localized around the atoms
  real_function_3d gaus = real_factory_3d(world);
  for (auto atom : molecule.get_atoms()) {
    real_function_3d x =
        real_factory_3d(world).functor(real_functor_3d(new GaussianGuess<3>(
            atom.get_coords(), 0.01, std::vector<int>{0, 0, 0})));
    gaus = gaus + x;
  }

  // Lambda function to add in noise
  auto lambda = [](const Key<3>& key, Tensor<double>& x) mutable {
    Tensor<double> y(x.size());
    y.fillrandom();
    y.scale(1e3);
    x = x + y;
  };

  // Go through each function in f_copy and add in random noise
  for (unsigned int i = 0; i < f.size(); i++) {
    // Add in random noise using rng and a the defined lambda function
    f[i].unaryop(lambda);

    // Apply mask to get boundary condition right
    f[i] = mask * f[i] * gaus;
  }

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, grounds);
  for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

  // Normalize
  for (unsigned int i = 0; i < f.size(); i++) {
    double norm = f[i].norm2();
    f[i].scale(1.0 / norm);
  }

  return f;
}

// Creates an initial guess function from nwchem output files
response_space TDDFT::create_nwchem_guess(World& world, size_t m) {
  // Basic output
  if (world.rank() == 0)
    print("   Creating an initial guess from NWChem file", Rparams.nwchem);

  // Create empty containers
  response_space f;

  // Create the nwchem reader
  slymer::NWChem_Interface nwchem(Rparams.nwchem, std::cout);

  // For parallel runs, silencing all but 1 slymer instance
  if (world.rank() != 0) {
    std::ostream dev_null(nullptr);
    nwchem.err = dev_null;
  }

  // Read in basis set
  nwchem.read(slymer::Properties::Basis);

  // Read in the molecular orbital coefficients, energies,
  // and occupancies
  nwchem.read(slymer::Properties::MOs | slymer::Properties::Energies |
              slymer::Properties::Occupancies);

  // Create the nwchem orbitals as madness functions
  std::vector<real_function_3d> temp1;
  for (auto basis :
       slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set)) {
    // Get the center of gaussian as its special point
    std::vector<coord_3d> centers;
    coord_3d r;
    r[0] = basis.get().center[0];
    r[1] = basis.get().center[1];
    r[2] = basis.get().center[2];
    centers.push_back(r);

    // Now make the function
    temp1.push_back(FunctionFactory<double, 3>(world).functor(
        std::shared_ptr<FunctionFunctorInterface<double, 3>>(
            new slymer::Gaussian_Functor(basis.get(), centers))));

    // Let user know something is going on
    if (temp1.size() % 10 == 0 and world.rank() == 0)
      print("Created", temp1.size(), "functions.");
  }
  if (world.rank() == 0) print("Finished creating", temp1.size(), "functions.");

  // Normalize ao's
  madness::normalize(world, temp1);

  // Transform ao's now
  std::vector<real_function_3d> temp = madness::transform(
      world, temp1, nwchem.MOs, FunctionDefaults<3>::get_thresh(), true);

  // Now save the unoccupied orbitals
  std::vector<real_function_3d> temp2;
  size_t num_virt = 0;
  for (size_t i = 0; i < temp1.size(); i++) {
    if (nwchem.occupancies[i] == 0) {
      temp2.push_back(temp[i]);
      num_virt++;
    }
  }

  // Create as many vectors of functions as we can from these nwchem
  // virtual orbitals putting 1 virtual orbital from nwchem per vector
  for (size_t i = 0; i < std::max(m, num_virt); i++) {
    // Create the vector to add the new function to
    std::vector<real_function_3d> v1 =
        zero_functions_compressed<double, 3>(world, Gparams.orbitals.size());

    // Put the "new" function into the vector
    v1[i % v1.size()] = temp2[i];

    // Add vector to return container
    f.push_back(v1);

    // See if we've made enough functions
    if (f.size() >= m) break;
  }
  if (world.rank() == 0)
    print("Created", f.size(), "guess functions from provided NWChem data.");

  // If not enough functions have been made, start adding symmetry adapted
  // functions
  size_t n = f.size();

  // If still not enough functions have been made, add in random guesses
  if (n < m) {
    // Tell user the bad news
    if (world.rank() == 0)
      print("\n   Only",
            n,
            "guess functions were provided by augmenting NWChem "
            "functions.\n   "
            "Augmenting with random functions.");

    // Create the random guess
    response_space rand = create_random_guess(
        world, m - n, Gparams.num_orbitals, Gparams.orbitals, Gparams.molecule);

    // Add to vector of functions
    for (unsigned int i = 0; i < rand.size(); i++) f.push_back(rand[i]);
  }

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, Gparams.orbitals);
  for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

  // Truncate and normalize
  f.truncate_rf();
  normalize(world, f);

  return f;
}

// Creates potentials using the ResponsePotential object
// Potentials are modified in place
void TDDFT::create_all_potentials(World& world,
                                  response_space& x,
                                  response_space& x_gamma,
                                  response_space& x_V0,
                                  ResponsePotential& potentials,
                                  size_t print_level) {
  // Intermediaries
  response_space gammaJ, gammaK, groundJ, groundK;

  // Calc. coulomb like terms
  if (print_level >= 1) molresponse::start_timer(world);
  potentials.coulomb_terms(x, gammaK, groundJ);
  if (print_level >= 1) molresponse::end_timer(world, "Coulomb terms:");

  // Calc. exchange like terms
  if (print_level >= 1) molresponse::start_timer(world);
  potentials.exchange_terms(x, gammaJ, groundK);
  if (print_level >= 1) molresponse::end_timer(world, "Exchange terms:");

  // Assemble pieces together
  x_gamma = gammaJ - gammaK;
  x_V0 = groundJ - groundK;

  // Debugging output
  if (print_level >= 2) {
    // Coulomb
    if (world.rank() == 0) printf("   Coulomb Deriv matrix:\n");
    Tensor<double> temp = expectation(world, x, gammaJ);
    if (world.rank() == 0) print(temp);

    // Exchange or VXC
    if (Rparams.xc == "hf" and world.rank() == 0)
      printf("   Exchange Deriv matrix:\n");
    if (Rparams.xc != "hf" and world.rank() == 0)
      printf("   Negative of XC Deriv matrix:\n");
    temp = expectation(world, x, gammaK);
    if (world.rank() == 0) print(temp);

    // Total Gamma
    if (world.rank() == 0) printf("   Gamma matrix:\n");
    temp = expectation(world, x, x_gamma);
    if (world.rank() == 0) print(temp);

    // Coulomb (ground)
    if (world.rank() == 0) printf("   Coulomb + Nuclear potential matrix:\n");
    temp = expectation(world, x, groundJ);
    if (world.rank() == 0) print(temp);

    // Exchange or VXC (ground)
    if (Rparams.xc == "hf" and world.rank() == 0)
      printf("   Exchange potential matrix:\n");
    if (Rparams.xc != "hf" and world.rank() == 0)
      printf("   XC potential matrix:\n");
    temp = expectation(world, x, groundK);
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0) printf("   Total Potential Energy matrix:\n");
    temp = expectation(world, x, x_V0);
    if (world.rank() == 0) print(temp);
  }
}

// Iterates the response functions until converged or out of iterations
void TDDFT::IteratePolarizability(World& world, response_space& dipoles) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
  QProjector<double, 3> projector(
      world, Gparams.orbitals);     // Projector to project out ground state
  size_t n = Gparams.num_orbitals;  // Number of ground state orbitals
  size_t m = Rparams.states;        // Number of excited states
  Tensor<double> x_norms(m);
  // Holds the norms of x function residuals (for convergence)
  Tensor<double> y_norms(
      m);  // Holds the norms of y function residuals (for convergence)
  Tensor<double> x_shifts(m);    // Holds the shifted energy values
  Tensor<double> y_shifts(m);    // Holds the shifted energy values
  response_space bsh_x_resp;     // Holds wave function corrections
  response_space bsh_y_resp;     // Holds wave function corrections
  response_space x_differences;  // Holds wave function corrections
  response_space y_differences;  // Holds wave function corrections
  response_space x_gamma;        // Holds the perturbed two electron piece
  response_space y_gamma;        // Holds the perturbed two electron piece
  response_space Hx;             // Holds the perturbed two electron piece
  response_space Hy;             // Holds the perturbed two electron piece
  response_space Gx;             // Holds the perturbed two electron piece
  response_space Gy;             // Holds the perturbed two electron piece
  response_space x_fe;  // Holds the ground state-fock and energy scaled x
  // response components
  response_space y_fe;  // Holds the ground state-fock and energy scaled y
  // response components
  response_space V_x_response;  // Holds V^0 applied to response functions
  response_space V_y_response;  // Holds V^0 applied to response functions
  response_space B_x;  // Holds the off diagonal perturbed piece of y equation
  response_space B_y;  // Holds the off diagonal perturbed piece of x equation
  response_space shifted_V_x_response;  // Holds the shifted V^0 applied to
  // response functions
  response_space shifted_V_y_response;  // Holds the shifted V^0 applied to
  // response functions
  response_space old_x_response;  // Holds the old x_response vector of vectors
  response_space old_y_response;  // Holds the old y_response vector of vectors
  real_function_3d v_xc;          // For TDDFT
  bool converged = false;         // Converged flag

  // If DFT, initialize the XCOperator
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  // The KAIN solver
  XNonlinearSolver<response_space, double, TDHF_allocator> kain(
      TDHF_allocator(world, (Rparams.omega != 0.0) ? 2 * m : m, n), false);

  // Setting max sub size for KAIN solver
  if (Rparams.kain) kain.set_maxsub(Rparams.maxsub);

  // Set omega (its constant here,
  // and has only 1 entry for each axis)
  omega = Tensor<double>(3);
  omega = Rparams.omega;

  // Verify if any shift is needed (NEEDS CHECKING)
  if ((Gparams.energies[n - 1] + Rparams.omega) > 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    // for all \eps, \omega
    x_shifts =
        create_shift(world, Gparams.energies, omega, Rparams.print_level, "x");
    y_shifts = Gparams.energies[n - 1] + Rparams.omega + 0.05;
  }

  // Construct BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
      bsh_x_operators = create_bsh_operators(world,
                                             x_shifts,
                                             Gparams.energies,
                                             omega,
                                             Rparams.small,
                                             FunctionDefaults<3>::get_thresh());
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
      bsh_y_operators;

  // Negate omega to make this next set of BSH operators \eps - omega
  if (Rparams.omega != 0.0) {
    omega = -omega;
    bsh_y_operators = create_bsh_operators(world,
                                           y_shifts,
                                           Gparams.energies,
                                           omega,
                                           Rparams.small,
                                           FunctionDefaults<3>::get_thresh());
  }

  // Now to iterate
  while (iteration < Rparams.max_iter and !converged) {
    // Start a timer for this iteration
    molresponse::start_timer(world);

    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n",
               static_cast<int>(iteration),
               wall_time());
      if (world.rank() == 0) print(" -------------------------------");
    }

    // If omega = 0.0, x = y
    if (Rparams.omega == 0.0) y_response = x_response.copy();

    // Save current to old
    old_x_response = x_response.copy();
    if (Rparams.omega != 0.0) old_y_response = y_response.copy();
    //      world.gop.fence(); // Norm calc. below sometimes hangs without
    //      this
    //      (?)

    // Get norms
    for (size_t i = 0; i < m; i++)
      x_norms[i] = sqrt(inner(x_response[i], x_response[i]) -
                        inner(y_response[i], y_response[i]));

    // Scale x and y
    Tensor<double> rec_norms(m);
    for (size_t i = 0; i < m; i++)
      rec_norms(i) = 1.0 / std::max(1.0, x_norms(i));
    x_response.scale(rec_norms);
    y_response.scale(rec_norms);

    // Here I will first compute Hx Hy Gx Gy
    // Create gamma
    // If TDA we only compute Hx and Gy
    Hx = ComputeHf(world,
                   x_response,
                   Gparams.orbitals,
                   Rparams.small,
                   FunctionDefaults<3>::get_thresh(),
                   Rparams.print_level,
                   "x");
    Gy = ComputeGf(world,
                   y_response,
                   Gparams.orbitals,
                   Rparams.small,
                   FunctionDefaults<3>::get_thresh(),
                   Rparams.print_level,
                   "y");
    // else Compute everything
    if (Rparams.omega != 0.0) {  // not sure why this is the condition
      Hy = ComputeHf(world,
                     y_response,
                     Gparams.orbitals,
                     Rparams.small,
                     FunctionDefaults<3>::get_thresh(),
                     Rparams.print_level,
                     "y");
      Gx = ComputeGf(world,
                     x_response,
                     Gparams.orbitals,
                     Rparams.small,
                     FunctionDefaults<3>::get_thresh(),
                     Rparams.print_level,
                     "x");
    }

    x_gamma = CreateGamma(world,
                          x_response,
                          y_response,
                          Gparams.orbitals,
                          Rparams.small,
                          FunctionDefaults<3>::get_thresh(),
                          Rparams.print_level,
                          "x");

    if (Rparams.omega != 0.0)  // what and why?
      y_gamma = CreateGamma(world,
                            y_response,
                            x_response,
                            Gparams.orbitals,
                            Rparams.small,
                            FunctionDefaults<3>::get_thresh(),
                            Rparams.print_level,
                            "y");

    // Create \hat{V}^0 applied to response functions
    V_x_response =
        CreatePotential(world, x_response, xc, Rparams.print_level, "x");
    if (Rparams.omega != 0.0)
      V_y_response =
          CreatePotential(world, y_response, xc, Rparams.print_level, "y");

    // Apply shift
    V_x_response = apply_shift(world, x_shifts, V_x_response, x_response);
    if (Rparams.omega != 0.0)
      V_y_response = apply_shift(world, y_shifts, V_y_response, y_response);

    // Create \epsilon applied to response functions
    x_fe =
        x_response * ham_no_diag;  // scale_2d(world, x_response, ham_no_diag);
    if (Rparams.omega != 0.0)
      y_fe = y_response *
             ham_no_diag;  // scale_2d(world, y_response, ham_no_diag);
    if (Rparams.print_level >= 2) {
      Tensor<double> t = expectation(world, x_response, x_fe);
      if (world.rank() == 0) {
        print("   Energy scaled response orbitals for x components:");
        print(t);
      }

      if (Rparams.omega != 0.0) {
        t = expectation(world, y_response, y_fe);
        if (world.rank() == 0) {
          print("   Energy scaled response orbitals for y components:");
          print(t);
        }
      }
    }

    // Load balance
    // Only balancing on x-components. Smart?
    // Only balance on first two iterations or every 5th iteration
    if (world.size() > 1 && (iteration < 2 or iteration % 5 == 0)) {
      // Start a timer
      if (Rparams.print_level >= 1) molresponse::start_timer(world);
      if (world.rank() == 0) print("");  // Makes it more legible

      LoadBalanceDeux<3> lb(world);
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < Rparams.states; k++) {
          lb.add_tree(x_response[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(V_x_response[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(x_gamma[k][j], lbcost<double, 3>(1.0, 8.0), true);
        }
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, "Load balancing:");
    }

    // Calculate coupling terms
    response_space B_y =
        createBf(world, Gx, Gparams.orbitals, Rparams.print_level);
    if (Rparams.omega != 0.0)
      response_space B_x =
          createBf(world, Gx, Gparams.orbitals, Rparams.print_level);

    // Scale dipoles by same value
    response_space dip_copy(dipoles);
    dip_copy.scale(rec_norms);

    // Construct RHS of equation
    response_space rhs_x, rhs_y;
    rhs_x = V_x_response - x_fe + dip_copy + x_gamma + B_y;
    if (Rparams.omega != 0.0)
      rhs_y = V_y_response - y_fe + dip_copy + y_gamma + B_x;

    // Project out ground state
    for (size_t i = 0; i < m; i++) rhs_x[i] = projector(rhs_x[i]);
    if (Rparams.omega != 0.0) {
      for (size_t i = 0; i < m; i++) rhs_y[i] = projector(rhs_y[i]);
    }

    // Debugging output
    if (Rparams.print_level >= 2) {
      if (world.rank() == 0)
        print("   Norms of RHS of main equation x components:");
      print_norms(world, rhs_x);

      if (Rparams.omega != 0.0) {
        if (world.rank() == 0)
          print("   Norms of RHS of main equation y components:");
        print_norms(world, rhs_y);
      }
    }

    // Apply BSH and get updated response components
    if (Rparams.print_level >= 1) molresponse::start_timer(world);
    bsh_x_resp = apply(world, bsh_x_operators, rhs_x);
    if (Rparams.omega != 0.0) bsh_y_resp = apply(world, bsh_y_operators, rhs_y);
    if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply BSH:");

    // Scale by -2.0 (coefficient in eq. 37 of reference paper)
    for (size_t i = 0; i < m; i++)
      bsh_x_resp[i] = bsh_x_resp[i] * (std::max(1.0, x_norms[i]) * -2.0);
    if (Rparams.omega != 0.0) {
      for (size_t i = 0; i < m; i++)
        bsh_y_resp[i] = bsh_y_resp[i] * (std::max(1.0, x_norms[i]) * -2.0);
    }

    // Debugging output
    if (Rparams.print_level >= 2) {
      if (world.rank() == 0)
        print("   Norms after application of BSH to x components:");
      print_norms(world, bsh_x_resp);

      if (Rparams.omega != 0.0) {
        if (world.rank() == 0)
          print("   Norms after application of BSH to y components:");
        print_norms(world, bsh_y_resp);
      }
    }

    // Update orbitals
    x_response = bsh_x_resp;
    if (Rparams.omega != 0.0) y_response = bsh_y_resp;

    // Get the difference between old and new
    x_differences = old_x_response - x_response;
    if (Rparams.omega != 0.0) y_differences = old_y_response - y_response;

    // Next calculate 2-norm of these vectors of differences
    // Remember: the entire vector is one state
    for (size_t i = 0; i < m; i++) x_norms(i) = norm2(world, x_differences[i]);
    if (Rparams.omega != 0.0) {
      for (size_t i = 0; i < m; i++)
        y_norms(i) = norm2(world, y_differences[i]);
    }

    // Basic output
    if (Rparams.print_level >= 1 and world.rank() == 0) {
      print("\n   2-norm of response function residuals of x components:");
      print(x_norms);

      if (Rparams.omega != 0.0) {
        print("   2-norm of response function residuals of y components:");
        print(y_norms);
      }
    }

    // Check convergence
    if (std::max(x_norms.absmax(), y_norms.absmax()) < Rparams.dconv and
        iteration > 0) {
      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, "This iteration:");
      if (world.rank() == 0) print("\n   Converged!");
      converged = true;
      break;
    }

    // KAIN solver update
    // Returns next set of components
    // If not kain, save the new components
    if (Rparams.kain) {
      if (Rparams.omega != 0.0) {
        // Add y functions to bottom of x functions
        // (for KAIN)
        for (size_t i = 0; i < m; i++) {
          x_response.push_back(y_response[i]);
          x_differences.push_back(y_differences[i]);
        }
      }

      molresponse::start_timer(world);
      x_response = kain.update(
          x_response, x_differences, FunctionDefaults<3>::get_thresh(), 3.0);
      molresponse::end_timer(world, " KAIN update:");

      if (Rparams.omega != 0.0) {
        // Add new functions back into y and
        // reduce x size back to original
        for (size_t i = 0; i < m; i++) y_response[i] = x_response[m + i];
        for (size_t i = 0; i < m; i++) {
          x_response.pop_back();
          x_differences.pop_back();
        }
      }
    }

    // Apply mask
    for (size_t i = 0; i < m; i++) x_response[i] = mask * x_response[i];
    if (Rparams.omega != 0.0) {
      for (size_t i = 0; i < m; i++) y_response[i] = mask * y_response[i];
    }

    // Update counter
    iteration += 1;

    // Done with the iteration.. truncate
    x_response.truncate_rf();
    if (Rparams.omega != 0.0) y_response.truncate_rf();

    // Save
    if (Rparams.save) {
      molresponse::start_timer(world);
      save(world, Rparams.save_file);
      if (Rparams.print_level >= 1) molresponse::end_timer(world, "Save:");
    }
    // Basic output
    if (Rparams.print_level >= 1)
      molresponse::end_timer(world, " This iteration:");
  }
}  // Done with iterate_polarizability

/*
print("x norms in iteration after truncation Plot: ", iteration);
print(x_response.norm2());

print("y norms in iteration after truncation Plot: ", iteration);
print(y_response.norm2());
*/

// Calculates polarizability according to
// alpha_ij(\omega) = -sum_{ directions } < x_j | r_i | 0 > + < 0 | r_i |
// y_j
// >
void TDDFT::polarizability(World& world, Tensor<double> polar) {
  // Get transition density
  // std::vector<real_function_3d> rhos = transition_density(world,
  // Gparams.orbitals, x_response, y_response);
  std::vector<real_function_3d> rhos;
  if (Rparams.omega == 0)
    rhos = transition_density(world, Gparams.orbitals, x_response, x_response);
  else
    rhos = transition_density(world, Gparams.orbitals, x_response, y_response);

  // For each r_axis
  for (size_t axis = 0; axis < 3; axis++) {
    real_function_3d drho = rhos[axis];

    // Run over axis and calc.
    // the polarizability
    for (size_t i = 0; i < 3; i++) {
      // Create dipole operator in the 'i' direction
      std::vector<int> f(3, 0);
      f[i] = 1;
      real_function_3d dip = real_factory_3d(world).functor(
          real_functor_3d(new BS_MomentFunctor(f)));
      polar(axis, i) = -2.0 * dip.inner(drho);
    }
  }
}

void TDDFT::PrintPolarizabilityAnalysis(World& world,
                                        const Tensor<double> polar_tensor,
                                        const Tensor<double> omega) {
  // Final polarizability analysis
  // diagonalize
  Tensor<double> V, epolar;
  syev(polar_tensor, V, epolar);
  double Dpolar_average = 0.0;
  double Dpolar_iso = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    Dpolar_average = Dpolar_average + epolar[i];
  Dpolar_average = Dpolar_average / 3.0;
  Dpolar_iso =
      sqrt(.5) * sqrt(std::pow(polar_tensor(0, 0) - polar_tensor(1, 1), 2) +
                      std::pow(polar_tensor(1, 1) - polar_tensor(2, 2), 2) +
                      std::pow(polar_tensor(2, 2) - polar_tensor(0, 0), 2));

  if (world.rank() == 0) {
    print("\nTotal Dynamic Polarizability Tensor");
    printf("\nFrequency  = %.6f a.u.\n\n", omega(0, 0));
    // printf("\nWavelength = %.6f a.u.\n\n", Rparams.omega * ???);
    print(polar_tensor);
    printf("\tEigenvalues = ");
    printf("\t %.6f \t %.6f \t %.6f \n", epolar[0], epolar[1], epolar[2]);
    printf("\tIsotropic   = \t %.6f \n", Dpolar_average);
    printf("\tAnisotropic = \t %.6f \n", Dpolar_iso);
    printf("\n");
  }
}
// TODO It makes sense to align the plots with the direction of the perturbation
// ??
// TODO In excited state calculations we are going to get various directions.
// We should plot the functions with respect to

/* @brief plot tranistions and ground orbitals in x y z direction
 *
 * @param world
 * @param iteration
 * @param x_response
 * @param y_response
 * @param Rparams
 * @param Gparams
 */
void TDDFT::PlotGroundandResponseOrbitals(World& world,
                                          size_t iteration,
                                          response_space& x_response,
                                          response_space& y_response,
                                          ResponseParameters const& Rparams,
                                          GroundParameters const& Gparams) {
  std::filesystem::create_directories("plots/xy");
  std::filesystem::create_directory("plots/ground");
  std::filesystem::create_directory("plots/transition_density");

  // TESTING
  // get transition density
  // num orbitals
  size_t n = x_response[0].size();
  size_t m = x_response.size();

  real_function_3d ground_density =
      dot(world, Gparams.orbitals, Gparams.orbitals);
  std::vector<real_function_3d> densities =
      transition_density(world, Gparams.orbitals, x_response, y_response);
  std::string dir("xyz");
  // for plotname size
  size_t buffSize = 500;
  char plotname[buffSize];
  double Lp = std::min(Gparams.L, 24.0);
  // Doing line plots along each axis
  for (int d = 0; d < 3; d++) {
    // print ground_state
    plotCoords plt(0, Lp);
    if (iteration == 1) {
      snprintf(
          plotname, buffSize, "plots/ground/ground_density_%c.plt", dir[d]);
      plot_line(plotname, 5001, plt.lo, plt.hi, ground_density);
    }

    for (int b = 0; b < static_cast<int>(m); b++) {
      for (int i = 0; i < static_cast<int>(n); i++) {
        // print ground_state
        snprintf(plotname,
                 buffSize,
                 "plots/ground/ground_%c_%d.plt",
                 dir[d],
                 static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, Gparams.orbitals[i]);
      }
      for (int i = 0; i < static_cast<int>(n); i++) {
        // print ground_state
        // plot x function  x_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/xy/x_K_%d_iter_%d_dir_%d_res_%d_orb_%d",
                 FunctionDefaults<3>::get_k(),
                 static_cast<int>(iteration - 1),
                 dir[d],
                 static_cast<int>(b),
                 static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, x_response[b][i]);

        // plot y functione  y_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/xy/y_K_%d_iter_%d_dir_%d_res_%d_orb_%d",
                 FunctionDefaults<3>::get_k(),
                 static_cast<int>(iteration - 1),
                 dir[d],
                 static_cast<int>(b),
                 static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, y_response[b][i]);
      }
    }
  }
  world.gop.fence();

  // END TESTING
}
// Main function, makes sure everything happens in correct order
// Solves for polarizability
void TDDFT::solve_polarizability(World& world, Property& p) {
  // Get start time
  molresponse::start_timer(world);

  // Warm and fuzzy
  if (world.rank() == 0) {
    print("\n\n    Response Calculation");
    print("   ------------------------");
  }

  // Create the polarizability tensor
  Tensor<double> polar_tensor(3, 3);

  // Keep a copy of dipoles * MO (needed explicitly in eq.)
  response_space dipoles;

  // For each protocol
  for (unsigned int proto = 0; proto < Rparams.protocol_data.size(); proto++) {
    // Set defaults inside here
    set_protocol<3>(world, Rparams.protocol_data[proto]);

    // Do something to ensure all functions have same k value
    check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());

    // Create guesses if no response functions
    // If restarting, load here
    if (proto == 0) {
      if (Rparams.restart) {
        if (world.rank() == 0)
          print("   Initial guess from file:", Rparams.restart_file);
        load(world, Rparams.restart_file);
        check_k(
            world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());
      } else {  // Dipole guesses

        x_response = PropertyRHS(world, p);
        y_response = x_response.copy();
      }
    }

    // Set the dipoles (ground orbitals are probably
    // more accurate now, so recalc the dipoles)
    dipoles = PropertyRHS(world, p);
    // why is it called dipole guess.
    // This is just orbitals times dipole operator

    // Now actually ready to iterate...
    IteratePolarizability(world, dipoles);
  }

  // Have response function, now calculate polarizability for this axis
  polarizability(world, polar_tensor);

  // Final polarizability analysis
  // diagonalize
  Tensor<double> V, epolar;
  syev(polar_tensor, V, epolar);
  double Dpolar_average = 0.0;
  double Dpolar_iso = 0.0;
  for (unsigned int i = 0; i < 3; ++i)
    Dpolar_average = Dpolar_average + epolar[i];
  Dpolar_average = Dpolar_average / 3.0;
  Dpolar_iso =
      sqrt(.5) * sqrt(std::pow(polar_tensor(0, 0) - polar_tensor(1, 1), 2) +
                      std::pow(polar_tensor(1, 1) - polar_tensor(2, 2), 2) +
                      std::pow(polar_tensor(2, 2) - polar_tensor(0, 0), 2));

  if (world.rank() == 0) {
    print("\nTotal Dynamic Polarizability Tensor");
    printf("\nFrequency  = %.6f a.u.\n\n", Rparams.omega);
    // printf("\nWavelength = %.6f a.u.\n\n", Rparams.omega * ???);
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
  molresponse::end_timer(world, "total:");
}
// compute the frequency response, Rparams sets the the calculation type.
// options are dipole,nuclear,order2, order3.  Computes the density respone
// No matter the calculation type we do the same iteration.
// The only difference is the number of response states as well as the
// number of right hand side vectors.
void TDDFT::ComputeFrequencyResponse(World& world,
                                     std::string property,
                                     response_space& x,
                                     response_space& y) {
  // Get start time
  this->x_response = x;
  this->y_response = y;

  molresponse::start_timer(world);

  // Warm and fuzzy
  if (world.rank() == 0) {
    print("\n\n    Response Calculation");
    print("   ------------------------");
  }

  // Keep a copy of dipoles * MO (needed explicitly in eq.)

  // For each protocol
  for (unsigned int proto = 0; proto < Rparams.protocol_data.size(); proto++) {
    // Set defaults inside here
    // default value of
    set_protocol<3>(world, Rparams.protocol_data[proto]);

    check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());
    // Do something to ensure all functions have same k value

    if (property.compare("dipole") == 0) {
      if (world.rank() == 0) print("creating dipole property operator");
      p = Property(world, "dipole");
    } else if (property.compare("nuclear") == 0) {
      if (world.rank() == 0) print("creating nuclear property operator");
      p = Property(world, "nuclear", Gparams.molecule);
    }

    // Create guesses if no response functions
    // If restarting, load here
    if (proto == 0) {
      if (Rparams.restart) {
        if (world.rank() == 0)
          print("   Initial guess from file:", Rparams.restart_file);
        load(world, Rparams.restart_file);
        check_k(
            world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());

        if (Rparams.dipole) {
          // set states
          this->P = PropertyRHS(world, p);
          this->Q = P.copy();

          // set RHS_Vector
        } else if (Rparams.nuclear) {
          // set guesses
          // print("Creating X for Nuclear Operators");

          // print("x norms:");
          // print(x_response.norm2());

          this->P = PropertyRHS(world, p);
          this->Q = P.copy();
        } else if (Rparams.order2) {
          //
        } else if (Rparams.order3) {
          //
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      } else {  // Dipole guesses

        if (Rparams.dipole) {
          // set states
          this->P = PropertyRHS(world, p);
          this->Q = P.copy();
          //
          print("okay this is not a good idea if it comes up more than once");
          // set RHS_Vector
        } else if (Rparams.nuclear) {
          // set guesses
          // print("Creating X for Nuclear Operators");
          // create zero guesses
          // print("x norms:");
          // print(x_response.norm2());

          this->P = PropertyRHS(world, p);
          this->Q = P.copy();
        } else if (Rparams.order2) {
          //
        } else if (Rparams.order3) {
          //
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      }
    }
    //
    // Here i should print some information about the calculation we are
    // about to do
    print("Preiteration Information");
    print("Property : ", property);
    print("Number of Response States: ", Rparams.states);
    print("Number of Ground States: ", Gparams.num_orbitals);
    print("k = ", FunctionDefaults<3>::get_k());
    print("protocol threshold = ", FunctionDefaults<3>::get_k());

    print("Property rhs func k = ", P[0][0].k());
    print("Property func k thresh= ", P[0][0].thresh());

    print("Property rhs func Q k = ", Q[0][0].k());
    print("Property func Q k thresh = ", Q[0][0].thresh());

    // Now actually ready to iterate...
    IterateFrequencyResponse(world, P, Q);
  }  // end for --finished reponse density

  // Have response function, now calculate polarizability for this axis
  // polarizability(world, polar_tensor);
  // PrintPolarizabilityAnalysis(world, polar_tensor, omega);

  // Print total time
  // Precision is set to 10 coming in, drop it to 2
  std::cout.precision(2);
  std::cout << std::fixed;

  // Get start time
  molresponse::end_timer(world, "total:");
}  // end compute frequency response

// Exactam eam
// Deuces