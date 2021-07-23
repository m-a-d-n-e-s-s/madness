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

#include <../chem/NWChem.h>  // For nwchem interface
#include <../chem/SCFOperators.h>
#include <../chem/molecule.h>
#include <chem/potentialmanager.h>
#include <chem/projector.h>  // For easy calculation of (1 - \hat{\rho}^0)
#include <madness/mra/funcdefaults.h>
#include <madness/world/worldmem.h>
#include <math.h>
#include <molresponse/Plot_VTK.h>
#include <molresponse/basic_operators.h>
#include <molresponse/density.h>
#include <molresponse/global_functions.h>
#include <molresponse/property.h>
#include <molresponse/response_functions.h>
#include <molresponse/timer.h>
#include <molresponse/x_space.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

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

TDDFT::TDDFT(World& world, density_vector& rho)
    : r_params(rho.r_params),
      g_params(rho.g_params),
      molecule(rho.molecule),
      omega(rho.omega),
      xcf(rho.xcf),
      rho(rho) {
  // Start the timer
  Chi = rho.Chi;
  PQ = rho.PQ;

  ground_orbitals = g_params.orbitals();
  ground_energies = g_params.get_energies();

  molresponse::start_timer(world);

  // Broadcast to all other nodes
  world.gop.broadcast_serializable(r_params, 0);
  world.gop.broadcast_serializable(molecule, 0);

  // Read in archive
  // Create the projector Qhat to be used in any calculation

  // Set some function defaults
  FunctionDefaults<3>::set_cubic_cell(-r_params.L(), r_params.L());
  FunctionDefaults<3>::set_truncate_mode(1);
  FunctionDefaults<3>::set_truncate_on_project(true);

  // Create the masking function
  mask = real_function_3d(
      real_factory_3d(world).f(mask3).initial_level(4).norefine());

  if (world.size() > 1) {
    // Start a timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);
    if (world.rank() == 0) print("");  // Makes it more legible

    LoadBalanceDeux<3> lb(world);
    for (unsigned int j = 0; j < r_params.num_orbitals(); j++) {
      lb.add_tree(ground_orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
    }
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

    if (r_params.print_level() >= 1)
      molresponse::end_timer(world, "Load balancing:");
  }
}
X_space& TDDFT::GetXspace() { return Chi; }
X_space& TDDFT::GetPQspace() { return PQ; }
// Get response parameters
ResponseParameters TDDFT::GetResponseParameters() { return r_params; }
GroundParameters TDDFT::GetGroundParameters() { return g_params; }
PropertyBase TDDFT::GetPropertyObject() { return p; }
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
  ar& r_params.archive();
  ar& r_params.tda();
  ar& r_params.num_orbitals();
  ar& r_params.n_states();
  ar& omega;

  for (size_t i = 0; i < r_params.n_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++) ar& Chi.X[i][j];
  if (not r_params.tda()) {
    for (size_t i = 0; i < r_params.n_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++) ar& Chi.Y[i][j];
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

  ar& r_params.archive();
  ar& r_params.tda();
  ar& r_params.num_orbitals();
  ar& r_params.n_states();
  ar& omega;

  Chi = X_space(world, r_params.n_states(), r_params.num_orbitals());

  for (size_t i = 0; i < r_params.n_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++) ar& Chi.X[i][j];
  world.gop.fence();

  if (not r_params.tda()) {
    for (size_t i = 0; i < r_params.n_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++) ar& Chi.Y[i][j];
    world.gop.fence();
  }
}

void TDDFT::initial_load_bal(World& world) {
  LoadBalanceDeux<3> lb(world);
  real_function_3d vnuc;
  vnuc = potentialmanager->vnuclear();
  lb.add_tree(vnuc,
              lbcost<double, 3>(r_params.vnucextra() * 1.0,
                                r_params.vnucextra() * 8.0));
  FunctionDefaults<3>::redistribute(world,
                                    lb.load_balance(r_params.loadbalparts()));
}

void TDDFT::loadbal(World& world,
                    vecfuncT rho_omega,
                    X_space Chi,
                    X_space Chi_old) {
  molresponse::start_timer(world);
  if (world.size() == 1) return;

  LoadBalanceDeux<3> lb(world);
  real_function_3d vnuc;
  vnuc = potentialmanager->vnuclear();
  lb.add_tree(
      vnuc,
      lbcost<double, 3>(r_params.vnucextra() * 1.0, r_params.vnucextra() * 8.0),
      false);
  for (size_t i = 0; i < Chi.X.size(); ++i) {
    lb.add_tree(rho_omega[i], lbcost<double, 3>(1.0, 8.0), false);
  }
  for (size_t i = 0; i < Chi.X.size(); ++i) {
    for (size_t j = 0; j < Chi.X.size_orbitals(); ++j) {
      lb.add_tree(Chi.X[i][j], lbcost<double, 3>(1.0, 8.0), false);
    }
  }
  if (r_params.omega() != 0) {
    for (size_t i = 0; i < Chi.X.size(); ++i) {
      for (size_t j = 0; j < Chi.X.size_orbitals(); ++j) {
        lb.add_tree(Chi.Y[i][j], lbcost<double, 3>(1.0, 8.0), false);
      }
    }
  }

  world.gop.fence();

  FunctionDefaults<3>::redistribute(
      world,
      lb.load_balance(r_params.loadbalparts()));  // 6.0 needs retuning after
                                                  // param.vnucextra

  world.gop.fence();
  molresponse::end_timer(world, "Load balancing");
  print_meminfo(world.rank(), "Load balancing");
}
// compute pmap based on ground and first order orbitals
// set default pmap to new pmap
// make orbital copies using new pmap
void TDDFT::orbital_load_balance(World& world,
                                 vecfuncT& psi0,
                                 vecfuncT& psi0_copy,
                                 X_space& Chi,
                                 X_space& Chi_copy) {
  size_t m = r_params.n_states();
  size_t n = r_params.num_orbitals();

  molresponse::start_timer(world);
  if (world.size() > 1) {
    LoadBalanceDeux<3> lb(world);
    for (unsigned int i = 0; i < m; ++i) {
      lb.add_tree(psi0[i], lbcost<double, 3>(1.0, 8.0), false);
      for (unsigned int j = 0; j < n; ++j) {
        // add a tree for orbitals
        lb.add_tree(Chi.X[i][j], lbcost<double, 3>(1.0, 8.0), false);
        lb.add_tree(Chi.Y[i][j], lbcost<double, 3>(1.0, 8.0), false);
      }
    }

    world.gop.fence();

    // newpamap is the new pmap just based on the orbitals
    std::shared_ptr<WorldDCPmapInterface<Key<3>>> newpmap =
        lb.load_balance(r_params.loadbalparts());
    molresponse::end_timer(world, "Gamma compute loadbal");
    // default process map
    // We set the newpmap
    molresponse::start_timer(world);
    FunctionDefaults<3>::set_pmap(newpmap);  // set default to be new

    world.gop.fence();
    // copy orbitals using new pmap
    Chi_copy = Chi.copy(newpmap, false);
    world.gop.fence();  // then fence

    psi0_copy = copy(world, ground_orbitals, newpmap, false);
    world.gop.fence();  // then fence
  }
  molresponse::end_timer(world, "Gamma redist");
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
void TDDFT::normalize(World& world, X_space& Chi) {
  // Run over rows
  for (size_t i = 0; i < size_states(Chi); i++) {
    // Get the normalization constant
    // (Sum included inside inner)
    double normf = inner(Chi.X[i], Chi.X[i]);
    double normg = inner(Chi.Y[i], Chi.Y[i]);
    double norm = sqrt(normf - normg);
    print("---------------------Normalize--------------");
    print(normf);
    print(normg);
    print(norm);

    // Doing this to deal with zero functions.
    // Maybe not smrt.
    if (norm == 0) continue;
    Chi.X[i] = Chi.X[i] * (1.0 / norm);
    Chi.Y[i] = Chi.Y[i] * (1.0 / norm);
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
X_space TDDFT::create_trial_functions(World& world,
                                      size_t k,
                                      std::vector<real_function_3d>& orbitals,
                                      size_t print_level) {
  // Get size
  print("In create trial functions");
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
  response_space trials_X;

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
      trials_X.push_back(temp);
      count++;
    }

    // Stop when we first get beyond k components
    if (count >= k) break;
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) print("   Norms of guess functions:");
    print_norms(world, trials_X);
  }

  // Truncate
  madness::truncate(
      world, trials_X, madness::FunctionDefaults<3>::get_thresh(), true);

  X_space trials(world, count, n);
  trials.X = trials_X.copy();
  trials_X.clear();

  // Done
  return trials;
}

// Returns initial guess functions as
// ground MO * <x,y,z>
X_space TDDFT::create_trial_functions2(World& world,
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

  X_space trials(world, 3 * n * n, n);
  for (size_t i = 0; i < n; i++) {
    for (size_t d = 0; d < directions; d++) {
      for (size_t o = 0; o < n; o++) {
        //        trials[i + j + o][o] = functions[i][j];
        trials.X[count][o] = copy(functions.at(d).at(o));
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
    print_norms(world, trials.X);
  }

  // Truncate
  madness::truncate(world, trials.X);

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

response_space TDDFT::PropertyRHS(World& world, PropertyBase& p) const {
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  response_space rhs(world, p.num_operators, r_params.num_orbitals());

  reconstruct(world, ground_orbitals);
  QProjector<double, 3> Qhat(world, ground_orbitals);
  // Set the dipoles (ground orbitals are probably
  // more accurate now, so recalc the dipoles)
  // why is it called dipole guess.
  // This is just orbitals times dipole operator
  std::vector<real_function_3d> orbitals = ground_orbitals;

  print("num operators ", p.num_operators);
  for (size_t i = 0; i < p.num_operators; i++) {
    // question here....MolecularDerivativeFunctor takes derivative with
    // respect to axis atom and axis
    // here we save
    // need to project

    rhs[i] =
        mul(world, p.operator_vector.at(i), ground_orbitals, r_params.lo());

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
// Calculates ground state coulomb potential
real_function_3d TDDFT::Coulomb(World& world) {
  // Coulomb operator
  real_convolution_3d op =
      CoulombOperator(world, r_params.lo(), FunctionDefaults<3>::get_thresh());

  // Get density
  std::vector<real_function_3d> vsq = square(world, ground_orbitals);
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
response_space TDDFT::exchange(World& world, response_space& f) {
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // Container for results and others
  response_space result(world, m, n);
  real_function_3d psif = real_function_3d(world);

  // Modified 'small memory' algorithm from SCF.cc
  f.reconstruct_rf();
  // K[rho0]xp=\sum_i |i><i|J|xp>
  // for each state
  // for each occ orbital
  //		- create the vector of orbital products
  //		- gaxpy to collect into result
  // compute product <i|xp>

  // Run over each excited state
  for (size_t k = 0; k < m; k++) {
    // And run over each occupied state
    for (size_t j = 0; j < n; j++) {
      // Get a vector of transition densities
      //      v[i]=<j|f[k][i]>  phi[i]*f[k][:]
      auto phix = mul_sparse(
          world, ground_orbitals[j], f[k], FunctionDefaults<3>::get_thresh());
      // Clean up
      truncate(world, phix);
      // Apply operator to each member of vector
      phix = apply(world, *coulop, phix);
      // Clean up
      truncate(world, phix);
      // Final multiplication of each member of vector by a single function
      phix = mul_sparse(
          world, ground_orbitals[j], phix, FunctionDefaults<3>::get_thresh());
      // Add the vector to result
      gaxpy(world, 1.0, result[k], 1.0, phix);
    }
  }

  // Truncate
  madness::truncate(world, result);

  // Done!
  return result;
}

void TDDFT::make_nuclear_potential(World& world) {
  molresponse::start_timer(world);
  potentialmanager = std::shared_ptr<PotentialManager>(
      new PotentialManager(molecule, r_params.core_type()));
  potentialmanager->make_nuclear_potential(world);
}

// Returns the ground state potential applied to functions f
// (V0 f) V0=(Vnuc+J0-K0+EXC0)
// J0=J[rho0]
// K0=K[rho0]f
// EXC0=EXC[rho0]

// result(i,j) = inner(a[i],b[j]).sum()
Tensor<double> TDDFT::expectation(World& world,
                                  const response_space& A,
                                  const response_space& B) {
  // Get sizes
  MADNESS_ASSERT(A.size() > 0);
  MADNESS_ASSERT(A.size() == B.size());
  MADNESS_ASSERT(A[0].size() > 0);
  MADNESS_ASSERT(A[0].size() == B[0].size());

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
// compute rms and maxabsval of vector of doubles
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
    if (anorm[i] > r_params.maxrotn()) {
      double s = r_params.maxrotn() / anorm[i];
      ++nres;
      if (world.rank() == 0) {
        if (nres == 1 and (r_params.print_level() > 1))
          printf("  restricting step for %s orbitals:", spin.c_str());
        printf(" %d", i);
      }
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
    }
  }
  if (nres > 0 && world.rank() == 0 and (r_params.print_level() > 1))
    printf("\n");

  world.gop.fence();
  double rms, maxval;
  vector_stats(anorm, rms, maxval);
  if (world.rank() == 0 and (r_params.print_level() > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}

double TDDFT::do_step_restriction(World& world,
                                  const vecfuncT& x,
                                  const vecfuncT& y,
                                  vecfuncT& x_new,
                                  vecfuncT& y_new,
                                  std::string spin) const {
  // sub(world, x, x_new)
  vecfuncT x_diff = sub(world, x, x_new);
  vecfuncT y_diff = sub(world, y, y_new);
  // sub(world, x, x_new)
  std::vector<double> anorm_x = norm2s(world, x_diff);
  std::vector<double> anorm_y = norm2s(world, y_diff);
  std::vector<double> anorm;
  for (unsigned int i = 0; i < x.size(); ++i) {
    anorm.push_back(std::sqrt(anorm_x.at(i) * anorm_x.at(i) +
                              anorm_y.at(i) * anorm_y.at(i)));
  }
  size_t nres = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    print("anorm ", i, " : ", anorm[i]);
    if (anorm[i] > r_params.maxrotn()) {
      double s = r_params.maxrotn() / anorm[i];
      ++nres;
      if (world.rank() == 0) {
        if (nres == 1 and (r_params.print_level() > 1))
          printf("  restricting step for %s orbitals:", spin.c_str());
        printf(" %d", i);
      }
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
      y_new[i].gaxpy(s, y[i], 1.0 - s, false);
    }
  }
  if (nres > 0 && world.rank() == 0 and (r_params.print_level() > 1))
    printf("\n");

  world.gop.fence();
  double rms, maxval;
  vector_stats(anorm, rms, maxval);
  if (world.rank() == 0 and (r_params.print_level() > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}
// Construct the Hamiltonian
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
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

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
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Apply shift:");

  // Done
  return shifted_V;
}

// Returns the given shift applied to the given potential
response_space TDDFT::apply_shift(World& world,
                                  double& shift,
                                  response_space& V,
                                  response_space& f) {
  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

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
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Apply shift:");

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
                            double lo,
                            double thresh) {
  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

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
      double mu = sqrt(-2.0 * (ground(p) + omega(k) + shift(k, p)));
      print("res state ", k, " orb ", p, " bsh exponent mu :", mu);
      temp[p] = std::shared_ptr<SeparatedConvolution<double, 3>>(
          BSHOperatorPtr3D(world, mu, lo, thresh));
    }

    // Add intermediary to return container
    operators.push_back(temp);
  }

  // End timer
  if (r_params.print_level() >= 1)
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
                                       double lo,
                                       double thresh) {
  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

  // Sizes inferred from ground and omega
  size_t n = ground.size();  // number of orbitals
  size_t num_states = r_params.n_states();
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
                             lo,
                             thresh));
      }
      operators.push_back(temp);
    }

    // Add intermediary to return container
  }

  // End timer
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Creating BSH ops:");

  // Done
  return operators;
}

std::vector<poperatorT> TDDFT::make_bsh_operators_response(
    World& world,
    double& shift,
    double& omega) const {
  if (r_params.print_level() >= 1) molresponse::start_timer(world);
  double tol = FunctionDefaults<3>::get_thresh();

  // Sizes inferred from ground and omega
  size_t num_orbitals = ground_energies.size();  // number of orbitals
  std::vector<poperatorT> ops(num_orbitals);
  // Run over occupied components
  for (size_t p = 0; p < num_orbitals; p++) {
    double mu = sqrt(-2.0 * (ground_energies(p) + omega + shift));
    ops[p] = poperatorT(BSHOperatorPtr3D(world, mu, r_params.lo(), tol));
  }
  return ops;
  // End timer
}
// shift
std::vector<std::shared_ptr<real_convolution_3d>>
TDDFT::CreateBSHOperatorPropertyVector(World& world,
                                       double& shift,
                                       Tensor<double>& ground,
                                       double& omega,
                                       double lo,
                                       double eps) {
  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

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
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Creating BSH ops:");

  // Done
  return ghat_operators;
}
void TDDFT::update_x_space_response(World& world,
                                    X_space& old_Chi,
                                    X_space& Chi,
                                    X_space& res,
                                    XCOperator<double, 3>& xc,
                                    std::vector<poperatorT>& bsh_x_ops,
                                    std::vector<poperatorT>& bsh_y_ops,
                                    QProjector<double, 3>& projector,
                                    double& x_shifts,
                                    double& omega_n,
                                    NonLinearXsolver kain_x_space,
                                    std::vector<X_vector> Xvector,
                                    std::vector<X_vector> Xresidual,
                                    Tensor<double>& bsh_residualsX,
                                    Tensor<double>& bsh_residualsY,
                                    size_t iteration) {
  size_t m = Chi.num_states();
  bool compute_y = r_params.omega() != 0.0;
  // size_t n = Chi.num_orbitals();

  Tensor<double> errX(m);
  Tensor<double> errY(m);

  X_space theta_X = Compute_Theta_X(world, Chi, xc, r_params.calc_type());
  // compute residual X_space
  // compute errX and errY which are max orbital residuals for each response
  print("BSH update iter = ", iteration);
  X_space temp = bsh_update_response(
      world, old_Chi, Chi, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);

  Tensor<double> G;
  G = -2 * inner(temp, PQ);
  print("Polarizability Tensor bsh update <temp|r|PQ>");
  print(G);

  // computes residual from old Chi and temp
  res = compute_residual(
      world, old_Chi, temp, bsh_residualsX, bsh_residualsY, compute_y);

  print("print residual norms from bsh update");
  print_residual_norms(world, res, compute_y, iteration);

  // kain update with temp adjusts temp
  // TODO test if default zero init guess
  if (r_params.kain() && iteration > 0) {
    print("Kain update iter = ", iteration);
    kain_x_space_update(world, temp, res, kain_x_space, Xvector, Xresidual);
    G = -2 * inner(temp, PQ);
    print("Polarizability Tensor kain update");
    print(G);
  }

  if (iteration > 0 && true) {
    print("Restrict Step update iter = ", iteration);
    x_space_step_restriction(world, old_Chi, temp, compute_y);
    G = -2 * inner(temp, PQ);
    print("Polarizability Tensor step restriction update");
    print(G);
  }

  // truncate x
  temp.X.truncate_rf();
  // truncate y if compute y
  if (compute_y) temp.Y.truncate_rf();
  //	if not compute y then copy x in to y
  if (!compute_y) temp.Y = temp.X.copy();

  Chi = temp.copy();
  G = -2 * inner(Chi, PQ);
  print("Polarizability Tensor  copy temp into Chi update");
  print(G);
  // print x norms

  if (r_params.print_level() >= 1) {
    print("Chi.x norms in iteration after truncate: ", iteration);
    print(Chi.X.norm2());

    print("Chi.y norms in iteration after truncate: ", iteration);
    print(Chi.Y.norm2());
  }
}
X_space TDDFT::compute_residual(World& world,
                                X_space& old_Chi,
                                X_space& temp,
                                Tensor<double>& bsh_residualsX,
                                Tensor<double>& bsh_residualsY,
                                bool compute_y) {
  size_t m = old_Chi.X.size();
  size_t n = old_Chi.X.size_orbitals();
  molresponse::start_timer(world);
  //	compute residual
  X_space res(world, m, n);
  res.X = old_Chi.X - temp.X;
  if (compute_y) {
    res.Y = old_Chi.Y - temp.Y;
  }
  // Basic outp:w
  //
  //*************************
  Tensor<double> errX(m);
  Tensor<double> errY(m);
  // rmsX and maxvalX for each m response states
  std::vector<double> rmsX(m), maxvalX(m);
  std::vector<std::vector<double>> rnormsX;
  std::vector<std::vector<double>> rnormsY;
  // find the norms of each of the residual response vectors
  for (size_t i = 0; i < m; i++) {
    // the 2norms of each of the orbitals in response vector
    rnormsX.push_back(norm2s(world, res.X[i]));
    if (world.rank() == 0 and (r_params.print_level() > 1))
      print("residuals X: state ", i, " : ", rnormsX[i]);
    // maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
    // maxvalX= largest abs(v[i])
    vector_stats(rnormsX[i], rmsX[i], maxvalX[i]);
    // errX[i] is the largest residual orbital value
    errX[i] = maxvalX[i];
    if (world.rank() == 0 and (r_params.print_level() > 1))
      print("BSH residual: rms", rmsX[i], "   max", maxvalX[i]);
  }
  if (compute_y) {
    std::vector<double> rmsY(m), maxvalY(m);
    for (size_t i = 0; i < m; i++) {
      rnormsY.push_back(norm2s(world, res.Y[i]));
      if (world.rank() == 0 and (r_params.print_level() > 1))
        print("residuals Y: state ", i, " : ", rnormsY[i]);
      vector_stats(rnormsY[i], rmsY[i], maxvalY[i]);
      errY[i] = maxvalY[i];
      if (world.rank() == 0 and (r_params.print_level() > 1))
        print("BSH residual: rms", rmsY[i], "   max", maxvalY[i]);
    }
  }
  molresponse::end_timer(world, "BSH residual");

  if (r_params.print_level() >= 1) {
    print("res.X norms in iteration after compute_residual function: ");
    print(res.X.norm2());

    print("res.Y norms in iteration after compute_residual function: ");
    print(res.Y.norm2());
  }

  bsh_residualsX = errX;
  bsh_residualsY = errY;
  // Apply shifts and rhs
  // Next calculate 2-norm of these vectors of differences
  return res;
}

void TDDFT::print_residual_norms(World& world,
                                 X_space& res,
                                 bool compute_y,
                                 size_t iteration) {
  size_t m = res.num_states();

  Tensor<double> x_norms(m);
  Tensor<double> y_norms(m);
  for (size_t i = 0; i < m; i++) x_norms(i) = norm2(world, res.X[i]);
  if (compute_y) {
    for (size_t i = 0; i < m; i++) y_norms(i) = norm2(world, res.Y[i]);
  }
  if (r_params.print_level() >= 0 and world.rank() == 0) {
    if (compute_y) {
      std::cout << "res " << iteration << " X :";
      for (size_t i(0); i < m; i++) {
        std::cout << x_norms[i] << "  ";
      }
      std::cout << " Y :";
      for (size_t i(0); i < m; i++) {
        std::cout << y_norms[i] << "  ";
      }
      std::cout << endl;
    } else {
      print("resX ", iteration, " :", x_norms);
    }
  }
}

X_space TDDFT::bsh_update_response(World& world,
                                   X_space& old_Chi,
                                   X_space& Chi,
                                   X_space& theta_X,
                                   std::vector<poperatorT>& bsh_x_ops,
                                   std::vector<poperatorT>& bsh_y_ops,
                                   QProjector<double, 3>& projector,
                                   double& x_shifts) {
  size_t m = Chi.X.size();
  size_t n = Chi.X.size_orbitals();
  bool compute_y = r_params.omega() != 0.0;

  molresponse::start_timer(world);

  theta_X.X += Chi.X * x_shifts;
  theta_X.X += PQ.X;
  theta_X.X = theta_X.X * -2;
  theta_X.X.truncate_rf();

  if (compute_y) {
    theta_X.Y += PQ.Y;
    theta_X.Y = theta_X.Y * -2;
    theta_X.Y.truncate_rf();
  }
  molresponse::end_timer(world, "Compute residual stuff theta_X");

  // apply bsh
  molresponse::start_timer(world);
  X_space bsh_X(world, m, n);

  bsh_X.X = apply(world, bsh_x_ops, theta_X.X);
  if (compute_y) {
    bsh_X.Y = apply(world, bsh_y_ops, theta_X.Y);
  }
  molresponse::end_timer(world, "Apply BSH to theta_X");

  molresponse::start_timer(world);
  // Project out ground state
  for (size_t i = 0; i < m; i++) bsh_X.X[i] = projector(bsh_X.X[i]);
  if (compute_y) {
    for (size_t i = 0; i < m; i++) bsh_X.Y[i] = projector(bsh_X.Y[i]);
  }
  molresponse::end_timer(world, "Project out BSH_X");

  molresponse::start_timer(world);
  X_space temp(world, m, n);
  temp.X = bsh_X.X.copy();
  if (compute_y) {
    temp.Y = bsh_X.Y.copy();
  }
  temp.truncate();
  molresponse::end_timer(world, "Trucate bsh_X");

  return temp;
}

void TDDFT::update_x_space_excited(World& world,
                                   X_space& old_Chi,
                                   X_space& Chi,
                                   X_space& old_Lambda_X,
                                   X_space& res,
                                   XCOperator<double, 3>& xc,
                                   QProjector<double, 3>& projector,
                                   Tensor<double>& omega,
                                   NonLinearXsolver kain_x_space,
                                   std::vector<X_vector> Xvector,
                                   std::vector<X_vector> Xresidual,
                                   Tensor<double>& energy_residuals,
                                   Tensor<double>& old_energy,
                                   Tensor<double>& bsh_residualsX,
                                   Tensor<double>& bsh_residualsY,
                                   Tensor<double>& S,
                                   Tensor<double>& old_S,
                                   Tensor<double>& A,
                                   Tensor<double>& old_A,
                                   std::vector<bool>& converged,
                                   size_t iter) {
  size_t m = Chi.num_states();
  bool compute_y = not r_params.tda();

  Tensor<double> errX(m);
  Tensor<double> errY(m);
  Tensor<double> x_shifts(m);
  Tensor<double> y_shifts(m);

  X_space Lambda_X = Compute_Lambda_X(world, Chi, xc, r_params.calc_type());
  // This diagonalizes XAX and computes new omegas
  // updates Chi
  compute_new_omegas_transform(world,
                               old_Chi,
                               Chi,
                               old_Lambda_X,
                               Lambda_X,
                               omega,
                               old_energy,
                               S,
                               old_S,
                               A,
                               old_A,
                               energy_residuals,
                               iter);

  old_Chi = Chi.copy();

  // Analysis gets messed up if BSH is last thing applied
  // so exit early if last iteration
  if (iter == r_params.maxiter() - 1) {
    print("Reached max iter");
  } else {
    X_space theta_X = Compute_Theta_X(world, Chi, xc, r_params.calc_type());
    //  Calculates shifts needed for potential / energies
    X_space temp =
        bsh_update_excited(world, old_Chi, Chi, theta_X, projector, converged);

    res = compute_residual(
        world, old_Chi, temp, bsh_residualsX, bsh_residualsY, compute_y);

    if (iter > 0) {
      x_space_step_restriction(world, old_Chi, temp, compute_y);
    }
    if (r_params.kain() && (iter > 0 || r_params.first_run())) {
      kain_x_space_update(world, temp, res, kain_x_space, Xvector, Xresidual);
    }
    temp.X.truncate_rf();
    if (!compute_y) temp.Y = temp.X.copy();
    if (compute_y) temp.Y.truncate_rf();
    if (r_params.print_level() >= 1) {
      print("Chi.x norms in iteration after truncate: ", iter);
      print(Chi.X.norm2());

      print("Chi.y norms in iteration after truncate: ", iter);
      print(Chi.Y.norm2());
    }
    // print x norms
    Chi = temp.copy();
  }

  // Apply mask
  /*
for (size_t i = 0; i < m; i++) Chi.X[i] = mask * Chi.X[i];
if (not r_params.tda()) {
for (size_t i = 0; i < m; i++) Chi.Y[i] = mask * Chi.Y[i];
}
  */
}

// Load Balancing
void TDDFT::compute_new_omegas_transform(World& world,
                                         X_space& old_Chi,
                                         X_space& Chi,
                                         X_space& old_Lambda_X,
                                         X_space& Lambda_X,
                                         Tensor<double>& omega,
                                         Tensor<double>& old_energy,
                                         Tensor<double>& S,
                                         Tensor<double>& old_S,
                                         Tensor<double>& A,
                                         Tensor<double>& old_A,
                                         Tensor<double>& energy_residuals,
                                         size_t iter) {
  size_t m = Chi.X.size();
  // Basic output
  if (r_params.print_level() >= 1 and world.rank() == 0) {
    print("Before Deflate");
    print("\n   Excitation Energies:");
    print("i=", iter, " roots: ", iter, omega);
  }
  if (r_params.tda()) {
    deflateTDA(world,
               Chi,
               old_Chi,
               Lambda_X,
               old_Lambda_X,
               S,
               old_S,
               old_A,
               omega,
               iter,
               m);
    // Constructing S
    // Full TDHF
  } else {
    deflateFull(world,
                Chi,
                old_Chi,
                Lambda_X,
                old_Lambda_X,
                S,
                old_S,
                old_A,
                omega,
                iter,
                m);
  }

  // Basic output
  if (r_params.print_level() >= 1 and world.rank() == 0) {
    print("After Deflate");
    print("\n   Excitation Energies:");
    print("i=", iter, " roots: ", iter, omega);
  }

  // Calculate energy residual and update old_energy
  energy_residuals = abs(omega - old_energy);
  old_energy = copy(omega);
}
X_space TDDFT::bsh_update_excited(World& world,
                                  X_space& old_Chi,
                                  X_space& Chi,
                                  X_space& theta_X,
                                  QProjector<double, 3>& projector,
                                  std::vector<bool>& converged) {
  size_t m = Chi.X.size();
  size_t n = Chi.X.size_orbitals();
  X_space res(world, m, n);
  Tensor<double> x_shifts(m);
  Tensor<double> y_shifts(m);
  x_shifts =
      create_shift(world, ground_energies, omega, r_params.print_level(), "x");
  if (not r_params.tda()) {
    omega = -omega;  // Negative here is so that these Greens functions are
    // (eps - omega)
    y_shifts = create_shift_target(world,
                                   ground_energies,
                                   omega,
                                   ground_energies[n - 1],
                                   r_params.print_level(),
                                   "y");
    omega = -omega;
  }
  // Compute Theta X
  // Apply the shifts
  theta_X.X = apply_shift(world, x_shifts, theta_X.X, Chi.X);
  theta_X.X = theta_X.X * -2;
  theta_X.X.truncate_rf();
  if (not r_params.tda()) {
    y_shifts = -y_shifts;
    theta_X.Y = apply_shift(world, y_shifts, theta_X.Y, Chi.Y);
    theta_X.Y = theta_X.Y * -2;
    theta_X.Y.truncate_rf();
  }
  // Construct BSH operators
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_ops =
      create_bsh_operators(world,
                           x_shifts,
                           ground_energies,
                           omega,
                           r_params.lo(),
                           FunctionDefaults<3>::get_thresh());
  std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_y_ops;
  if (not r_params.tda()) {
    omega = -omega;
    bsh_y_ops = create_bsh_operators(world,
                                     y_shifts,
                                     ground_energies,
                                     omega,
                                     r_params.lo(),
                                     FunctionDefaults<3>::get_thresh());
    omega = -omega;
  }
  X_space bsh_update;
  // Apply BSH and get updated response components
  if (r_params.print_level() >= 1) molresponse::start_timer(world);
  bsh_update.X = apply(world, bsh_x_ops, theta_X.X);
  if (not r_params.tda()) bsh_update.Y = apply(world, bsh_y_ops, theta_X.Y);
  if (r_params.print_level() >= 1) molresponse::end_timer(world, "Apply BSH:");

  // Project out ground state
  for (size_t i = 0; i < m; i++) bsh_update.X[i] = projector(bsh_update.X[i]);
  if (not r_params.tda()) {
    for (size_t i = 0; i < m; i++) bsh_update.Y[i] = projector(bsh_update.Y[i]);
  }

  X_space temp = Chi.copy();
  // Only update non-converged components
  for (size_t i = 0; i < m; i++) {
    if (not converged[i]) {
      temp.X[i] = bsh_update.X[i];
      temp.X[i] = mask * temp.X[i];
      if (not r_params.tda()) {
        temp.Y[i] = bsh_update.Y[i];
        temp.Y[i] = mask * temp.Y[i];
      }
    }
  }
  temp.truncate();

  return temp;
}

void TDDFT::kain_x_space_update(World& world,
                                X_space& temp,
                                X_space& res,
                                NonLinearXsolver& kain_x_space,
                                std::vector<X_vector>& Xvector,
                                std::vector<X_vector>& Xresidual) {
  size_t m = Chi.X.size();
  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    Xvector[b] = (X_vector(temp, b));
    Xresidual[b] = (X_vector(res, b));
  }

  for (size_t b = 0; b < m; b++) {
    X_vector kain_X = kain_x_space[b].update(
        Xvector[b], Xresidual[b], FunctionDefaults<3>::get_thresh(), 3.0);
    temp.X[b].assign(kain_X.X[0].begin(), kain_X.X[0].end());
    temp.Y[b].assign(kain_X.Y[0].begin(), kain_X.Y[0].end());
  }
  molresponse::end_timer(world, " KAIN update:");
}

void TDDFT::x_space_step_restriction(World& world,
                                     X_space& old_Chi,
                                     X_space& temp,
                                     bool restrict_y) {
  size_t m = Chi.X.size();
  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    if (restrict_y) {
      do_step_restriction(world,
                          old_Chi.X[b],
                          temp.X[b],
                          old_Chi.Y[b],
                          temp.Y[b],
                          "x and y_response");
    } else {
      do_step_restriction(world, old_Chi.X[b], temp.X[b], "x_response");
    }
  }
  molresponse::end_timer(world, " Step Restriction:");
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

vecfuncT TDDFT::make_density(World& world, X_space& Chi, bool compute_y) {
  molresponse::start_timer(world);
  vecfuncT rho_omega;
  if (compute_y) {
    rho_omega = transition_density(world, ground_orbitals, Chi.X, Chi.Y);
  } else {
    rho_omega = transition_density(world, ground_orbitals, Chi.X, Chi.X);
  }

  molresponse::end_timer(world, "Make density omega");
  print_meminfo(world.rank(), "Make density omega");
  return rho_omega;
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
// knowledge of ground_orbitals and g_params.energies. Function sets
// act_orbitals and num_act_orbitals.
void TDDFT::select_active_subspace(World& world) {
  // Default output
  if (r_params.print_level() >= 0) {
    // Set print output to something reasonable
    std::cout.precision(2);
    std::cout << std::fixed;

    if (world.rank() == 0)
      print(
          "   Selecting ground state subspace to excite from for "
          "components.");
    if (world.rank() == 0)
      print("   This is all orbitals between",
            r_params.e_range_lo(),
            "and",
            r_params.e_range_hi(),
            "\n");

    // Reset precision
    std::cout.precision(10);
    std::cout << std::scientific;
  }

  // Determine active orbitals based on energy differences
  // from HOMO
  for (unsigned int i = 0; i < r_params.num_orbitals(); i++) {
    if (r_params.e_range_lo() < g_params.get_energies()(i) and
        g_params.get_energies()(i) < r_params.e_range_hi()) {
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
    act_orbitals.push_back(ground_orbitals[active[i]]);
    act_ground_energies(i) =
        g_params.get_energies()(active[i]);  // Put energies on diagonal
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
  if (r_params.print_level() >= 2 and world.rank() == 0) {
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
    if (r_params.print_level() >= 2 and world.rank() == 0) {
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
                                            X_space& Chi,
                                            X_space& Lambda_X,
                                            Tensor<double>& evals,
                                            Tensor<double>& A,
                                            Tensor<double>& S,
                                            const double thresh) {
  // compute the unitary transformation matrix U that diagonalizes
  // the fock matrix
  Tensor<double> U = get_fock_transformation(world, S, A, evals, thresh);

  // Sort into ascending order
  Tensor<int> selected = sort_eigenvalues(world, evals, U);

  // Debugging output
  if (r_params.print_level() >= 2 and world.rank() == 0) {
    print("   U:");
    print(U);
  }

  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

  // transform the orbitals and the potential
  // Truncate happens inside here
  Chi.X = transform(world, Chi.X, U);
  Lambda_X.X = transform(world, Lambda_X.X, U);

  // End timer
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Transform orbs.:");

  // Normalize x
  normalize(world, Chi.X);

  // Debugging output
  if (r_params.print_level() >= 2 and world.rank() == 0) {
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
                    X_space& Chi,
                    X_space& old_Chi,
                    X_space& Lambda_X,
                    X_space& last_Lambda_X,
                    Tensor<double>& S,
                    Tensor<double>& A,
                    Tensor<double>& old_S,
                    Tensor<double>& old_A,
                    size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = Chi.X.size();
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
  // Calculate correct inner products of upper off diagonal

  Tensor<double> off = response_space_inner(Chi.X, last_Lambda_X.X);
  temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);  // top right
  // Now for lower off diagonal
  off = response_space_inner(old_Chi.X, Lambda_X.X);
  temp_A(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);  // bottom left
  temp_A(Slice(0, m - 1), Slice(0, m - 1)) = copy(A);        // xAx top left
  temp_A(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) =
      copy(old_A);  // xoldAxold bottom right
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
  off = expectation(world, Chi.X, old_Chi.X);
  // Use slicing to put in correct spot
  temp_S(Slice(0, m - 1), Slice(m, 2 * m - 1)) =
      copy(off);  // top right <x|xold>
  // Now the lower off diagonal block
  // (Go ahead and cheat and use the transpose...)
  off = transpose(off);  // just transpose <xold|x>
  // Use slicing to put in correct spot
  temp_S(Slice(m, 2 * m - 1), Slice(0, m - 1)) =
      copy(off);  // bottom right <xold|x>
  // Put together the rest of S
  temp_S(Slice(0, m - 1), Slice(0, m - 1)) = copy(S);  // top left <x|x>
  temp_S(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_S);  //<xold|xold>
  // Save temp_S as S_x
  S = copy(temp_S);
  // Add in old vectors to current vectors for the appropriate ones
  // Augment the vectors step
  for (size_t i = 0; i < m; i++) {
    Chi.X.push_back(old_Chi.X[i]);
    Lambda_X.X.push_back(last_Lambda_X.X[i]);
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
                         X_space& Chi,
                         X_space& old_Chi,
                         X_space& Lambda_X,
                         X_space& last_Lambda_X,
                         Tensor<double>& S,
                         Tensor<double>& A,
                         Tensor<double>& old_S,
                         Tensor<double>& old_A,
                         size_t print_level) {
  // Basic output
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = Chi.X.size();
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
  // Calculate correct inner products of upper off diagonal

  Tensor<double> off = inner(Chi, last_Lambda_X);
  temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);  // top right
  // Now for lower off diagonal
  off = inner(old_Chi, Lambda_X);
  temp_A(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);  // bottom left
  temp_A(Slice(0, m - 1), Slice(0, m - 1)) = copy(A);        // xAx top left
  temp_A(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) =
      copy(old_A);  // xoldAxold bottom right
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
  off = response_space_inner(Chi.X, old_Chi.X) -
        response_space_inner(Chi.Y, old_Chi.Y);
  // Use slicing to put in correct spot
  temp_S(Slice(0, m - 1), Slice(m, 2 * m - 1)) =
      copy(off);  // top right <x|xold>
  // Now the lower off diagonal block
  // (Go ahead and cheat and use the transpose...)
  off = transpose(off);  // just transpose <xold|x>
  // Use slicing to put in correct spot
  temp_S(Slice(m, 2 * m - 1), Slice(0, m - 1)) =
      copy(off);  // bottom right <xold|x>
  // Put together the rest of S
  temp_S(Slice(0, m - 1), Slice(0, m - 1)) = copy(S);  // top left <x|x>
  temp_S(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_S);  //<xold|xold>
  // Save temp_S as S_x
  S = copy(temp_S);
  // Add in old vectors to current vectors for the appropriate ones
  // Augment the vectors step
  for (size_t i = 0; i < m; i++) {
    Chi.X.push_back(old_Chi.X[i]);
    Chi.Y.push_back(old_Chi.X[i]);
    Lambda_X.X.push_back(last_Lambda_X.X[i]);
    Lambda_X.Y.push_back(last_Lambda_X.X[i]);
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
// If using a larger subspace to diagonalize in, after diagonalization this
// will put everything in the right spot

void TDDFT::unaugment(World& world,
                      X_space& Chi,
                      X_space& old_Chi,
                      X_space& Lambda_X,
                      X_space& last_Lambda_X,
                      Tensor<double>& omega,
                      Tensor<double>& S_x,
                      Tensor<double>& A_x,
                      Tensor<double>& old_S,
                      Tensor<double>& old_A,
                      size_t num_states,
                      size_t iter,
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
      Chi.X.pop_back();
      Lambda_X.X.pop_back();
      //  x_fe.pop_back();
      // V_x_response.pop_back();
      // x_gamma.pop_back();
      // x_response.pop_back();
    }
  }

  // Save the "current" into "old"
  // old_x_fe = x_fe.copy();
  // old_x_gamma = x_gamma.copy();
  // old_V_x_response = V_x_response.copy();

  // Project out ground state
  // QProjector<double, 3> projector(world, ground_orbitals);
  // for(size_t i = 0; i < m; i++) x_response[i] = projector(x_response[i]);
  old_Chi.X = Chi.X.copy();
  last_Lambda_X.X = Lambda_X.X.copy();

  // Copy S into old_S
  // S is identity
  old_S = response_space_inner(Chi.X, Chi.X);

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
  // QProjector<double, 3> projector(world, ground_orbitals);
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

void TDDFT::unaugment_full(World& world,
                           X_space& Chi,
                           X_space& old_Chi,
                           X_space& Lambda_X,
                           X_space& last_Lambda_X,
                           Tensor<double>& omega,
                           Tensor<double>& S_x,
                           Tensor<double>& A_x,
                           Tensor<double>& old_S,
                           Tensor<double>& old_A,
                           size_t num_states,
                           size_t iter,
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
      Chi.X.pop_back();
      Lambda_X.X.pop_back();
      Chi.Y.pop_back();
      Lambda_X.Y.pop_back();
    }
  }
  old_Chi = Chi.copy();
  last_Lambda_X = Lambda_X.copy();

  old_S =
      response_space_inner(Chi.X, Chi.X) - response_space_inner(Chi.Y, Chi.Y);

  old_A = Tensor<double>(num_states, num_states);
  for (size_t i = 0; i < num_states; i++) old_A(i, i) = omega(i);

  // End the timer
  if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}
// Diagonalize the full response matrix, taking care of degenerate components
// Why diagonalization and then transform the x_fe vectors

Tensor<double> TDDFT::diagonalizeFullResponseMatrix(World& world,
                                                    X_space& Chi,
                                                    X_space& Lambda_X,
                                                    Tensor<double>& omega,
                                                    Tensor<double>& S,
                                                    Tensor<double>& A,
                                                    const double thresh,
                                                    size_t print_level) {
  // compute the unitary transformation matrix U that diagonalizes
  // the response matrix
  Tensor<double> U = GetFullResponseTransformation(world, S, A, omega, thresh);

  // Sort into ascending order
  Tensor<int> selected = sort_eigenvalues(world, omega, U);

  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

  Chi.X = transform(world, Chi.X, U);
  Chi.Y = transform(world, Chi.Y, U);

  Lambda_X.X = transform(world, Lambda_X.X, U);
  Lambda_X.Y = transform(world, Lambda_X.Y, U);
  // Transform the vectors of functions
  // Truncate happens in here
  // we do transform here
  // End timer
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "Transform orbs.:");

  // Normalize x and y
  normalize(world, Chi);

  // Debugging output
  if (world.rank() == 0 and print_level >= 2) {
    print("   Eigenvector coefficients from diagonalization:");
    print(U);
  }

  // Return the selected functions
  return U;
}
// Similar to what robert did above in "get_fock_transformation"
Tensor<double> TDDFT::GetFullResponseTransformation(
    World& world,
    Tensor<double>& S,
    Tensor<double>& A,
    Tensor<double>& evals,
    const double thresh_degenerate) {
  // Start timer
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

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
  if (r_params.print_level() >= 2 and world.rank() == 0) {
    print("\n   Singular values of overlap matrix:");
    print(s_vals);
    print("   Left singular vectors of overlap matrix:");
    print(l_vecs);
  }

  // Check how many singular values are less than 10*thresh_degen
  size_t num_zero = 0;
  for (int64_t i = 0; i < s_vals.dim(0); i++) {
    if (s_vals(i) < 10 * thresh_degenerate) {
      if (world.rank() == 0 and num_zero == 0) print("");
      if (world.rank() == 0)
        printf(
            "   Detected singular value (%.8f) below threshold (%.8f). "
            "Reducing subspace size.\n",
            s_vals(i),
            10 * thresh_degenerate);
      num_zero++;
    }
    if (world.rank() == 0 and i == s_vals.dim(0) - 1 and num_zero > 0)
      print("");
  }

  // Going to use these a lot here, so just calculate them
  size_t size_l = s_vals.dim(0);      // number of singular values
  size_t size_s = size_l - num_zero;  // smaller subspace size
  /**
   * @brief l_vecs_s(m,1)
   *
   * @return Tensor<double>
   */
  Tensor<double> l_vecs_s(
      size_l,
      num_zero);                   // number of sv by number smaller than thress
  Tensor<double> copyA = copy(A);  // we copy xAx

  // Transform into this smaller space if necessary
  if (num_zero > 0) {
    // Cut out the singular values that are small
    // (singular values come out in descending order)

    // S(m-sl,m-sl)
    S = Tensor<double>(size_s, size_s);  // create size of new size
    for (size_t i = 0; i < size_s; i++) S(i, i) = s_vals(i);
    // Copy the active vectors to a smaller container
    // left vectors [m,m-sl]
    l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
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
    if (r_params.print_level() >= 2 and world.rank() == 0) {
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
  if (world.rank() == 0 and r_params.print_level() >= 2)
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
  if (num_zero > 0) {
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
  if (r_params.print_level() >= 1)
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

void TDDFT::deflateGuesses(World& world,
                           X_space& Chi,
                           X_space& Lambda_X,
                           Tensor<double>& S,
                           Tensor<double>& omega,
                           size_t& iteration,
                           size_t& m) {
  // XX =Omega XAX
  S = response_space_inner(Chi.X, Chi.X);
  Tensor<double> XAX = response_space_inner(Chi.X, Lambda_X.X);

  // Debugging output
  if (r_params.print_level() >= 2 and world.rank() == 0) {
    print("   Overlap matrix:");
    print(S);
  }
  // Just to be sure dimensions work out, clear omega
  omega.clear();
  diagonalizeFockMatrix(
      world, Chi, Lambda_X, omega, XAX, S, FunctionDefaults<3>::get_thresh());
}
void TDDFT::deflateTDA(World& world,
                       X_space& Chi,
                       X_space& old_Chi,
                       X_space& Lambda_X,
                       X_space& old_Lambda_X,
                       Tensor<double>& S,
                       Tensor<double> old_S,
                       Tensor<double> old_A,
                       Tensor<double>& omega,
                       size_t& iteration,
                       size_t& m) {
  S = response_space_inner(Chi.X, Chi.X);
  Tensor<double> XAX = response_space_inner(Chi.X, Lambda_X.X);

  // Debugging output
  if (r_params.print_level() >= 2 and world.rank() == 0) {
    print("   Overlap matrix:");
    print(S);
  }

  // Augment S_x, A_x, x_gamma, x_response, V_x_response and x_gamma
  // if using a larger subspace and not iteration zero (TODO ---Gotta
  // look at this and make sure it uses my new functions molresponse )
  // by default r_params.larger_subspace() = 0 therefore never uses this
  if (iteration < r_params.larger_subspace() and iteration > 0) {
    print("Using augmented subspace");
    augment(world,
            Chi,
            old_Chi,
            Lambda_X,
            old_Lambda_X,
            S,
            XAX,
            old_S,
            old_A,
            r_params.print_level());
  }

  // Solve Ax = Sxw
  // Just to be sure dimensions work out, clear omega
  omega.clear();
  diagonalizeFockMatrix(
      world, Chi, Lambda_X, omega, XAX, S, FunctionDefaults<3>::get_thresh());

  // If larger subspace, need to "un-augment" everything
  if (iteration < r_params.larger_subspace()) {
    print("Unaugmenting subspace");
    unaugment(world,
              Chi,
              old_Chi,
              Lambda_X,
              old_Lambda_X,
              omega,
              S,
              XAX,
              old_S,
              old_A,
              r_params.n_states(),
              iteration,
              r_params.print_level());
  }
}
void TDDFT::deflateFull(World& world,
                        X_space& Chi,
                        X_space& old_Chi,
                        X_space& Lambda_X,
                        X_space& old_Lambda_X,
                        Tensor<double>& S,
                        Tensor<double> old_S,
                        Tensor<double> old_A,
                        Tensor<double>& omega,
                        size_t& iteration,
                        size_t& m) {
  // Debugging output
  S = response_space_inner(Chi.X, Chi.X) - response_space_inner(Chi.Y, Chi.Y);
  if (world.rank() == 0 and r_params.print_level() >= 2) {
    print("\n   Overlap Matrix:");
    print(S);
  }

  Tensor<double> A = inner(Chi, Lambda_X);
  if (world.rank() == 0 and r_params.print_level() >= 2) {
    print("\n   Lambda Matrix:");
    print(A);
  }

  // Larger subspace augmentation BROKEN!!!!!
  // if(iteration < r_params.larger_subspace() and iteration > 0)
  //{
  //   augment_full(world, S, A,
  //                B_x, x_gamma, x_response, V_x_response, x_fe,
  //                B_y, y_gamma, y_response, V_y_response, y_fe,
  //                old_S, old_A,
  //                old_B_x, old_x_gamma, old_x_response,
  //                old_V_x_response, old_x_fe, old_B_y, old_y_gamma,
  //                old_y_response, old_V_y_response, old_y_fe,
  //                r_params.print_level());
  //}

  // Diagonalize
  // Just to be sure dimensions work out, clear omega
  omega.clear();

  Tensor<double> U =
      diagonalizeFullResponseMatrix(world,
                                    Chi,
                                    Lambda_X,
                                    omega,
                                    S,
                                    A,
                                    FunctionDefaults<3>::get_thresh(),
                                    r_params.print_level());
  // Larger subspace un-augmentation BROKEN!!!!
  // if(iteration < r_params.larger_subspace())
  //{
  //   unaugment_full(world, m, iteration, U, omega, S, A,
  //                  x_gamma, x_response, V_x_response, x_fe, B_x,
  //                  y_gamma, y_response, V_y_response, y_fe, B_y,
  //                  old_S, old_A,
  //                  old_x_gamma, old_x_response, old_V_x_response,
  //                  old_x_fe, old_B_x, old_y_gamma, old_y_response,
  //                  old_V_y_response, old_y_fe, old_B_y,/
  //                  r_params.print_level());
  //}
}
// const double thresh, int print_level) {
// Creates the XCOperator<double,3>  object and initializes it with correct
// parameters
XCOperator<double, 3> TDDFT::create_XCOperator(
    World& world,
    std::vector<real_function_3d> orbitals,
    std::string xc) {
  // First calculate the ground state density
  std::vector<real_function_3d> vsq =
      square(world, ground_orbitals);  // we square each orbital
  compress(world, vsq);  // compress into multiwavelet representation
  real_function_3d rho = real_factory_3d(world);  // create function rho
  rho.compress();
  for (unsigned int i = 0; i < vsq.size(); ++i) {
    rho.gaxpy(1.0, vsq[i], 1.0, false);
  }
  world.gop.fence();

  // And create the object using r_params.xc()
  XCOperator<double, 3> xcop(
      world,
      xc,
      false,
      rho,
      rho);  // world,which xc, spin_polarized? ,spinup, spindown

  return xcop;
}

// Uses an XCOperator<double,3>  to construct v_xc for the ground state
// density Returns d^2/d rho^2 E_xc[rho]
std::vector<real_function_3d> TDDFT::create_fxc(
    World& world,
    std::vector<real_function_3d>& orbitals,
    response_space& f,
    response_space& g) {
  // Create the xcop
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());

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

// Uses an XCOperator<double,3>  to construct v_xc for the ground state
// density Returns d^2/d rho^2 E_xc[rho]
std::vector<real_function_3d> TDDFT::GetWxcOnFDensities(
    World& world,
    const std::vector<real_function_3d>& orbitals,
    const response_space& f) {
  // Create the xcop
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());
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
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());

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
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());
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
// Create and diagonalize the CIS matrix for improved initial guess
response_space TDDFT::diagonalize_CIS_guess(
    World& world,
    std::vector<real_function_3d>& virtuals,
    Tensor<double>& omega,
    std::vector<real_function_3d>& orbitals,
    Tensor<double>& energies,
    double lo,
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
  if (r_params.print_level() >= 1) molresponse::start_timer(world);

  size_t I = -1;  // combined index from i and a, start is -1 so that initial
                  // value
  // is 0
  size_t J = -1;  // combined index from j and b, start is -1 so that initial
                  // value
  // is 0

  const size_t m = virtuals.size();
  const size_t n = orbitals.size();

  Tensor<double> MCIS(m * n, m * n);
  real_convolution_3d op = CoulombOperator(world, lo, thresh);

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

// Simplified iterate scheme for guesses
void TDDFT::iterate_guess(World& world, X_space& guesses) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
  QProjector<double, 3> projector(
      world,
      ground_orbitals);                // Projector to project out ground state
  size_t m = r_params.n_states();      // Number of excited states
  size_t n = r_params.num_orbitals();  // Number of ground state orbitals
  Tensor<double> x_shifts;             // Holds the shifted energy values
  response_space bsh_resp(world, m, n);  // Holds wave function corrections
  response_space V;  // Holds V^0 applied to response functions
  response_space
      shifted_V;          // Holds the shifted V^0 applied to response functions
  Tensor<double> S;       // Overlap matrix of response components for x states
  real_function_3d v_xc;  // For TDDFT

  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());
  // Useful to have
  response_space zeros(world, m, n);

  // Now to iterate
  while (iteration < r_params.guess_max_iter()) {
    // Start a timer for this iteration
    molresponse::start_timer(world);
    //
    size_t N0 = guesses.X.size();
    size_t Np = 2 * m;
    size_t p = r_params.guess_max_iter() - 1;
    size_t Ni = N0;

    // No*exp(log(Np/N0)/p*t) to exponential decay
    // the number of states down to 2*r_params.states
    if (iteration > 1 && r_params.guess_xyz() && Ni > Np) {
      Ni = std::ceil(N0 * std::exp(std::log(static_cast<double>(Np) /
                                            static_cast<double>(N0)) /
                                   static_cast<double>(p) * iteration));
      sort(world, omega, guesses.X);
      print(omega);
      for (size_t i = 0; i < Ni; i++) {
        Ni = 0;
        if (omega[i] < 1) {
          Ni++;
        }
      }
      // this function selects k functions
      guesses.X =
          select_functions(world, guesses.X, omega, Ni, r_params.print_level());
    }

    // Basic output
    if (r_params.print_level() >= 1) {
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
      if (r_params.print_level() >= 1) molresponse::start_timer(world);
      if (world.rank() == 0) print("");  // Makes it more legible

      LoadBalanceDeux<3> lb(world);
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < r_params.n_states(); k++) {
          lb.add_tree(guesses.X[k][j], lbcost<double, 3>(1.0, 8.0), true);
        }
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

      if (r_params.print_level() >= 1)
        molresponse::end_timer(world, "Load balancing:");
    }

    // compute rho_omega
    rho_omega = transition_densityTDA(world, ground_orbitals, guesses.X);
    // Project out ground state
    for (size_t i = 0; i < Ni; i++) guesses.X[i] = projector(guesses.X[i]);

    // Truncate before doing expensive things
    guesses.X.truncate_rf();

    // Normalize after projection
    if (r_params.tda()) normalize(world, guesses.X);
    // (TODO why not normalize if not tda)
    // compute Y = false
    X_space Lambda_X = Compute_Lambda_X(world, guesses, xc, "tda");
    deflateGuesses(world, guesses, Lambda_X, S, omega, iteration, m);
    // Debugging output

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
    if (r_params.print_level() >= 1 and world.rank() == 0) {
      print("\n   Excitation Energies:");
      print("gi=", iteration, " roots: ", omega);
    }

    // Only do BSH if not the last iteration
    if (iteration + 1 < r_params.guess_max_iter()) {
      //  Calculates shifts needed for potential / energies
      //  If none needed, the zero tensor is returned
      x_shifts = create_shift(
          world, ground_energies, omega, r_params.print_level(), "x");

      X_space theta_X = Compute_Theta_X(world, guesses, xc, "tda");
      theta_X.X = apply_shift(world, x_shifts, theta_X.X, guesses.X);
      theta_X.X = theta_X.X * -2;
      theta_X.X.truncate_rf();

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
          bsh_x_operators =
              create_bsh_operators(world,
                                   x_shifts,
                                   ground_energies,
                                   omega,
                                   r_params.lo(),
                                   FunctionDefaults<3>::get_thresh());

      // Apply BSH and get updated components
      if (r_params.print_level() >= 1) molresponse::start_timer(world);
      bsh_resp = apply(world, bsh_x_operators, theta_X.X);
      if (r_params.print_level() >= 1)
        molresponse::end_timer(world, "Apply BSH:");

      // Project out ground state
      for (size_t i = 0; i < Ni; i++) bsh_resp[i] = projector(bsh_resp[i]);
      // Save new components
      guesses.X = bsh_resp;
      // Apply mask
      for (size_t i = 0; i < Ni; i++) guesses.X[i] = mask * guesses.X[i];
    }
    // Update counter
    iteration += 1;
    // Done with the iteration.. truncate
    guesses.X.truncate_rf();

    // Basic output
    if (r_params.print_level() >= 1) {  //
      molresponse::end_timer(world, " This iteration:");
    }
  }
}  // Done with iterate gues
// Simplified iterate scheme for guesses

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
  PotentialManager manager(g_params.molecule(), "a");
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
  if (r_params.store_potential()) {
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
      CoulombOperator(world, r_params.lo(), FunctionDefaults<3>::get_thresh());
  std::vector<real_function_3d> Kf =
      zero_functions_compressed<double, 3>(world, m);
  for (size_t i = 0; i < m; ++i) {
    std::vector<real_function_3d> psif =
        mul_sparse(world, f[i], f, FunctionDefaults<3>::get_thresh());
    truncate(world, psif);
    psif = apply(world, op, psif);
    truncate(world, psif);

    // Save the potential here if we are saving it
    if (r_params.store_potential()) {
      stored_potential.push_back(psif);
    }

    psif = mul_sparse(world, f[i], psif, FunctionDefaults<3>::get_thresh());
    gaxpy(world, 1.0, Kf, 1.0, psif);
  }

  // Only use the exchange above if HF:
  Tensor<double> V;
  real_function_3d v_xc;

  if (r_params.xc() == "hf") {
    // Construct V
    V = matrix_inner(world, f, vf) - matrix_inner(world, f, Kf);
  } else {  // DFT

    XCOperator<double, 3> xcop = create_XCOperator(world, f, r_params.xc());

    real_function_3d v_xc = xcop.make_xc_potential();
    v = v + v_xc;
    std::vector<real_function_3d> vf = v * f;
    if ((*xcop.xc).hf_exchange_coefficient() > 0.0) {
      // XCOperator<double,3>  has member variable xc, which is an
      // xcfunctional which has the hf_exchange_coeff we need here
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
functionT TDDFT::make_ground_density(World& world, const vecfuncT& v) {
  tensorT occ = g_params.get_occ();
  vecfuncT vsq = square(world, v);
  compress(world, vsq);
  functionT rho = factoryT(world);
  rho.compress();
  for (unsigned int i = 0; i < vsq.size(); ++i) {
    rho.gaxpy(1.0, vsq[i], double(1.0), false);
  }
  world.gop.fence();
  vsq.clear();
  return rho;
}

// Creates the transition densities
std::vector<real_function_3d> TDDFT::transition_density(
    World& world,
    std::vector<real_function_3d>& orbitals,
    response_space& x,
    response_space& y) {
  // Get sizes
  size_t m = x.size();

  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);
  x.truncate_rf();
  y.truncate_rf();
  truncate(world, orbitals);
  for (size_t b = 0; b < m; b++) {
    // Run over occupied...
    // y functions are zero if TDA is active
    densities[b] = dot(world, x[b], orbitals);
    densities[b] += dot(world, orbitals, y[b]);
  }

  truncate(world, densities);
  world.gop.fence();
  // Done!
  return densities;
}
std::vector<real_function_3d> TDDFT::transition_densityTDA(
    World& world,
    std::vector<real_function_3d> const& orbitals,
    response_space& x) {
  // Get sizes
  size_t m = x.size();
  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);
  x.truncate_rf();
  truncate(world, ground_orbitals);
  for (size_t b = 0; b < m; b++) {
    // y functions are zero if TDA is active
    densities[b] = dot(world, x[b], ground_orbitals);
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

Tensor<double> TDDFT::polarizability() { return -2 * inner(Chi, PQ); }

void TDDFT::PrintPolarizabilityAnalysis(World& world,
                                        const Tensor<double> polar_tensor) {
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
    // printf("\nWavelength = %.6f a.u.\n\n", r_params.omega() * ???);
    print(polar_tensor);
    printf("\tEigenvalues = ");
    printf("\t %.6f \t %.6f \t %.6f \n", epolar[0], epolar[1], epolar[2]);
    printf("\tIsotropic   = \t %.6f \n", Dpolar_average);
    printf("\tAnisotropic = \t %.6f \n", Dpolar_iso);
    printf("\n");
  }
}
// TODO It makes sense to align the plots with the direction of the
// perturbation
// ??
// TODO In excited state calculations we are going to get various directions.
// We should plot the functions with respect to

/* @brief plot tranistions and ground orbitals in x y z direction
 *
 * @param world
 * @param iteration
 * @param x_response
 * @param y_response
 * @param r_params
 * @param g_params
 */
void TDDFT::PlotGroundandResponseOrbitals(World& world,
                                          size_t iteration,
                                          response_space& x_response,
                                          response_space& y_response,
                                          ResponseParameters const& r_params,
                                          GroundParameters const& g_params) {
  std::filesystem::create_directories("plots/xy");
  std::filesystem::create_directory("plots/ground");
  std::filesystem::create_directory("plots/transition_density");

  // TESTING
  // get transition density
  // num orbitals
  size_t n = x_response[0].size();
  size_t m = x_response.size();

  real_function_3d ground_density =
      dot(world, ground_orbitals, ground_orbitals);
  std::vector<real_function_3d> densities =
      transition_density(world, ground_orbitals, x_response, y_response);
  std::string dir("xyz");
  // for plotname size
  size_t buffSize = 500;
  char plotname[buffSize];
  double Lp = std::min(r_params.L(), 24.0);
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
        plot_line(plotname, 5001, plt.lo, plt.hi, ground_orbitals[i]);
      }
      for (int i = 0; i < static_cast<int>(n); i++) {
        // print ground_state
        // plot x function  x_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/xy/x_direction_%c_res_%d_orb_%d",
                 dir[d],
                 static_cast<int>(b),
                 static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, x_response[b][i]);

        // plot y functione  y_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/xy/y_direction_%c_res_%d_orb_%d",
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

void TDDFT::plot_excited_states(World& world,
                                size_t iteration,
                                response_space& x_response,
                                response_space& y_response,
                                ResponseParameters const& r_params,
                                GroundParameters const& g_params) {
  std::filesystem::create_directories("plots/virtual");
  // num orbitals
  size_t n = x_response[0].size();
  size_t m = x_response.size();

  std::string dir("xyz");
  // for plotname size
  size_t buffSize = 500;
  char plotname[buffSize];
  double Lp = std::min(r_params.L(), 24.0);
  // Doing line plots along each axis
  for (int d = 0; d < 3; d++) {
    // print ground_state
    plotCoords plt(0, Lp);
    for (int b = 0; b < static_cast<int>(m); b++) {
      for (int i = 0; i < static_cast<int>(n); i++) {
        // print ground_state
        // plot x function  x_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/virtual/x_direction_%c_res_%d_orb_%d",
                 dir[d],
                 static_cast<int>(b),
                 static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, x_response[b][i]);

        // plot y functione  y_dir_b_i__k_iter
        snprintf(plotname,
                 buffSize,
                 "plots/xy/y_direction_%c_res_%d_orb_%d",
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
// compute the frequency response, r_params sets the the calculation type.
// options are dipole,nuclear,order2, order3.  Computes the density respone
// No matter the calculation type we do the same iteration.
// The only difference is the number of response states as well as the
// number of right hand side vectors.

// Exactam eam
// Deuces
