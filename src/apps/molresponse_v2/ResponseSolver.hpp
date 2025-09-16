#pragma once
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseVector.hpp"
#include <madness/mra/nonlinsol.h>

struct response_vector_allocator {
  World &world;
  const size_t n_orbtials;
  response_vector_allocator(World &world, size_t n_orbtials)
      : world(world), n_orbtials(n_orbtials) {}
  vector_real_function_3d operator()() {
    return zero_functions<double, 3>(world, static_cast<int>(n_orbtials));
  }
};

using response_solver = XNonlinearSolver<vector_real_function_3d, double,
                                         response_vector_allocator>;

// template <typename ResponseType>
// class ResponseSolverPolicy; // forward declaration

// // Primary template left undefined; we’ll specialize it:
// template <typename ResponseType> class ResponseSolverPolicy {
//   // nothing here
// };


// alpha factor as a constexpr function (nice for inlining/constants)
constexpr double alpha_factor(const StaticRestrictedResponse&)     { return -4.0; }
constexpr double alpha_factor(const DynamicRestrictedResponse&)    { return -2.0; }
constexpr double alpha_factor(const StaticUnrestrictedResponse&)   { return -2.0; }
constexpr double alpha_factor(const DynamicUnrestrictedResponse&)  { return -2.0; }

// compute_density
real_function_3d compute_density(World&, const StaticRestrictedResponse&,
                                 const vector_real_function_3d&);

real_function_3d compute_density(World&, const DynamicRestrictedResponse&,
                                 const vector_real_function_3d&);

real_function_3d compute_density(World&, const StaticUnrestrictedResponse&,
                                 const vector_real_function_3d&);

real_function_3d compute_density(World&, const DynamicUnrestrictedResponse&,
                                 const vector_real_function_3d&);

// make_bsh_operators
std::vector<poperatorT> make_bsh_operators(World&, const ResponseManager&, double freq,
                                           const Tensor<double>& orbital_energies, int n,
                                           ResponseDebugLogger&, const StaticRestrictedResponse&);

std::vector<poperatorT> make_bsh_operators(World&, const ResponseManager&, double freq,
                                           const Tensor<double>& orbital_energies, int n,
                                           ResponseDebugLogger&, const DynamicRestrictedResponse&);

std::vector<poperatorT> make_bsh_operators(World&, const ResponseManager&, double freq,
                                           const Tensor<double>& orbital_energies, int n,
                                           ResponseDebugLogger&, const StaticUnrestrictedResponse&);

std::vector<poperatorT> make_bsh_operators(World&, const ResponseManager&, double freq,
                                           const Tensor<double>& orbital_energies, int n,
                                           ResponseDebugLogger&, const DynamicUnrestrictedResponse&);

// CoupledResponseEquations (one per variant)
vector_real_function_3d CoupledResponseEquations(World&, const GroundStateData&,
  const StaticRestrictedResponse&, const vector_real_function_3d& vp,
  const std::vector<poperatorT>& bsh_x, const ResponseManager&, ResponseDebugLogger&);

vector_real_function_3d CoupledResponseEquations(World&, const GroundStateData&,
  const DynamicRestrictedResponse&, const vector_real_function_3d& vp,
  const std::vector<poperatorT>& bsh_x, const ResponseManager&, ResponseDebugLogger&);

vector_real_function_3d CoupledResponseEquations(World&, const GroundStateData&,
  const StaticUnrestrictedResponse&, const vector_real_function_3d& vp,
  const std::vector<poperatorT>& bsh_x, const ResponseManager&, ResponseDebugLogger&);

vector_real_function_3d CoupledResponseEquations(World&, const GroundStateData&,
  const DynamicUnrestrictedResponse&, const vector_real_function_3d& vp,
  const std::vector<poperatorT>& bsh_x, const ResponseManager&, ResponseDebugLogger&);


// template <> class ResponseSolverPolicy<StaticRestrictedResponse> {
// public:
//   static constexpr double alpha_factor = -4.0;

//   static real_function_3d compute_density(World &world,
//                                           const StaticRestrictedResponse &rvec,
//                                           const vector_real_function_3d &phi0) {
//     auto &x = rvec.x_alpha;
//     auto xphi = mul(world, x, phi0, true);
//     return 2.0 * sum(world, xphi, true);
//   }

//   static vector_real_function_3d CoupledResponseEquations(
//       World &world, const GroundStateData &gs,
//       const StaticRestrictedResponse &vecs, const vector_real_function_3d &vp,
//       const std::vector<poperatorT> &bsh_x, const ResponseManager &rm,
//       ResponseDebugLogger &logger);

//   // Helper functions
//   static std::vector<poperatorT>
//   make_bsh_operators(World &world, const ResponseManager &rm, double freq,
//                      const Tensor<double> &orbital_energies, int n,
//                      ResponseDebugLogger &logger);
// };



// template <> class ResponseSolverPolicy<DynamicRestrictedResponse> {
// public:
//   static constexpr double alpha_factor = -2.0;
//   static std::vector<poperatorT>
//   make_bsh_operators(World &world, const ResponseManager &rm, double freq,
//                      const Tensor<double> &orbital_energies, int n,
//                      ResponseDebugLogger &logger);
//   static real_function_3d compute_density(World &world,
//                                           const DynamicRestrictedResponse &rvec,
//                                           const vector_real_function_3d &phi0) {
//     auto phi_phi = phi0;
//     phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
//     auto &x = rvec.flat;

//     auto xphi = mul(world, x, phi_phi, true);
//     return 2.0 * sum(world, xphi, true);
//   }

//   static vector_real_function_3d CoupledResponseEquations(
//       World &world, const GroundStateData &gs,
//       const DynamicRestrictedResponse &vecs, const vector_real_function_3d &vp,
//       const std::vector<poperatorT> &bsh_x, const ResponseManager &rm,
//       ResponseDebugLogger &logger);
// };

// template <> class ResponseSolverPolicy<StaticUnrestrictedResponse> {
// public:
//   static constexpr double alpha_factor = -2.0;
//   static std::vector<poperatorT>
//   make_bsh_operators(World &world, const ResponseManager &rm, double freq,
//                      const Tensor<double> &orbital_energies, int n,
//                      ResponseDebugLogger &logger);
//   static real_function_3d
//   compute_density(World &world, const StaticUnrestrictedResponse &rvec,
//                   const vector_real_function_3d &phi0) {
//     auto phi_phi = phi0;
//     phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
//     auto &x = rvec.flat;

//     auto xphi = mul(world, x, phi_phi, true);
//     return 2.0 * sum(world, xphi, true);
//   }
//   static vector_real_function_3d CoupledResponseEquations(
//       World &world, const GroundStateData &gs, const ResponseVector &vecs,
//       const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
//       const ResponseManager &rm, ResponseDebugLogger &logger);
// };

// template <> class ResponseSolverPolicy<DynamicUnrestrictedResponse> {
// public:
//   static constexpr double alpha_factor = -2.0;
//   static std::vector<poperatorT>
//   make_bsh_operators(World &world, const ResponseManager &rm, double freq,
//                      const Tensor<double> &orbital_energies, int n,
//                      ResponseDebugLogger &logger);
//   static real_function_3d
//   compute_density(World &world, const DynamicUnrestrictedResponse &rvec,
//                   const vector_real_function_3d &phi0) {
//     auto phi_phi = phi0;
//     phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
//     auto &x = rvec.flat;

//     auto xphi = mul(world, x, phi_phi, true);
//     return 2.0 * sum(world, xphi, true);
//   }
//   static vector_real_function_3d CoupledResponseEquations(
//       World &world, const GroundStateData &gs, const ResponseVector &vecs,
//       const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
//       const ResponseManager &rm, ResponseDebugLogger &logger);
// };
