#pragma once
#include "GroundStateData.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseState.hpp"
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

/*struct ResponseSolver {*/
/*  virtual ~ResponseSolver() = default;*/
/**/
/*  bool iterate(World &world, const ResponseManager &rm,*/
/*               const GroundStateData &gs, const ResponseState &state,*/
/*               ResponseVector &response, ResponseDebugLogger &logger,*/
/*               size_t max_iter, double conv);*/
/*};*/

// I don't need to make a base class, just keep the same interface for each
// specilized solver
//
// Each solver class will be passed to a standard iterate function,
//
// Each solver implements a set of specialized functions for each ResponseVector
// variant
//
// Interface for each solver
// --
// -- compute_density
//

class StaticRestrictedSolver {
public:
  std::vector<poperatorT>
  make_bsh_operators(World &world, const ResponseManager &rm, double freq,
                     const Tensor<double> &orbital_energies, int n,
                     ResponseDebugLogger &logger);

  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, ResponseDebugLogger &logger,
                      size_t max_iter, double conv_thresh);

  static real_function_3d compute_density(World &world,
                                          const vector_real_function_3d &x,
                                          const vector_real_function_3d &phi0) {

    auto xphi = mul(world, x, phi0, true);
    return 2.0 * sum(world, xphi, true);
  }

  static vector_real_function_3d CoupledResponseEquations(
      World &world, const GroundStateData &gs, const ResponseVector &vecs,
      const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
      const ResponseManager &rm, ResponseDebugLogger &logger);
};

class DynamicRestrictedSolver {
public:
  static std::vector<poperatorT>
  make_bsh_operators(World &world, const ResponseManager &rm, double freq,
                     const Tensor<double> &orbital_energies, int n,
                     ResponseDebugLogger &logger);
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, ResponseDebugLogger &logger,
                      size_t max_iter, double conv_thresh);
  static real_function_3d compute_density(World &world,
                                          const vector_real_function_3d &x,
                                          const vector_real_function_3d &phi0) {
    auto phi_phi = phi0;
    phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());

    auto xphi = mul(world, x, phi_phi, true);
    return 2.0 * sum(world, xphi, true);
  }

  static vector_real_function_3d CoupledResponseEquations(
      World &world, const GroundStateData &gs, const ResponseVector &vecs,
      const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
      const ResponseManager &rm, ResponseDebugLogger &logger);
};

class StaticUnrestrictedSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, ResponseDebugLogger &logger,
                      size_t max_iter, double conv_thresh) {
    throw std::runtime_error("StaticUnrestrictedSolver not implemented");
  }
};
class DynamicUnrestrictedSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, ResponseDebugLogger &logger,
                      size_t max_iter, double conv_thresh) {
    throw std::runtime_error("DynamicUnrestrictedSolver not implemented");
  }
};
