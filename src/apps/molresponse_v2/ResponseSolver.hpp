#pragma once
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include <madness/mra/nonlinsol.h>

struct response_vector_allocator {
  World &world;
  const size_t n_orbtials;
  response_vector_allocator(World &world, size_t n_orbtials)
      : world(world), n_orbtials(n_orbtials) {}
  // overloading the default constructor () operator
  vector_real_function_3d operator()() {
    // print("allocator called with ", int(n_orbtials), " orbitals");
    //  returning constructor of x_vector
    return zero_functions<double, 3>(world, n_orbtials);
  }
};

typedef XNonlinearSolver<vector_real_function_3d, double,
                         response_vector_allocator>
    response_solver;

struct ResponseSolver {
  virtual ~ResponseSolver() = default;

  bool iterate(World &world, const ResponseManager &rm,
               const GroundStateData &gs, const ResponseState &state,
               ResponseVector &response, size_t max_iter, double conv);
};

class StaticRestrictedSolver : public ResponseSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, size_t max_iter,
                      double conv_thresh);

  static real_function_3d compute_density(World &world,
                                          const vector_real_function_3d &x,
                                          const vector_real_function_3d &phi0) {

    auto xphi = mul(world, x, phi0, true);
    return 2.0 * sum(world, xphi, true);
  }

  static vector_real_function_3d
  ComputeRSH(World &world, const GroundStateData &gs,
             const ResponseVector &vecs, const vector_real_function_3d &vp,
             const std::vector<poperatorT> &bsh_x, const ResponseManager &rm);
};

class DynamicRestrictedSolver : public ResponseSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, size_t max_iter,
                      double conv_thresh);
  static real_function_3d compute_density(World &world,
                                          const vector_real_function_3d &x,
                                          const vector_real_function_3d &phi0) {
    auto xphi = mul(world, x, phi0, true);
    return 2.0 * sum(world, xphi, true);
  }

  static vector_real_function_3d
  ComputeRSH(World &world, const GroundStateData &gs,
             const ResponseVector &vecs, const vector_real_function_3d &vp,
             const std::vector<poperatorT> &bsh_x, const ResponseManager &rm);
};

class StaticUnrestrictedSolver : public ResponseSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, size_t max_iter,
                      double conv_thresh) {
    throw std::runtime_error("StaticUnrestrictedSolver not implemented");
  }
};
class DynamicUnrestrictedSolver : public ResponseSolver {
public:
  static bool iterate(World &world, const ResponseManager &rm,
                      const GroundStateData &gs, const ResponseState &state,
                      ResponseVector &response, size_t max_iter,
                      double conv_thresh) {
    throw std::runtime_error("DynamicUnrestrictedSolver not implemented");
  }
};
