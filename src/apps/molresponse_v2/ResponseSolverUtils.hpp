#pragma once
#include <madness/chem/SCF.h>
#include <madness/mra/nonlinsol.h>

#include <madness/external/nlohmann_json/json.hpp>

using json = nlohmann::json;

namespace ResponseSolverUtils {

using namespace madness;

inline void print_iteration_line(int iter, double residual, double deltaE, double xVp, double density_target, double x_residual_target) {
  // 1) Print header once
  static bool printed_header = false;
  if (!printed_header) {
    std::cout << std::setw(6) << "Iter" << " │ " << std::setw(12) << "Residual" << " │ " << std::setw(12) << "Δρ" << " │ " << std::setw(12) << "<x|Vp>" << " │ "
              << "Targets\n";
    std::cout << std::string(6, '-') << "─┼─" << std::string(12, '-') << "─┼─" << std::string(12, '-') << "─┼─" << std::string(12, '-') << "─┼─"
              << std::string(19, '-') << "\n";
    // printed_header = true;
  }

  // 2) Compose the target string
  std::ostringstream tgt;
  tgt << "||x||≤" << std::scientific << std::setprecision(3) << x_residual_target << ",||Δρ||≤" << std::scientific << std::setprecision(3) << density_target;

  // 3) Print the iteration line
  std::cout << std::setw(6) << iter << " │ " << std::setw(12) << std::scientific << std::setprecision(3) << residual << " │ " << std::setw(12)
            << std::scientific << std::setprecision(3) << deltaE << " │ " << std::setw(12) << std::scientific << std::setprecision(3) << xVp << " │ "
            << tgt.str() << "\n";
}

inline std::vector<poperatorT> make_bsh_operators_response(World &world, const double shift, const double omega, const Tensor<double> &ground_energies,
                                                           const double &lo) {
  double tol = FunctionDefaults<3>::get_thresh();
  // Sizes inferred from ground and omega
  size_t num_orbitals = ground_energies.size();  // number of orbitals
  std::vector<poperatorT> ops(num_orbitals);
  // Run over occupied components
  int p = 0;
  std::for_each(ops.begin(), ops.end(), [&](auto &operator_p) {
    double mu = sqrt(-2.0 * (ground_energies(p++) + omega + shift));
    operator_p = poperatorT(BSHOperatorPtr3D(world, mu, lo, tol));
  });
  return ops;
  // End timer
}

inline double inner(World &world, const vector_real_function_3d &x, const vector_real_function_3d &y) {
  double result = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    result += x[i].inner(y[i]);
  }
  return result;
}

inline void do_step_restriction(World &world, const vecfuncT &x, vecfuncT &x_new, const double &anorm, const std::string &spin, const double &maxrotn) {
  int nres = 0;
  if (anorm > maxrotn) {
    if (world.rank() == 0) print("Doing step restriction, norm of change: ", anorm);
    double s = maxrotn / anorm;
    gaxpy(s, x_new, 1.0 - s, x, true);
  }

  world.gop.fence();
  if (world.rank() == 0) print("Norm of vector changes", spin, ": ", anorm);
}
}  // namespace ResponseSolverUtils
