#ifndef MADNESS_APPS_DIRAC_DFCONVERGENCE_H_INCLUDED
#define MADNESS_APPS_DIRAC_DFCONVERGENCE_H_INCLUDED

#include <madness/world/madness_exception.h>
#include <cmath>
#include <string>
#include <string_view>

namespace madness {

enum class DFConvergenceCriterion { bsh_residual, energy_density_residual };

/// Map an input keyword onto the criterion. Returns false for an unknown keyword
/// instead of throwing. MADNESS_EXCEPTION stores the \c const \c char* it is given
/// and does not copy it. A message built here from the keyword dangles after the
/// throw. The caller prints the keyword, then throws.
inline bool df_convergence_criterion_from_string(const std::string& keyword,
                                                 DFConvergenceCriterion& criterion) {
    if (keyword == "bsh_residual") {
        criterion = DFConvergenceCriterion::bsh_residual;
        return true;
    }
    if (keyword == "energy_density_residual") {
        criterion = DFConvergenceCriterion::energy_density_residual;
        return true;
    }
    return false;
}

inline const char* df_convergence_criterion_name(const DFConvergenceCriterion criterion) {
    switch (criterion) {
        case DFConvergenceCriterion::bsh_residual:
            return "bsh_residual";
        case DFConvergenceCriterion::energy_density_residual:
            return "energy_density_residual";
    }
    MADNESS_EXCEPTION("invalid Dirac-Fock convergence criterion", 0);
}

struct DFConvergenceMetrics {
    double current_total_energy;
    double previous_total_energy;
    double energy_tolerance;
    double density_residual;
    double density_tolerance;
    double max_bsh_residual;
    double bsh_tolerance;
};

inline bool df_iteration_converged(const DFConvergenceCriterion criterion,
                                   const DFConvergenceMetrics& m) {
    switch (criterion) {
        case DFConvergenceCriterion::bsh_residual:
            return m.max_bsh_residual <= m.bsh_tolerance;
        case DFConvergenceCriterion::energy_density_residual: {
            const double relative_energy_change =
                m.current_total_energy == 0.0
                    ? std::abs(m.current_total_energy - m.previous_total_energy)
                    : std::abs((m.current_total_energy - m.previous_total_energy) /
                               m.current_total_energy);
            return relative_energy_change <= m.energy_tolerance &&
                   m.density_residual <= m.density_tolerance &&
                   m.max_bsh_residual <= 100.0 * m.bsh_tolerance;
        }
    }
    MADNESS_EXCEPTION("invalid Dirac-Fock convergence criterion", 0);
}

/// Why the SCF loop stopped. The message text is separate from the print, so a
/// unit test can assert the exact wording.
enum class DFStopReason { converged_bsh, converged_combined, max_iterations };

inline std::string_view df_stop_message(const DFStopReason reason) {
    switch (reason) {
        case DFStopReason::converged_bsh:
            return "Converged due to residuals";
        case DFStopReason::converged_combined:
            return "Converged due to energy, density, and residuals";
        case DFStopReason::max_iterations:
            return "WARNING: maximum iterations reached without convergence";
    }
    MADNESS_EXCEPTION("invalid Dirac-Fock stop reason", 0);
}

} // namespace madness

#endif
