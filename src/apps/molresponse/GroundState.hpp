#ifndef MOLRESPONSE_V3_GROUNDSTATE_HPP
#define MOLRESPONSE_V3_GROUNDSTATE_HPP

#include <madness/chem/SCF.h>
#include <madness/chem/projector.h>
#include <madness/world/world.h>
#include <memory>
#include <string>

namespace molresponse_v3 {

using namespace madness;

/// Thin wrapper around SCF for response property calculations.
///
/// Delegates all ground-state data access to a shared SCF object.
/// Extends with response-specific cached operators (V_local, Fock matrices,
/// QProjector) that are rebuilt per protocol step via prepare().
///
/// Supports both closed-shell (restricted) and open-shell (unrestricted)
/// through SCF's native alpha/beta orbital storage. Fock matrices are
/// computed using SCF::apply_potential() and SCF::make_fock_matrix()
/// directly — no reimplementation of the Fock construction pipeline.
class GroundState {
public:
    /// Factory: load from a moldft restart archive.
    /// Reads the archive header to extract L, xc, etc., constructs an SCF
    /// with matching parameters, then delegates to SCF::load_mos().
    static GroundState from_archive(World& world,
                                     const std::string& archive_path,
                                     const Molecule& molecule);

    // -----------------------------------------------------------------
    // Delegated accessors — no local copies, all point to scf_->
    // -----------------------------------------------------------------

    const Molecule& molecule() const { return scf_->molecule; }

    /// Alpha orbitals (convenience alias for restricted case)
    const vecfuncT& orbitals() const { return scf_->get_amo(); }
    const vecfuncT& orbitals_alpha() const { return scf_->get_amo(); }
    const vecfuncT& orbitals_beta() const { return scf_->get_bmo(); }

    /// Orbital energies
    const tensorT& energies() const { return scf_->aeps; }
    const tensorT& energies_alpha() const { return scf_->aeps; }
    const tensorT& energies_beta() const { return scf_->beps; }

    /// Occupation numbers
    const tensorT& occupations_alpha() const { return scf_->get_aocc(); }
    const tensorT& occupations_beta() const { return scf_->get_bocc(); }

    /// Orbital counts
    long num_orbitals() const { return static_cast<long>(scf_->get_amo().size()); }
    long num_alpha() const { return static_cast<long>(scf_->get_amo().size()); }
    long num_beta() const { return static_cast<long>(scf_->get_bmo().size()); }

    bool is_spin_restricted() const { return scf_->is_spin_restricted(); }
    double L() const { return scf_->param.L(); }
    int k() const { return FunctionDefaults<3>::get_k(); }
    const CalculationParameters& params() const { return scf_->param; }
    double hf_exchange_coefficient() const { return scf_->xc.hf_exchange_coefficient(); }

    // -----------------------------------------------------------------
    // Response-specific cached operators (built by prepare())
    // -----------------------------------------------------------------

    /// Prepare ground-state operators for the current protocol.
    /// Must be called after FunctionDefaults<3> are set for the current
    /// protocol step. Reprojects orbitals if k changed, then builds
    /// V_local and Fock matrices using SCF::apply_potential() and
    /// SCF::make_fock_matrix().
    void prepare(World& world, double vtol,
                 const poperatorT& coulop,
                 const std::string& fock_json_file = "",
                 bool verbose = true);   // F2d-ii: false in subworlds (suppress print_info)

    /// Local potential: V_nuc + V_coul + V_xc (for response multiplicative term)
    const real_function_3d& V_local() const;

    /// Alpha Fock matrix (also the only Fock matrix for restricted)
    const tensorT& focka() const;
    /// Beta Fock matrix (only meaningful for unrestricted)
    const tensorT& fockb() const;

    /// Alpha Fock with diagonal zeroed (off-diagonal coupling)
    const tensorT& focka_no_diag() const;
    /// Beta Fock with diagonal zeroed
    const tensorT& fockb_no_diag() const;

    /// Convenience: alpha Fock (matches v2 interface name)
    const tensorT& hamiltonian() const { return focka(); }
    const tensorT& hamiltonian_no_diag() const { return focka_no_diag(); }

    /// Alpha projector onto virtual space: Q_alpha = 1 - |phi_alpha><phi_alpha|
    const QProjector<double, 3>& Q() const;
    const QProjector<double, 3>& Q_alpha() const;

    /// Beta projector (only for unrestricted): Q_beta = 1 - |phi_beta><phi_beta|
    const QProjector<double, 3>& Q_beta() const;

    // -----------------------------------------------------------------
    // Direct SCF access for advanced use
    // -----------------------------------------------------------------

    SCF& scf() { return *scf_; }
    const SCF& scf() const { return *scf_; }

    void print_info() const;

private:
    /// Construct a shell around an SCF (shared ownership). Private: GroundState
    /// is only built via from_archive, so prepare() can always reload pristine
    /// MOs from the checkpoint on each protocol climb — the restart-safe
    /// invariant the madqc and CLI paths both rely on.
    explicit GroundState(World& world, std::shared_ptr<SCF> scf);

    std::shared_ptr<SCF> scf_;
    int original_k_;
    int current_k_ = 0;  // k that orbitals are currently projected to

    // Cached response-specific data (rebuilt per protocol step)
    real_function_3d v_local_;
    tensorT focka_;
    tensorT fockb_;
    tensorT focka_no_diag_;
    tensorT fockb_no_diag_;
    QProjector<double, 3> q_alpha_;
    QProjector<double, 3> q_beta_;
    bool prepared_ = false;

public:
    /// Read just the archive header to extract L, k, xc, etc.
    /// Used by from_archive to set up CalculationParameters before load_mos,
    /// and by callers that need to set FunctionDefaults<3> via
    /// `set_response_protocol(...)` *before* the archive is opened so that
    /// orbitals load directly into the intended (k, thresh) representation.
    struct ArchiveHeader {
        unsigned int version = 0;
        double energy = 0.0;
        bool spin_restricted = true;
        double L = 0.0;
        int k = 0;
        std::string xc;
        std::string localize_method;
        double converged_for_thresh = 0.0;
        unsigned int nmo_alpha = 0;
    };

    static ArchiveHeader read_archive_header(World& world,
                                              const std::string& archive_path);

private:

    /// Build V_local (V_nuc + V_coul + V_xc) for the response solver.
    /// This is the multiplicative local potential used in the response
    /// equations, separate from the Fock matrix construction.
    void build_v_local(World& world, double vtol, const poperatorT& coulop);

    /// Build Fock matrices using SCF::apply_potential + SCF::make_fock_matrix.
    void build_fock_matrices(World& world, double vtol,
                              const std::string& fock_json_file);
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_GROUNDSTATE_HPP
