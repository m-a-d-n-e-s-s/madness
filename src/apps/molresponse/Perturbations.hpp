#ifndef MOLRESPONSE_V3_PERTURBATIONS_HPP
#define MOLRESPONSE_V3_PERTURBATIONS_HPP

#include "GroundState.hpp"
#include <madness/chem/SCF.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/molecular_functors.h>
#include <madness/mra/mra.h>

namespace molresponse_v3 {

using namespace madness;

/// Stable identity for an FD perturbation source (doc 13). `description()`
/// is the canonical string used in archive filenames and as the
/// `<pert>` key in response_metadata.json — it must round-trip across runs
/// and contain only filename-safe characters.
///
/// This is identity only; the existing dipole_perturbation / magnetic /
/// nuclear free functions below still build the RHS. Callers dispatch on
/// `kind` and pass `axis`/`atom` to whichever builder applies.
struct Perturbation {
    enum class Kind { Dipole, NuclearDisplacement, Magnetic };

    Kind kind = Kind::Dipole;
    int  axis = 0;        // 0=x, 1=y, 2=z
    int  atom = -1;       // for NuclearDisplacement; <0 = all-atoms sentinel ("*")

    /// Canonical filename/metadata-key fragment:
    ///   dipole_x   magnetic_y   nuc_3_z   nuc_*_z
    /// A negative `atom` is the all-atoms sentinel — the planner emits one
    /// such record and the calc manager expands it to one per atom against
    /// the concrete molecule (mirrors the ES-root "*" sentinel).
    std::string description() const {
        static const char *ax[] = {"x", "y", "z"};
        const char *a = (axis >= 0 && axis <= 2) ? ax[axis] : "?";
        switch (kind) {
            case Kind::Dipole:               return std::string("dipole_") + a;
            case Kind::Magnetic:             return std::string("magnetic_") + a;
            case Kind::NuclearDisplacement: {
                const std::string atom_tag = (atom < 0) ? "*" : std::to_string(atom);
                return std::string("nuc_") + atom_tag + "_" + a;
            }
        }
        return "unknown";
    }

    bool is_all_atoms() const {
        return kind == Kind::NuclearDisplacement && atom < 0;
    }

    static Perturbation dipole(int axis)            { return {Kind::Dipole,   axis, -1}; }
    static Perturbation magnetic(int axis)          { return {Kind::Magnetic, axis, -1}; }
    static Perturbation nuclear(int atom, int axis) { return {Kind::NuclearDisplacement, axis, atom}; }
    /// All-atoms sentinel — calc manager expands per molecule.
    static Perturbation nuclear_all(int axis)       { return {Kind::NuclearDisplacement, axis, -1}; }
};

/// Construct the dipole perturbation right-hand side for one direction.
///
/// Returns Q * (mu_dir * phi_i) for each occupied orbital phi_i.
/// Uses MomentFunctor from SCF.h for the coordinate function.
///
/// @param world  MADNESS world
/// @param gs     ground state (provides orbitals and Q projector)
/// @param axis   0=x, 1=y, 2=z
/// @return       Q-projected perturbation vector (vecfuncT)
inline vector_real_function_3d
dipole_perturbation(World& world, const GroundState& gs, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);

    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());

    auto rhs = mul(world, mu, gs.orbitals_alpha());
    rhs = gs.Q()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the dipole perturbation for beta orbitals (open-shell).
inline vector_real_function_3d
dipole_perturbation_beta(World& world, const GroundState& gs, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    MADNESS_CHECK(!gs.is_spin_restricted());

    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());

    auto rhs = mul(world, mu, gs.orbitals_beta());
    rhs = gs.Q_beta()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the nuclear displacement perturbation right-hand side.
///
/// Uses DNuclear from SCFOperators.h which handles nuclear correlation
/// factors and core potentials correctly.
///
/// @param world  MADNESS world
/// @param gs     ground state
/// @param atom   atom index
/// @param axis   0=x, 1=y, 2=z
/// @return       Q-projected perturbation vector
inline vector_real_function_3d
nuclear_perturbation(World& world, const GroundState& gs, int atom, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    // FD source from the SAME operator the VBC/beta contraction uses:
    // madchem::MolecularDerivativeFunctor (the molecule's smoothed nuclear-
    // attraction derivative, with special_points AT the nucleus so the MRA
    // refines the sharp feature). Matches molresponse_v2 perturbation_vector
    // (nuclear). The old DNuclear(SCF/NCF) path was a DIFFERENT, under-resolved
    // operator -> source/operator mismatch (wrong Raman) + a peaked source that
    // blew the FD past the explosion guard.
    real_function_3d dV =
        real_factory_3d(world)
            .functor(real_functor_3d(
                new madchem::MolecularDerivativeFunctor(gs.molecule(), atom, axis)))
            .truncate_on_project()
            .truncate_mode(0);
    dV.get_impl()->set_truncate_mode(1);
    auto rhs = mul(world, dV, gs.orbitals_alpha());
    rhs = gs.Q()(rhs);
    truncate(world, rhs, FunctionDefaults<3>::get_thresh());
    return rhs;
}

/// Construct the magnetic (Lz) perturbation right-hand side.
///
/// Uses Lz from SCFOperators.h. Returns complex functions since
/// Lz = -i(x d/dy - y d/dx) maps real orbitals to complex.
///
/// NOTE: The Lz<T,NDIM> template in SCFOperators.h has a known issue
/// when T=double — its operator() returns real but internal code expects
/// complex. Must be instantiated as Lz<std::complex<double>, 3> with
/// complex input orbitals, or the SCFOperators template needs fixing.
/// For now, this is left as future work for magnetic response.
///
// inline std::vector<Function<std::complex<double>, 3>>
// magnetic_perturbation_lz(World& world, const GroundState& gs);

/// Construct the raw dipole operator function (without applying to orbitals).
/// Useful for property computation: alpha_ij = <x_i | mu | phi_j>
inline real_function_3d
dipole_operator(World& world, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    std::vector<int> dir(3, 0);
    dir[axis] = 1;
    real_function_3d mu = real_factory_3d(world).functor(
        real_functor_3d(new MomentFunctor(dir)));
    mu.truncate(FunctionDefaults<3>::get_thresh());
    return mu;
}

/// Construct all three dipole operators (x, y, z).
inline std::array<real_function_3d, 3>
dipole_operators(World& world) {
    return {dipole_operator(world, 0),
            dipole_operator(world, 1),
            dipole_operator(world, 2)};
}

/// Raw nuclear-displacement operator dV_nuc/dQ(atom, axis) as a real_function_3d,
/// for the VBC/beta contraction (the nuclear analog of dipole_operator; used as
/// VC_op in Raman beta(dipole; dipole, nuclear)). Matches molresponse_v2
/// make_perturbation_operator(NuclearDisplacement): MolecularDerivativeFunctor on
/// gs.molecule(), truncate-on-project then truncate_mode 1.
inline real_function_3d
nuclear_operator(World& world, const GroundState& gs, int atom, int axis) {
    MADNESS_CHECK(axis >= 0 && axis <= 2);
    real_function_3d dvdx =
        real_factory_3d(world)
            .functor(real_functor_3d(
                new madchem::MolecularDerivativeFunctor(gs.molecule(), atom, axis)))
            .truncate_on_project()
            .truncate_mode(0);
    dvdx.get_impl()->set_truncate_mode(1);
    return dvdx;
}

/// Dispatch the raw perturbation operator for the VBC/beta contraction by kind:
/// dipole -> r (dipole_operator); nuclear -> dV_nuc/dQ (nuclear_operator). This
/// is what lets a VBC pair mix perturbation types (e.g. Raman = dipole x nuclear).
inline real_function_3d
perturbation_operator(World& world, const GroundState& gs, const Perturbation& p) {
    switch (p.kind) {
        case Perturbation::Kind::Dipole:
            return dipole_operator(world, p.axis);
        case Perturbation::Kind::NuclearDisplacement:
            return nuclear_operator(world, gs, p.atom, p.axis);
        case Perturbation::Kind::Magnetic:
            throw std::runtime_error(
                "perturbation_operator: magnetic operator not implemented");
    }
    throw std::runtime_error("perturbation_operator: unknown perturbation kind");
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_PERTURBATIONS_HPP
