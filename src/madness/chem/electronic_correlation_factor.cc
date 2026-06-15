// electronic_correlation_factor.cc — diagnostic machinery of CorrelationFactor.
//
// The production methods (apply_U_local/semilocal/mixed, apply_U, f(), …) are
// header-only in electronic_correlation_factor.h; the diagnostic bodies live
// here to keep the widely-included header lean.  See operator_diagnostics.md
// for the derivation of the 3D-only references (Schwinger fit of G on the
// bra, Q12 expanded on the ket as scalar contractions).

#include <madness/chem/electronic_correlation_factor.h>

namespace madness {

// ---------------------------------------------------------------------------
//  apply_U diagnostics block
// ---------------------------------------------------------------------------

void CorrelationFactor::run_apply_U_diagnostics(
        const real_function_6d& local,
        const real_function_6d& semilocal,
        const real_function_6d& mixed,
        const real_function_3d& phi_i, const real_function_3d& phi_j,
        const real_convolution_6d& op_mod,
        const std::vector<real_function_3d>& U1nuc,
        const std::vector<real_function_3d>& ao,
        const std::vector<real_function_3d>& occ_ket,
        const std::vector<real_function_3d>& occ_bra,
        const real_function_3d& phi_i_bra,
        const real_function_3d& phi_j_bra) const
{
    std::cout << std::scientific << std::setprecision(8);
    const double mu     = op_mod.mu();
    const double energy = -0.5 * mu * mu;
    real_convolution_6d G = BSHOperator<6>(world, mu, lo, 1.e-6);

    // Ue diagnostics (occ/energy ignored by diagnose_Ue)
    auto dm = diagnose_Ue(local, semilocal, phi_i, phi_j, ao,
                          occ_ket, occ_bra, energy, mixed, U1nuc);
    print("diagnosis for Ue term:");
    dm.print_report("Ue");

    // GUe diagnostics: apply G to each Ue piece (occ ignored by diagnose_GUe)
    real_function_6d Glocal     = apply(G, copy(local)).truncate();
    real_function_6d Gsemilocal = apply(G, copy(semilocal)).truncate();
    real_function_6d Gmixed     = apply(G, copy(mixed)).truncate();
    auto dm1 = diagnose_GUe(Glocal, Gsemilocal, phi_i, phi_j, ao,
                            occ_ket, occ_bra, energy, Gmixed, U1nuc);
    print("diagnosis for GUe term:");
    dm1.print_report("GUe");

    // GQUe diagnostics: Q12 from occ_bra/occ_ket
    if (occ_ket.size() > 0) {
        MADNESS_CHECK_THROW(occ_bra.size()>0,"no bra in run_apply_U_diagnostics");

        StrongOrthogonalityProjector<double,3> Q12(world);
        Q12.set_spaces(occ_bra, occ_ket, occ_bra, occ_ket);

        real_function_6d GQlocal     = apply(G, Q12(local)).truncate();
        real_function_6d GQsemilocal = apply(G, Q12(semilocal)).truncate();
        real_function_6d GQmixed     = apply(G, Q12(mixed)).truncate();

        auto dm2 = diagnose_GQUe(GQlocal, GQsemilocal, phi_i, phi_j, ao,
                                 occ_ket, occ_bra, energy, GQmixed, U1nuc);
        print("diagnosis for GQUe term:");
        dm2.print_report("GQUe");

        // MP2 pair-energy diagnostics: raw pair-bra basis {i_bra, j_bra};
        // the (0,1) element is <i_bra j_bra|G Q12 Ue_X|ij>
        if (phi_i_bra.is_initialized() && phi_j_bra.is_initialized()) {
            std::vector<real_function_3d> pair_bra = {phi_i_bra, phi_j_bra};
            auto dm3 = diagnose_GUe(Glocal, Gsemilocal, phi_i, phi_j, pair_bra,
                                     occ_ket, occ_bra, energy, Gmixed, U1nuc,
                                     /*orthonormalize_basis=*/false);
            print("diagnosis for GUe term in the raw pair-bra basis:");
            dm3.print_report("GUe-pair");

            print_pair_energy_report(dm3, {"Glocal", "Gsemilocal", "Gmixed"},
                                     {1.0, 1.0, -1.0},
                                     "MP2 energy contribution <i_bra j_bra|G Ue|ij>"
                                     " (Ue = local + semilocal - mixed):");

            auto dm4 = diagnose_GQUe(GQlocal, GQsemilocal, phi_i, phi_j, pair_bra,
                                     occ_ket, occ_bra, energy, GQmixed, U1nuc,
                                     /*orthonormalize_basis=*/false);
            print("diagnosis for GQUe term in the raw pair-bra basis:");
            dm4.print_report("GQUe-pair");

            print_pair_energy_report(dm4, {"GQlocal", "GQsemilocal", "GQmixed"},
                                     {1.0, 1.0, -1.0},
                                     "MP2 energy contribution <i_bra j_bra|G Q12 Ue|ij>"
                                     " (Ue = local + semilocal - mixed):");

            // RI MP2 pair-energy of the Ue term: contract <ab|G Q12 Ue|ij> with
            // the strong-orthogonality-projected g12 weights W^Q on the SAME
            // orthonormal observer basis dm2 uses.  Expr3 (6d ket, .result) vs
            // Expr2 (3d Schwinger ket, .ref).  Ue enters g~ with +, so the
            // constant-part contribution carries factor -2*facE.
            const double thresh3 = FunctionDefaults<3>::get_thresh();
            std::shared_ptr<SeparatedConvolution<double,3>> g12coul(CoulombOperatorPtr(world, lo, thresh3));
            const bool symmetric = ((phi_i - phi_j).norm2() < 1.e-8);
            const double facE = symmetric ? 1.0 : 2.0;
            Tensor<double> Wq = pair_energy_weights<double,6>(world, dm2.ao_basis,
                    phi_i_bra, phi_j_bra, occ_ket, occ_bra, *g12coul);
            print_pair_energy_report_RI(dm2, {"GQlocal","GQsemilocal","GQmixed"}, {1.0,1.0,-1.0},
                    Wq, -2.0*facE,
                    "RI MP2 pair-energy contribution of Ue, Expr2 3d-ket vs Expr3 6d-ket:");
        }
    }
}

// ---------------------------------------------------------------------------
//  diagnose_GQUe
// ---------------------------------------------------------------------------

DiagnosticMatrix<> CorrelationFactor::diagnose_GQUe(
        const real_function_6d& GQ12Uphi_local,
        const real_function_6d& GQ12Uphi_semilocal,
        const real_function_3d& phi_i,
        const real_function_3d& phi_j,
        const std::vector<real_function_3d>& aobasis,
        const std::vector<real_function_3d>& occ_ket,
        const std::vector<real_function_3d>& occ_bra,
        const double energy,
        const real_function_6d& GQ12Uphi_mixed,
        const std::vector<real_function_3d>& U1nuc,
        const bool orthonormalize_basis) const
{
    MADNESS_CHECK_THROW(energy < 0.0, "diagnose_GQUe: energy must be negative");
    const bool do_mixed = GQ12Uphi_mixed.is_initialized() && !U1nuc.empty();
    const double wall0 = wall_time();
    timer t(world);

    DiagnosticMatrix<> dm(world, aobasis, orthonormalize_basis);
    dm.init("GQlocal");
    dm.init("GQsemilocal");
    if (do_mixed) dm.init("GQmixed");

    // Result: project simple <ab| onto the caller-supplied G Q12 Ue|ij>
    dm.entries["GQlocal"].result     = dm.project_ab(GQ12Uphi_local);
    dm.entries["GQsemilocal"].result = dm.project_ab(GQ12Uphi_semilocal);
    if (do_mixed) dm.entries["GQmixed"].result = dm.project_ab(GQ12Uphi_mixed);
    t.tag("compute 6D projections");

    // Ref: Schwinger fit of G on the bra + Q12 ket-expansion (ref_GQab)
    dm.entries["GQlocal"].ref = dm.ref_GQab(
        ue_local_provider(phi_i, phi_j), occ_ket, occ_bra, energy, lo);
    dm.entries["GQsemilocal"].ref = dm.ref_GQab(
        ue_semilocal_provider(phi_i, phi_j), occ_ket, occ_bra, energy, lo);
    if (do_mixed)
        dm.entries["GQmixed"].ref = dm.ref_GQab(
            ue_mixed_provider(phi_i, phi_j, U1nuc), occ_ket, occ_bra, energy, lo);
    t.tag("3D ref with Q12 ket expansion");

    dm.compute_errors();
    dm.time = wall_time() - wall0;
    return dm;
}

// ---------------------------------------------------------------------------
//  diagnose_Ue_impl — unified core for diagnose_Ue and diagnose_GUe
// ---------------------------------------------------------------------------

DiagnosticMatrix<> CorrelationFactor::diagnose_Ue_impl(
        DiagnosticMatrix<> dm,
        const std::vector<CCPairFunction<double,6>>& bra,
        const real_function_6d& Uphi_local,
        const real_function_6d& Uphi_semilocal,
        const real_function_3d& phi_i, const real_function_3d& phi_j,
        const std::string& name_local, const std::string& name_semilocal,
        const real_function_6d& Uphi_mixed,
        const std::vector<real_function_3d>& U1nuc,
        const std::string& name_mixed) const
{
    const bool do_mixed = Uphi_mixed.is_initialized() && !U1nuc.empty();
    const double wall0  = wall_time();
    timer t(world);

    dm.init(name_local);
    dm.init(name_semilocal);
    if (do_mixed) dm.init(name_mixed);

    // Result: project B-applied 6D kets with simple bra <ab|
    dm.entries[name_local].result     = dm.project_ab(Uphi_local);
    dm.entries[name_semilocal].result = dm.project_ab(Uphi_semilocal);
    if (do_mixed) dm.entries[name_mixed].result = dm.project_ab(Uphi_mixed);
    t.tag("compute 6D projections");

    {
        auto local_el = ue_local_provider(phi_i, phi_j);
        for (const auto& bk : bra)
            dm.entries[name_local].ref += local_el(bk.get_a(), bk.get_b());
    }
    t.tag("3D ref for local");

    {
        auto semilocal_el = ue_semilocal_provider(phi_i, phi_j);
        for (const auto& bk : bra)
            dm.entries[name_semilocal].ref += semilocal_el(bk.get_a(), bk.get_b());
    }
    t.tag("3D ref for semilocal");

    if (do_mixed) {
        auto mixed_el = ue_mixed_provider(phi_i, phi_j, U1nuc);
        for (const auto& bk : bra)
            dm.entries[name_mixed].ref += mixed_el(bk.get_a(), bk.get_b());
        t.tag("3D ref for mixed");
    }

    dm.compute_errors();
    dm.time = wall_time() - wall0;
    return dm;
}

// ---------------------------------------------------------------------------
//  Ue element providers — the operator-specific unit plugged into
//  DiagnosticMatrix::ref_Gab / ref_GQab
// ---------------------------------------------------------------------------

DiagnosticMatrix<>::ElementProvider CorrelationFactor::ue_local_provider(
        const real_function_3d& phi_i, const real_function_3d& phi_j) const
{
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto fg_op     = std::make_shared<SeparatedConvolution<double,3>>(
            world, OperatorInfo(_gamma, dcut, thresh, OT_FG12));
    auto slater_op = std::make_shared<SeparatedConvolution<double,3>>(
            world, OperatorInfo(_gamma, dcut, thresh, OT_SLATER));
    const double gamma = _gamma;

    return [fg_op, slater_op, gamma, phi_i, phi_j](
            const std::vector<real_function_3d>& p1,
            const std::vector<real_function_3d>& p2) {
        Tensor<double> M(static_cast<long>(p1.size()), static_cast<long>(p2.size()));
        for (size_t x = 0; x < p1.size(); ++x) {
            real_function_3d ket_x = 2.0*gamma*(*fg_op)(p1[x]*phi_i)
                                   + 0.5*gamma*(*slater_op)(p1[x]*phi_i);
            for (size_t y = 0; y < p2.size(); ++y)
                M(x,y) = inner(phi_j * p2[y], ket_x);
        }
        return M;
    };
}

DiagnosticMatrix<>::ElementProvider CorrelationFactor::ue_semilocal_provider(
        const real_function_3d& phi_i, const real_function_3d& phi_j) const
{
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto bsh_op = std::make_shared<SeparatedConvolution<double,3>>(
            world, OperatorInfo(_gamma, dcut, thresh, OT_BSH));

    std::vector<real_function_3d> r(3);
    for (int ax = 0; ax < 3; ++ax)
        r[ax] = real_factory_3d(world).functor(
            [ax](const coord_3d& xyz){ return xyz[ax]; });
    auto dphi_i = grad(phi_i);
    auto dphi_j = grad(phi_j);
    auto rphi_i = phi_i * r;
    auto rphi_j = phi_j * r;
    real_function_3d itilde = dot(world, r, dphi_i);
    real_function_3d jtilde = dot(world, r, dphi_j);

    return [bsh_op, phi_i, phi_j, dphi_i, dphi_j, rphi_i, rphi_j, itilde, jtilde](
            const std::vector<real_function_3d>& p1,
            const std::vector<real_function_3d>& p2) {
        Tensor<double> M(static_cast<long>(p1.size()), static_cast<long>(p2.size()));
        for (size_t x = 0; x < p1.size(); ++x) {
            auto term1 = (*bsh_op)(p1[x] * itilde);   // A
            auto term2 = (*bsh_op)(p1[x] * rphi_i);   // B
            auto term3 = (*bsh_op)(p1[x] * dphi_i);   // C
            auto term4 = (*bsh_op)(p1[x] * phi_i);    // D
            for (size_t y = 0; y < p2.size(); ++y) {
                double tmp = 0.0;
                tmp += inner(p2[y] * phi_j,  term1);   // +A
                tmp -= inner(p2[y] * dphi_j, term2);   // -B
                tmp -= inner(p2[y] * rphi_j, term3);   // -C
                tmp += inner(p2[y] * jtilde, term4);   // +D
                M(x,y) = -2.0 * constants::pi * tmp;
            }
        }
        return M;
    };
}

DiagnosticMatrix<>::ElementProvider CorrelationFactor::ue_mixed_provider(
        const real_function_3d& phi_i, const real_function_3d& phi_j,
        const std::vector<real_function_3d>& U1nuc) const
{
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto bsh_op = std::make_shared<SeparatedConvolution<double,3>>(
            world, OperatorInfo(_gamma, dcut, thresh, OT_BSH));

    std::vector<real_function_3d> r(3);
    for (int ax = 0; ax < 3; ++ax)
        r[ax] = real_factory_3d(world).functor(
            [ax](const coord_3d& xyz){ return xyz[ax]; });
    auto Wphii  = U1nuc * phi_i;
    auto Wphij  = U1nuc * phi_j;
    auto rphi_i = phi_i * r;
    auto rphi_j = phi_j * r;
    real_function_3d itilde_nuc = dot(world, r, Wphii);
    real_function_3d jtilde_nuc = dot(world, r, Wphij);

    return [bsh_op, phi_i, phi_j, Wphii, Wphij, rphi_i, rphi_j, itilde_nuc, jtilde_nuc](
            const std::vector<real_function_3d>& p1,
            const std::vector<real_function_3d>& p2) {
        Tensor<double> M(static_cast<long>(p1.size()), static_cast<long>(p2.size()));
        for (size_t x = 0; x < p1.size(); ++x) {
            auto term1 = (*bsh_op)(p1[x] * itilde_nuc);
            auto term2 = (*bsh_op)(p1[x] * rphi_i);
            auto term3 = (*bsh_op)(p1[x] * Wphii);
            auto term4 = (*bsh_op)(p1[x] * phi_i);
            for (size_t y = 0; y < p2.size(); ++y) {
                double tmp = 0.0;
                tmp += inner(p2[y] * phi_j,      term1);
                tmp -= inner(p2[y] * Wphij,      term2);
                tmp -= inner(p2[y] * rphi_j,     term3);
                tmp += inner(p2[y] * jtilde_nuc, term4);
                M(x,y) = -2.0 * constants::pi * tmp;
            }
        }
        return M;
    };
}

} // namespace madness
