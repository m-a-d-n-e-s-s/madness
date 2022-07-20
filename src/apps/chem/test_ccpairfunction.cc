//
// Created by Florian Bischoff on 6/27/22.
//


#include<madness/mra/mra.h>
#include<chem/ccpairfunction.h>
#include<chem/correlationfactor.h>
#include<chem/electronic_correlation_factor.h>
#include<chem/CCStructures.h>

#include<madness/world/test_utilities.h>

using namespace madness;


int test_overlap(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                  const CCParameters& parameters) {
    CCTimer timer(world, "testing");
    auto R2 = ncf->square();

    CCConvolutionOperator g12(world, OT_G12, parameters);
    CCConvolutionOperator f12(world, OT_F12, parameters);

    CorrelationFactor corrfac(world, 1.0, 1.e-7, molecule);

    auto g = [](const coord_3d& r) { return exp(-1.0 * r.normf()); };
    real_function_3d f1 = real_factory_3d(world).f(g);
    real_function_3d R2f = R2 * f1;
    double norm = inner(f1, R2f);
    f1.scale(1.0 / sqrt(norm));
    R2f.scale(1.0 / sqrt(norm));
    CC_vecfunction mo_ket_(std::vector<real_function_3d>(1, f1), HOLE);
    CC_vecfunction mo_bra_(std::vector<real_function_3d>(1, R2f), HOLE);

    bool passed_lo = true;
    bool passed_hi = true;
    bool success=true;
    int isuccess=0;
    const double lo = parameters.thresh_6D();
    const double hi = parameters.thresh_3D();
    const double hi_loose = parameters.thresh_3D()*5.0;

    print("threshold 3D / 6D", hi,lo);
    test_output t1("CCPairFunction::inner");

    // f^2 = 1/(4y^2)(1 - 2*f2(y) + f2(2y)) , f2(2y) =f2(y)^2
    const double y = parameters.gamma();
    SeparatedConvolution<double, 3> f = SlaterF12Operator(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, 3> f2 = SlaterOperator(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, 3> ff = SlaterOperator(world, 2.0 * y, parameters.lo(), parameters.thresh_bsh_3D());

    real_function_3d a = mo_ket_(0).function;
    real_function_3d b = mo_ket_(0).function;
    real_function_6d fab_6d = CompositeFactory<double, 6, 3>(world).g12(corrfac.f()).particle1(copy(a)).particle2(
            copy(b));
    fab_6d.fill_tree().truncate().reduce_rank();
    fab_6d.print_size("fab_6d");

    real_function_3d aR = mo_bra_(0).function;
    real_function_3d bR = mo_bra_(0).function;

    const real_function_3d aa = (aR * a).truncate();
    const real_function_3d bb = (bR * b).truncate();
    const real_function_3d af2a = f2(aa);
    const real_function_3d affa = ff(aa);
    const real_function_3d afa = f(aa);

    const double prefactor = 1.0 / (4 * y * y);
    const double ab_f2_ab = prefactor * (aR.inner(a) * bR.inner(b) - 2.0 * bb.inner(af2a) + bb.inner(affa));
    const double ab_f_ab = inner(f(aa), bb);

    const long rank = world.rank();
    auto printer = [&rank](std::string msg, const double result, const double reference, const double time) {
        if (rank == 0) {
            long len = std::max(0l, long(40 - msg.size()));
            msg += std::string(len, ' ');
            std:
            cout << msg << std::fixed << std::setprecision(8) << "result, ref, diff "
                 << result << " " << reference << " " << std::setw(9) << std::setprecision(2) << std::scientific
                 << result - reference
                 << ", elapsed time " << std::fixed << std::setw(5) << std::setprecision(2) << time << std::endl;
        }
    };
    t1.checkpoint(true, "initialization");

    double time_init = timer.reset();
    if (world.rank() == 0) print("time spent in initializing ", time_init);
    {
        CCPairFunction fab(world, &f12, a, b);
        const double test1 = inner(fab, fab, R2);
        const double diff = test1 - ab_f2_ab;
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;

        printer(" <a b |f f | a b> op_dec/op_dec  : ", test1, ab_f2_ab, timer.reset());
        success=(fabs(diff) < hi);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "op_dec/op_dec");
    {
        CCPairFunction ab(world, mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
        const double test1 = inner(ab, ab, R2);
        const double test2 = double(mo_ket_.size()); // mos are normed
        const double diff = test1 - test2;
        printer(" <a b | a b> dec/dec  : ", test1, double(mo_ket_.size()), timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;
        success=(fabs(diff) < hi);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "dec/dec");
    {
        CCPairFunction fab(world, fab_6d);
        const double test1 = inner(fab, fab, R2);
        const double diff = test1 - ab_f2_ab;
        printer(" <abf | abf> pure/pure : ", test1, ab_f2_ab, timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;
        success=(fabs(diff) < hi_loose);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "pure/pure");

    {
        CCPairFunction ab1(world, mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
        CCPairFunction fab(world, &f12, a, b);
        timer.reset();
        const double test1 = inner(ab1, fab, R2);
        const double test2 = ab_f_ab;
        const double diff = test1 - test2;
        printer(" <a b | fab > dec/op_dec :  ", test1, ab_f_ab, timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;
        success=(fabs(diff) < hi);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "dec/op_dec");

    {
        // the next tests evaulate <ab|f|ab> in different ways
        CCPairFunction fab(world, fab_6d);
        CCPairFunction ab2(world, mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
        const double test1 = inner(fab, ab2, R2);
        const double test2 = bb.inner(afa);
        const double diff = test1 - test2;
        printer(" <a b | fab > dec/pure", test1, test2, timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;

        success=(fabs(diff) < hi_loose);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "dec/pure");
//        {
//            CCPairFunction fab(world, fab_6d);
//            CCPairFunction ab2(world, ab_6d);
//            const double test1 = overlap(fab, ab2);
//            const double test2 = bb.inner(afa);
//            const double diff = test1 - test2;
//            if (world.rank() == 0)
//                std::cout << "Overlap Test 6 : " << std::fixed << std::setprecision(10) << "result=" << test1
//                          << ", test=" << test2 << ", diff=" << diff << "\n";
//            if (fabs(diff) > lo) passed_lo = false;
//            if (fabs(diff) > hi) passed_hi = false;
//        }
    {
        CCPairFunction fab(world, &f12, a, b);
        CCPairFunction ab2(world, fab_6d);
        const double test1 = inner(fab, ab2, R2);
        const double test2 = ab_f2_ab;
        const double diff = test1 - test2;
        printer(" <a b f | fab > op_dec/pure  : ", test1, test2, timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;
        success=(fabs(diff) < hi_loose);                      // be a bit loose here ..
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "op_dec/pure");

    {
        CCPairFunction fab(world, &f12, a, b);
        CCPairFunction ab2(world, mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
        timer.reset();
        const double test1 = inner(fab, ab2, R2);
//            const double test2 = bb.inner(afa);
        const double test2 = ab_f_ab;
        const double diff = test1 - test2;
        printer(" <a b f| a b > op_dec/dec : ", test1, test2, timer.reset());
        if (fabs(diff) > lo) passed_lo = false;
        if (fabs(diff) > hi) passed_hi = false;
        success=(fabs(diff) < hi);
    }
    if (not success) isuccess++;
    t1.checkpoint(success, "op_dec/dec");

    return isuccess;
}

int test_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;

    return success;
}

int test_scale(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;

    return success;
}

int test_swap_particles(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                        const CCParameters& parameters) {
    int success = 0;
    test_output t1("CCPairFunction::swap_particles");

    CCTimer timer(world, "testing swap_particles");

    CCConvolutionOperator f12(world, OT_F12, parameters);

    auto one = [](const coord_3d& r) { return 1.0; };
    real_function_3d R2 = real_factory_3d(world).f(one);

    // prepare
    auto g1 = [](const coord_3d& r) { return exp(-1.0 * inner(r, r)); };
    auto g2 = [](const coord_3d& r) { return exp(-2.0 * inner(r, r)); };
    auto g3 = [](const coord_3d& r) { return exp(-3.0 * inner(r, r)); };
    real_function_3d f1 = real_factory_3d(world).f(g1);
    real_function_3d f2 = real_factory_3d(world).f(g2);
    real_function_3d f3 = real_factory_3d(world).f(g3);
    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f3, f1};

    // test decomposed
    {
        CCPairFunction p1(world, a, b);
        CCPairFunction p2(world, b, a);

        double norm1 = inner(p1, p2.swap_particles(), R2);
        double norm1a = inner(p1, p1, R2);
        // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
        double norm2 = matrix_inner(world, a, a).emul(matrix_inner(world, b, b)).sum();
        print("norm1 ", norm1);
        print("norm1a", norm1a);
        print("norm2 ", norm2);
        t1.checkpoint(std::abs(norm1 - norm2) < FunctionDefaults<3>::get_thresh(), "swap_particles a,b");
    }

    // test pure
    {
        auto g = [](const coord_6d& r) {
            double r1=r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            double r2=r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
            return exp(-1.0*r1 - 2.0*r2);
        };
        real_function_6d f = real_factory_6d(world).f(g);
        CCPairFunction p(world, f);
        CCPairFunction p_swapped=p.swap_particles();
        double pnorm=p.get_function().norm2();
        double psnorm=p_swapped.get_function().norm2();
        print("p/s norm",pnorm,psnorm);

        CCPairFunction p1(world, {f1}, {f2});
        CCPairFunction p2(world, {f2}, {f1});
        double ref1=inner(f1,f1)*inner(f2,f2);
        double ref2=inner(f1,f2)*inner(f2,f1);
        print("ref1/2",ref1,ref2);
        print("pref1/2",inner(p1,p1),inner(p1,p2));

        double norm12_12=inner(p,p1);
        double norm12_21=inner(p,p1.swap_particles());
        double norm12_12_again=inner(p,p1);
        double norm21_12=inner(p_swapped,p1);
        double norm21_21=inner(p_swapped,p1.swap_particles());

        print("norms in exp(-12 - 12):",norm12_12);
        print("norms in exp(-12 - 21):",norm12_21);
        print("norms in exp(-12 - 12) again:",norm12_12_again);
        print("norms in exp(-21 - 12):",norm21_12);
        print("norms in exp(-21 - 21):",norm21_21);

        double ref_12_12=inner(p1,p1);
        double ref_12_21=inner(p1,p2);

        print("difference norms in exp(-12 - 12):",norm12_12,ref_12_12);
        print("difference norms in exp(-12 - 21):",norm12_21,ref_12_21);
        print("difference norms in exp(-21 - 12):",norm21_12,ref_12_21);
        print("difference norms in exp(-21 - 21):",norm21_21,ref_12_12);

        double total_error= fabs(norm12_12-ref_12_12)+ fabs(norm12_21-ref_12_21)
                + fabs(norm21_12-ref_12_21)+ fabs(norm21_21-ref_12_12);


        t1.checkpoint(total_error < FunctionDefaults<3>::get_thresh(), "swap_particles u");
    };

    return success;
}

int test_dirac_convolution(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                           const CCParameters& parameter) {
    int success=0;

    return success;
}

int test_partial_inner(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                       const CCParameters& parameters) {
    int success=0;

    return success;
}

/** functionality
 *
 *  - ctor
 *  - assignment
 *  - add
 *  - scalar multiplication
 *  - inner
 *  - inner_partial
 *  - swap_particles
 *  - apply
 *  - apply_partial (i.e. exchange)
 *  - serialize
 *  - callapse_to_pure (excl g!)
 *  - mul_partial
 */

int main(int argc, char **argv) {


    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    int isuccess=0;
#ifdef USE_GENTENSOR
    {
        parser.set_keyval("geometry", "source=library,he");
        parser.print_map();
        Molecule mol(world, parser);
        mol.print();
        CCParameters ccparam(world, parser);

        std::shared_ptr<NuclearCorrelationFactor> ncf = create_nuclear_correlation_factor(world,
                         mol, nullptr, std::make_pair("slater", 2.0));

        isuccess+=test_overlap(world, ncf, mol, ccparam);
        isuccess+=test_swap_particles(world, ncf, mol, ccparam);
    }
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}