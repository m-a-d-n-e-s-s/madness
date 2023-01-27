//
// Created by Florian Bischoff on 6/27/22.
//


#include<madness/mra/mra.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/projector.h>

#include<madness/world/test_utilities.h>

using namespace madness;


struct data {
    real_function_3d f1,f2,f3,f4,f5;
    real_function_6d f;

    std::shared_ptr<CCConvolutionOperator> f12;
    data() {}

    data(World& world, CCParameters& parameters) {
        f12.reset(new CCConvolutionOperator(world,OT_F12,parameters));
        if (not f1.is_initialized()) initialize(world);
    }

    void initialize(World& world) {
        auto g1 = [](const coord_3d& r) { return exp(-1.0 * inner(r, r)); };
        auto g2 = [](const coord_3d& r) { return exp(-2.0 * inner(r, r)); };
        auto g3 = [](const coord_3d& r) { return exp(-3.0 * inner(r, r)); };
        auto g4 = [](const coord_3d& r) { return exp(-4.0 * inner(r, r)); };
        auto g5 = [](const coord_3d& r) { return exp(-5.0 * inner(r, r)); };
        f1=real_factory_3d(world).f(g1);
        f2=real_factory_3d(world).f(g2);
        f3=real_factory_3d(world).f(g3);
        f4=real_factory_3d(world).f(g4);
        f5=real_factory_3d(world).f(g5);

        auto g = [](const coord_6d& r) {
            double r1=r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            double r2=r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
            return exp(-1.0*r1 - 2.0*r2);
        };
        f = real_factory_6d(world).f(g);

    }
    void clear() {
        f12.reset();
        f1.clear();
        f2.clear();
        f3.clear();
        f4.clear();
        f5.clear();
        f.clear();
    }

    auto get_functions() const {
        return std::make_tuple(f1,f2,f3,f4,f5,f);
    }

    auto get_ccpairfunctions() const {
        CCPairFunction p1;
        CCPairFunction p2(f);
        CCPairFunction p3({f1,f2},{f1,f3});
        CCPairFunction p4(f12.get(),{f1,f2},{f1,f3});
        return std::make_tuple(p1,p2,p3);
    }

};

data data1;

int test_constructor(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;
    test_output t1("CCPairFunction::constructor");

    real_function_6d f=real_factory_6d(world);
    auto [f1,f2,f3,f4,f5,ff]=data1.get_functions();

    vector_real_function_3d  a= zero_functions<double,3>(world,3);
    vector_real_function_3d  b= zero_functions<double,3>(world,3);
    CCConvolutionOperator f12(world, OT_F12, parameter);
    t1.checkpoint(true,"preparation");

    CCPairFunction p1;
    CCPairFunction p2(f);
    CCPairFunction p3({f1,f2},{f1,f3});
    CCPairFunction p4(&f12,{f1,f2},{f2,f3});
    t1.checkpoint(true,"construction");

    {
        MADNESS_CHECK(p2.is_pure());
        MADNESS_CHECK(!p2.is_decomposed());
        MADNESS_CHECK(!p2.is_op_decomposed());
        auto p = p2.pure();
        auto ff = p2.pure().get_function();
        MADNESS_CHECK(!(ff.get_impl()==f.get_impl())); // deep copy of f
    }
    t1.checkpoint(true,"checks on pure");

    {
        MADNESS_CHECK(!p3.is_pure());
        MADNESS_CHECK(p3.is_decomposed());
        MADNESS_CHECK(!p3.is_op_decomposed());
        auto a1=p3.get_a();
        auto b1=p3.get_b();
    }
    t1.checkpoint(true,"checks on decomposed");

    {
        MADNESS_CHECK(!p4.is_pure());
        MADNESS_CHECK(p4.is_decomposed());
        MADNESS_CHECK(p4.is_op_decomposed());
        auto a1=p4.get_a();
        auto b1=p4.get_b();
        auto op=p4.get_operator();
    }
    t1.checkpoint(true,"checks on op_decomposed");

    {
        CCPairFunction tmp;
        tmp=p2;
        MADNESS_CHECK(tmp.is_pure());
        MADNESS_CHECK(tmp.component==p2.component);

        tmp=p3;
        MADNESS_CHECK(!tmp.is_pure());
        MADNESS_CHECK(tmp.component==p3.component);

    }
    t1.checkpoint(true,"checks on assignment");
    t1.end();

    return  (t1.get_final_success()) ? 0 : 1;
}

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
    t1.set_cout_to_terminal();
    // f^2 = 1/(4y^2)(1 - 2*f2(y) + f2(2y)) , f2(2y) =f2(y)^2
    // operator apply of SlaterF12Operator includes a factor of 1/(2 gamma) and the identity
    // operator apply of SlaterOperator has no further terms
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
    const double term1= prefactor * (aR.inner(a) * bR.inner(b));
    const double term2= prefactor * 2.0 * bb.inner(af2a) ;
    const double term3= prefactor *  bb.inner(affa);
    const double ab_f2_ab = prefactor * (aR.inner(a) * bR.inner(b) - 2.0 * bb.inner(af2a) + bb.inner(affa));
    const double ab_f_ab = inner(f(aa), bb);

    const long rank = world.rank();
    auto printer = [&rank](std::string msg, const double result, const double reference, const double time) {
        if (rank == 0) {
            long len = std::max(0l, long(40 - msg.size()));
            msg += std::string(len, ' ');
            std::cout << msg << std::fixed << std::setprecision(8) << "result, ref, diff "
                 << result << " " << reference << " " << std::setw(9) << std::setprecision(2) << std::scientific
                 << result - reference
                 << ", elapsed time " << std::fixed << std::setw(5) << std::setprecision(2) << time << std::endl;
        }
    };
    t1.checkpoint(true, "initialization");

    double time_init = timer.reset();
    if (world.rank() == 0) print("time spent in initializing ", time_init);
    {
        CCPairFunction fab(&f12, a, b);
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
        CCPairFunction ab(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
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
        CCPairFunction fab(fab_6d);
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
        CCPairFunction ab1(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
        CCPairFunction fab(&f12, a, b);
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
        CCPairFunction fab(fab_6d);
        CCPairFunction ab2(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
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
//            CCPairFunction fab(fab_6d);
//            CCPairFunction ab2(ab_6d);
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
        CCPairFunction fab(&f12, a, b);
        CCPairFunction ab2(fab_6d);
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
        CCPairFunction fab(&f12, a, b);
        CCPairFunction ab2(mo_ket_.get_vecfunction(), mo_ket_.get_vecfunction());
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

    return  (t1.get_final_success()) ? 0 : 1;
}



int test_partial_inner(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;

    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction::test_partial_inner");
    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();

    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f2, f3};

    CCConvolutionOperator f12(world, OT_F12, parameter);

    CCPairFunction p1(f);
    CCPairFunction p2(a,b);
    CCPairFunction p3(&f12,a,b);

    double g11=inner(f1,f1);
    double g22=inner(f2,f2);
    double g12=inner(f1,f2);
    double g13=inner(f1,f3);
    double g23=inner(f2,f3);
    real_function_3d gf11=f12(f1*f1);
    real_function_3d gf12=f12(f1*f2);
    real_function_3d gf13=f12(f1*f3);
    print("g11, g22",g11,g22);

    print("time in preparation",timer.reset());
    t1.checkpoint(true,"prep");

    // test pure
    {
        real_function_3d r=inner(p1,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g11 * g22;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"pure -- 1");
    }
    // test pure
    {
        real_function_3d r=inner(p1,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g12 * g12;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"pure -- 2");
    }
    // test decomposed
    {
        real_function_3d r=inner(p2,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g11 * g22 + g12 * g23;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"decomposed -- 1");
    }
    // test decomposed
    {
        real_function_3d r=inner(p2,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g12 * g12 + g13 * g22;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"decomposed -- 2");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f1 f2 > + < f1 f2 | f | f2 f3>
        real_function_3d r=inner(p3,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=inner(gf11*f2,f2) + inner(gf12*f2,f3);
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"op_decomposed -- 1");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f2 f1 > + < f1 f2 | f | f3 f2>
        real_function_3d r=inner(p3,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=inner(gf12*f1,f2) + inner(gf13*f2,f2);
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm<FunctionDefaults<3>::get_thresh());
        t1.checkpoint(good,"op_decomposed -- 2");
    }


    return  (t1.get_final_success()) ? 0 : 1;
}

int test_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;
    test_output t1("CCPairFunction::test_apply");

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_scalar_multiplication(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    int success=0;
    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction::test_scalar_multiplication");

    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();

    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f3, f1};

    print("time in preparation",timer.reset());
    t1.checkpoint(true,"prep");

    CCPairFunction p(f);
    CCPairFunction p1(a,b);
    double norm1=inner(p,p1);
    double pnorm=inner(p,p);
    double p1norm=inner(p1,p1);
    print("pnorm,p1norm",pnorm,p1norm);
    double anorm=norm2(world,p1.get_a());
    double bnorm=norm2(world,p1.get_b());
    print("anorm,bornm",anorm,bnorm);

    p*=2.0;
    p1*=2.0;
    anorm=norm2(world,p1.get_a());
    bnorm=norm2(world,p1.get_b());
    print("anorm,bornm",anorm,bnorm);
    pnorm=inner(p,p);
    p1norm=inner(p1,p1);
    print("pnorm,p1norm",pnorm,p1norm);
    double norm2=inner(p,p1);
    print("norm1,norm2",norm1,norm2);

    print("time in scaling and inner",timer.reset());
    bool bsuccess=fabs(4.0*norm1-norm2)<FunctionDefaults<3>::get_thresh();
    t1.checkpoint(bsuccess,"scaling");
    if (bsuccess) success++;

    t1.end();
    return  (t1.get_final_success()) ? 0 : 1;
}

int test_swap_particles(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                        const CCParameters& parameters) {
    int success = 0;
    test_output t1("CCPairFunction::swap_particles");
    CCTimer timer(world, "testing swap_particles");

    // prepare
    auto one = [](const coord_3d& r) { return 1.0; };
    real_function_3d R2 = real_factory_3d(world).f(one);

    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();
    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f3, f1};

    // test decomposed
    {
        CCPairFunction p1(a, b);
        CCPairFunction p2(b, a);

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
        CCPairFunction p(f);
        CCPairFunction p_swapped=p.swap_particles();
        double pnorm=p.get_function().norm2();
        double psnorm=p_swapped.get_function().norm2();
        print("p/s norm",pnorm,psnorm);

        CCPairFunction p1({f1}, {f2});
        CCPairFunction p2({f2}, {f1});
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
    t1.end();

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_projector(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                           const CCParameters& parameter) {
    int success=0;
    test_output t1("CCPairFunction::test_projector");
    CCTimer timer(world, "testing");

    t1.set_cout_to_logger();
    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();
    double nf1=f1.norm2();
    double nf2=f2.norm2();
    f1.scale(1.0/nf1);
    f2.scale(1.0/nf1);
    std::vector<real_function_3d> a = {f1+f2, f2+f3, f3+f4};
    std::vector<real_function_3d> b = {f3, f1, f2};
    std::vector<real_function_3d> o = orthonormalize_cd<double,3>({f1-f3, f5});  // projects on an orthonormal basis
    a = orthonormalize_cd<double,3>(a);  // projects on an orthonormal basis
    CCConvolutionOperator f12(world, OT_F12, parameter);

    {
        double n1=inner(f1,o[0]);
        double n2=inner(f2,o[0]);
        auto ovlp=matrix_inner(world,o,o);
        print("<o | o>", ovlp);
    }

    CCPairFunction p1(a,b);
    CCPairFunction p2(&f12,a,b);
    CCPairFunction p3(f); // outer (f1,f2)

    std::vector<CCPairFunction> vp1({p1});
    std::vector<CCPairFunction> vp2({p2});
    std::vector<CCPairFunction> vp3({p3});

    Projector<double,3> O(o,o);
    QProjector<double,3> Q(world,o,o);
    StrongOrthogonalityProjector<double,3> Q12(world);
    Q12.set_spaces(o);

    // some hardwire test
    {
        // <ab|k ><k |ab>  = <a|k><k|a><b|b>
        auto ak_tensor=matrix_inner(world,a,o);
        double ak=ak_tensor.trace(ak_tensor);
        auto aa=matrix_inner(world,a,a);
        auto bb=matrix_inner(world,b,b);
        Projector<double,3> O1(a,a);
        O.set_particle(0);
        O1.set_particle(0);
        Projector<double,3> O2(a,a);
        O2.set_particle(0);
        double n1=inner(vp1,O1(vp2));
        double n1a=inner(O1(vp1),vp2);
        print("n1,n1a",n1,n1a);
        double n2=inner(vp1,O2(vp2));
        double n2a=inner(O2(vp1),vp2);
        print("n2,n2a",n1,n1a);

    }
    timer.reset();


    // test {O,Q} with particles 1,2 on all three CCPairFunction types
    {
        int i=0;
        for (auto p : {p1,p2,p3}) {
            std::string s="CCPairFunction p"+std::to_string(++i);
            std::vector<CCPairFunction> vp({p});
            for (int particle : {0,1}) {
                O.set_particle(particle);
                Q.set_particle(particle);
                auto Op = O(vp);
                auto Qp = Q(vp);

                double n1 = inner(Op, vp3) + inner(Qp, vp3);
                double n2 = inner(vp, vp3);
                print("n1 (RI), n2", n1, n2, fabs(n1 - n2));
                t1.checkpoint(fabs(n1 - n2) < FunctionDefaults<3>::get_thresh(), "RI with particle "+std::to_string(particle)+" on "+s ,timer.reset());

                auto Op1 = O(Op);
                double n3=inner(Op, vp3);
                double n4=inner(Op1, vp3);
                print("n3, n4", n3, n4);
                t1.checkpoint(fabs(n3-n4) < FunctionDefaults<3>::get_thresh(), "idempotency with particle "+std::to_string(particle)+" on "+s,timer.reset());

            }
            // testing SO projector
            auto SOp=Q12(vp);
            double n1,n2,n3;
            {
                Q.set_particle(0);
                auto tmp = Q(vp);
                Q.set_particle(1);
                auto tmp1 = Q(tmp);
                n2=inner(vp3,tmp1);
            }
            {
                Q.set_particle(1);
                auto tmp = Q(vp);
                Q.set_particle(0);
                auto tmp1 = Q(tmp);
                n3=inner(vp3,tmp1);
            }
            n1 = inner(SOp, vp3);
            print("SO: n1, n2, n3", n1, n2, n3, fabs(n1 - n2));
            double zero=fabs(n1-n2) + fabs(n1-n3) + fabs(n2-n3);
            t1.checkpoint(zero < FunctionDefaults<3>::get_thresh(), "SO operator on "+s,timer.reset() );
        }
    }
    t1.end();

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_dirac_convolution(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                           const CCParameters& parameter) {
    int success=0;
    test_output t1("CCPairFunction::test_dirac_convolution");

    return  (t1.get_final_success()) ? 0 : 1;
}

/// testing <ij | g Q f | ij> = 0.032 mEh
int test_helium(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                   const CCParameters& parameters) {

    test_output t1("CCPairFunction::test_helium");
    CCTimer timer(world, "testing");

    real_function_3d Vnuc = real_factory_3d(world).f([](const coord_3d& r) {return -2.0/(r.normf()+1.e-8);});
    real_function_3d psi  = real_factory_3d(world).f([](const coord_3d& r) {return exp(-r.normf());});

    auto iterate=[&world](const real_function_3d& potential, real_function_3d& psi, double& eps) {
        real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
        real_function_3d Vpsi = (potential*psi);
        Vpsi.scale(-2.0).truncate();
        real_function_3d tmp = op(Vpsi).truncate();
        double norm = tmp.norm2();
        real_function_3d r = tmp-psi;
        double rnorm = r.norm2();
        double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
        if (world.rank() == 0) {
            print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        }
        psi = tmp.scale(1.0/norm);
        eps = eps_new;
    };
    psi.truncate();
    psi.scale(1.0/psi.norm2());
    double eps = -0.6;
    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    for (int iter=0; iter<10; iter++) {
        real_function_3d rho = square(psi).truncate();
        real_function_3d potential = Vnuc + op(rho).truncate();
        iterate(potential, psi, eps);
    }
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }
    real_function_3d rho = square(psi);
    double two_electron_energy = inner(op(rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho);
    double total_energy = kinetic_energy + two_electron_energy + nuclear_attraction_energy;

    t1.checkpoint(fabs(total_energy+2.8616533)<1.e-4,"helium iterations",timer.reset());
    print("ke, total", kinetic_energy, total_energy);


    CCConvolutionOperator::Parameters param;
    CCConvolutionOperator f(world,OT_F12,param);
    CCConvolutionOperator g(world,OT_G12,param);

    CCPairFunction fij(&f,psi,psi);
    CCPairFunction gij(&g,psi,psi);
    CCPairFunction ij({psi},{psi});
    std::vector<CCPairFunction> vfij={fij};
    std::vector<CCPairFunction> vgij={gij};
    std::vector<CCPairFunction> vij={ij};

    StrongOrthogonalityProjector<double,3> SO(world);
    SO.set_spaces({psi},{psi},{psi},{psi});
    std::vector<CCPairFunction> Qfij=SO(vfij);
    std::vector<CCPairFunction> Qgij=SO(vgij);

    double result1=inner(vgij,Qfij);
    print("<ij | g (Q f | ij>)",result1);
    double result2=inner(vfij,Qgij);
    print("(<ij | g Q) f | ij>",result2);
    bool good=fabs(result1-result2)<1.e-5;
    good=good and fabs(result1+3.2624783e-02)<1.e-5;
    t1.checkpoint(good,"V matrix element",timer.reset());

    return  (t1.get_final_success()) ? 0 : 1;

}

/** functionality
 *
 *  - ctor                      OK
 *  - assignment                OK
 *  - add
 *  - scalar multiplication     OK
 *  - inner                     OK
 *  - inner_partial             OK
 *  - swap_particles            OK
 *  - apply
 *  - apply_partial (i.e. exchange)
 *  - serialize
 *  - collapse_to_pure (excl g!)
 *  - mul_partial
 */

int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<3>::set_thresh(1.e-5);
    FunctionDefaults<3>::set_cubic_cell(-1.0,1.0);
    FunctionDefaults<6>::set_cubic_cell(-1.0,1.0);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    int isuccess=0;
#ifdef USE_GENTENSOR

    try {
        parser.set_keyval("geometry", "he");
        parser.print_map();
        Molecule mol(world, parser);
        mol.print();
        CCParameters ccparam(world, parser);

        data1=data(world,ccparam);

        std::shared_ptr<NuclearCorrelationFactor> ncf = create_nuclear_correlation_factor(world,
                         mol, nullptr, std::make_pair("slater", 2.0));

        isuccess+=test_constructor(world, ncf, mol, ccparam);
//        isuccess+=test_overlap(world, ncf, mol, ccparam);
//        isuccess+=test_swap_particles(world, ncf, mol, ccparam);
//        isuccess+=test_scalar_multiplication(world, ncf, mol, ccparam);
//        isuccess+=test_partial_inner(world, ncf, mol, ccparam);
//        isuccess+=test_projector(world, ncf, mol, ccparam);
        FunctionDefaults<3>::set_cubic_cell(-10,10);
        isuccess+=test_helium(world,ncf,mol,ccparam);
        data1.clear();
    } catch (std::exception& e) {
        madness::print("an error occured");
        madness::print(e.what());
        data1.clear();
    }
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}