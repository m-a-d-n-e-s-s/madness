//
// Created by Florian Bischoff on 6/27/22.
//

#include<chem/ccpairfunction.h>
#include<chem/CCStructures.h>
#include<chem/projector.h>

using namespace madness;

namespace madness {


madness::CCPairFunction
CCPairFunction::invert_sign() {
    (*this)*=(-1.0);
    return *this;
}

double
CCPairFunction::make_xy_u(const CCFunction& xx, const CCFunction& yy) const {
    double result = 0.0;
    if (is_pure()) {
        World& world=xx.function.world();
        real_function_6d ij = CompositeFactory<double, 6, 3>(world).particle1(madness::copy(xx.function))
                .particle2( madness::copy(yy.function));
        result = inner(pure().get_function(), ij);
    } else if (is_decomposed_no_op()) {
        for (size_t i = 0; i < get_a().size(); i++)
            result += (xx.function.inner(get_a()[i])) * (yy.function.inner(get_b()[i]));
    } else if (is_op_decomposed()) {
        const CCConvolutionOperator& op=*decomposed().get_operator_ptr();
        result = yy.function.inner(op(xx, get_a()[0]) * get_b()[0]);
    }
    return result;
}

real_function_3d CCPairFunction::project_out(const CCFunction& f, const size_t particle) const {
    MADNESS_ASSERT(particle == 1 or particle == 2);
    real_function_3d result;
    if (is_pure()) {
        result = pure().get_function().project_out(f.function, particle - 1); // this needs 0 or 1 for particle but we give 1 or 2
    } else if (is_decomposed_no_op()) {
        result = project_out_decomposed(f.function, particle);
    } else if (is_op_decomposed()) {
        result = project_out_op_decomposed(f, particle);
    }
    if (not result.is_initialized()) MADNESS_EXCEPTION("Result of project out on CCPairFunction was not initialized",
                                                       1);
    return result;
}

// result is: <x|op12|f>_particle
real_function_3d
CCPairFunction::dirac_convolution(const CCFunction& x, const CCConvolutionOperator& op, const size_t particle) const {
    real_function_3d result;
    if (is_pure()) {
        result = op(x, pure().get_function(), particle);
    } else if (is_decomposed_no_op()) {
        result = dirac_convolution_decomposed(x, op, particle);
    } else {
        MADNESS_EXCEPTION("op_decomposed dirac convolution not yet implemented", 1);
    }
    return result;
}

// CCPairFunction CCPairFunction::swap_particles() const {
//     switch (type) {
//         default: MADNESS_EXCEPTION("Undefined enum", 1);
//         case PT_FULL:
//             return swap_particles_pure();
//             break;
//         case PT_DECOMPOSED:
//             return swap_particles_decomposed();
//             break;
//         case PT_OP_DECOMPOSED:
//             return swap_particles_op_decomposed();
//             break;
//     }
//     MADNESS_EXCEPTION("swap_particles in CCPairFunction: we should not end up here", 1);
// }



real_function_3d inner(const CCPairFunction& c, const real_function_3d& f,
                       const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    auto v11=std::array<int,3>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)});
    auto v22=std::array<int,3>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)});

    return c.partial_inner(f,v11,v22);
}


real_function_3d CCPairFunction::partial_inner(const real_function_3d& f,
                                               const std::array<int, 3>& v1,
                                               const std::array<int, 3>& v2) const {
    auto a012=std::array<int,3>{0,1,2};
    auto a345=std::array<int,3>{3,4,5};
    MADNESS_CHECK(v2==a012 ); // only 3 dimension in f
    MADNESS_CHECK(v1==a012 or v1== a345); // 6 dimension in f
    int particle=-1;
    if (v1== a012) particle=0;
    if (v1== a345) particle=1;

    real_function_3d result;

    if (is_pure()) {
        result = pure().get_function().project_out(f, particle);
    } else if (is_decomposed_no_op()) {
        result = project_out_decomposed(f, particle+1);
    } else if (is_op_decomposed()) {
        result = project_out_op_decomposed(f, particle+1);
    } else {
        MADNESS_EXCEPTION("confused state in CCPairFunction::partial_inner",1);
    }
    return result;
}

real_function_3d CCPairFunction::project_out_decomposed(const real_function_3d& f, const size_t particle) const {
    World& world=f.world();
    real_function_3d result = real_factory_3d(world);
    const std::pair<vector_real_function_3d, vector_real_function_3d> decompf = assign_particles(particle);
    Tensor<double> c = inner(world, f, decompf.first);
    for (size_t i = 0; i < get_a().size(); i++) result += c(i) * decompf.second[i];
    return result;
}

real_function_3d CCPairFunction::project_out_op_decomposed(const CCFunction& f, const size_t particle) const {
    World& world=f.get().world();
    const CCConvolutionOperator& op=*decomposed().get_operator_ptr();
    if (particle == 1) {
//        return op(f, get_a()[0]) * get_b()[0];
        // result(2) = < f(1) | op(1,2) | a_i(1) b_i(2) >
        return sum(world,mul(world,op(f.f()* get_a()),get_b()));
    } else if (particle == 2) {
//        return op(f, get_b()[0]) * get_a()[0];
        return sum(world,mul(world,op(f.f()* get_b()),get_a()));
    } else {
        MADNESS_EXCEPTION("project_out_op_decomposed: particle must be 1 or 2", 1);
        return real_factory_3d(world);
    }
}

real_function_3d CCPairFunction::dirac_convolution_decomposed(const CCFunction& bra, const CCConvolutionOperator& op,
                                                              const size_t particle) const {
    World& world=bra.function.world();
    const std::pair<vector_real_function_3d, vector_real_function_3d> f = assign_particles(particle);
    const vector_real_function_3d braa = mul(world, bra.function, f.first);
    const vector_real_function_3d braga = op(braa);
    real_function_3d result = real_factory_3d(world);
    for (size_t i = 0; i < braga.size(); i++) result += braga[i] * f.second[i];
    return result;
}


const std::pair<vector_real_function_3d, vector_real_function_3d>
CCPairFunction::assign_particles(const size_t particle) const {
    if (particle == 1) {
        return std::make_pair(get_a(), get_b());
    } else if (particle == 2) {
        return std::make_pair(get_b(), get_a());
    } else {
        MADNESS_EXCEPTION("project_out_decomposed: Particle is neither 1 nor 2", 1);
        return std::make_pair(get_a(), get_b());
    }
}


double CCPairFunction::inner_internal(const CCPairFunction& other, const real_function_3d& R2) const {
    const CCPairFunction& f1=*this;
    const CCPairFunction& f2=other;

    double thresh=FunctionDefaults<6>::get_thresh();

    double result = 0.0;
    if (f1.is_pure() and f2.is_pure()) {
        CCTimer tmp(world(), "making R1R2|u>");
        if (R2.is_initialized()) {
            real_function_6d R1u = multiply<double, 6, 3>(::copy(f1.pure().get_function()), ::copy(R2), 1);
            real_function_6d R1R2u = multiply<double, 6, 3>(R1u, ::copy(R2), 2);     // R1u function now broken
            tmp.info();
            result = f2.pure().get_function().inner(R1R2u);
        } else {
            result = f2.pure().get_function().inner(f1.pure().get_function());
        }
    } else if (f1.is_pure() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = R2.is_initialized() ?  R2* f2.get_a() : f2.get_a();
        vector_real_function_3d b = R2.is_initialized() ?  R2* f2.get_b() : f2.get_b();
        for (size_t i = 0; i < a.size(); i++) {
            real_function_6d ab = CompositeFactory<double, 6, 3>(world()).particle1(::copy(a[i])).particle2(::copy(b[i]));
            result += f1.pure().get_function().inner(ab);
        }
    } else if (f1.is_pure() and f2.is_op_decomposed()) {
        real_function_3d x = R2.is_initialized() ? R2 * f2.get_a()[0] : f2.get_a()[0];
        real_function_3d y = R2.is_initialized() ? R2 * f2.get_b()[0] : f2.get_a()[0];
        real_function_6d op;
        if (f2.decomposed().get_operator_ptr()->type() == OT_F12)
            op = TwoElectronFactory(world()).dcut(f2.get_operator().parameters.lo).gamma(f2.get_operator().parameters.gamma).f12()
                    .thresh(thresh);
        else if (f2.get_operator().type() == OT_G12)
            op = TwoElectronFactory(world()).dcut(f2.get_operator().parameters.lo).thresh(thresh);
        else MADNESS_EXCEPTION(
                ("6D Overlap with operatortype " + assign_name(f2.get_operator().type()) + " not supported").c_str(), 1);
        real_function_6d opxy = CompositeFactory<double, 6, 3>(world()).g12(op).particle1(::copy(x)).particle2(::copy(y));
        return f1.pure().get_function().inner(opxy);
    } else if (f1.is_decomposed_no_op() and f2.is_pure()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        vector_real_function_3d a = R2.is_initialized() ?  R2* f1.get_a() : f1.get_a();
        vector_real_function_3d b = R2.is_initialized() ?  R2* f1.get_b() : f1.get_b();
        for (size_t i = 0; i < a.size(); i++) {
            real_function_6d ab = CompositeFactory<double, 6, 3>(world()).particle1(::copy(a[i])).particle2(::copy(b[i]));
            result += f2.pure().get_function().inner(ab);
        }
    } else if (f1.is_decomposed_no_op() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = R2.is_initialized() ?  R2* f2.get_a() : f2.get_a();
        vector_real_function_3d b = R2.is_initialized() ?  R2* f2.get_b() : f2.get_b();
        // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
        result = (matrix_inner(world(), a, f1.get_a()).emul(matrix_inner(world(), b, f1.get_b()))).sum();

//         for (size_t i1 = 0; i1 < f1.get_a().size(); i1++) {
//             for (size_t i2 = 0; i2 < f2.get_a().size(); i2++) {
//                 const double aa = a[i1].inner(f2.get_a()[i2]);
//                 const double bb = b[i1].inner(f2.get_b()[i2]);
//                 result += aa * bb;
//             }
//         }
    } else if (f1.is_decomposed_no_op() and f2.is_op_decomposed()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        vector_real_function_3d a = f1.get_a();
        vector_real_function_3d b = f1.get_b();
        real_function_3d x = R2.is_initialized() ? R2 * f2.get_a()[0] : f2.get_a()[0];
        real_function_3d y = R2.is_initialized() ? R2 * f2.get_b()[0] : f2.get_a()[0];
        for (size_t i = 0; i < a.size(); i++) {
            real_function_3d ax = (a[i] * x).truncate();
            real_function_3d aopx = f2.get_operator()(ax);
            result += b[i].inner(aopx * y);
        }
    } else if (f1.is_op_decomposed() and f2.is_pure()) {
        real_function_3d x = R2.is_initialized() ? R2 * f1.get_a()[0] : f1.get_a()[0];
        real_function_3d y = R2.is_initialized() ? R2 * f1.get_b()[0] : f1.get_a()[0];
        real_function_6d op;
        if (f1.get_operator().type() == OT_F12)
            op = TwoElectronFactory(world()).dcut(f1.get_operator().parameters.lo).gamma(f1.get_operator().parameters.gamma).f12()
                    .thresh(thresh);
        else if (f1.get_operator().type() == OT_G12)
            op = TwoElectronFactory(world()).dcut(f1.get_operator().parameters.lo).thresh(thresh);
        else MADNESS_EXCEPTION(
                ("6D Overlap with operatortype " + assign_name(f1.get_operator().type()) + " not supported").c_str(), 1);
        real_function_6d opxy = CompositeFactory<double, 6, 3>(world()).g12(op).particle1(::copy(x)).particle2(::copy(y));
        return f2.pure().get_function().inner(opxy);
    } else if (f1.is_op_decomposed() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = f2.get_a();
        vector_real_function_3d b = f2.get_b();
        real_function_3d x = R2.is_initialized() ? R2 * f1.get_a()[0] : f1.get_a()[0];
        real_function_3d y = R2.is_initialized() ? R2 * f1.get_b()[0] : f1.get_a()[0];
        for (size_t i = 0; i < a.size(); i++) {
            real_function_3d ax = (a[i] * x).truncate();
            real_function_3d aopx = f1.get_operator()(ax);
            result += b[i].inner(aopx * y);
        }
    } else if (f1.is_op_decomposed() and f2.is_op_decomposed()) {
        MADNESS_ASSERT(f1.get_operator().type() == OT_F12 and f2.get_operator().type() == OT_F12);     // in principle we have everything for f12g12 but we will not need it
        // we use the slater operator which is S = e^(-y*r12), y=gamma
        // the f12 operator is: 1/2y*(1-e^(-y*r12)) = 1/2y*(1-S)
        // so the squared f12 operator is: f*f = 1/(4*y*y)(1-2S+S*S), S*S = S(2y) = e(-2y*r12)
        // we have then: <xy|f*f|xy> = 1/(4*y*y)*(<xy|xy> - 2*<xy|S|xy> + <xy|SS|xy>)
        const double gamma =f1.get_operator().parameters.gamma;
        const double prefactor = 1.0 / (4.0 * gamma * gamma);
        SeparatedConvolution<double, 3> S = SlaterOperator(world(), gamma, f1.get_operator().parameters.lo, f1.get_operator().parameters.thresh_op);
        SeparatedConvolution<double, 3> S2 = SlaterOperator(world(), 2.0 * gamma, f1.get_operator().parameters.lo,
                                                            f1.get_operator().parameters.thresh_op);
        real_function_3d x1 = R2.is_initialized() ? R2 * f1.get_a()[0] : f2.get_a()[0];
        real_function_3d y1 = R2.is_initialized() ? R2 * f1.get_b()[0] : f1.get_a()[0];
        real_function_3d xx = f2.get_a()[0] * x1;
        real_function_3d yy = f2.get_b()[0] * y1;
        real_function_3d xSx = S(xx);
        real_function_3d xSSx = S2(xx);
        result = prefactor * (x1.inner(f2.get_a()[0]) * y1.inner(f2.get_b()[0]) - 2.0 * yy.inner(xSx) + yy.inner(xSSx));
    } else MADNESS_EXCEPTION(
            ("CCPairFunction Overlap not supported for combination " + f1.name() + " and " + f2.name()).c_str(), 1)

    ;
    return result;
}

CCPairFunction apply(const ProjectorBase& projector, const CCPairFunction& argument) {
    auto result=madness::apply(projector,std::vector<CCPairFunction> (1,argument));
    MADNESS_CHECK(result.size()==1);
    return result[0];
}

std::vector<CCPairFunction> apply(const ProjectorBase& projector, const std::vector<CCPairFunction>& argument) {
    if (argument.size()==0) return argument;
    World& world=argument.front().world();
    if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) MADNESS_CHECK(P->get_particle()==0 or P->get_particle()==1);
    if (auto Q=dynamic_cast<const Projector<double,3>*>(&projector)) MADNESS_CHECK(Q->get_particle()==0 or Q->get_particle()==1);
    std::vector<CCPairFunction> result;
    for (const auto& pf : argument) {
        if (pf.is_pure()) {
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,3>*>(&projector)) {
                auto tmp=(*SO)(pf.get_function());
                auto tmp2=CCPairFunction(tmp);
                result.push_back(tmp2);
            } else if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) {
                result.push_back(CCPairFunction((*P)(pf.get_function(),P->get_particle()+1)));

            } else if (auto Q=dynamic_cast<const QProjector<double,3>*>(&projector)) {
                result.push_back(CCPairFunction((*Q)(pf.get_function(),Q->get_particle()+1)));

            } else {
                MADNESS_EXCEPTION("CCPairFunction: unknown projector type",1);
            }
        } else if (pf.is_decomposed_no_op()) {  // pair function is sum_i | a_i b_i >
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,3>*>(&projector)) {
                // Q12 | kl > = (1-O1)(1-O2) |kl> = |(1-O1)k (1-O2)l>
                QProjector<double,3> Q1(world,SO->bra1(),SO->ket1());
                QProjector<double,3> Q2(world,SO->bra2(),SO->ket2());
                result.push_back(CCPairFunction(Q1(pf.get_a()),Q2(pf.get_b())));

            } else if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) {
                // P1 | kl > = P1 |kl> = |P1 k l>
                if (P->get_particle()==0) result.push_back(CCPairFunction((*P)(pf.get_a()),pf.get_b()));
                // P2 | kl > = P2 |kl> = |k P2 l>
                if (P->get_particle()==1) result.push_back(CCPairFunction(pf.get_a(),(*P)(pf.get_b())));

            } else if (auto Q=dynamic_cast<const QProjector<double,3>*>(&projector)) {
                // Q1 | kl > = Q1 |kl> = |Q1 k l>
                if (Q->get_particle()==0) result.push_back(CCPairFunction((*Q)(pf.get_a()),pf.get_b()));
                // P2 | kl > = Q2 |kl> = |k Q2 l>
                if (Q->get_particle()==1) result.push_back(CCPairFunction(pf.get_a(),(*Q)(pf.get_b())));
            } else {
                MADNESS_EXCEPTION("CCPairFunction: unknown projector type",1);
            }
        } else if (pf.is_op_decomposed()) {
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,3>*>(&projector)) {
                // Q12 = 1 - O1 (1 - 1/2 O2) - O2 (1 - 1/2 O1)
                QProjector<double,3> Q1(world,SO->bra1(),SO->ket1());
                Q1.set_particle(0);
                QProjector<double,3> Q2(world,SO->bra2(),SO->ket2());
                Q2.set_particle(1);
                auto tmp=Q1(Q2(std::vector<CCPairFunction>({pf})));
                for (auto& t: tmp) result.push_back(t);

            } else if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) {
                std::vector<real_function_3d> tmp= zero_functions_compressed<double,3>(world,P->get_ket_vector().size());

                for (std::size_t i=0; i<pf.get_a().size(); ++i) {
                    real_function_3d a=pf.get_a()[i];
                    real_function_3d b=pf.get_b()[i];
                    if (P->get_particle()==1) std::swap(a,b);

                    std::vector<real_function_3d> kb=mul(world,b,P->get_bra_vector());
                    real_convolution_3d& op=*(pf.get_operator().get_op());
                    std::vector<real_function_3d> f_kb=apply(world,op,kb);
                    std::vector<real_function_3d> a_f_kb=mul(world,a,f_kb);
                    tmp+=a_f_kb;
                }
                truncate(world,tmp);
                if (P->get_particle()==0) result.push_back(CCPairFunction(P->get_ket_vector(),tmp));
                if (P->get_particle()==1) result.push_back(CCPairFunction(tmp,P->get_ket_vector()));

            } else if (auto Q=dynamic_cast<const QProjector<double,3>*>(&projector)) {
                // Q1 f12 |a_i b_i> = f12 |a_i b_i> - \sum_k |k(1) a_i(2)*f_(kb_i)(2) >
                result.push_back(pf);
                // reuse the projector code above
                std::vector<CCPairFunction> tmp=madness::apply(Q->get_P_projector(),std::vector<CCPairFunction>(1,pf));
                for (auto& t : tmp) {
                    t*=-1.0;
                    result.push_back(t);
                }

            } else {
                MADNESS_EXCEPTION("CCPairFunction: unknown projector type",1);
            }

        } else {
            MADNESS_EXCEPTION("confused type in CCPairFunction",1);
        }

    }
//    print("working on ",argument[0].name(),"with ",(&projector)->type(),": result has",result.size(),"components");
    return result;
};


} // namespace madness