//
// Created by Florian Bischoff on 6/27/22.
//

#include<chem/ccpairfunction.h>
#include <chem/CCStructures.h>

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
        return op(f, get_a()[0]) * get_b()[0];
    } else if (particle == 2) {
        return op(f, get_b()[0]) * get_a()[0];
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
        real_function_6d R1u = multiply<double, 6, 3>(::copy(f1.pure().get_function()), copy(R2), 1);
        real_function_6d R1R2u = multiply<double, 6, 3>(R1u, copy(R2), 2);     // R1u function now broken
        tmp.info();
        result = f2.pure().get_function().inner(R1R2u);
    } else if (f1.is_pure() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = R2* f2.get_a();
        vector_real_function_3d b = R2* f2.get_b();
        for (size_t i = 0; i < a.size(); i++) {
            real_function_6d ab = CompositeFactory<double, 6, 3>(world()).particle1(copy(a[i])).particle2(copy(b[i]));
            result += f1.pure().get_function().inner(ab);
        }
    } else if (f1.is_pure() and f2.is_op_decomposed()) {
        real_function_3d x = R2 * f2.get_a()[0];
        real_function_3d y = R2 * f2.get_b()[0];
        real_function_6d op;
        if (f2.decomposed().get_operator_ptr()->type() == OT_F12)
            op = TwoElectronFactory(world()).dcut(f2.get_operator().parameters.lo).gamma(f2.get_operator().parameters.gamma).f12()
                    .thresh(thresh);
        else if (f2.get_operator().type() == OT_G12)
            op = TwoElectronFactory(world()).dcut(f2.get_operator().parameters.lo).thresh(thresh);
        else MADNESS_EXCEPTION(
                ("6D Overlap with operatortype " + assign_name(f2.get_operator().type()) + " not supported").c_str(), 1);
        real_function_6d opxy = CompositeFactory<double, 6, 3>(world()).g12(op).particle1(copy(x)).particle2(copy(y));
        return f1.pure().get_function().inner(opxy);
    } else if (f1.is_decomposed_no_op() and f2.is_pure()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        vector_real_function_3d a = R2* f1.get_a();
        vector_real_function_3d b = R2* f1.get_b();
        for (size_t i = 0; i < a.size(); i++) {
            real_function_6d ab = CompositeFactory<double, 6, 3>(world()).particle1(copy(a[i])).particle2(copy(b[i]));
            result += f2.pure().get_function().inner(ab);
        }
    } else if (f1.is_decomposed_no_op() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = R2* f1.get_a();
        vector_real_function_3d b = R2* f1.get_b();
        for (size_t i1 = 0; i1 < f1.get_a().size(); i1++) {
            for (size_t i2 = 0; i2 < f2.get_a().size(); i2++) {
                const double aa = a[i1].inner(f2.get_a()[i2]);
                const double bb = b[i1].inner(f2.get_b()[i2]);
                result += aa * bb;
            }
        }
    } else if (f1.is_decomposed_no_op() and f2.is_op_decomposed()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        vector_real_function_3d a = f1.get_a();
        vector_real_function_3d b = f1.get_b();
        real_function_3d x = (R2 * f2.get_a()[0]).truncate();
        real_function_3d y = (R2 * f2.get_b()[0]).truncate();
        for (size_t i = 0; i < a.size(); i++) {
            real_function_3d ax = (a[i] * x).truncate();
            real_function_3d aopx = f2.get_operator()(ax);
            result += b[i].inner(aopx * y);
        }
    } else if (f1.is_op_decomposed() and f2.is_pure()) {
        real_function_3d x = R2 * f1.get_a()[0];
        real_function_3d y = R2 * f1.get_b()[0];
        real_function_6d op;
        if (f1.get_operator().type() == OT_F12)
            op = TwoElectronFactory(world()).dcut(f1.get_operator().parameters.lo).gamma(f1.get_operator().parameters.gamma).f12()
                    .thresh(thresh);
        else if (f1.get_operator().type() == OT_G12)
            op = TwoElectronFactory(world()).dcut(f1.get_operator().parameters.lo).thresh(thresh);
        else MADNESS_EXCEPTION(
                ("6D Overlap with operatortype " + assign_name(f1.get_operator().type()) + " not supported").c_str(), 1);
        real_function_6d opxy = CompositeFactory<double, 6, 3>(world()).g12(op).particle1(copy(x)).particle2(copy(y));
        return f2.pure().get_function().inner(opxy);
    } else if (f1.is_op_decomposed() and f2.is_decomposed_no_op()) {
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());
        vector_real_function_3d a = f2.get_a();
        vector_real_function_3d b = f2.get_b();
        real_function_3d x = (R2 * f1.get_a()[0]).truncate();
        real_function_3d y = (R2 * f1.get_b()[0]).truncate();
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
        real_function_3d x1 = f1.get_a()[0] * R2;
        real_function_3d y1 = f1.get_b()[0] * R2;
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


} // namespace madness