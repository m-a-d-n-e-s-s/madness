//
// Created by Florian Bischoff on 6/27/22.
//

#include<chem/ccpairfunction.h>
#include <chem/CCStructures.h>


namespace madness {
madness::CCPairFunction
CCPairFunction::operator=(const CCPairFunction& other) {
    MADNESS_ASSERT(type == other.type);
    a = other.a;
    b = other.b;
    op = other.op;
    x = other.x;
    y = other.y;
    u = other.u;
    return *this;
}

madness::CCPairFunction
CCPairFunction::copy() const {
    if (type == PT_FULL) {
        return CCPairFunction(world, madness::copy(u));
    } else if (type == PT_DECOMPOSED) {
        return CCPairFunction(world, madness::copy(world, a), madness::copy(world, b));
    } else if (type == PT_OP_DECOMPOSED) {
        return CCPairFunction(world, op, CCFunction(madness::copy(x.function), x.i, x.type),
                              CCFunction(madness::copy(y.function), y.i, y.type));
    } else MADNESS_EXCEPTION("Unknown type", 1)

    ;
}

madness::CCPairFunction
CCPairFunction::invert_sign() {
    if (type == PT_FULL) {
        real_function_6d uc = madness::copy(u);
        uc.scale(-1.0);
        u = uc;
    } else if (type == PT_DECOMPOSED) {
        vector_real_function_3d ac = madness::copy(world, a);
        scale(world, ac, -1.0);
        a = ac;
    } else if (type == PT_OP_DECOMPOSED) {
        real_function_3d tmp = madness::copy(x.function);
        tmp.scale(-1.0);
        x.function = tmp;
    } else MADNESS_EXCEPTION("wrong type in CCPairFunction invert_sign", 1)

    ;
    return *this;
}

void
CCPairFunction::print_size() const {
    if (type == PT_FULL) {
        u.print_size(name());
    } else if (type == PT_DECOMPOSED) {
        madness::print_size(world, a, "a from " + name());
        madness::print_size(world, b, "b from " + name());
    } else if (type == PT_OP_DECOMPOSED) {
        x.function.print_size(x.name() + " from " + name());
        y.function.print_size(y.name() + " from " + name());
    } else MADNESS_EXCEPTION("Unknown type in CCPairFunction, print_size", 1)

    ;
}


double
CCPairFunction::make_xy_u(const CCFunction& xx, const CCFunction& yy) const {
    double result = 0.0;
    switch (type) {
        default: MADNESS_EXCEPTION("Undefined enum", 1);
        case PT_FULL: {
            real_function_6d ij = CompositeFactory<double, 6, 3>(world).particle1(madness::copy(xx.function)).particle2(
                    madness::copy(yy.function));
            result = inner(u, ij);
        }
            break;
        case PT_DECOMPOSED: {
            for (size_t i = 0; i < a.size(); i++)
                result += (xx.function.inner(a[i])) * (yy.function.inner(b[i]));
        }
            break;
        case PT_OP_DECOMPOSED: {
            result = yy.function.inner((*op)(xx, x) * y.function);
        }
            break;
    }
    return result;
}

std::string CCPairFunction::name() const {
    if (type == PT_FULL) return "|u>";
    else if (type == PT_DECOMPOSED) return "|ab>";
    else if (type == PT_OP_DECOMPOSED) return op->name() + "|xy>";
    return "???";
}

real_function_3d CCPairFunction::project_out(const CCFunction& f, const size_t particle) const {
    MADNESS_ASSERT(particle == 1 or particle == 2);
    real_function_3d result;
    switch (type) {
        default: MADNESS_EXCEPTION("Undefined enum", 1);
        case PT_FULL :
            result = u.project_out(f.function, particle - 1); // this needs 0 or 1 for particle but we give 1 or 2
            break;
        case PT_DECOMPOSED :
            result = project_out_decomposed(f.function, particle);
            break;
        case PT_OP_DECOMPOSED:
            result = project_out_op_decomposed(f, particle);
            break;
    }
    if (not result.is_initialized()) MADNESS_EXCEPTION("Result of project out on CCPairFunction was not initialized",
                                                       1);
    return result;
}

// result is: <x|op12|f>_particle
real_function_3d
CCPairFunction::dirac_convolution(const CCFunction& x, const CCConvolutionOperator& op, const size_t particle) const {
    real_function_3d result;
    switch (type) {
        default: MADNESS_EXCEPTION("Undefined enum", 1);
        case PT_FULL:
            result = op(x, u, particle);
            break;
        case PT_DECOMPOSED :
            result = dirac_convolution_decomposed(x, op, particle);
            break;
        case PT_OP_DECOMPOSED: MADNESS_EXCEPTION("op_decomposed dirac convolution not yet implemented", 1);
    }
    return result;
}

CCPairFunction CCPairFunction::swap_particles() const {
    switch (type) {
        default: MADNESS_EXCEPTION("Undefined enum", 1);
        case PT_FULL:
            return swap_particles_pure();
            break;
        case PT_DECOMPOSED:
            return swap_particles_decomposed();
            break;
        case PT_OP_DECOMPOSED:
            return swap_particles_op_decomposed();
            break;
    }
    MADNESS_EXCEPTION("swap_particles in CCPairFunction: we should not end up here", 1);
}

real_function_3d CCPairFunction::project_out_decomposed(const real_function_3d& f, const size_t particle) const {
    real_function_3d result = real_factory_3d(world);
    const std::pair<vector_real_function_3d, vector_real_function_3d> decompf = assign_particles(particle);
    Tensor<double> c = inner(world, f, decompf.first);
    for (size_t i = 0; i < a.size(); i++) result += c(i) * decompf.second[i];
    return result;
}

real_function_3d CCPairFunction::project_out_op_decomposed(const CCFunction& f, const size_t particle) const {
    if (particle == 1) {
        return (*op)(f, x) * y.function;
    } else if (particle == 2) {
        return (*op)(f, y) * x.function;
    } else {
        MADNESS_EXCEPTION("project_out_op_decomposed: particle must be 1 or 2", 1);
        return real_factory_3d(world);
    }
}

real_function_3d CCPairFunction::dirac_convolution_decomposed(const CCFunction& bra, const CCConvolutionOperator& op,
                                                              const size_t particle) const {
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
        return std::make_pair(a, b);
    } else if (particle == 2) {
        return std::make_pair(b, a);
    } else {
        MADNESS_EXCEPTION("project_out_decomposed: Particle is neither 1 nor 2", 1);
        return std::make_pair(a, b);
    }
}

CCPairFunction CCPairFunction::swap_particles_pure() const {
    // CC_Timer timer_swap(world,"swap particles");
    // this could be done more efficiently for SVD, but it works decently
    std::vector<long> map(6);
    map[0] = 3;
    map[1] = 4;
    map[2] = 5;     // 2 -> 1
    map[3] = 0;
    map[4] = 1;
    map[5] = 2;     // 1 -> 2
    // timer_swap.info();
    real_function_6d swapped_u = mapdim(u, map);
    return CCPairFunction(world, swapped_u);
}

CCPairFunction CCPairFunction::swap_particles_decomposed() const {
    return CCPairFunction(world, b, a);
}

CCPairFunction CCPairFunction::swap_particles_op_decomposed() const {
    return CCPairFunction(world, op, y, x);
}


} // namespace madness