//
// Created by Florian Bischoff on 6/27/22.
//

#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/projector.h>
#include<madness/mra/operator.h>

using namespace madness;

namespace madness {


madness::CCPairFunction
CCPairFunction::invert_sign() {
    (*this)*=(-1.0);
    return *this;
}

bool CCPairFunction::is_convertible_to_pure_no_op() const {
    if (has_operator()) {
        const auto type=get_operator().type();
        if (not (type==OpType::OT_SLATER or type==OpType::OT_F12)) return false;
    }
    if (is_decomposed() and (get_a().size()>2)) return false;
    return true;
};


void CCPairFunction::convert_to_pure_no_op_inplace() {
    pureT tmp;
    pureT result=real_factory_6d(world());
    if (is_pure_no_op()) {
        return;
    } else if (is_pure()) {
        tmp= CompositeFactory<double, 6, 3>(world())
                .g12(get_operator().get_kernel())
                .ket(get_function());
        tmp.fill_tree();
        result=tmp;
    } else if (is_decomposed()) {
        MADNESS_CHECK_THROW(get_a().size()<3,"a.size not <3 in convert_to_pure_no_op_inplace");
        for (int i=0; i<get_a().size(); ++i) {
            if (is_op_decomposed()) {
                tmp= CompositeFactory<double, 6, 3>(world())
                        .g12(get_operator().get_kernel())
                        .particle1(get_a()[i])
                        .particle2(get_b()[i]);
            } else if (is_decomposed_no_op()) {
                tmp= CompositeFactory<double, 6, 3>(world())
                        .particle1(get_a()[i])
                        .particle2(get_b()[i]);
            }
            tmp.fill_tree();
            result+=tmp;
        }
    }
    component.reset(new TwoBodyFunctionPureComponent<T>(result));
};

std::vector<CCPairFunction> consolidate(const std::vector<CCPairFunction>& other) {

    std::vector<CCPairFunction> result;
    std::vector<real_function_6d> all_pure;
    for (auto& c : other) {
        if (c.is_pure_no_op()) all_pure.push_back(c.get_function());
        else result.push_back(c);
    }
    if (not all_pure.empty()) {
        for (std::size_t i=1; i<all_pure.size(); ++i) all_pure.front()+=all_pure[i];
        all_pure.front().truncate();
        result.push_back(CCPairFunction(all_pure.front()));
    }

    return result;
}

CCPairFunction multiply(const CCPairFunction& other, const real_function_3d& f, const std::array<int, 3>& v1) {
    auto a012=std::array<int,3>{0,1,2};
    auto a345=std::array<int,3>{3,4,5};
    int particle=-1;
    if (v1== a012) particle=0;
    if (v1== a345) particle=1;
    MADNESS_CHECK(particle==0 or particle==1);
    World& world=other.world();

    if (other.is_decomposed()) {
        if (particle == 0) {
            return CCPairFunction(other.get_operator_ptr(), f * other.get_a(), copy(world, other.get_b()));
        } else {
            return CCPairFunction(other.get_operator_ptr(), copy(world, other.get_a()), f * other.get_b());
        }
    } else if (other.is_pure()) {
        auto tmp=multiply(other.get_function(),f,particle+1);
        return CCPairFunction(other.get_operator_ptr(),tmp);
    } else  {
        MADNESS_EXCEPTION("confused CCPairFunction in multiply",1);
    }
};


/// multiplication with a 2-particle function
CCPairFunction& CCPairFunction::multiply_with_op_inplace(const std::shared_ptr<CCConvolutionOperator> op) {
    if (has_operator()) {
        auto newop=combine(get_operator_ptr(),op);
        reset_operator(newop);
    } else {
        reset_operator(op);
    }
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
    if (is_pure_no_op()) {
        result = pure().get_function().project_out(f.function,
                                                   particle - 1); // this needs 0 or 1 for particle but we give 1 or 2
    } else if (is_op_pure()) {
        MADNESS_EXCEPTION("implement CCPairFunction::project_out for op_pure",1);
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

CCPairFunction inner(const CCPairFunction& c1, const CCPairFunction& c2,
                       const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    auto v11=std::array<int,3>({std::get<0>(v1),std::get<1>(v1),std::get<2>(v1)});
    auto v22=std::array<int,3>({std::get<0>(v2),std::get<1>(v2),std::get<2>(v2)});

    return c1.partial_inner(c2,v11,v22);
}

std::vector<CCPairFunction> inner(const std::vector<CCPairFunction>& c1,
                                  const std::vector<CCPairFunction>& c2,
                                  const std::tuple<int,int,int> v1, const std::tuple<int,int,int> v2) {
    std::vector<CCPairFunction> result;
    for (const auto& cc1 : c1) {
        for (const auto& cc2 : c2) {
            print("inner of ",cc1.name(), cc2.name());
            result.push_back(inner(cc1,cc2,v1,v2));
        }
    }
    return result;
}

CCPairFunction CCPairFunction::partial_inner(const CCPairFunction& other,
                                               const std::array<int, 3>& v1,
                                               const std::array<int, 3>& v2) const {
    auto a012=std::array<int,3>{0,1,2};
    auto a345=std::array<int,3>{3,4,5};
    MADNESS_CHECK(v1==a012 or v1== a345);
    MADNESS_CHECK(v2==a012 or v2== a345);
    MADNESS_CHECK(not this->is_op_pure()); // not implemented yet
    MADNESS_CHECK(not other.is_op_pure()); // not implemented yet

    auto integration_index=[&a012,&a345](auto v) {return (v==a012) ? 0l : 1l;};
    auto remaining_index=[&integration_index](auto v) {return (integration_index(v)+1)%2;};

    CCPairFunction result;
    if (this->is_pure()) {
        if (other.is_pure()) {
            real_function_6d tmp=madness::innerXX<6>(this->get_function(),other.get_function(),v1,v2);
            return CCPairFunction(tmp);

        } else if (other.is_decomposed_no_op()) {
            // \int \sum_i f(1,2) a_i(1) b_i(3) d1  =  \sum_i b_i(3) \int a_i(1) f(1,2) d1
            vector_real_function_3d tmp;
            for (auto& a : other.get_vector(integration_index(v2))) {
//                tmp.push_back(innerXX<3>(this->get_function(),a,v1,a012));  // a012 is correct, referring to 3D function
                tmp.push_back(this->get_function().project_out(a,integration_index(v1)));
            }
            return CCPairFunction(tmp,other.get_vector(remaining_index(v2)));

        } else if (other.is_op_decomposed()) {

            // \int \sum_i h(1,2) f(1,3) c_i(1) d_i(3) d1
            //  = \sum_i d_i(3) \int h_c_i(1,2) f(1,3) d1
            //  = \sum_i d_i(3) H_i(3,2)
            const auto& h=this->pure().get_function();
            const auto& c=other.get_vector(integration_index(v2));
            const auto& d=other.get_vector(remaining_index(v2));
            auto& op=*(other.get_operator().get_op());
            op.particle()=integration_index(v1)+1;

            const vector_real_function_6d tmp=partial_mul(h,c,integration_index(v1)+1);
            auto H=apply(world(),op,tmp);
            real_function_6d result=real_factory_6d(world());
//            const vector_real_function_6d result=partial_mul(H,d,integration_index(v1)+1);
            for (int i=0; i<H.size(); ++i) {
                result+=multiply(H[i],d[i],integration_index(v1)+1);
            }
            return CCPairFunction(result);
        } else {
            MADNESS_EXCEPTION("confused CCPairfunction",1);
        }

    } else if (this->is_decomposed_no_op()) {
        if (other.is_pure()) {
            return other.partial_inner(*this,v2,v1);
        } else if (other.is_decomposed_no_op()) {
            // \int \sum_i a_i(1) b_i(2) \sum_j c_j(1) d_j(3) d1
            //    =  \sum_ij  <a_i|c_j>  b_i(2) d_j(3)
            //    =  \sum_i  b~_i(2) d~_i(3)        // SVD decomposition of S_ac
            Tensor<T> ovlp=matrix_inner(world(),this->get_vector(integration_index(v1)),other.get_vector(integration_index(v2)));
            Tensor< typename Tensor<T>::scalar_type > s;
            Tensor<T> U,VT;
            svd(ovlp,U,s,VT);
            for (int i=0; i<s.size(); ++i) U(_,i)*=s(i);
            auto left=transform(world(),this->get_vector(remaining_index(v1)),U);
            auto right=transform(world(),other.get_vector(remaining_index(v2)),transpose(VT));
            return CCPairFunction(left,right);

        } else if (other.is_op_decomposed()) {
            // \int \sum_ij a_i(1) b_i(2) f(1,3) c_j(1) d_j(3) d1
            //   = \sum_ij b_i(2) d_j(3) \int ac_ij(1) f(1,3) d1
            //   = \sum_i b_i(2) \sum_j d_j(3) g_ij(3)
            //   = \sum_i b_i(2) h_i(3)
            const auto& a=this->get_vector(integration_index(v1));
            const auto& b=this->get_vector(remaining_index(v1));
            const auto& c=other.get_vector(integration_index(v2));
            const auto& d=other.get_vector(remaining_index(v2));
            const auto& op=*(other.get_operator().get_op());
            std::decay_t<decltype(a)> h(a.size());        // /same type as a, without reference&
            for (int i=0; i<a.size(); ++i) {
                const auto ac=a[i]*c;
                const auto g=madness::apply(world(),op,ac);
                h[i]=dot(world(),d,g);
            }
            return CCPairFunction(b,h);
        } else {
            MADNESS_EXCEPTION("confused CCPairfunction",1);
        }

    } else if (this->is_op_decomposed()) {
        if (other.is_pure()) {
            return other.partial_inner(*this,v2,v1);
        } else if (other.is_decomposed_no_op()) {
            return other.partial_inner(*this,v2,v1);
        } else if (other.is_op_decomposed()) {
            if (this->is_convertible_to_pure_no_op()) {
                CCPairFunction tmp=copy(*this);
                tmp.convert_to_pure_no_op_inplace();
                return tmp.partial_inner(other,v1,v2);
            } else if (other.is_convertible_to_pure_no_op()) {
                CCPairFunction tmp=copy(other);
                tmp.convert_to_pure_no_op_inplace();
                return this->partial_inner(tmp,v1,v2);
            } else {
                MADNESS_EXCEPTION("no partial_inner for this combination: <op_decomposed|op_decomposed>",1);
            }
        } else {
            MADNESS_EXCEPTION("confused CCPairfunction",1);
        }
    } else {
        MADNESS_EXCEPTION("confused CCPairfunction",1);
    }
    return result;

}

real_function_3d CCPairFunction::partial_inner(const real_function_3d& f,
                                               const std::array<int, 3>& v1,
                                               const std::array<int, 3>& v2) const {
    auto a012=std::array<int,3>{0,1,2};
    auto a345=std::array<int,3>{3,4,5};
    MADNESS_CHECK(v2==a012 ); // only 3 dimension in f
    MADNESS_CHECK(v1==a012 or v1== a345); // 6 dimension in f
    MADNESS_CHECK(not this->is_op_pure()); // not implemented yet
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

/// compute the inner product of this and other

/// there are 4 possible components: pure/decomposed with and without operator, gives us 16 pair combinations..
double CCPairFunction::inner_internal(const CCPairFunction& other, const real_function_3d& R2) const {
    const CCPairFunction& f1=*this;
    const CCPairFunction& f2=other;

    double thresh=FunctionDefaults<6>::get_thresh()*0.1;

    double result = 0.0;
    if (f1.is_pure() and f2.is_pure()) {        // these are 4 combinations pure/pure
        pureT bra=f1.get_function();
        pureT ket=f2.get_function();
        // include the operator(s), if any
        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        real_function_6d tmp1;
        if (op) {
            if (R2.is_initialized()) {
                tmp1 = CompositeFactory<double, 6, 3>(world()).g12(op->get_kernel()).ket(ket).particle1(R2).particle2(R2);
            } else {
                tmp1 = CompositeFactory<double, 6, 3>(world()).g12(op->get_kernel()).ket(ket);
            }
        } else {
            if (R2.is_initialized()) {
                tmp1 = CompositeFactory<double, 6, 3>(world()).ket(ket).particle1(R2).particle2(R2);
            } else {
                tmp1 = CompositeFactory<double, 6, 3>(world()).ket(ket);
            }
        }
        result=inner(bra,tmp1);
    } else if (f1.is_pure() and f2.is_decomposed()) {       // with or without operator
        const vector_real_function_3d a = R2.is_initialized() ? R2 * f2.get_a() : copy(world(), f2.get_a());
        const vector_real_function_3d b = R2.is_initialized() ? R2 * f2.get_b() : copy(world(), f2.get_b());
        const pureT& bra=f1.get_function();

        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        if (op) {
            double bla=0.0;
            for (int i=0; i<a.size(); ++i) {
                real_function_6d tmp = CompositeFactory<double, 6, 3>(world()).g12(op->get_kernel()).particle1(a[i]).particle2(b[i]);
                bla += inner(bra, tmp);
            }
            result+=bla;
        } else { // no operators
            for (int i=0; i<a.size(); ++i) {
                real_function_6d tmp = CompositeFactory<double, 6, 3>(world()).particle1(a[i]).particle2(b[i]);
                result+=inner(bra,tmp);
            }
        }
    } else if (f1.is_decomposed() and f2.is_pure()) {     // with or without op
        result= f2.inner_internal(f1,R2);

    } else if (f1.is_decomposed() and f2.is_decomposed()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());

        const vector_real_function_3d& a1 = f1.get_a();
        const vector_real_function_3d& b1 = f1.get_b();
        const vector_real_function_3d a2 = R2.is_initialized() ?  R2* f2.get_a() : f2.get_a();
        const vector_real_function_3d b2 = R2.is_initialized() ?  R2* f2.get_b() : f2.get_b();


//        MADNESS_EXCEPTION("still to debug",1);
        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        if (not op) {
            // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
            result = (matrix_inner(world(), a1, a2)).trace(matrix_inner(world(),b1,b2));
        } else {
            // <a_i b_i | op | a_j b_j>  =  <a_i * a_j | op(b_i*b_j) >
            result=0.0;
            for (size_t i = 0; i < a1.size(); i++) {
                vector_real_function_3d aa = truncate(a1[i] * a2);
                vector_real_function_3d bb = truncate(b1[i] * b2);
                vector_real_function_3d aopx = (*op)(aa);
                result +=  inner(bb, aopx);
            }
        }
    } else MADNESS_EXCEPTION(
            ("CCPairFunction Overlap not supported for combination " + f1.name() + " and " + f2.name()).c_str(), 1) ;
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
//    print("apply projector on argument with terms",argument.size());
    if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) {
//        print("P->get_particle()",P->get_particle());
        MADNESS_CHECK_THROW(P->get_particle()==0 or P->get_particle()==1,"P Projector particle must be 0 or 1 in CCPairFunction");
    }
    if (auto Q=dynamic_cast<const QProjector<double,3>*>(&projector)) {
//        print("Q->get_particle()",Q->get_particle());
        MADNESS_CHECK_THROW(Q->get_particle()==0 or Q->get_particle()==1,"Q Projector particle must be 0 or 1 in CCPairFunction");
    }
    std::vector<CCPairFunction> result;
    for (const auto& pf : argument) {
        if (pf.is_pure()) {
            MADNESS_CHECK(not pf.has_operator()); // not yet implemented
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
//                CCTimer t(world,"SO block");
                // Q12 = 1 - O1 (1 - 1/2 O2) - O2 (1 - 1/2 O1)
//                print("entering SO block");
                QProjector<double,3> Q1(world,SO->bra1(),SO->ket1());
                Q1.set_particle(0);
                QProjector<double,3> Q2(world,SO->bra2(),SO->ket2());
                Q2.set_particle(1);

                Projector<double,3> O1(SO->bra1(),SO->ket1());
                O1.set_particle(0);
                Projector<double,3> O2(SO->bra2(),SO->ket2());
                O2.set_particle(1);

//                auto arg=std::vector<CCPairFunction>({pf});
//                auto o1arg=O1(arg);
//                auto o2arg=O2(arg);
//                auto o1o2arg=O1(o2arg);
//
//                result.push_back(pf);
//                for (auto& t: o1arg) result.push_back(-1.0*t);
//                for (auto& t: o2arg) result.push_back(-1.0*t);
//                for (auto& t: o1o2arg) result.push_back(t);


                auto tmp=Q1(Q2(std::vector<CCPairFunction>({pf})));
//                auto tmp=Q2(Q1(std::vector<CCPairFunction>({pf})));
//                print("result of SO");
//                for (auto& t: tmp) t.print_size();
                for (auto& t: tmp) result.push_back(t);
//                for (auto& t: result) t.print_size();

            } else if (auto P=dynamic_cast<const Projector<double,3>*>(&projector)) {
//                CCTimer t(world,"P block");
//                print("entering P block");
                std::vector<real_function_3d> tmp= zero_functions_compressed<double,3>(world,P->get_ket_vector().size());

                // per term a_i b_i:
                // P1 f |a b> = \sum_k |k(1)> |f_ak(2)*b(2)>
                for (std::size_t i=0; i<pf.get_a().size(); ++i) {
                    real_function_3d a=pf.get_a()[i];
                    real_function_3d b=pf.get_b()[i];
                    if (P->get_particle()==1) std::swap(a,b);

                    std::vector<real_function_3d> ka=a*P->get_bra_vector();
                    real_convolution_3d& op=*(pf.get_operator().get_op());
                    std::vector<real_function_3d> f_ka=apply(world,op,ka);
                    std::vector<real_function_3d> b_f_ka=f_ka*b;
                    tmp+=b_f_ka;
                }
                truncate(world,tmp);
//                print("size of tmp",tmp.size());

                if (P->get_particle()==0) result.push_back(CCPairFunction(P->get_ket_vector(),tmp));
                if (P->get_particle()==1) result.push_back(CCPairFunction(tmp,P->get_ket_vector()));
//                t.print();

            } else if (auto Q=dynamic_cast<const QProjector<double,3>*>(&projector)) {
//                CCTimer t(world,"Q block");
                // Q1 f12 |a_i b_i> = f12 |a_i b_i> - \sum_k |k(1) a_i(2)*f_(kb_i)(2) >
                result.push_back(pf);
//                print("entering Q block");
                // reuse the projector code above
                std::vector<CCPairFunction> tmp=madness::apply(Q->get_P_projector(),std::vector<CCPairFunction>(1,pf));
//                for (auto& t : tmp) t.print_size();
                for (auto& t : tmp) {
                    t*=-1.0;
                    result.push_back(t);
                }
//                t.print();

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

/// apply the operator to the argument

/// the operator is applied to one particle only, the other one is left untouched
/// note the ordering of the particles, cf the corresponding comment in mra.h
/// op.particle==1 :  op(f(x,y)) = op(x,x') f(x',y) = result(x,y);
/// op.particle==2 :  op(f(x,y)) = op(y,y') f(x,y') = result(y,x);
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction> apply(const SeparatedConvolution<T,NDIM>& op, const std::vector<CCPairFunction>& argument) {
    if (argument.size()==0) return argument;
    World& world=argument.front().world();
    std::vector<CCPairFunction> result;
    timer t(world);
    for (const auto& arg : argument) {
        bool convert_to_pure=(arg.has_operator() or arg.is_pure());

        if (convert_to_pure) {
            auto tmp=arg.to_pure().get_function();
            tmp=op(tmp);

            // !! confusing ordering of the result variables!!
            if (op.particle()==2) tmp=swap_particles(tmp);
            result.push_back(CCPairFunction(tmp));

        } else if (arg.is_decomposed_no_op()) {
            MADNESS_CHECK(op.particle()==1 or op.particle()==2);
            if (op.particle()==1) {
                auto tmp= madness::apply(world,op,arg.get_a());
                result.push_back(CCPairFunction(tmp,arg.get_b()));
            } else if (op.particle()==2) {
                auto tmp= madness::apply(world,op,arg.get_b());
                result.push_back(CCPairFunction(tmp,arg.get_a()));
            }

        } else {
            MADNESS_CHECK_THROW(false,"confused type in apply(CCPairFunction)");
        }
//        t.tag("applying op to CCPairFunction "+arg.name());
    }

    return result;
}


template std::vector<CCPairFunction> apply(const SeparatedConvolution<double,3>& op, const std::vector<CCPairFunction>& argument);

} // namespace madness