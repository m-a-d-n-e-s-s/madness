//
// Created by Florian Bischoff on 6/27/22.
//

#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/projector.h>
#include<madness/chem/lowrankfunction.h>
#include<madness/mra/operator.h>

using namespace madness;

namespace madness {

template<typename T, std::size_t NDIM>
madness::CCPairFunction<T,NDIM>
CCPairFunction<T,NDIM>::invert_sign() {
    (*this)*=(-1.0);
    return *this;
}

template<typename T, std::size_t NDIM>
bool CCPairFunction<T,NDIM>::is_convertible_to_pure_no_op() const {
    if (has_operator()) {
        const auto type=get_operator().type();
        if (not (type==OpType::OT_SLATER or type==OpType::OT_F12)) return false;
    }
    return true;
};


template<typename T, std::size_t NDIM>
void CCPairFunction<T,NDIM>::convert_to_pure_no_op_inplace() {
    pureT result;
    if (is_pure_no_op()) {
        return;
    } else if (is_pure()) {
        result= CompositeFactory<T, NDIM, LDIM>(world())
                .g12(get_operator().get_kernel())
                .ket(get_function());
    } else if (is_decomposed_no_op()) {
        result= CompositeFactory<T, NDIM, LDIM>(world())
                .particle1(get_a())
                .particle2(get_b());
    } else if (is_op_decomposed()) {
        result= CompositeFactory<T, NDIM, LDIM>(world())
                .g12(get_operator().get_kernel())
                .particle1(get_a())
                .particle2(get_b());
    } else {
        MADNESS_EXCEPTION("error in convert_to_pure_no_op_inplace",1);
    }
    result.fill_tree();
    result.truncate(FunctionDefaults<NDIM>::get_thresh()*0.1);
    component.reset(new TwoBodyFunctionPureComponent<T,NDIM>(result));
};

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::op_pure_to_pure(const std::vector<CCPairFunction<T,NDIM>>& other) {
    World& world=other.front().world();
    std::vector<CCPairFunction<T,NDIM>> result;
    Function<T,NDIM> pure=FunctionFactory<T,NDIM>(world);
    for (const auto& c : other) {
        if (c.is_pure_no_op()) {
            pure+=c.get_function();
        } else if (c.is_op_pure()) {
            Function<T,NDIM> tmp=CompositeFactory<T,NDIM,LDIM>(world).g12(c.get_operator().get_kernel()).ket(c.get_function());
            tmp.fill_tree();
            pure+=tmp;
        } else if (c.is_decomposed()) {
            result.push_back(c);
        }
    }
    if (pure.is_initialized()) {
        pure.truncate(FunctionDefaults<NDIM>::get_thresh()*0.1);
        result.push_back(CCPairFunction<T,NDIM>(pure));
    }
    return result;
}

/// turn decomposed functions with operator into decomposed functions using LowRankFunction
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::op_dec_to_dec(const std::vector<CCPairFunction<T,NDIM>>& other,
                                                                          const std::vector<Vector<double,CCPairFunction<T,NDIM>::LDIM>>& centers) {
    LowRankFunctionParameters lrparameters;
    lrparameters.set_derived_value("tol",1.e-10);
    auto builder = LowRankFunctionFactory<T,NDIM>(lrparameters,centers);
//    builder.set_volume_element(3.e-2);
    if (other.front().world().rank()==0) {
        builder.parameters.print("lrparameters");
        print("centers",centers);
    }
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& c : other) {
        if (c.is_op_decomposed()) {
            LRFunctorF12<T,NDIM> functor(c.get_operator_ptr()->get_op(),c.get_a(),c.get_b());
            LowRankFunction<T,NDIM> tmp=builder.project(functor);
//            double l2error=tmp.l2error(functor);
            tmp.optimize(functor);
            result.push_back(CCPairFunction<T,NDIM>(tmp.get_g(),tmp.get_h()));
        } else {
            result.push_back(c);
        }
    }
    return result;
}

/// turn decomposed functions with operator into pure functions
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::dec_to_pure(const std::vector<CCPairFunction<T,NDIM>>& other) {
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& c : other) {
        if (c.is_decomposed_no_op()) {
            CCPairFunction<T,NDIM> tmp=copy(c);
            tmp.convert_to_pure_no_op_inplace();
            result.push_back(tmp);
        } else {
            result.push_back(c);
        }
    }
    return result;
}


/// turn decomposed functions with operator into pure functions
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::op_dec_to_pure(const std::vector<CCPairFunction<T,NDIM>>& other) {
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& c : other) {
        if (c.is_op_decomposed()) {
            CCPairFunction<T,NDIM> tmp=copy(c);
            tmp.convert_to_pure_no_op_inplace();
            result.push_back(tmp);
        } else {
            result.push_back(c);
        }
    }
    return result;

}

/// turn decomposed functions with operator into decomposed functions using LowRankFunction
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::remove_linearly_dependent_terms(const std::vector<CCPairFunction<T,NDIM>>& other,
    double thresh) {
    if (thresh<0.0) thresh=FunctionDefaults<3>::get_thresh()*0.1;
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& c : other) {
        if (c.is_pure()) result.push_back(c);
        else if (c.is_decomposed()) {

            LowRankFunction<T,NDIM> lrf(c.get_a(),c.get_b(),thresh,"canonical");
            lrf.reorthonormalize();
            result.push_back(CCPairFunction<T,NDIM>(c.get_operator_ptr(),lrf.get_g(),lrf.get_h()));

        } else {
            MADNESS_EXCEPTION("you should not be here",1);
        }
    }

    return result;
}

/// collect all terms with of similiar type: pure, op_pure, decomposed, op_decomposed
template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::collect_same_types(const std::vector<CCPairFunction<T, NDIM>>& other) {

    if (other.size()==0) return other;
    if (is_collected(other)) return other;

    World& world=other.front().world();

    /// vector includes OT_ONE, meaning no operator
    std::vector<std::vector<Function<T,NDIM>>> op_pure(OT_SIZE);
    std::vector<std::vector<Function<T,LDIM>>> op_decomposed_a(OT_SIZE);
    std::vector<std::vector<Function<T,LDIM>>> op_decomposed_b(OT_SIZE);
    std::vector<std::shared_ptr<CCConvolutionOperator<T,LDIM>>> ops(OT_SIZE);

    // collect terms of the same type
    for (const auto& c : other) {
        int iop= (c.has_operator()) ? int(c.get_operator().type()) : OT_ONE;
        ops[iop]=c.get_operator_ptr();
        if (c.is_decomposed()) {
            op_decomposed_a[iop]=append(op_decomposed_a[iop],c.get_a());
            op_decomposed_b[iop]=append(op_decomposed_b[iop],c.get_b());
        } else if (c.is_pure()) {
            op_pure[iop].push_back(c.get_function());
        }
    }

    std::vector<CCPairFunction<T,NDIM>> result;

    // accumulate all terms of the same type
    for (int opint=OT_ONE; opint<OT_SIZE; ++opint) {
        if (op_pure[opint].size()>0) {
            auto op=ops[opint];
            if (op_pure[opint].size()>1) {
                // Function<T,NDIM> tmp=CompositeFactory<T,NDIM,LDIM>(world).ket(op_pure[opint]);
                // tmp.fill_tree();
                Tensor<double> c(op_pure[opint].size(),1);
                c=1.0;
                Function<T,NDIM> tmp=transform_reconstructed(world, op_pure[opint],c,true)[0];
                result.push_back(CCPairFunction<T,NDIM>(op,tmp));
            } else {
                MADNESS_CHECK_THROW(op_pure[opint].size()==1,"op_pure[opint].size()!=1");
                result.push_back(CCPairFunction<T,NDIM>(op,copy(op_pure[opint].front())));
            }
        }
        if (op_decomposed_a[opint].size()>0) {
            result.push_back(CCPairFunction<T,NDIM>(ops[opint],op_decomposed_a[opint],op_decomposed_b[opint]));
        }
    }
    return result;

}
template<typename T, std::size_t NDIM>
bool CCPairFunction<T,NDIM>::is_collected(const std::vector<CCPairFunction<T,NDIM>>& other) {

    // simply count the occurence of each term
    std::map<int,int> counter;
    for (const auto& c : other) {
        int index=0;
        if (c.has_operator()) index=int(c.get_operator().type()) *100;
        if (c.is_op_pure()) index+=1;
        if (c.is_op_decomposed()) index+=2;
        if (c.is_pure_no_op()) index+=3;
        if (c.is_decomposed_no_op()) index+=4;
        counter[index]++;
    }
    for (const auto& c : counter) if (c.second>1) return false;
    return true;
}

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::consolidate(const std::vector<CCPairFunction<T,NDIM>>& other,
                                                const std::vector<std::string>& options,
                                                const std::vector<Vector<double,CCPairFunction<T,NDIM>::LDIM>>& centers) const {

    // convert op_pure functions to pure
    bool op_pure_to_pure=find(options.begin(),options.end(),"op_pure_to_pure")!=options.end();
    // convert op_dec functions to dec (via LowRankFunctions
    bool op_dec_to_dec=find(options.begin(),options.end(),"op_dec_to_dec")!=options.end();
    // convert op_dec functions to pure (via fill_tree)
    bool op_dec_to_pure=find(options.begin(),options.end(),"op_dec_to_pure")!=options.end();
    // convert dec functions to pure (via hartree product)
    bool dec_to_pure=find(options.begin(),options.end(),"dec_to_pure")!=options.end();
    // reorthogonalize decomposed functions and op_decomposed functions
    bool lindep=find(options.begin(),options.end(),"remove_lindep")!=options.end();

    // always collect all terms of the same type
    auto result= is_collected(other) ? other : collect_same_types(other);
    if (lindep) result=CCPairFunction<T,NDIM>::remove_linearly_dependent_terms(result);

    if (op_dec_to_dec) result=CCPairFunction<T,NDIM>::op_dec_to_dec(result,centers);
    if (op_dec_to_pure) result=CCPairFunction<T,NDIM>::op_dec_to_pure(result);
    if (dec_to_pure) result=CCPairFunction<T,NDIM>::dec_to_pure(result);
    if (op_pure_to_pure) result=CCPairFunction<T,NDIM>::op_pure_to_pure(result);

    if (not is_collected(result)) result=collect_same_types(result);


    return result;
}

/// multiplication with a 2-particle function
template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM>& CCPairFunction<T,NDIM>::multiply_with_op_inplace(const std::shared_ptr<CCConvolutionOperator<T,CCPairFunction<T,NDIM>::LDIM>> op) {
    if (has_operator()) {
        auto newop=combine(get_operator_ptr(),op);
        reset_operator(newop);
    } else {
        reset_operator(op);
    }
    return *this;
}

template<typename T, std::size_t NDIM>
double
CCPairFunction<T,NDIM>::make_xy_u(const CCFunction<T,LDIM>& xx, const CCFunction<T,LDIM>& yy) const {
    CCPairFunction<T,NDIM> bra(xx.function,yy.function);
    return inner(bra,*this);
    T result = 0.0;
    if (is_pure()) {
        World& world=xx.function.world();
        Function<T,NDIM> ij = CompositeFactory<double, NDIM, LDIM>(world).particle1(madness::copy(xx.function))
                .particle2( madness::copy(yy.function));
        result = inner(pure().get_function(), ij);
    } else if (is_decomposed_no_op()) {
        for (size_t i = 0; i < get_a().size(); i++)
            result += (xx.function.inner(get_a()[i])) * (yy.function.inner(get_b()[i]));
    } else if (is_op_decomposed()) {
        const CCConvolutionOperator<T,LDIM>& op=*decomposed().get_operator_ptr();
        result = yy.function.inner(op(xx, get_a()[0]) * get_b()[0]);
    }
    return result;
}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM>
CCPairFunction<T,NDIM>::project_out(const CCFunction<T,LDIM>& f, const size_t particle) const {
    MADNESS_ASSERT(particle == 1 or particle == 2);
    Function<T,LDIM> result;
    if (is_pure_no_op()) {
        result = pure().get_function().project_out(f.function,
                                                   particle - 1); // this needs 0 or 1 for particle but we give 1 or 2
    } else if (is_op_pure()) {
        MADNESS_EXCEPTION("implement CCPairFunction<T,NDIM>::project_out for op_pure",1);
    } else if (is_decomposed_no_op()) {
        result = project_out_decomposed(f.function, particle);
    } else if (is_op_decomposed()) {
        result = project_out_op_decomposed(f, particle);
    }
    if (not result.is_initialized()) MADNESS_EXCEPTION("Result of project out on CCPairFunction<T,NDIM> was not initialized",
                                                       1);
    return result;
}

// result is: <x|op12|f>_particle
template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM>
CCPairFunction<T,NDIM>::dirac_convolution(const CCFunction<T,LDIM>& x, const CCConvolutionOperator<T,CCPairFunction<T,NDIM>::LDIM>& op, const size_t particle) const {
    Function<T,CCPairFunction<T,NDIM>::LDIM> result;
    if (is_pure()) {
        result = op(x, pure().get_function(), particle);
    } else if (is_decomposed_no_op()) {
        result = dirac_convolution_decomposed(x, op, particle);
    } else {
        MADNESS_EXCEPTION("op_decomposed dirac convolution not yet implemented", 1);
    }
    return result;
}

template<typename T, std::size_t NDIM>
CCPairFunction<T,NDIM> CCPairFunction<T,NDIM>::partial_inner(const CCPairFunction<T,NDIM>& other,
                                               const std::array<int, CCPairFunction<T,NDIM>::LDIM>& v1,
                                               const std::array<int, CCPairFunction<T,NDIM>::LDIM>& v2) const {
    auto a012=std::array<int,LDIM>();
    auto a345=std::array<int,LDIM>();
    for (size_t i=0; i<LDIM; ++i) {
        a012[i]=i;
        a345[i]=i+LDIM;
    }
    MADNESS_CHECK(v1==a012 or v1== a345);
    MADNESS_CHECK(v2==a012 or v2== a345);
    MADNESS_CHECK(not this->is_op_pure()); // not implemented yet
    MADNESS_CHECK(not other.is_op_pure()); // not implemented yet

    auto integration_index=[&a012](auto v) {return (v==a012) ? 0l : 1l;};
    auto remaining_index=[&integration_index](auto v) {return (integration_index(v)+1)%2;};

    CCPairFunction<T,NDIM> result;
    if (this->is_pure()) {
        if (other.is_pure()) {
            Function<T,NDIM> tmp=madness::innerXX<NDIM>(this->get_function(),other.get_function(),v1,v2);
            return CCPairFunction<T,NDIM>(tmp);

        } else if (other.is_decomposed_no_op()) {
            // \int \sum_i f(1,2) a_i(1) b_i(3) d1  =  \sum_i b_i(3) \int a_i(1) f(1,2) d1
            std::vector<Function<T,LDIM>> tmp;
            auto avec=other.get_vector(integration_index(v2));
            change_tree_state(avec,redundant);
            for (auto& a : other.get_vector(integration_index(v2))) {
//                tmp.push_back(innerXX<3>(this->get_function(),a,v1,a012));  // a012 is correct, referring to 3D function
                tmp.push_back(this->get_function().project_out(a,integration_index(v1)));
            }
            return CCPairFunction<T,NDIM>(tmp,other.get_vector(remaining_index(v2)));

        } else if (other.is_op_decomposed()) {

            // \int \sum_i h(1,2) f(1,3) c_i(1) d_i(3) d1
            //  = \sum_i d_i(3) \int h_c_i(1,2) f(1,3) d1
            //  = \sum_i d_i(3) H_i(3,2)
            const auto& h=this->pure().get_function();
            const auto& c=other.get_vector(integration_index(v2));
            const auto& d=other.get_vector(remaining_index(v2));
            auto& op=*(other.get_operator().get_op());
            op.particle()=integration_index(v1)+1;

            const std::vector<Function<T,NDIM>> tmp=partial_mul(h,c,integration_index(v1)+1);
            auto H=madness::apply(world(),op,tmp);
            Function<T,NDIM> result=FunctionFactory<T,NDIM>(world());
//            const vector_real_function_6d result=partial_mul(H,d,integration_index(v1)+1);
            for (size_t i=0; i<H.size(); ++i) {
                result+=multiply(H[i],d[i],integration_index(v1)+1);
            }
            return CCPairFunction<T,NDIM>(result);
        } else {
            MADNESS_EXCEPTION("confused CCPairFunction<T,NDIM>",1);
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
            return CCPairFunction<T,NDIM>(left,right);

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
            for (size_t i=0; i<a.size(); ++i) {
                const auto ac=a[i]*c;
                const auto g=madness::apply(world(),op,ac);
                h[i]=dot(world(),d,g);
            }
            return CCPairFunction<T,NDIM>(b,h);
        } else {
            MADNESS_EXCEPTION("confused CCPairFunction<T,NDIM>",1);
        }

    } else if (this->is_op_decomposed()) {
        if (other.is_pure()) {
            return other.partial_inner(*this,v2,v1);
        } else if (other.is_decomposed_no_op()) {
            return other.partial_inner(*this,v2,v1);
        } else if (other.is_op_decomposed()) {
            if (this->is_convertible_to_pure_no_op()) {
                CCPairFunction<T,NDIM> tmp=copy(*this);
                tmp.convert_to_pure_no_op_inplace();
                return tmp.partial_inner(other,v1,v2);
            } else if (other.is_convertible_to_pure_no_op()) {
                CCPairFunction<T,NDIM> tmp=copy(other);
                tmp.convert_to_pure_no_op_inplace();
                return this->partial_inner(tmp,v1,v2);
            } else {
                MADNESS_EXCEPTION("no partial_inner for this combination: <op_decomposed|op_decomposed>",1);
            }
        } else {
            MADNESS_EXCEPTION("confused CCPairFunction<T,NDIM>",1);
        }
    } else {
        MADNESS_EXCEPTION("confused CCPairFunction<T,NDIM>",1);
    }
    return result;

}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM> CCPairFunction<T,NDIM>::partial_inner(
        const Function<T,CCPairFunction<T,NDIM>::LDIM>& f,
        const std::array<int, CCPairFunction<T,NDIM>::LDIM>& v1,
        const std::array<int, CCPairFunction<T,NDIM>::LDIM >& v2) const {
//    auto a012=std::array<int,LDIM>{0,1,2};
//    auto a345=std::array<int,LDIM>{3,4,5};
    auto a012=std::array<int,LDIM>();
    auto a345=std::array<int,LDIM>();
    for (size_t i=0; i<LDIM; ++i) {
        a012[i]=i;
        a345[i]=i+LDIM;
    }
    MADNESS_CHECK(v2==a012 ); // only 3 dimension in f
    MADNESS_CHECK(v1==a012 or v1== a345); // 6 dimension in f
    MADNESS_CHECK(not this->is_op_pure()); // not implemented yet
    int particle=-1;
    if (v1== a012) particle=0;
    if (v1== a345) particle=1;

    Function<T,CCPairFunction<T,NDIM>::LDIM> result;

    if (is_pure()) {
        result = pure().get_function().project_out(f, particle);
    } else if (is_decomposed_no_op()) {
        result = project_out_decomposed(f, particle+1);
    } else if (is_op_decomposed()) {
        result = project_out_op_decomposed(f, particle+1);
    } else {
        MADNESS_EXCEPTION("confused state in CCPairFunction<T,NDIM>::partial_inner",1);
    }
    return result;
}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM> CCPairFunction<T,NDIM>::project_out_decomposed(
        const Function<T,LDIM>& f, const size_t particle) const {
    World& world=f.world();
    Function<T,LDIM> result = FunctionFactory<T,LDIM>(world);
    const std::pair<std::vector<Function<T,LDIM>>, std::vector<Function<T,LDIM>>> decompf = assign_particles(particle);
    Tensor<double> c = inner(world, f, decompf.first);
    for (size_t i = 0; i < get_a().size(); i++) result += c(i) * decompf.second[i];
    return result;
}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM> CCPairFunction<T,NDIM>::project_out_op_decomposed(const CCFunction<T,LDIM>& f, const size_t particle) const {
    World& world=f.get().world();
    const CCConvolutionOperator<T,LDIM>& op=*decomposed().get_operator_ptr();
    if (particle == 1) {
//        return op(f, get_a()[0]) * get_b()[0];
        // result(2) = < f(1) | op(1,2) | a_i(1) b_i(2) >
        return sum(world,mul(world,op(f.f()* get_a()),get_b()));
    } else if (particle == 2) {
//        return op(f, get_b()[0]) * get_a()[0];
        return sum(world,mul(world,op(f.f()* get_b()),get_a()));
    } else {
        MADNESS_EXCEPTION("project_out_op_decomposed: particle must be 1 or 2", 1);
        return FunctionFactory<T,LDIM>(world);
    }
}

template<typename T, std::size_t NDIM>
Function<T,CCPairFunction<T,NDIM>::LDIM> CCPairFunction<T,NDIM>::dirac_convolution_decomposed(const CCFunction<T,LDIM>& bra,
                                                                      const CCConvolutionOperator<T,LDIM>& op,
                                                              const size_t particle) const {
    World& world=bra.function.world();
    const std::pair<std::vector<Function<T,LDIM>>, std::vector<Function<T,LDIM>>> f = assign_particles(particle);
    const std::vector<Function<T,LDIM>> braa = mul(world, bra.function, f.first);
    const std::vector<Function<T,LDIM>> braga = op(braa);
    Function<T,LDIM> result = FunctionFactory<T,LDIM>(world);
    for (size_t i = 0; i < braga.size(); i++) result += braga[i] * f.second[i];
    return result;
}


template<typename T, std::size_t NDIM>
const std::pair<std::vector<Function<T,CCPairFunction<T,NDIM>::LDIM>>, std::vector<Function<T,CCPairFunction<T,NDIM>::LDIM>>>
CCPairFunction<T,NDIM>::assign_particles(const size_t particle) const {
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
template<typename T, std::size_t NDIM>
double CCPairFunction<T,NDIM>::inner_internal(const CCPairFunction<T,NDIM>& other, const Function<T,LDIM>& R2) const {
    const CCPairFunction<T,NDIM>& f1=*this;
    const CCPairFunction<T,NDIM>& f2=other;

    double result = 0.0;
    if (f1.is_pure() and f2.is_pure()) {        // these are 4 combinations pure/pure
        pureT bra=f1.get_function();
        pureT ket=f2.get_function();
        // include the operator(s), if any
        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        Function<T,NDIM> tmp1;
        if (op) {
            if (R2.is_initialized()) {
                tmp1 = CompositeFactory<T,NDIM,LDIM>(world()).g12(op->get_kernel()).ket(ket).particle1(R2).particle2(R2);
            } else {
                tmp1 = CompositeFactory<T,NDIM,LDIM>(world()).g12(op->get_kernel()).ket(ket);
            }
        } else {
            if (R2.is_initialized()) {
                tmp1 = CompositeFactory<T,NDIM,LDIM>(world()).ket(ket).particle1(R2).particle2(R2);
            } else {
                tmp1 = CompositeFactory<T,NDIM,LDIM>(world()).ket(ket);
            }
        }
        result=inner(bra,tmp1);
    } else if (f1.is_pure() and f2.is_decomposed()) {       // with or without operator
        const std::vector<Function<T,LDIM>> a = R2.is_initialized() ? R2 * f2.get_a() : copy(world(), f2.get_a());
        const std::vector<Function<T,LDIM>> b = R2.is_initialized() ? R2 * f2.get_b() : copy(world(), f2.get_b());
        const pureT& bra=f1.get_function();

        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        if (op) {
            double bla=0.0;
            for (size_t i=0; i<a.size(); ++i) {
                Function<T,NDIM> tmp = CompositeFactory<T,NDIM,LDIM>(world()).g12(op->get_kernel()).particle1(a[i]).particle2(b[i]);
                bla += inner(bra, tmp);
            }
            result+=bla;
        } else { // no operators
            for (size_t i=0; i<a.size(); ++i) {
                Function<T,NDIM> tmp = CompositeFactory<T,NDIM,LDIM>(world()).particle1(a[i]).particle2(b[i]);
                result+=inner(bra,tmp);
            }
        }
    } else if (f1.is_decomposed() and f2.is_pure()) {     // with or without op
        result= f2.inner_internal(f1,R2);

    } else if (f1.is_decomposed() and f2.is_decomposed()) {
        MADNESS_ASSERT(f1.get_a().size() == f1.get_b().size());
        MADNESS_ASSERT(f2.get_a().size() == f2.get_b().size());

        const std::vector<Function<T,LDIM>>& a1 = f1.get_a();
        const std::vector<Function<T,LDIM>>& b1 = f1.get_b();
        const std::vector<Function<T,LDIM>> a2 = R2.is_initialized() ?  R2* f2.get_a() : f2.get_a();
        const std::vector<Function<T,LDIM>> b2 = R2.is_initialized() ?  R2* f2.get_b() : f2.get_b();


//        MADNESS_EXCEPTION("still to debug",1);
        auto op=combine(f1.get_operator_ptr(),f2.get_operator_ptr());
        if (not op) {
            // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
            result = (matrix_inner(world(), a1, a2)).trace(matrix_inner(world(),b1,b2));
        } else {
            // <a_i b_i | op | a_j b_j>  =  <a_i * a_j | op(b_i*b_j) >
            result=0.0;
            for (size_t i = 0; i < a1.size(); i++) {
                std::vector<Function<T,LDIM>> aa = truncate(a1[i] * a2);
                std::vector<Function<T,LDIM>> bb = truncate(b1[i] * b2);
                std::vector<Function<T,LDIM>> aopx = (*op)(aa);
                result +=  inner(bb, aopx);
            }
        }
    } else MADNESS_EXCEPTION(
            ("CCPairFunction<T,NDIM> Overlap not supported for combination " + f1.name() + " and " + f2.name()).c_str(), 1) ;
    return result;
}

template<typename T, std::size_t NDIM>
std::vector<CCPairFunction<T,NDIM>> CCPairFunction<T,NDIM>::apply(const ProjectorBase& projector, const std::vector<CCPairFunction<T,NDIM>>& argument) {
    if (argument.size()==0) return argument;
    World& world=argument.front().world();
    constexpr std::size_t LDIM=CCPairFunction<T,NDIM>::LDIM;
//    print("apply projector on argument with terms",argument.size());
    if (auto P=dynamic_cast<const Projector<double,LDIM>*>(&projector)) {
        MADNESS_CHECK_THROW(P->get_particle()==0 or P->get_particle()==1,"P Projector particle must be 0 or 1 in CCPairFunction<T,NDIM>");
    }
    if (auto Q=dynamic_cast<const QProjector<double,LDIM>*>(&projector)) {
        MADNESS_CHECK_THROW(Q->get_particle()==0 or Q->get_particle()==1,"Q Projector particle must be 0 or 1 in CCPairFunction<T,NDIM>");
    }
    std::vector<CCPairFunction<T,NDIM>> result;
    for (const auto& pf : argument) {
        if (pf.is_pure()) {
            MADNESS_CHECK(not pf.has_operator()); // not yet implemented
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,LDIM>*>(&projector)) {
                auto tmp=(*SO)(pf.get_function());
                auto tmp2=CCPairFunction<T,NDIM>(tmp);
                result.push_back(tmp2);
            } else if (auto P=dynamic_cast<const Projector<double,LDIM>*>(&projector)) {
                // result.push_back(CCPairFunction<T,NDIM>((*P)(pf.get_function())));
                auto [left,right]=P->get_vectors_for_outer_product(pf.get_function());
                result.push_back(CCPairFunction<T,NDIM>(left,right));


            } else if (auto Q=dynamic_cast<const QProjector<double,LDIM>*>(&projector)) {
                // result.push_back(CCPairFunction<T,NDIM>((*Q)(pf.get_function())));
                result.push_back(pf);
                result.push_back(-1.0*Q->get_P_projector()(pf));

            } else {
                MADNESS_EXCEPTION("CCPairFunction<T,NDIM>: unknown projector type",1);
            }
        } else if (pf.is_decomposed_no_op()) {  // pair function is sum_i | a_i b_i >
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,LDIM>*>(&projector)) {
                // Q12 | kl > = (1-O1)(1-O2) |kl> = |(1-O1)k (1-O2)l>
                QProjector<double,LDIM> Q1(SO->bra1(),SO->ket1());
                QProjector<double,LDIM> Q2(SO->bra2(),SO->ket2());
                result.push_back(CCPairFunction<T,NDIM>(Q1(pf.get_a()),Q2(pf.get_b())));

            } else if (auto P=dynamic_cast<const Projector<double,LDIM>*>(&projector)) {
                // P1 | kl > = P1 |kl> = |P1 k l>
                if (P->get_particle()==0) result.push_back(CCPairFunction<T,NDIM>((*P)(pf.get_a()),pf.get_b()));
                // P2 | kl > = P2 |kl> = |k P2 l>
                if (P->get_particle()==1) result.push_back(CCPairFunction<T,NDIM>(pf.get_a(),(*P)(pf.get_b())));

            } else if (auto Q=dynamic_cast<const QProjector<double,LDIM>*>(&projector)) {
                // Q1 | kl > = Q1 |kl> = |Q1 k l>
                if (Q->get_particle()==0) result.push_back(CCPairFunction<T,NDIM>((*Q)(pf.get_a()),pf.get_b()));
                // P2 | kl > = Q2 |kl> = |k Q2 l>
                if (Q->get_particle()==1) result.push_back(CCPairFunction<T,NDIM>(pf.get_a(),(*Q)(pf.get_b())));
            } else {
                MADNESS_EXCEPTION("CCPairFunction<T,NDIM>: unknown projector type",1);
            }
        } else if (pf.is_op_decomposed()) {
            if (auto SO=dynamic_cast<const StrongOrthogonalityProjector<double,LDIM>*>(&projector)) {
//                CCTimer t(world,"SO block");
                // Q12 = 1 - O1 (1 - 1/2 O2) - O2 (1 - 1/2 O1)
//                print("entering SO block");
                QProjector<double,LDIM> Q1(SO->bra1(),SO->ket1());
                Q1.set_particle(0);
                QProjector<double,LDIM> Q2(SO->bra2(),SO->ket2());
                Q2.set_particle(1);

                Projector<double,LDIM> O1(SO->bra1(),SO->ket1());
                O1.set_particle(0);
                Projector<double,LDIM> O2(SO->bra2(),SO->ket2());
                O2.set_particle(1);

//                auto arg=std::vector<CCPairFunction<T,NDIM>>({pf});
//                auto o1arg=O1(arg);
//                auto o2arg=O2(arg);
//                auto o1o2arg=O1(o2arg);
//
//                result.push_back(pf);
//                for (auto& t: o1arg) result.push_back(-1.0*t);
//                for (auto& t: o2arg) result.push_back(-1.0*t);
//                for (auto& t: o1o2arg) result.push_back(t);


                auto tmp=Q1(Q2(std::vector<CCPairFunction<T,NDIM>>({pf})));
//                auto tmp=Q2(Q1(std::vector<CCPairFunction<T,NDIM>>({pf})));
//                print("result of SO");
//                for (auto& t: tmp) t.print_size();
                for (auto& t: tmp) result.push_back(t);
//                for (auto& t: result) t.print_size();

            } else if (auto P=dynamic_cast<const Projector<double,LDIM>*>(&projector)) {
//                CCTimer t(world,"P block");
//                print("entering P block");
                std::vector<Function<T,LDIM>> tmp= zero_functions_compressed<double,LDIM>(world,P->get_ket_vector().size());

                // per term a_i b_i:
                // P1 f |a b> = \sum_k |k(1)> |f_ak(2)*b(2)>
                for (std::size_t i=0; i<pf.get_a().size(); ++i) {
                    Function<T,LDIM> a=pf.get_a()[i];
                    Function<T,LDIM> b=pf.get_b()[i];
                    if (P->get_particle()==1) std::swap(a,b);

                    std::vector<Function<T,LDIM>> ka=a*P->get_bra_vector();
                    SeparatedConvolution<T,LDIM>& op=*(pf.get_operator().get_op());
                    std::vector<Function<T,LDIM>> f_ka=madness::apply(world,op,ka);
                    std::vector<Function<T,LDIM>> b_f_ka=f_ka*b;
                    tmp+=b_f_ka;
                }
                truncate(world,tmp);
//                print("size of tmp",tmp.size());

                if (P->get_particle()==0) result.push_back(CCPairFunction<T,NDIM>(P->get_ket_vector(),tmp));
                if (P->get_particle()==1) result.push_back(CCPairFunction<T,NDIM>(tmp,P->get_ket_vector()));
//                t.print();

            } else if (auto Q=dynamic_cast<const QProjector<double,LDIM>*>(&projector)) {
//                CCTimer t(world,"Q block");
                // Q1 f12 |a_i b_i> = f12 |a_i b_i> - \sum_k |k(1) a_i(2)*f_(kb_i)(2) >
                result.push_back(pf);
//                print("entering Q block");
                // reuse the projector code above
                std::vector<CCPairFunction<T,NDIM>> tmp=madness::apply(Q->get_P_projector(),std::vector<CCPairFunction<T,NDIM>>(1,pf));
//                for (auto& t : tmp) t.print_size();
                for (auto& t : tmp) {
                    t*=-1.0;
                    result.push_back(t);
                }
//                t.print();

            } else {
                MADNESS_EXCEPTION("CCPairFunction<T,NDIM>: unknown projector type",1);
            }

        } else {
            MADNESS_EXCEPTION("confused type in CCPairFunction<T,NDIM>",1);
        }

    }
//    print("working on ",argument[0].name(),"with ",(&projector)->type(),": result has",result.size(),"components");
    return result;
};


template class CCPairFunction<double,6>;
template class CCPairFunction<double,4>;
template class CCPairFunction<double,2>;

} // namespace madness
