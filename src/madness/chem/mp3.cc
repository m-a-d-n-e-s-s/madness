//
// Created by Florian Bischoff on 2/15/24.
//

#include "mp3.h"

namespace madness {
double MP3::compute_mp3_cd(const Pairs<CCPair>& mp2pairs, const Info& info) const {
    print_header3("compute term_CD of the MP3 energy with R2_bra");
    // compute the MP3 energy
    std::size_t nocc=info.mo_ket.size();
    print("freeze, nocc",info.parameters.freeze(),nocc);
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }
    const auto& R2=info.R_square;
    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    double result = 0.0;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            auto bra = clusterfunctions(i, j);
            double tmp1 = inner(bra, g12 * clusterfunctions(i, j), R2);
            double tmp2 = inner(bra, g12 * clusterfunctions(j, i), R2);
            double fac = (i == j) ? 0.5 : 1.0;
            double tmp = fac * (4.0 * tmp1 - 2.0 * tmp2);
            printf("mp3 energy: term_CD %2ld %2ld: %12.8f\n", i, j, tmp);
            result+= tmp;
        }
    }
    printf("MP3 energy: term_CD %12.8f\n", result);
    return result;
};

double MP3::compute_mp3_ef(const Pairs<CCPair>& mp2pairs, const Info& info) const {

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }



    const auto& R2=info.R_square;
    double result=0.0;
    const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);



    print_header3("computing term EF of the MP3 energy with R2_bra");
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            double tmp=0.0;
            for (std::size_t k=info.parameters.freeze(); k< nocc; ++k) {
                for (std::size_t l=info.parameters.freeze(); l<nocc; ++l) {
                    auto bra_ik = clusterfunctions(i,k);
                    auto bra_ki = clusterfunctions(k,i);
                    double ovlp_E=inner(bra_ik,clusterfunctions(j,l),R2);
                    double ovlp_F=inner(bra_ki,clusterfunctions(j,l),R2);
                    auto ket_i=nemo_orbital[i];
                    auto ket_k=nemo_orbital[k];
                    auto bra_j=R2_orbital[j];
                    auto bra_l=R2_orbital[l];

                    double g_jlik=inner(bra_j*ket_i, (*g12)(bra_l*ket_k));
                    //                    print("<jl | g | ik>",g_jlik);
                    tmp+=(2.0*ovlp_E - ovlp_F)*g_jlik;
                }
            }
            printf("mp3 energy: term_EF %2ld %2ld %12.8f\n",i,j,tmp);
            result+=tmp;
        }
    }
    printf("MP3 energy: term_EF %12.8f\n",result);
    return result;
};

/// permutations in physisists notation: <ij | kl>  = (ik | jl), loop over ik<jl

/// std::vector<permutation> all_permutations;
/// for (int ik=0; ik<npair(); ++ik) {
///     for (int jl=ik; jl<npair(); ++jl) {
///         auto [i,k]=ij_to_i_and_j(ik);
///         auto [j,l]=ij_to_i_and_j(jl);
///         permutation p(i,j,k,l);
///         auto perms=p.make_all_permutations();
///         for (const auto& p:perms) all_permutations.push_back(p);
///         int number_of_unique_permutations=perms.size();
///         print("ij, kl, i, j, k, l",ik,jl," -- ", i,j,k,l," - ",perms.size(),perms);
///     }
/// }
struct permutation {
    int i,j,k,l;
    permutation(int i, int j, int k, int l) : i(i), j(j), k(k), l(l) {}
    std::vector<permutation> make_all_eri_permutations() const{
        permutation p1(i,j,k,l);
        permutation p2(k,j,i,l);
        permutation p3(i,l,k,j);
        permutation p4(k,l,i,j);
        permutation p5(j,i,l,k);
        permutation p6(l,i,j,k);
        permutation p7(j,k,l,i);
        permutation p8(l,k,j,i);
        std::vector<permutation> result({p1,p2,p3,p4,p5,p6,p7,p8});
        return remove_duplicates(result);
    }

    /// <\tau_{ij} | \tau_{kl}> = <\tau_{ji} | \tau_{lk}> = <\tau_{kl} | \tau_{ij}>= <\tau_{lk} | \tau_{ji}>
    std::vector<permutation> make_all_tau_permutations() const{
        permutation p1(i,j,k,l);
        permutation p2(j,i,l,k);
        permutation p3(k,l,i,j);
        permutation p4(l,k,j,i);
        return remove_duplicates({p1,p2,p3,p4});
    }

    bool operator<(const permutation& other) const {
        return std::tie(i,j,k,l) < std::tie(other.i,other.j,other.k,other.l);
    }
    bool operator==(const permutation& other) const {
        return std::tie(i,j,k,l) == std::tie(other.i,other.j,other.k,other.l);
    }
    bool operator!=(const permutation& other) const {
        return (not (*this==other));
    }
    static std::vector<permutation> remove_duplicates(std::vector<permutation> v) {
        // remove duplicates
        std::sort(v.begin(), v.end(),[](permutation a, permutation b) { return a < b; });
        auto last = std::unique(v.begin(), v.end());
        v.erase(last, v.end());
        return v;
    }
};
std::ostream& operator<<(std::ostream& os, const permutation& p) {
    os << "(" << p.i << p.j << p.k << p.l << ")";
    return os;
}
/// compute the EF term of the MP3 energy with permutational symmetry

/// there is the usual 8-fold permutational symmetry for 2-electron integrals
/// <ij | kl >      = <ji | lk > = <kj | il >  = < il | kj >
///    = <kl | ij>  = <lk | ji > = <il | kj >  = < kj | il >
/// implemented as loop over (ij) = \sum_i<j  and (ij) < (kl)
double MP3::compute_mp3_ef_with_permutational_symmetry(const Pairs<CCPair>& mp2pairs, const Info& info) const {
    print_header3("computing term EF of the MP3 energy with R2_bra");
    World& world=mp2pairs.allpairs.begin()->second.function().world();
    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    std::size_t nfrozen=info.parameters.freeze();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }

    const auto& R2=info.R_square;
    double result=0.0;
    const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    // number of pairs
    auto npair = [&nocc, &nfrozen]() { return (nocc - nfrozen) * (nocc - nfrozen + 1) / 2; };

    // turn composite index ij into i and j, taking care of frozen orbitals
    PairVectorMap map=PairVectorMap::triangular_map(nfrozen,nocc);
    auto ij_to_i_and_j = [&map](const int ij) { return map.map[ij]; };

    /// <ij | kl >  = (ik | jl)
    std::vector<permutation> all_tau_permutations;
    // loop over unique pairs (ij)
    for (size_t ij=0; ij<npair(); ++ij) {
        auto [i,j]=ij_to_i_and_j(ij);
        double tmp=0;
        std::vector<CCPairFunction<double,6>> tmp_tau;
        // loop over all k and l
        for (std::size_t k=nfrozen; k<nocc; ++k) {
            for (std::size_t l=nfrozen; l<nocc; ++l) {

                // make all possible permutations of the 4 indices i,j,k,l
                permutation p0(i,j,k,l);
                auto perms=p0.make_all_tau_permutations();  // permutations are sorted
                // continue only if this permutation is the canonical one
                if (p0!=perms.front()) continue;
                for (const auto& p:perms) all_tau_permutations.push_back(p);
                const double weight=perms.size();

                // terms C+D = <tau_ij | tau_kl> (2*<ij|g|kl> - <ji|g|kl>)
                //           = <tau_ij | tau_kl> (2*(ik|jl) - (jk|il))
                const auto& ket_i=nemo_orbital[i];
                const auto& ket_k=nemo_orbital[k];
                const auto& ket_l=nemo_orbital[l];
                const auto& bra_i=R2_orbital[i];
                const auto& bra_j=R2_orbital[j];
                const auto& bra_l=R2_orbital[l];
                double g_ikjl=inner(bra_i*ket_k, (*g12)(bra_j*ket_l));
                double g_jkil=inner(bra_j*ket_k, (*g12)(bra_l*ket_i));
                // double ovlp=inner(clusterfunctions(i,j),clusterfunctions(k,l),R2);
                tmp_tau+=weight*(2.0* g_ikjl - g_jkil)*clusterfunctions(k,l);
                // tmp+=weight*ovlp*(2.0* g_ikjl - g_jkil);
            }
        }
        tmp_tau=consolidate(tmp_tau,{"remove_lindep"});
        tmp=inner(clusterfunctions(i,j),tmp_tau,R2);
        printf("mp3 energy: term_EF %2d %2d %12.8f\n",i,j,tmp);
        result+=tmp;
    }
    printf("MP3 energy: term_EF %12.8f\n",result);

    // sanity check
    int npermutations=all_tau_permutations.size();
    all_tau_permutations=permutation::remove_duplicates(all_tau_permutations);
    int nuniquepermutations=all_tau_permutations.size();
    int ntotalpermutations=std::pow(nocc-nfrozen,4);
    MADNESS_CHECK_THROW(npermutations==nuniquepermutations,"incorrect number of unique permutations");
    MADNESS_CHECK_THROW(npermutations==ntotalpermutations,"incorrect number of unique permutations");

    return result;
};

double MP3::compute_mp3_ef_low_scaling(const Pairs<CCPair>& mp2pairs,
    const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const {
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    print_header3("computing term EF of the MP3 energy with R2_bra, low-scaling version");

    std::size_t nocc=info.mo_ket.size();
    std::size_t nfrozen=info.parameters.freeze();

    const auto& R2=info.R_square;
    double result=0.0;
    const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    timer timer_sum(world);
    timer_sum.interrupt();
    timer timer_inner(world);
    timer_inner.interrupt();
    timer timer_consolidate(world);
    timer_consolidate.interrupt();

    /// <ij | kl >  = (ik | jl)
    std::vector<permutation> all_tau_permutations;
    // loop over unique pairs (ij)
    for (std::size_t i=nfrozen; i<nocc; ++i) {
        for (std::size_t j=nfrozen; j<nocc; ++j) {
            timer_sum.resume();
            std::vector<CCPairFunction<double,6>> sigma;
            for (std::size_t k=nfrozen; k<nocc; ++k) {
                for (std::size_t l=nfrozen; l<nocc; ++l) {
                    const auto& ket_i=nemo_orbital[i];
                    const auto& ket_k=nemo_orbital[k];
                    const auto& ket_l=nemo_orbital[l];
                    const auto& bra_i=R2_orbital[i];
                    const auto& bra_j=R2_orbital[j];
                    const auto& bra_l=R2_orbital[l];
                    double g_ikjl=inner(bra_i*ket_k, (*g12)(bra_j*ket_l));
                    double g_jkil=inner(bra_j*ket_k, (*g12)(bra_l*ket_i));
                    double g_ijkl=(2.0* g_ikjl - g_jkil);
                    sigma+=g_ijkl*clusterfunctions(k,l);
                }
            }
            timer_sum.interrupt();
            timer_consolidate.resume();
            sigma=consolidate(sigma,{"remove_lindep"});
            timer_consolidate.interrupt();
            timer_inner.resume();
            double tmp=inner(clusterfunctions(i,j),sigma,R2);
            printf("mp3 energy: term_EF %2ld %2ld %12.8f\n",i,j,tmp);
            result+=tmp;
            timer_inner.interrupt();
        }
    }

    timer_sum.print("summation in EF term");
    timer_consolidate.print("consolidation in EF term");
    timer_inner.print("inner in EF term");
    printf("MP3 energy: term_EF %12.8f\n",result);
    return result;

}

double MP3::compute_mp3_ghij(const Pairs<CCPair>& mp2pairs,
    const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const {
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    double result=0.0;

    const auto& R2=info.R_square;
    //const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    timer t2(world);


    // compute intermediates for terms G, I, H, and J

    // \sum_j tau_ij(1,2) * phi_j(2)
    std::vector<ClusterFunction> tau_kk_i(nocc);
    // \sum_j tau_ij(1,2) * phi_j(1)
    std::vector<ClusterFunction> tau_ij_j(nocc);
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            auto tmp2 = multiply(clusterfunctions(i, j), R2_orbital[i], {0, 1, 2});
            for (auto& t: tmp2) tau_kk_i[j].push_back(t);

            auto tmp3 = multiply(clusterfunctions(i, j), R2_orbital[j], {0, 1, 2});
            for (auto& t: tmp3) tau_ij_j[i].push_back(t);
        }
    }
    std::vector<std::string> consolidation={"op_pure_to_pure","remove_lindep"};
    print("consolidating with ",consolidation);
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        tau_kk_i[i] = consolidate(tau_kk_i[i],consolidation);
        tau_ij_j[i] = consolidate(tau_ij_j[i], consolidation);
        // for (auto& c: tau_kk_i[i]) c.info();
    }

    t2.tag("GHIJ term prep");

    // terms G, I, H, J of Bartlett/Silver 1975
    real_convolution_3d& g = *(g12->get_op());
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        // tmp(1,2) = g(1,1') | tau_ij(1',2) j(2) >
        timer t4(world, "gtau");
        g.set_particle(1);
        auto gtau_same = g(tau_kk_i[i]);
        t4.tag("compute gtau_same");

        // tmp(1',2) = g(1',1) | tau_ij(1,2) j(1) >
        g.set_particle(1);
        auto gtau_other = g(tau_ij_j[i]); // < tau_ij(1,2) j(1) | g(1,1') |
        t4.tag("compute gtau_other");

        auto bra_kk_i = multiply(tau_kk_i[i], R2, {3, 4, 5});
        auto bra_ij_j = multiply(tau_ij_j[i], R2, {3, 4, 5});

        double G = inner(bra_kk_i, gtau_same);
        printf("G     %12.8f\n", G);

        double H = inner(bra_ij_j, gtau_other);
        printf("H     %12.8f\n", H);

        double I = inner(bra_kk_i, gtau_other);
        printf("I     %12.8f\n", I);

        double J = inner(bra_ij_j, gtau_same);
        printf("J     %12.8f\n", J);

        t4.tag("compute inner products");
        double tmp = (8.0 * G - 4.0 * I + 2.0 * H - 4.0 * J);
        printf("mp3 energy: term_GHIJ  %2ld %12.8f\n", i, tmp);
        result += tmp;
    }
    printf("MP3 energy: term_GHIJ %12.8f\n", result);
    t2.tag("GHIJ term");
    return result;
};


double MP3::compute_mp3_ghij_fast(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const {
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    print_header3("entering compute_mp3_ghij_fast");

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    double result=0.0;

    const auto& R2=info.R_square;
    //const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    timer t2(world);

    // compute intermediates for terms G, I, H, and J

    // \sum_j tau_ij(1,2) * phi_j(2)
    std::vector<ClusterFunction> tau_kk_i(nocc);
    // \sum_j tau_ij(1,2) * phi_j(1)
    std::vector<ClusterFunction> tau_ij_j(nocc);

    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            auto tmp2 = multiply(clusterfunctions(i, j), R2_orbital[i], {0, 1, 2});
            for (auto& t: tmp2) tau_kk_i[j].push_back(t);

            auto tmp3 = multiply(clusterfunctions(i, j), R2_orbital[j], {0, 1, 2});
            for (auto& t: tmp3) tau_ij_j[i].push_back(t);
        }
    }

    std::vector<std::string> consolidation={"op_pure_to_pure","remove_lindep","op_dec_to_dec"};
    print("consolidating with ",consolidation);
    std::vector<ClusterFunction> intermediate(nocc);
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        intermediate[i]=2.0*tau_kk_i[i];
        intermediate[i]-=tau_ij_j[i];
        intermediate[i]=consolidate(intermediate[i],consolidation,info.molecular_coordinates);
    }

    t2.tag("GHIJ term prep");

    // terms G, I, H, J of Bartlett/Silver 1975
    real_convolution_3d& g = *(g12->get_op());
    g.set_particle(1);
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        // tmp(1,2) = g(1,1') | tau_ij(1',2) j(2) >
        timer t4(world, "gtau");
        auto gintermediate = g(intermediate[i]);
        t4.tag("compute gintermediate");
        auto bra_intermediate = multiply(intermediate[i], R2, {3, 4, 5});
        t4.tag("multiply");
        double tmp = 2.0*inner(bra_intermediate, gintermediate);
        printf("mp3 energy: term_GHIJ  %2ld %12.8f\n", i, tmp);
        t4.tag("inner");
        result += tmp;
    }
    printf("MP3 energy: term_GHIJ %12.8f\n", result);
    t2.tag("GHIJ term");
    return result;
};

double MP3::compute_mp3_klmn_fast(const Pairs<CCPair>& mp2pairs, const Info& info) const {
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }
    double result=0.0;

    const auto& R2=info.R_square;
    const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    // compute the term <i(1) |g(1,2) | j(1)>(2)
    madness::Pairs<real_function_3d> gij;

    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            gij.insert(i,j,(*g12)(nemo_orbital[i]*R2_orbital[j]));
        }
    }


    timer multiply_KLMN(world, "multiplication in KLMN term");
    multiply_KLMN.interrupt();
    timer inner_KLMN(world, "inner in KLMN term");
    inner_KLMN.interrupt();
    // prepare intermediates for terms K, L, M, N of Bartlett/Silver 1975
    // tau_g_ij(1,2) = \sum_k tau_ik(1,1') g_jk(2)
    Pairs<ClusterFunction> tau_ik_g_kj, tau_kj_g_ki;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            multiply_KLMN.resume();
            std::vector<CCPairFunction<double, 6>> rhs;
            for (std::size_t k = info.parameters.freeze(); k < nocc; ++k) {
                rhs += +2.0 * multiply(clusterfunctions(j, k), gij(k, i), {0, 1, 2});   // M
                rhs += +2.0 * multiply(clusterfunctions(k, i), gij(k, j), {0, 1, 2});    //  N
                rhs += -4.0 * multiply(clusterfunctions(i, k), gij(k, j), {3, 4, 5});   // K
                rhs += -4.0 * multiply(clusterfunctions(k, j), gij(k, i), {3, 4, 5});   // L
            }
            rhs = consolidate(rhs, {});
            multiply_KLMN.interrupt();

            inner_KLMN.resume();
            double tmp = inner(clusterfunctions(i, j), rhs, R2);
            inner_KLMN.interrupt();

            printf("mp3 energy: term_KLMN with particle=1 %2ld %2ld %12.8f\n", i, j, tmp);
            result += tmp;
        }
    }
    printf("MP3 energy: term_KLMN (KLMN) %12.8f\n", result);
    multiply_KLMN.print("multiplication in KLMN term");
    inner_KLMN.print("inner in KLMN term");

    return result;

};
double MP3::compute_mp3_klmn(const Pairs<CCPair>& mp2pairs, const Info& info) const {
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }
    double result=0.0;

    const auto& R2=info.R_square;
    const std::vector<real_function_3d>& nemo_orbital=info.mo_ket;
    const std::vector<real_function_3d>& R2_orbital=info.mo_bra;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    // compute the term <i(1) |g(1,2) | j(1)>(2)
    madness::Pairs<real_function_3d> gij;

    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            gij.insert(i,j,(*g12)(nemo_orbital[i]*R2_orbital[j]));
        }
    }


    timer multiply_KLMN(world, "multiplication in KLMN term");
    multiply_KLMN.interrupt();
    timer inner_KLMN(world, "inner in KLMN term");
    inner_KLMN.interrupt();
    // prepare intermediates for terms K, L, M, N of Bartlett/Silver 1975
    // tau_g_ij(1,2) = \sum_k tau_ik(1,1') g_jk(2)
    Pairs<ClusterFunction> tau_ik_g_kj, tau_kj_g_ki;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            multiply_KLMN.resume();
            std::vector<CCPairFunction<double, 6>> tmp1, tmp2;
            for (std::size_t k = info.parameters.freeze(); k < nocc; ++k) {
                tmp1 += multiply(clusterfunctions(i, k), gij(k, j), {3, 4, 5});
                tmp2 += multiply(clusterfunctions(k, j), gij(k, i), {3, 4, 5});
            }
            tmp1 = consolidate(tmp1, {});
            tmp2 = consolidate(tmp2, {});
            tau_ik_g_kj(i, j) = tmp1;
            tau_kj_g_ki(i, j) = tmp2;
            multiply_KLMN.interrupt();

            inner_KLMN.resume();
            double K = inner(clusterfunctions(i, j), tau_ik_g_kj(i, j), R2);
            double L = inner(clusterfunctions(i, j), tau_kj_g_ki(i, j), R2);
            double M = inner(clusterfunctions(j, i), tau_kj_g_ki(i, j), R2);
            double N = inner(clusterfunctions(j, i), tau_ik_g_kj(i, j), R2);
            inner_KLMN.interrupt();

            double tmp = -4 * K - 4 * L + 2 * M + 2 * N;
            printf("mp3 energy: term_KLMN with particle=1 %2ld %2ld %12.8f\n", i, j, tmp);
            result += tmp;
        }
    }
    printf("MP3 energy: term_KLMN (KLMN) %12.8f\n", result);
    multiply_KLMN.print("multiplication in KLMN term");
    inner_KLMN.print("inner in KLMN term");

    return result;

};

double MP3::mp3_test(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions, const Info& info) const {
    print_header2("entering mp3 test");
    World& world=mp2pairs.allpairs.begin()->second.function().world();

    auto R2 = info.R_square;
    const std::size_t nocc=info.mo_ket.size();
    std::vector<real_function_3d> nemo_orbital=info.mo_ket;
    std::vector<real_function_3d> R2_orbital=info.mo_bra;
    //    std::vector<real_function_3d> nemo_orbital=mo_ket().get_vecfunction();
    std::vector<CCPairFunction<double,6>> ij;
    std::size_t i=1, j=1;
    ij.push_back(CCPairFunction<double,6>(nemo_orbital[i],nemo_orbital[j]));
    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);
    double eri=0.0;
    for (std::size_t i=0; i<nocc; ++i) {
        for (std::size_t j=0; j<nocc; ++j) {
            std::vector<CCPairFunction<double, 6>> ii;
            ii.push_back(CCPairFunction<double, 6>(nemo_orbital[i], nemo_orbital[j]));
            double tmp = inner(ii, g12 * ii, R2);
            print("eri for pair", i,j, tmp);
            eri += tmp;
        }
    }
    print("total eri",eri);

    print("eri(1,1)",eri);
    double mp2_energy=inner(clusterfunctions(i,j),ij,R2);
    print("mp2 energy for pair ", i, j, mp2_energy);

    double mp3_contrib=inner(clusterfunctions(i,j),clusterfunctions(i,j),R2);
    print("<tau_ij | tau_ij>",mp3_contrib);

    return 0.0;
}

double MP3::compute_mp3_cd(World& world,
                           const long i, const long j,
                           const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                           const Info& info,
                           // const std::vector<Function<double,3>>& mo_ket,
                           // const std::vector<Function<double,3>>& mo_bra,
                           // const CCParameters& parameters,
                           // const Molecule& molecule,
                           // const Function<double,3>& Rsquare,
                           const std::vector<std::string>& argument) {

    CCConvolutionOperator<double,3>::Parameters cparam(info.parameters);
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    auto bra = pair_square(i, j);
    double tmp1 = inner(bra, g12 * pair_square(i, j), info.R_square);
    double tmp2 = inner(bra, g12 * pair_square(j, i), info.R_square);
    double fac = (i == j) ? 0.5 : 1.0;
    double result = fac * (4.0 * tmp1 - 2.0 * tmp2);
    constexpr std::size_t bufsize=256;
    char buf[bufsize];
    snprintf(buf,bufsize,"mp3 energy: term_CD %2ld %2ld: %12.8f\n", i, j, result);
    print(std::string(buf));
    return result;
}

double MP3::compute_mp3_ef(World& world,
                           const long i, const long j,
                           const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                           const Info& info,
                           // const std::vector<Function<double,3>>& mo_ket,
                           // const std::vector<Function<double,3>>& mo_bra,
                           // const CCParameters& parameters,
                           // const Molecule& molecule,
                           // const Function<double,3>& Rsquare,
                           const std::vector<std::string>& argument) {

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    std::size_t nfrozen=info.parameters.freeze();
    std::size_t nocc=info.mo_ket.size();

    /// <ij | kl >  = (ik | jl)
    std::vector<permutation> all_tau_permutations;
    std::vector<CCPairFunction<double,6>> tmp_tau;

    long counter=0;
    timer timer_prep(world);
    // loop over all k and l
    for (std::size_t k=nfrozen; k<nocc; ++k) {
        for (std::size_t l=nfrozen; l<nocc; ++l) {
            counter++;

            // make all possible permutations of the 4 indices i,j,k,l
            permutation p0(i,j,k,l);
            auto perms=p0.make_all_tau_permutations();  // permutations are sorted
            // continue only if this permutation is the canonical one
            if (p0!=perms.front()) continue;
            // for (const auto& p:perms) all_tau_permutations.push_back(p);
            const double weight=perms.size();

            // terms C+D = <tau_ij | tau_kl> (2*<ij|g|kl> - <ji|g|kl>)
            //           = <tau_ij | tau_kl> (2*(ik|jl) - (jk|il))
            const auto& ket_i=info.mo_ket[i];
            const auto& ket_k=info.mo_ket[k];
            const auto& ket_l=info.mo_ket[l];
            const auto& bra_i=info.mo_bra[i];
            const auto& bra_j=info.mo_bra[j];
            const auto& bra_l=info.mo_bra[l];
            double g_ikjl=inner(bra_i*ket_k, (*g12)(bra_j*ket_l));
            double g_jkil=inner(bra_j*ket_k, (*g12)(bra_l*ket_i));
            tmp_tau+=weight*(2.0* g_ikjl - g_jkil)*pair_square(k,l);

            // save some memory
            // if (counter%10) tmp_tau=consolidate(tmp_tau,{"remove_lindep"});
        }
    }
    timer_prep.tag("preparing in EF term");
    tmp_tau=consolidate(tmp_tau,{"remove_lindep"});
    timer_prep.tag("consolidation in EF term");
    timer timer_inner(world);
    double result=inner(pair_square(i,j),tmp_tau,info.R_square);
    timer_prep.tag("inner in EF term");

    constexpr std::size_t bufsize=256;
    char buf[bufsize];
    snprintf(buf, bufsize,"mp3 energy: term_EF %2ld %2ld %12.8f\n",i,j,result);
    print(std::string(buf));

    // can't do sanity check, because we are working on a single pair
//    int npermutations=all_tau_permutations.size();
//    all_tau_permutations=permutation::remove_duplicates(all_tau_permutations);
//    int nuniquepermutations=all_tau_permutations.size();
//    int ntotalpermutations=std::pow(nocc-nfrozen,4);
//    MADNESS_CHECK_THROW(npermutations==nuniquepermutations,"incorrect number of unique permutations");
//    MADNESS_CHECK_THROW(npermutations==ntotalpermutations,"incorrect number of unique permutations");

    return result;
}

double MP3::compute_mp3_ghij(World& world,
                           const long i, const long jj,
                           const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                           const Info& info,
                           // const std::vector<Function<double,3>>& mo_ket,
                           // const std::vector<Function<double,3>>& mo_bra,
                           // const CCParameters& parameters,
                           // const Molecule& molecule,
                           // const Function<double,3>& Rsquare,
                           const std::vector<std::string>& argument) {

    // this scales linearly, use only the diagaonal contributions
    if (i!=jj) return 0.0;

    print_header3("entering compute_mp3_ghij_fast");

    // prepare cluster functions
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    double result=0.0;

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    timer t2(world);

    // compute intermediates for terms G, I, H, and J

    // \sum_j tau_ij(1,2) * phi_j(2)
    ClusterFunction tau_kk_i;
    // \sum_j tau_ij(1,2) * phi_j(1)
    ClusterFunction tau_ij_j;

    for (std::size_t ii = info.parameters.freeze(); ii < nocc; ++ii) {
        // for (std::size_t j = parameters.freeze(); j < nocc; ++j) {
            auto tmp2 = multiply(pair_square(ii, i), info.mo_bra[ii], {0, 1, 2});
            for (auto& t: tmp2) tau_kk_i.push_back(t);
        // }
    }

    // for (std::size_t i = parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = info.parameters.freeze(); j < nocc; ++j) {
            auto tmp3 = multiply(pair_square(i, j), info.mo_bra[j], {0, 1, 2});
            for (auto& t: tmp3) tau_ij_j.push_back(t);
        }
    // }

    std::vector<std::string> consolidation={"op_pure_to_pure","remove_lindep","op_dec_to_dec"};
    // std::vector<std::string> consolidation={"op_pure_to_pure","remove_lindep","op_dec_to_pure"};
    print("consolidating with ",consolidation);
    ClusterFunction intermediate;
    // for (std::size_t i = parameters.freeze(); i < nocc; ++i) {
        intermediate=2.0*tau_kk_i;
        intermediate-=tau_ij_j;
        intermediate=consolidate(intermediate,consolidation,info.molecular_coordinates);
    // }

    t2.tag("GHIJ term prep");

    // terms G, I, H, J of Bartlett/Silver 1975
    real_convolution_3d& g = *(g12->get_op());
    g.set_particle(1);
    // for (std::size_t i = parameters.freeze(); i < nocc; ++i) {
        // tmp(1,2) = g(1,1') | tau_ij(1',2) j(2) >
        timer t4(world, "gtau");
        auto gintermediate = g(intermediate);
        t4.tag("compute gintermediate");
        auto bra_intermediate = multiply(intermediate, info.R_square, {3, 4, 5});
        t4.tag("multiply");
        double tmp = 2.0*inner(bra_intermediate, gintermediate);
        constexpr std::size_t bufsize=256;
        char buf[bufsize];
        snprintf(buf, bufsize,"mp3 energy: term_GHIJ  %2ld %12.8f\n", i, tmp);
        print(std::string(buf));
        t4.tag("inner");
        result += tmp;
    // }
    t2.tag("GHIJ term");
    return result;
}


double MP3::compute_mp3_klmn(World& world,
                           const long i, const long j,
                           const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                           const Info& info,
                           // const std::vector<Function<double,3>>& mo_ket,
                           // const std::vector<Function<double,3>>& mo_bra,
                           // const CCParameters& parameters,
                           // const Molecule& molecule,
                           // const Function<double,3>& Rsquare,
                           const std::vector<std::string>& argument) {

    print_header3("starting term KLMN with particles "+std::to_string(i)+" "+std::to_string(j));
    std::size_t nocc=info.mo_ket.size();

    CCConvolutionOperator<double,3>::Parameters cparam;
    auto g12=CCConvolutionOperatorPtr<double,3>(world,OpType::OT_G12,cparam);

    // compute the term <i(1) |g(1,2) | j(1)>(2)
    madness::Pairs<real_function_3d> gij;

    for (std::size_t k = info.parameters.freeze(); k < nocc; ++k) {
        for (std::size_t l = info.parameters.freeze(); l < nocc; ++l) {
            gij.insert(k,l,(*g12)(info.mo_ket[k]*info.mo_bra[l]));
        }
    }


    timer multiply_KLMN(world, "multiplication in KLMN term");
    multiply_KLMN.interrupt();
    timer inner_KLMN(world, "inner in KLMN term");
    inner_KLMN.interrupt();

    multiply_KLMN.resume();
    std::vector<CCPairFunction<double, 6>> rhs;
    for (std::size_t k = info.parameters.freeze(); k < nocc; ++k) {
        // rhs += +2.0 * multiply(pair_square(j, k), gij(k, i), {0, 1, 2});   // M
        // rhs += +2.0 * multiply(pair_square(k, i), gij(k, j), {0, 1, 2});    //  N
        // rhs += -4.0 * multiply(pair_square(i, k), gij(k, j), {3, 4, 5});   // K
        // rhs += -4.0 * multiply(pair_square(k, j), gij(k, i), {3, 4, 5});   // L
        rhs += multiply(pair_square(k, j), gij(k, i), {0, 1, 2});
        rhs += multiply(pair_square(k, j), gij(k, i), {3, 4, 5});
    }
    rhs = consolidate(rhs, {});

    std::vector<CCPairFunction<double, 6>> lhs;
    lhs+=2.0*pair_square(i,j);
    lhs-= pair_square(j,i);
    lhs=consolidate(lhs,{});
    multiply_KLMN.interrupt();

    inner_KLMN.resume();
    // double result = inner(pair_square(i, j), rhs, Rsquare);
    double result = -2.0*inner(lhs, rhs, info.R_square);
    inner_KLMN.interrupt();

    constexpr std::size_t bufsize=256;
    char buf[bufsize];
    snprintf(buf,bufsize,"mp3 energy: term_KLMN %2ld %2ld %12.8f\n", i, j, result);
    print(std::string(buf));

    multiply_KLMN.print("multiplication in KLMN term with remove_lindep");
    inner_KLMN.print("inner in KLMN term");

    return result;
}

double MP3::mp3_energy_contribution(const Pairs<CCPair>& mp2pairs) const {

    print_header2("computing the MP3 correlation energy");
    print("mp2pairs.size()",mp2pairs.allpairs.size());

    timer t2(world);
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }

    t2.tag("make cluster functions");
//    mp3_test(mp2pairs,clusterfunctions);
    double term_CD=0.0, term_EF=0.0, term_GHIJ=0.0, term_KLMN=0.0;
    term_CD=compute_mp3_cd(mp2pairs,info);
    t2.tag("CD term");

    term_GHIJ=compute_mp3_ghij_fast(mp2pairs,clusterfunctions,info);
    t2.tag("GHIJ term");
    term_KLMN=compute_mp3_klmn_fast(mp2pairs,info);
    t2.tag("KLMN term fast");
    term_EF=compute_mp3_ef_with_permutational_symmetry(mp2pairs,info);
    t2.tag("EF term, permutational symmetry");

    printf("term_CD    %12.8f\n",term_CD);
    printf("term_GHIJ  %12.8f\n",term_GHIJ);
    printf("term_KLMN  %12.8f\n",term_KLMN);
    printf("term_EF    %12.8f\n",term_EF);
    double mp3_energy=term_CD+term_GHIJ+term_KLMN+term_EF;
    printf("MP3 energy contribution  %12.8f\n",mp3_energy);
    return mp3_energy;
}

double MP3::mp3_energy_contribution_macrotask_driver(const Pairs<CCPair>& mp2pairs) const {

    if (world.rank()==0) {
        print_header2("computing the MP3 correlation energy, macrotask version");
        print("mp2pairs.size()",mp2pairs.allpairs.size());
    }

    // compute all ij pairs
    timer t2(world);
    std::size_t nocc=info.mo_ket.size();
    typedef std::vector<CCPairFunction<double,6>> ClusterFunction;
    Pairs<ClusterFunction> clusterfunctions;
    for (std::size_t i = info.parameters.freeze(); i < nocc; ++i) {
        for (std::size_t j = i; j < nocc; ++j) {
//            clusterfunctions(i,j)=mp2pairs(i,j).functions;
            clusterfunctions(i,j)=CCPotentials::make_pair_mp2(world,mp2pairs(i,j).function(),i,j,info,true).functions;
            if (i!=j) {
                for (const auto& t : clusterfunctions(i,j)) {
                    clusterfunctions(j, i).push_back(t.swap_particles());
                }
            }
        }
    }
    // turn Pair into vector for cloud and stuff -- will be reversed later on
    PairVectorMap square_map=PairVectorMap::quadratic_map(info.parameters.freeze(),info.mo_ket.size());
    auto clusterfunc_vec=Pairs<std::vector<CCPairFunction<double,6>>>::pairs2vector(clusterfunctions,square_map);

    t2.tag("make cluster functions");

    // create dummy scheduling vector of length npair=nocc*(nocc+1)/2, for the macrotask scheduler, triangular form
    std::vector<int> ij_triangular(mp2pairs.allpairs.size());
    // create dummy scheduling vector of length nact for the macrotask scheduler, square form
    std::vector<int> nact(nocc-info.parameters.freeze());
    std::vector<int> dummy;


    // get mos
    //const std::vector<real_function_3d>& ket=info.mo_ket;
    //const std::vector<real_function_3d>& bra=info.mo_bra;

    auto taskq=std::shared_ptr<MacroTaskQ>(new MacroTaskQ(MacroTaskQFactory(world)));
    // taskq->set_printlevel(20);
    // taskq->cloud.set_debug(true);
    // MacroTaskMP3 task_triangular("triangular");
    // MacroTaskMP3 task_square("square");
    MacroTaskMP3 ghij_task("ghij");
    MacroTaskMP3 klmn_task("klmn");
    MacroTaskMP3 cd_task("cd");
    MacroTaskMP3 ef_task("ef");
    MacroTask macrotask_ghij(world,ghij_task,taskq);
    MacroTask macrotask_klmn(world,klmn_task,taskq);
    MacroTask macrotask_cd(world,cd_task,taskq);
    MacroTask macrotask_ef(world,ef_task,taskq);
    auto ghij_future=macrotask_ghij(ij_triangular, dummy, clusterfunc_vec,info, std::vector<std::string>());
    auto klmn_future=macrotask_klmn(nact, nact, clusterfunc_vec, info, std::vector<std::string>());
    auto cd_future=macrotask_cd(ij_triangular, dummy, clusterfunc_vec, info, std::vector<std::string>());
    auto ef_future=macrotask_ef(ij_triangular, dummy, clusterfunc_vec, info, std::vector<std::string>());
    taskq->print_taskq();
    taskq->run_all();

    double term_CD=cd_future.get();
    double term_EF=ef_future.get();
    double term_GHIJ=ghij_future.get();
    double term_KLMN=klmn_future.get();
    double mp3_energy=term_CD+term_GHIJ+term_KLMN+term_EF;
    if (world.rank()==0) {
        printf("term_CD    %12.8f\n",term_CD);
        printf("term_GHIJ  %12.8f\n",term_GHIJ);
        printf("term_KLMN  %12.8f\n",term_KLMN);
        printf("term_EF    %12.8f\n",term_EF);
        printf("MP3 energy contribution  %12.8f\n",mp3_energy);
    }
    return mp3_energy;
    // return 0.0;
}
}
