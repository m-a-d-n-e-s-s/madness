#include<madness/chem/exchangeoperator.h>

#include<madness/chem/SCF.h>
#include<madness/chem/nemo.h>
#include <madness/mra/macrotaskpartitioner.h>

using namespace madness;

namespace madness {


template<typename T, std::size_t NDIM>
Exchange<T, NDIM>::ExchangeImpl::ExchangeImpl(World& world, const SCF *calc, const int ispin)
        : world(world), symmetric_(false), lo(calc->param.lo()) {

    if constexpr (std::is_same_v<T,double_complex>) {
        if (ispin == 0) { // alpha spin
            mo_ket = convert<double, T, NDIM>(world, calc->amo);        // deep copy necessary if T==double_complex
        } else if (ispin == 1) {  // beta spin
            mo_ket = convert<double, T, NDIM>(world, calc->bmo);
        }
        mo_bra = conj(world, mo_ket);
    } else {
        if (ispin == 0) { // alpha spin
            mo_ket = calc->amo;        // deep copy necessary if T==double_complex
        } else if (ispin == 1) {  // beta spin
            mo_ket = calc->bmo;
        }
        mo_bra = mo_ket;
    }

}

template<typename T, std::size_t NDIM>
Exchange<T, NDIM>::ExchangeImpl::ExchangeImpl(World& world, const Nemo *nemo,
                            const int ispin) // @suppress("Class members should be properly initialized")
        : ExchangeImpl(world, nemo->get_calc().get(), ispin) {

    if (ispin == 0) { // alpha spin
        mo_ket = convert<double, T, NDIM>(world,
                                          nemo->get_calc()->amo);        // deep copy necessary if T==double_complex
    } else if (ispin == 1) {  // beta spin
        mo_ket = convert<double, T, NDIM>(world,
                                          nemo->get_calc()->bmo);        // deep copy necessary if T==double_complex
    }

    mo_bra = mul(world, nemo->ncf->square(), mo_ket);
    truncate(world, mo_bra);
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Exchange<T, NDIM>::ExchangeImpl::operator()(
        const std::vector<Function<T, NDIM> >& vket) const {

    reconstruct(world, mo_bra, false);
    reconstruct(world, mo_ket, false);
    world.gop.fence();
    norm_tree(world, mo_bra, false);
    norm_tree(world, mo_ket, false);
    world.gop.fence();

    reconstruct(world, vket);
    norm_tree(world, vket);

    // pick your algorithm.
    // Note that the macrotask algorithm partitions the exchange matrix into tiles. The final truncation
    // takes place for each sum over tile elements, instead of the sum over all matrix elements.
    // Truncation and addition doesn't commute, so truncation is done after the final accumulation.
    // Other truncations are elementwise and are not affected.
    reset_timer();
    vecfuncT Kf;
    if (algorithm_ == multiworld_efficient) {
        Kf = K_macrotask_efficient(vket, mul_tol);
    } else if (algorithm_ == multiworld_efficient_row) {
        Kf = K_macrotask_efficient_row(vket, mul_tol);
    } else if (algorithm_ == small_memory) {
        Kf = K_small_memory(vket, mul_tol);     // Smaller memory algorithm ... possible 2x saving using i-j sym
    } else if (algorithm_ == large_memory) {
        Kf = K_large_memory(vket, mul_tol);
    } else {
        MADNESS_EXCEPTION("unknown algorithm in exchangeoperator", 1);
    }
    if (printdebug()) {
        auto size = get_size(world, Kf);
        if (world.rank() == 0) print("total size of Kf before truncation", size);
    }
    truncate(world, Kf);
    if (printdebug()) {
        auto size=get_size(world,Kf);
        if (world.rank()==0) print("total size of Kf after truncation",size);
    }
    if (printlevel >= 3) print_timer(world);
    return Kf;
}

/// apply the exchange operator by tiling the exchange matrix

/// compute the matrix N_{ik} = N \phi_i \phi_k by tiles, with i,k \in batches A,B,
/// do a local reduce within the tiles: K_{iB} = \sum_{k \in batch B} \phi_k N_{ik}
/// and a universe-wide reduce of the tiles: K\phi_i = \sum_{batches B} K_{iB}
/// saving up to half of the cpu time compared to the naive algorithm
/// \tparam T       number type
/// \tparam NDIM    physical dimension of the argument vket
/// \param vf     argument of the exchange operator
/// \param mul_tol  cutoff parameter for sparse multiplication
/// \return         the exchange operator applied on vket
template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::ExchangeImpl::K_macrotask_efficient(const vecfuncT& vf, const double mul_tol) const {

    if (world.rank()==0 and printdebug()) print("\nentering macrotask_efficient version:");

    // the result is a vector of functions living in the universe
    const long nresult = vf.size();
    MacroTaskExchangeSimple xtask(nresult, lo, mul_tol, is_symmetric());
    vecfuncT Kf;

    // deferred execution if a taskq is provided by the user
    if (taskq) {
        taskq->set_printlevel(printlevel);
        MacroTask mtask(world, xtask, taskq);
        Kf = mtask(vf, mo_bra, mo_ket);
    } else {
        auto taskq_ptr = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(world, world.size()));
        taskq_ptr->set_printlevel(printlevel);
        MacroTask mtask(world, xtask, taskq_ptr);
        Kf = mtask(vf, mo_bra, mo_ket);
        taskq_ptr->run_all();
        if (printdebug()) taskq_ptr->cloud.print_timings(world);
        taskq_ptr->cloud.clear_timings();
        world.gop.fence();
    }
    return Kf;
}

/// compute each row of the exchange matrix in different subworlds
template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::ExchangeImpl::K_macrotask_efficient_row(const vecfuncT& vf, const double mul_tol) const {

    if (world.rank()==0 and printdebug()) print("\nentering macrotask_efficient_row version:");

    // the result is a vector of functions living in the universe
    const long nresult = vf.size();
    MacroTaskExchangeRow xtask(nresult, lo, mul_tol);
    vecfuncT Kf;

    // print the size of the amos
    if (printdebug()) {
        auto size=get_size(world,vf);
        if (world.rank()==0) print("total size of vf before iteration",size);
    }

    // deferred execution if a taskq is provided by the user
    if (taskq) {
        taskq->set_printlevel(printlevel);
        MacroTask mtask(world, xtask, taskq);
        mtask.set_name("K_macrotask_efficient_row");        
        Kf = mtask(vf, mo_bra, mo_ket);
    } else {
        auto taskq_ptr = std::shared_ptr<MacroTaskQ>(new MacroTaskQ(world, world.size()));
        taskq_ptr->set_printlevel(printlevel);
        MacroTask mtask(world, xtask, taskq_ptr);
        mtask.set_name("K_macrotask_efficient_row");        
        Kf = mtask(vf, mo_bra, mo_ket);
        taskq_ptr->run_all();
        if (printdebug()) taskq_ptr->cloud.print_timings(world);
        taskq_ptr->cloud.clear_timings();
        world.gop.fence();
    }
    return Kf;
}


template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Exchange<T, NDIM>::ExchangeImpl::K_small_memory(const vecfuncT& vket, const double mul_tol) const {

    const long nocc = mo_ket.size();
    const long nf = vket.size();
    vecfuncT Kf = zero_functions_compressed<T, NDIM>(world, nf);
    auto poisson = set_poisson(world, lo);

    for (int i = 0; i < nocc; ++i) {
        vecfuncT psif = mul_sparse(world, mo_bra[i], vket, mul_tol); /// was vtol
        truncate(world, psif);
        psif = apply(world, *poisson.get(), psif);
        truncate(world, psif);
        psif = mul_sparse(world, mo_ket[i], psif, mul_tol); /// was vtol
        gaxpy(world, 1.0, Kf, 1.0, psif);
    }
    truncate(world, Kf);
    return Kf;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Exchange<T, NDIM>::ExchangeImpl::K_large_memory(const vecfuncT& vket,
                                                                  const double mul_tol) const {    // Larger memory algorithm ... use i-j sym if psi==f

    auto poisson = set_poisson(world, lo);
    vecfuncT result = compute_K_tile(world, mo_bra, mo_ket, vket, poisson, is_symmetric(), mul_tol);
    truncate(world, result);
    return result;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::ExchangeImpl::compute_K_tile(World& world, const vecfuncT& mo_bra, const vecfuncT& mo_ket,
                                  const vecfuncT& vket, std::shared_ptr<real_convolution_3d> poisson,
                                  const bool symmetric, const double mul_tol) {

    double cpu0 = cpu_time();
    const long nf = vket.size();
    const long nocc = mo_ket.size();
    vecfuncT Kf = zero_functions_compressed<T, NDIM>(world, nf);

    vecfuncT psif;
    for (int i = 0; i < nocc; ++i) {
        int jtop = nf;
        if (symmetric)
            jtop = i + 1;
        for (int j = 0; j < jtop; ++j) {
            psif.push_back(mul_sparse(mo_bra[i], vket[j], mul_tol, false));
        }
    }

    world.gop.fence();
    truncate(world, psif);
    double cpu1 = cpu_time();
    mul1_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = cpu_time();
    psif = apply(world, *poisson.get(), psif);
    truncate(world, psif);
    cpu1 = cpu_time();
    apply_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = cpu_time();
    reconstruct(world, psif);
    norm_tree(world, psif);
    vecfuncT psipsif = zero_functions<T, NDIM>(world, nf * nocc);
    int ij = 0;
    for (int i = 0; i < nocc; ++i) {
        int jtop = nf;
        if (symmetric)
            jtop = i + 1;
        for (int j = 0; j < jtop; ++j, ++ij) {
            psipsif[i * nf + j] = mul_sparse(psif[ij], mo_ket[i], mul_tol, false);
            if (symmetric && i != j) {
                psipsif[j * nf + i] = mul_sparse(psif[ij], mo_ket[j], mul_tol, false);
            }
        }
    }

    world.gop.fence();
    cpu1 = cpu_time();
    mul2_timer += long((cpu1 - cpu0) * 1000l);
    psif.clear();
    world.gop.fence();
    compress(world, psipsif);
    for (int i = 0; i < nocc; ++i) {
        for (int j = 0; j < nf; ++j) {
            Kf[j].gaxpy(1.0, psipsif[i * nf + j], 1.0, false);
        }
    }
    // !! NO TRUNCATION AT THIS POINT !!
    world.gop.fence();
    psipsif.clear();
    world.gop.fence();
    return Kf;

}


/// compute a batch of the exchange matrix, with non-identical ranges

/// \param subworld     the world we're computing in
/// \param cloud        where to store the results
/// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
/// \param ket_batch    the ket batch of orbitals, also the orbitals to premultiply with
/// \param vf_batch     the argument of the exchange operator
template<typename T, std::size_t NDIM>
std::pair<std::vector<Function<T, NDIM>>, std::vector<Function<T, NDIM>>>
Exchange<T, NDIM>::ExchangeImpl::MacroTaskExchangeSimple::compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                          const vecfuncT& mo_ket,      // not batched
                                                                                          const vecfuncT& bra_batch,   // batched
                                                                                          const vecfuncT& vf_batch) const { // batched
    // orbital_product is a vector of vectors
    double cpu0 = cpu_time();
    std::vector<vecfuncT> orbital_product = matrix_mul_sparse<T, T, NDIM>(subworld, bra_batch, vf_batch, mul_tol);
    vecfuncT orbital_product_flat = flatten(orbital_product); // convert into a flattened vector
    truncate(subworld, orbital_product_flat);
    double cpu1 = cpu_time();
    mul1_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = cpu_time();
    auto poisson = set_poisson(subworld, lo);
    vecfuncT Nij = apply(subworld, *poisson.get(), orbital_product_flat);
    truncate(subworld, Nij);
    cpu1 = cpu_time();
    apply_timer += long((cpu1 - cpu0) * 1000l);

    // accumulate columns:      resultrow(i)=\sum_j j N_ij
    // accumulate rows:      resultcolumn(j)=\sum_i i N_ij
    cpu0 = cpu_time();

    // some helper functions
    std::size_t nrow = bra_batch.size();
    std::size_t ncolumn = vf_batch.size();
    auto ij = [&ncolumn](const int i, const int j) { return i * ncolumn + j; };

    auto Nslice = [&Nij, &ij, &ncolumn](const long irow, const Slice s) {
        vecfuncT result;
        MADNESS_CHECK(s.start == 0 && s.end == -1 && s.step == 1);
        for (std::size_t i = s.start; i <= s.end + ncolumn; ++i) {
            result.push_back(Nij[ij(irow, i)]);
        }
        return result;
    };
    auto Nslice1 = [&Nij, &ij, &nrow](const Slice s, const long jcolumn) {
        vecfuncT result;
        MADNESS_CHECK(s.start == 0 && s.end == -1 && s.step == 1);
        for (std::size_t i = s.start; i <= s.end + nrow; ++i) {
            result.push_back(Nij[ij(i, jcolumn)]);
        }
        return result;
    };

    // corresponds to bra_batch and ket_batch, but without the ncf R^2
    // result[i]        =                      sum_{k}                ket[k] \int bra[k] vf[i]
    // result[rowbatch] = \sum_{columnbatches} sum_{k in columnbatch} ket[k] \int bra[k] vf[rowbatch]
    MADNESS_CHECK(batch.input.size() == 2);
    auto row_range = batch.input[0];            // corresponds to bra_batch
    auto column_range = batch.input[1];         // corresponds to f_batch
    vecfuncT to_dot_with_bra = column_range.copy_batch(mo_ket);
    vecfuncT to_dot_with_vf = row_range.copy_batch(mo_ket);

    vecfuncT resultcolumn(nrow);
    for (std::size_t irow = 0; irow < nrow; ++irow) {
        resultcolumn[irow] = dot(subworld, to_dot_with_vf,
                                 Nslice(irow, _));  // sum over columns result=sum_j ket[j] N[j,i]
    }
    vecfuncT resultrow(ncolumn);
    for (std::size_t icolumn = 0; icolumn < ncolumn; ++icolumn) {
        resultrow[icolumn] = dot(subworld, to_dot_with_bra,
                                 Nslice1(_, icolumn));  // sum over rows result=sum_i ket[i] N[j,i]
    }

    // !! NO TRUNCATION AT THIS POINT !!
    subworld.gop.fence();
    cpu1 = cpu_time();
    mul2_timer += long((cpu1 - cpu0) * 1000l);

    return std::make_pair(resultcolumn, resultrow);
}


template 
class Exchange<double_complex, 3>::ExchangeImpl;

template
class Exchange<double, 3>::ExchangeImpl;

template<> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

template<> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(
        0);

Exchange<double,3>::ExchangeImpl junkjunkjunk(World& world, const SCF *calc, const int ispin) {return Exchange<double,3>::ExchangeImpl(world, calc, ispin);}

} /* namespace madness */
