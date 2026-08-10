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
    statistics=gather_statistics();
    double cpu0=wall_time();
    vecfuncT Kf;
    if (algorithm_ == multiworld_efficient) {
        Kf = K_macrotask_efficient(vket, mul_tol);
    } else if (algorithm_ == multiworld_efficient_row or algorithm_ == fetch_compute) {
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
    double cpu1=wall_time();
    elapsed_time=cpu1-cpu0;

    if (printtimings_detail()) print_timer(world);
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
    if (taskq) taskq->set_printlevel(printlevel);

    // construct MacroTask with or without user-provided taskq -> deferred execution or immediate execution
    auto mtask = (taskq) ? MacroTask(world, xtask, taskq)
                 : MacroTask(world, xtask, MacroTaskQFactory(world).set_printlevel(printlevel).set_policy(macro_task_info));

    // deferred execution if a taskq is provided by the user
    vecfuncT Kf = mtask(vf, mo_bra, mo_ket);
    world.gop.fence();
    statistics=gather_statistics();
    statistics["cloud"]=mtask.get_taskq()->get_cloud_statistics();
    statistics["macrotaskq"]=mtask.get_taskq()->get_taskq_statistics();

    return Kf;
}

/// compute each row of the exchange matrix in different subworlds
template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::ExchangeImpl::K_macrotask_efficient_row(const vecfuncT& vf, const double mul_tol) const {

    if (world.rank()==0 and printdebug()) print("\nentering macrotask_efficient_row version:");

    // the result is a vector of functions living in the universe
    const long nresult = vf.size();
    MacroTaskExchangeRow xtask(nresult, lo, mul_tol, algorithm_);

    // print the size of the amos
    if (printdebug()) {
        auto size=get_size(world,vf);
        if (world.rank()==0) print("total size of vf before iteration",size);
    }

    if (taskq) taskq->set_printlevel(printlevel);

    // construct MacroTask with or without user-provided taskq -> deferred execution or immediate execution
    auto mtask = (taskq) ? MacroTask(world, xtask, taskq)
                 : MacroTask(world, xtask, MacroTaskQFactory(world).set_printlevel(printlevel)
                     .set_policy(macro_task_info));

    // deferred execution if a taskq is provided by the user
    vecfuncT Kf = mtask(vf, mo_bra, mo_ket);
    world.gop.fence();
    statistics=gather_statistics();
    statistics["cloud"]=mtask.get_taskq()->get_cloud_statistics();
    statistics["macrotaskq"]=mtask.get_taskq()->get_taskq_statistics();

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

/// Streams the tile one bra row at a time, so only one row of intermediates is live at
/// once where building the tile in one go holds all nrow*ncolumn of them. Row `irow`
/// builds N_ij = P(bra[irow] vf[j]) over the whole column range and contributes to both
/// result ranges: ket[irow] N_ij to column j, and the sum over j of ket[j] N_ij to row irow.
/// \param subworld     the world we're computing in
/// \param mo_ket       the orbitals to premultiply with, not batched
/// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
/// \param vf_batch     the argument of the exchange operator
template<typename T, std::size_t NDIM>
std::pair<std::vector<Function<T, NDIM>>, std::vector<Function<T, NDIM>>>
Exchange<T, NDIM>::ExchangeImpl::MacroTaskExchangeSimple::compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                          const vecfuncT& mo_ket,      // not batched
                                                                                          const vecfuncT& bra_batch,   // batched
                                                                                          const vecfuncT& vf_batch) const { // batched
    MADNESS_CHECK(batch.input.size() == 2);
    // input[1] is the bra/row range and input[0] the vf/column range -- the same labelling
    // ExchangeImpl::operator() uses when it accumulates the two results back into Kf
    const auto row_range = batch.input[1];
    const auto column_range = batch.input[0];
    MADNESS_CHECK_THROW(row_range.size() == long(bra_batch.size()),
                        "symmetric offdiagonal tile: row range does not match the bra batch");
    MADNESS_CHECK_THROW(column_range.size() == long(vf_batch.size()),
                        "symmetric offdiagonal tile: column range does not match the vf batch");

    // result[i] = sum_k ket[k] \int bra[k] vf[i], so the tile needs the ket over both of
    // its ranges: over the rows to accumulate into the columns, and vice versa
    const vecfuncT ket_rows = row_range.copy_batch(mo_ket);
    const vecfuncT ket_columns = column_range.copy_batch(mo_ket);

    const long nrow = long(bra_batch.size());
    const long ncolumn = long(vf_batch.size());
    vecfuncT resultcolumn(nrow);                                                   // maps to the row range
    vecfuncT resultrow = zero_functions_compressed<T, NDIM>(subworld, ncolumn);    // maps to the column range
    auto poisson = set_poisson(subworld, lo);

    for (long irow = 0; irow < nrow; ++irow) {
        double cpu0 = cpu_time();
        vecfuncT Nij = mul_sparse(subworld, bra_batch[irow], vf_batch, mul_tol);
        truncate(subworld, Nij);
        double cpu1 = cpu_time();
        mul1_timer += long((cpu1 - cpu0) * 1000l);

        cpu0 = cpu_time();
        Nij = apply(subworld, *poisson.get(), Nij);
        truncate(subworld, Nij);
        cpu1 = cpu_time();
        apply_timer += long((cpu1 - cpu0) * 1000l);

        cpu0 = cpu_time();
        // every row contributes to every column, so the column result accumulates ...
        vecfuncT row_update = mul_sparse(subworld, ket_rows[irow], Nij, mul_tol);
        compress(subworld, row_update);
        gaxpy(subworld, 1.0, resultrow, 1.0, row_update);
        // ... while each row's own entry is written exactly once
        resultcolumn[irow] = dot(subworld, ket_columns, Nij, true, true, mul_tol);
        cpu1 = cpu_time();
        mul2_timer += long((cpu1 - cpu0) * 1000l);
    }

    // !! NO TRUNCATION AT THIS POINT !!
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
