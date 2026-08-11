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
    // Owner-pinned placement applies to the symmetric case only: it needs bra == ket so
    // that one set of batches serves every operand role and both of a task's batches are
    // owned by some rank. The asymmetric case keeps the size-driven partition.
    MacroTaskExchangeSimple xtask(nresult, lo, mul_tol, is_symmetric(),
                                  /*owner_pinned=*/is_symmetric(), batch_granularity_,
                                  world.rank(), accumulation_mode_, cost_aware_assign_);
    // The owner-pinned path stores the orbitals as batches itself and fetches the two it
    // needs per task, so the cloud must hold pointers rather than copy every operand into
    // every subworld. The algorithm therefore fixes the storage policy.
    const MacroTaskInfo policy = is_symmetric() ? MacroTaskInfo::preset("small_memory_owner")
                                                : macro_task_info;
    if (taskq) taskq->set_printlevel(printlevel);

    // construct MacroTask with or without user-provided taskq -> deferred execution or immediate execution
    auto mtask = (taskq) ? MacroTask(world, xtask, taskq)
                 : MacroTask(world, xtask, MacroTaskQFactory(world).set_printlevel(printlevel).set_policy(policy));

    // deferred execution if a taskq is provided by the user
    vecfuncT Kf = mtask(vf, mo_bra, mo_ket);
    world.gop.fence();

    // Every tile ran on exactly one rank, so summing the per-rank cost buffers unions them and
    // leaves every rank holding the same matrix -- which is what lets the next call's placement
    // be computed independently everywhere and still agree.
    if (is_symmetric() and cost_aware_assign_) {
        auto& costs = MacroTaskExchangeSimple::cost_this_call();
        if (not costs.empty()) {
            world.gop.sum(costs.data(), costs.size());
            // dump before committing: this is the matrix just measured, which is what the next
            // application will place against
            if (world.rank() == 0 and exch_task_profile_enabled())
                exch_write_cost_matrix(MacroTaskExchangeSimple::exchange_call_index(),
                                       FunctionDefaults<NDIM>::get_k(),
                                       MacroTaskExchangeSimple::cost_matrix_dimension(), costs);
            MacroTaskExchangeSimple::commit_cost_reference();
        }
    }

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
/// \param ket_rows     the orbitals to premultiply with, over the bra/row range
/// \param ket_columns  the orbitals to premultiply with, over the vf/column range
/// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
/// \param vf_batch     the argument of the exchange operator
template<typename T, std::size_t NDIM>
std::pair<std::vector<Function<T, NDIM>>, std::vector<Function<T, NDIM>>>
Exchange<T, NDIM>::ExchangeImpl::MacroTaskExchangeSimple::compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                                                                                          const vecfuncT& ket_rows,
                                                                                          const vecfuncT& ket_columns,
                                                                                          const vecfuncT& bra_batch,   // batched
                                                                                          const vecfuncT& vf_batch) const { // batched
    // result[i] = sum_k ket[k] \int bra[k] vf[i], so the tile contributes to both of its
    // ranges: the ket over the rows accumulates into the columns, and the other way round
    MADNESS_CHECK_THROW(ket_rows.size() == bra_batch.size(),
                        "symmetric offdiagonal tile: ket/bra row size mismatch");
    MADNESS_CHECK_THROW(ket_columns.size() == vf_batch.size(),
                        "symmetric offdiagonal tile: ket/vf column size mismatch");

    const long nrow = long(bra_batch.size());
    const long ncolumn = long(vf_batch.size());
    vecfuncT resultcolumn(nrow);                                                   // maps to the row range
    vecfuncT resultrow = zero_functions_compressed<T, NDIM>(subworld, ncolumn);    // maps to the column range
    auto poisson = set_poisson(subworld, lo);

    // per-stage wall, for the profiler only; `tick` is a no-op when it is off
    const bool prof_on = profile_active();
    auto tick = [&](double& acc, const double t0) { if (prof_on) acc += wall_time() - t0; };

    for (long irow = 0; irow < nrow; ++irow) {
        double cpu0 = cpu_time();
        double w0 = prof_on ? wall_time() : 0.0;
        vecfuncT Nij = mul_sparse(subworld, bra_batch[irow], vf_batch, mul_tol);
        tick(prof_.mul1_wall, w0);
        w0 = prof_on ? wall_time() : 0.0;
        truncate(subworld, Nij);
        tick(prof_.truncate_wall, w0);
        double cpu1 = cpu_time();
        mul1_timer += long((cpu1 - cpu0) * 1000l);

        cpu0 = cpu_time();
        w0 = prof_on ? wall_time() : 0.0;
        Nij = apply(subworld, *poisson.get(), Nij);
        tick(prof_.apply_wall, w0);
        w0 = prof_on ? wall_time() : 0.0;
        truncate(subworld, Nij);
        tick(prof_.truncate_wall, w0);
        cpu1 = cpu_time();
        apply_timer += long((cpu1 - cpu0) * 1000l);

        cpu0 = cpu_time();
        // every row contributes to every column, so the column result accumulates ...
        w0 = prof_on ? wall_time() : 0.0;
        vecfuncT row_update = mul_sparse(subworld, ket_rows[irow], Nij, mul_tol);
        tick(prof_.mul2_wall, w0);
        compress(subworld, row_update);
        gaxpy(subworld, 1.0, resultrow, 1.0, row_update);
        // ... while each row's own entry is written exactly once
        w0 = prof_on ? wall_time() : 0.0;
        resultcolumn[irow] = dot(subworld, ket_columns, Nij, true, true, mul_tol);
        tick(prof_.mul2_wall, w0);
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
