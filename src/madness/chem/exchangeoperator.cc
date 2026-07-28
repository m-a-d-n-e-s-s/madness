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
        bra_equals_ket_ = false;   // bra = conj(ket) != ket in general (complex)
    } else {
        if (ispin == 0) { // alpha spin
            mo_ket = calc->amo;        // deep copy necessary if T==double_complex
        } else if (ispin == 1) {  // beta spin
            mo_ket = calc->bmo;
        }
        mo_bra = mo_ket;
        bra_equals_ket_ = true;    // real, no nuclear correlation factor: bra == ket (moldft)
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
    bra_equals_ket_ = false;   // nemo: bra = R^2 * ket != ket
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Exchange<T, NDIM>::ExchangeImpl::operator()(
        const std::vector<Function<T, NDIM> >& vket) const {

    const bool needs_reconstructed_inputs =
            (algorithm_ == multiworld_efficient) ||
            (algorithm_ == large_memory) ||
            (algorithm_ == fetch_compute);

    // EXCH_OP split (printlevel>=4): attribute the "HF exchange" timer's exchange-operator part into
    // make_redundant (input prep + fences) / K_macrotask (mtask) / truncate (result). model §14.3.
    const double t_op_a = wall_time();
    if (needs_reconstructed_inputs) {
        reconstruct(world, mo_bra, false);
        reconstruct(world, mo_ket, false);
        world.gop.fence();
        reconstruct(world, vket);
        norm_tree(world, mo_bra, false);
        norm_tree(world, mo_ket, false);
        world.gop.fence();
        norm_tree(world, vket);
    } else {
        // Algorithms here drive multiplication via the redundant-tree mw-screening
        // kernel (mul_sparse + vmulXX). Pre-stage the inputs once here so subworld
        // kernels can skip make_redundant/fence around every mul_sparse call.
        // make_redundant routes through compress(redundant, ...), which populates
        // the per-node sum-coefficient + tree-norm structure that mw-screening uses
        // — no separate norm_tree pass is needed (the reconstructed branch above
        // does norm_tree only because reconstruct doesn't compute it).
        make_redundant(world, mo_bra, false);
        make_redundant(world, mo_ket, false);
        world.gop.fence();
        make_redundant(world, vket);
    }

    // pick your algorithm.
    // Note that the macrotask algorithm partitions the exchange matrix into tiles. The final truncation
    // takes place for each sum over tile elements, instead of the sum over all matrix elements.
    // Truncation and addition doesn't commute, so truncation is done after the final accumulation.
    // Other truncations are elementwise and are not affected.
    reset_timer();
    statistics=gather_statistics();
    double wall0=wall_time();
    double process_cpu0=process_cpu_time();
    vecfuncT Kf;
    if (algorithm_ == multiworld_efficient) {
        Kf = K_macrotask_efficient(vket, mul_tol);
    } else if (algorithm_ == small_memory_symmetric_mt
            or algorithm_ == small_memory_symmetric_mt_owner
            or algorithm_ == small_memory_symmetric_p2p_owner
            or algorithm_ == small_memory_mt_owner) {
        Kf = K_macrotask_efficient(vket, mul_tol);
    } else if (algorithm_ == multiworld_efficient_row or algorithm_ == fetch_compute
               ) {
        Kf = K_macrotask_efficient_row(vket, mul_tol);
    } else if (algorithm_ == small_memory) {
        Kf = K_small_memory(vket, mul_tol);     // Smaller memory algorithm ... possible 2x saving using i-j sym
    } else if (algorithm_ == small_memory_symmetric) {
        Kf = K_small_memory_symmetric(vket, mul_tol);
    } else if (algorithm_ == large_memory) {
        Kf = K_large_memory(vket, mul_tol);
    } else {
        MADNESS_EXCEPTION("unknown algorithm in exchangeoperator", 1);
    }
    if (printdebug()) {
        auto size = get_size(world, Kf);
        if (world.rank() == 0) print("total size of Kf before truncation", size);
    }
    const double t_op_c = wall_time();   // before truncate
    truncate(world, Kf);
    if (printdebug()) {
        auto size=get_size(world,Kf);
        if (world.rank()==0) print("total size of Kf after truncation",size);
    }
    double process_cpu1=process_cpu_time();
    double wall1=wall_time();
    elapsed_process_cpu_time=process_cpu1-process_cpu0;
    elapsed_time=wall1-wall0;
    if (printtimings_detail() and world.rank()==0) {
        printf("EXCH_OP k=%ld nmo=%lu | make_redundant=%.3f K_macrotask=%.3f truncate=%.3f (op total %.3f s)\n",
               long(FunctionDefaults<NDIM>::get_k()), (unsigned long)vket.size(),
               wall0 - t_op_a, t_op_c - wall0, wall1 - t_op_c, wall1 - t_op_a);
        fflush(stdout);
    }

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
    MacroTaskExchangeSimple xtask(nresult, lo, mul_tol, is_symmetric(),
                                  min_batch_size_, max_batch_size_, algorithm_, batch_granularity_);
    xtask.replicate_for_debug_ = replicate_for_debug_;
    xtask.accumulation_mode_ = accumulation_mode_;
    xtask.use_mflex_ = use_mflex_;
    xtask.mflex_max_exhaustive_ = mflex_max_exhaustive_;
    xtask.smallmem_mul_tol_ = smallmem_mul_tol_;  // TEMP debug knob
    // Graceful fallback: the cloud-batch fetch path pins one vf column-batch per
    // subworld via the row-owner partition (nbatch = min(nsubworld, n_MO)). When
    // nproc > n_MO that partition is degenerate — some ranks own no batch — which
    // the cloud path cannot fetch (it would deref an empty/missing record and
    // crash). Fall back to copy-via-universe in that case; the copy path and the
    // result reduction both handle the degenerate partition correctly.
    bool cloud_batch_active = use_cloud_batch_fetch_ and not replicate_for_debug_
                              and (algorithm_ == small_memory_mt_owner
                                   or algorithm_ == small_memory_symmetric_p2p_owner);
    // The asymmetric row-owner partition crashes when nproc > n_MO (empty batch
    // deref) → fall back to copy. The symmetric p2p-owner path handles nproc > n_MO
    // fine on the CLOUD path (M = min((base+level-1)*nproc, n_MO) < nproc batches,
    // extra ranks idle) now that the result-accumulation bad_weak_ptr is fixed — and
    // it MUST stay on the cloud path because its copy fallback is incorrect (uses the
    // held-vf accessor, wrong for the rotating triangular vf). So this guard is
    // asymmetric-only.
    if (cloud_batch_active and algorithm_ == small_memory_mt_owner and world.size() > nresult) {
        cloud_batch_active = false;
        if (world.rank()==0)
            print("WARNING: nproc (", world.size(), ") > n_MO (", nresult,
                  "); the row-owner cloud-batch fetch partition is degenerate, "
                  "falling back to copy-via-universe (cloud-batch fetch requires n_MO >= nproc)");
    }
    // smallmem_sym_p2p_owner stores ONE set and requires bra==ket (moldft). Nemo
    // (bra=R^2*ket != ket) is unsupported: the one-set store would be wrong and the
    // copy fallback is also incorrect for this algorithm. Fail loud — use the
    // asymmetric algorithm for nemo.
    if (algorithm_ == small_memory_symmetric_p2p_owner) {
        MADNESS_CHECK_THROW(bra_equals_ket_,
            "smallmem_sym_p2p_owner requires bra==ket (moldft / no nuclear correlation "
            "factor); use the asymmetric algorithm for nemo");
        // The triangular one-set algorithm ONLY works on the cloud-batch fetch path:
        // operator() fetches the rotating bra/vf batches by record key. With the cloud
        // path off (hfex_use_cloud_batch_fetch=false, or replicate_for_debug_),
        // cloud_batch_active is false and operator() would fall into the owner-aware
        // copy branch (acquire_current_vf), whose held-vf accessor is WRONG for the
        // rotating triangular vf → silently wrong energy. Fail loud instead. (The
        // nproc>n_MO fallback above is asymmetric-only, so for this algorithm
        // cloud_batch_active is false only when the user disabled the cloud path.)
        MADNESS_CHECK_THROW(cloud_batch_active,
            "smallmem_sym_p2p_owner requires the cloud-batch fetch path: set "
            "hfex_use_cloud_batch_fetch=true and do not use replicate_for_debug "
            "(the copy/replicate fallback is incorrect for the rotating triangular vf)");
    }
    xtask.use_cloud_batch_fetch_ = cloud_batch_active;
    xtask.log_diagnostics_ = printdebug();          // owner_map + batch record dump on universe rank 0
    xtask.universe_rank_ = world.rank();            // keys the per-task profile file (MAD_EXCH_TASK_PROFILE)
    xtask.cost_aware_assign_ = cost_aware_assign_;  // cost-aware owner assignment (from hfex_cost_aware_assign)
    if (taskq) taskq->set_printlevel(printlevel);
    auto effective_policy = macro_task_info;
    if (replicate_for_debug_) {
        effective_policy.storage_policy = MacroTaskInfo::StoreFunction;
        if (world.rank()==0) print("DEBUG: using StoreFunction policy to pre-replicate all data (zero communication during tasks)");
    }
    // cloud-batch fetch path (small_memory_mt_owner only): input batches live as
    // owner-pinned serialized records in the cloud's batch container. The argtuple
    // is still stored as pointers (shape); batches are stored/fetched separately.
    if (cloud_batch_active) {
        effective_policy.storage_policy = MacroTaskInfo::StoreFunctionBatched;
        effective_policy.ptr_target_distribution_policy = DistributionType::Distributed;
        if (world.rank()==0) print("using StoreFunctionBatched policy: input batches fetched owner-to-owner from the cloud");
    }
    auto taskq_factory = MacroTaskQFactory(world).set_printlevel(printlevel).set_policy(effective_policy);
    if (algorithm_ == small_memory_symmetric_mt_owner
     or algorithm_ == small_memory_symmetric_p2p_owner
     or algorithm_ == small_memory_mt_owner) {
        taskq_factory.set_nworld(world.size());
    }
    // TODO(option-ii): once small_memory_mt_owner is validated, evaluate setting
    //   effective_policy.ptr_target_distribution_policy = DistributionType::NodeReplicated;
    // for this algorithm so cross-rank ket fetches become single-owner reads via is_replicated.

    // construct MacroTask with or without user-provided taskq -> deferred execution or immediate execution
    // KMACRO split (printlevel>=4): ctor(taskq build) / mtask(run: store+span+finalize) / postfence /
    // costsync / stats. mtask − span (profiler) − finalize (RUN_PHASES) = store + framework. model §14.3.
    const double t_km_ctor = wall_time();
    auto mtask = (taskq) ? MacroTask(world, xtask, taskq)
                 : MacroTask(world, xtask, taskq_factory);

    // wire per-call cloud-batch debug logging (BATCH_REQ / BATCH_DESER / ...)
    // through to the cloud held by the taskq. Narrower than cloud.set_debug.
    if (mtask.get_taskq()) mtask.get_taskq()->cloud.set_batch_debug(printdebug());

    // deferred execution if a taskq is provided by the user
    const double t_km0 = wall_time();
    vecfuncT Kf = mtask(vf, mo_bra, mo_ket);
    const double t_km1 = wall_time();
    world.gop.fence();
    const double t_km2 = wall_time();

    // cost-aware load balancing (opt-in, sym_p2p): synchronize this call's per-task costs across
    // ranks so every rank holds an identical cost matrix for the NEXT call's deterministic
    // assignment. Each task ran on exactly one rank, so the sum unions the per-rank contributions.
    if ((MacroTaskExchangeSimple::cost_aware_resolve(cost_aware_assign_) or printdebug())
        and algorithm_ == small_memory_symmetric_p2p_owner) {
        auto& local = MacroTaskExchangeSimple::sym_cost_local_ref();
        if (not local.empty()) {
            world.gop.sum(local.data(), local.size());          // union per-rank task costs
            // debug (print_level>=10): dump the synchronized per-call cost matrix (rank 0) BEFORE the EMA blend
            if (world.rank() == 0 and printdebug())
                MacroTaskExchangeSimple::sym_dump_cost(FunctionDefaults<NDIM>::get_k());
            MacroTaskExchangeSimple::sym_commit_cost_global();  // EMA-blend into the next-call reference
        }
    }
    const double t_km3 = wall_time();

    statistics=gather_statistics();
    statistics["cloud"]=mtask.get_taskq()->get_cloud_statistics();
    statistics["macrotaskq"]=mtask.get_taskq()->get_taskq_statistics();
    if (printtimings_detail() and world.rank()==0) {
        printf("KMACRO k=%ld | ctor=%.3f mtask=%.3f postfence=%.3f costsync=%.3f stats=%.3f\n",
               long(FunctionDefaults<NDIM>::get_k()),
               t_km0 - t_km_ctor, t_km1 - t_km0, t_km2 - t_km1, t_km3 - t_km2, wall_time() - t_km3);
        fflush(stdout);
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
        vecfuncT psif = mul_sparse(world, mo_bra[i], vket, mul_tol*0.1, true, false); /// was vtol
        truncate(world, psif);
        psif = apply(world, *poisson.get(), psif);
        truncate(world, psif);
        make_redundant(world, psif, true);
        psif = mul_sparse(world, mo_ket[i], psif, mul_tol*0.1, true, false); /// was vtol
        gaxpy(world, 1.0, Kf, 1.0, psif);
    }
    truncate(world, Kf);
    return Kf;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> > Exchange<T, NDIM>::ExchangeImpl::K_small_memory_symmetric(const vecfuncT& vket,
                                                                                            const double mul_tol) const {

    if (!is_symmetric()) return K_small_memory(vket, mul_tol);

    const long nocc = mo_ket.size();
    const long nf = vket.size();
    if (nocc != nf) return K_small_memory(vket, mul_tol);

    vecfuncT Kf = zero_functions_compressed<T, NDIM>(world, nf);
    auto poisson = set_poisson(world, lo);

    for (int i = 0; i < nocc; ++i) {
        vecfuncT vket_subset(vket.begin(), vket.begin() + i + 1);
        vecfuncT psif = mul_sparse(world, mo_bra[i], vket_subset, mul_tol*0.1, true, false);

        truncate(world, psif);
        psif = apply(world, *poisson.get(), psif);
        truncate(world, psif);
        make_redundant(world, psif, true);

        // Build one full update vector and accumulate with vector gaxpy
        vecfuncT update_i = zero_functions_compressed<T, NDIM>(world, nf);
        compress(world, update_i);

        // Row contribution: update_i[j] += ket[i] * N_ij for j <= i
        vecfuncT row_contrib = mul_sparse(world, mo_ket[i], psif, mul_tol*0.1, true, false);
        compress(world, row_contrib);
        for (int j = 0; j <= i; ++j) {
            update_i[j] += row_contrib[j];
        }

        // Mirrored contribution: update_i[i] += ket[j] * N_ij for j < i
        for (int j = 0; j < i; ++j) {
            vecfuncT psif_single(1, psif[j]);
            vecfuncT mirrored = mul_sparse(world, mo_ket[j], psif_single, mul_tol*0.1, true, false);
            compress(world, mirrored);
            update_i[i] += mirrored[0];
        }

        truncate(world, update_i);
        gaxpy(world, 1.0, Kf, 1.0, update_i);
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

    double cpu0 = process_cpu_time();
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
    double cpu1 = process_cpu_time();
    mul1_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = process_cpu_time();
    psif = apply(world, *poisson.get(), psif);
    truncate(world, psif);
    cpu1 = process_cpu_time();
    apply_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = process_cpu_time();
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
    cpu1 = process_cpu_time();
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
    double cpu0 = process_cpu_time();
    std::vector<vecfuncT> orbital_product = matrix_mul_sparse<T, T, NDIM>(subworld, bra_batch, vf_batch, mul_tol);
    vecfuncT orbital_product_flat = flatten(orbital_product); // convert into a flattened vector
    truncate(subworld, orbital_product_flat);
    double cpu1 = process_cpu_time();
    mul1_timer += long((cpu1 - cpu0) * 1000l);

    cpu0 = process_cpu_time();
    auto poisson = set_poisson(subworld, lo);
    vecfuncT Nij = apply(subworld, *poisson.get(), orbital_product_flat);
    truncate(subworld, Nij);
    cpu1 = process_cpu_time();
    apply_timer += long((cpu1 - cpu0) * 1000l);

    // accumulate columns:      resultrow(i)=\sum_j j N_ij
    // accumulate rows:      resultcolumn(j)=\sum_i i N_ij
    cpu0 = process_cpu_time();

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
    cpu1 = process_cpu_time();
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
