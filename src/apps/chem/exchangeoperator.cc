//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <chem/exchangeoperator.h>

#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/orbital_partitioner.h>

using namespace madness;

namespace madness {


template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const SCF* calc, const int ispin)
        : world(world), symmetric_(false), lo(calc->param.lo()) {
    if (ispin==0) { // alpha spin
        mo_ket=convert<double,T,NDIM>(world,calc->amo);		// deep copy necessary if T==double_complex
    } else if (ispin==1) {  // beta spin
        mo_ket=convert<double,T,NDIM>(world,calc->bmo);
    }
    mo_bra=conj(world,mo_ket);
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const Nemo* nemo, const int ispin) // @suppress("Class members should be properly initialized")
    : Exchange<T,NDIM>(world,nemo->get_calc().get(),ispin) {

    if (ispin==0) { // alpha spin
        mo_ket=convert<double,T,NDIM>(world,nemo->get_calc()->amo);        // deep copy necessary if T==double_complex
    } else if (ispin==1) {  // beta spin
        mo_ket=convert<double,T,NDIM>(world,nemo->get_calc()->bmo);        // deep copy necessary if T==double_complex
    }

    mo_bra=mul(world,nemo->ncf->square(),mo_ket);
    truncate(world,mo_bra);
}

template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::operator()(
		const std::vector<Function<T,NDIM> >& vket, const double mul_tol) const {

    reconstruct(world, mo_bra,false);
    reconstruct(world, mo_ket,false);
    world.gop.fence();
    norm_tree(world, mo_bra,false);
    norm_tree(world, mo_ket,false);
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
    if (algorithm_==multiworld_efficient) {
        Kf=K_macrotask_efficient(vket,mul_tol);
    } else if (algorithm_==small_memory) {
        Kf=K_small_memory(vket,mul_tol);     // Smaller memory algorithm ... possible 2x saving using i-j sym
    } else if (algorithm_==large_memory) {
    	Kf=K_large_memory(vket,mul_tol);
    } else {
        MADNESS_EXCEPTION("unknown algorithm in exchangeoperator",1);
    }
    truncate(world,Kf);
    if (printlevel>=3) print_timer(world);
    return Kf;
}

    /// apply the exchange operator by tiling the exchange matrix

    /// compute the matrix N_{ik} = N \phi_i \phi_k by tiles, with i,k \in batches A,B,
    /// do a local reduce within the tiles: K_{iB} = \sum_{k \in batch B} \phi_k N_{ik}
    /// and a universe-wide reduce of the tiles: K\phi_i = \sum_{batches B} K_{iB}
    /// saving up to half of the cpu time compared to the naive algorithm
    /// \tparam T       number type
    /// \tparam NDIM    physical dimension of the argument vket
    /// \param vket     argument of the exchange operator
    /// \param mul_tol  cutoff parameter for sparse multiplication
    /// \return         the exchange operator applied on vket
    template<typename T, std::size_t NDIM>
    std::vector<Function<T,NDIM> > Exchange<T, NDIM>::K_macrotask_efficient(const vecfuncT &vket, const double mul_tol) const {

        if (printdebug()) print("\nentering macrotask_efficient version:");

        const long nocc=mo_ket.size();
        const long nsubworld=world.size();

        // the result is a vector of functions living in the universe
        const long nf=vket.size();
        vecfuncT Kf=zero_functions_compressed<T,NDIM>(world,nf);
        {

            double wall0 = wall_time();
            MacroTaskQ taskq(world, nsubworld);
            long start = 0;
//            for (int i = 0; i < nocc; ++i) {
//                taskq.cloud.store(world, mo_bra[i], i);
//                taskq.cloud.store(world, mo_ket[i], i + nocc);
//            }
//            start += 2 * nocc;
//            for (int i = 0; i < nf; ++i) {
//                taskq.cloud.store(world, vket[i], i + start);
//                taskq.cloud.store(world, Kf[i].get_impl(), start + i + nf); // store pointer to FunctionImpl
//            }
            double wall1 = wall_time();
            if (do_print_timings()) printf("wall time for storing %4.1fs\n", wall1 - wall0);

            // partition the exchange matrix into tiles
            // result[rowbatch] = \sum_{columnbatches} sum_{k in columnbatch} ket[k] \int bra[k] vf[rowbatch]
            std::vector<std::pair<long, long> > rowranges, columnranges;

            // symmetric exchange matrix, symmetric tiles
            if (is_symmetric()) {
                rowranges = OrbitalPartitioner::partition_for_exchange(5, nsubworld, nocc);
//            rowranges=OrbitalPartitioner::partition_for_exchange(5,nsubworld,nocc);
                columnranges = rowranges;
            } else {
                rowranges = OrbitalPartitioner::partition_for_exchange(5, nsubworld, nocc);
                columnranges = OrbitalPartitioner::partition_for_exchange(5, nsubworld, nf);
            }
            if (printdebug()) {
                print("symmetric", symmetric_);
                print("splitting nf into", rowranges.size(), "ranges");
                for (auto &r : rowranges) printf(" %2ld to %2ld", r.first, r.second);
                print("\nsplitting nocc into", columnranges.size(), "ranges");
                for (auto &r : columnranges) printf(" %2ld to %2ld", r.first, r.second);
                print("");
            }

            MacroTaskBase::taskqT vtask;
            for (auto &row_range : rowranges) {
                for (auto &column_range : columnranges) {
                    // if the exchange matrix is symmetric compute only the the upper triangular matrix
                    if (is_symmetric() and row_range.first > column_range.first) continue;

                    MacroTaskExchange task(row_range, column_range, nocc, nf, lo, mul_tol, symmetric_);
                    vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskExchange(task)));
                }
            }
            world.gop.fence();
            // sort the tasks according to their tile size, assuming that's a proxy for how much time they take
            std::sort(vtask.begin(), vtask.end(), [](auto a, auto b) { return a->get_priority() > b->get_priority(); });

            auto current_pmap = FunctionDefaults<3>::get_pmap();
            FunctionDefaults<3>::set_default_pmap(taskq.get_subworld());
            taskq.run_all(vtask);
            FunctionDefaults<3>::set_pmap(current_pmap);

            truncate(world, Kf);
            if (printlevel >= 3) taskq.cloud.print_timings(world);
            taskq.get_subworld().gop.fence();
        }
        world.gop.fence();
        return Kf;
    }


template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::K_small_memory(const vecfuncT& vket, const double mul_tol) const {

	const long nocc=mo_ket.size();
	const long nf=vket.size();
	vecfuncT Kf=zero_functions_compressed<T,NDIM>(world,nf);
	auto poisson=set_poisson(world,lo);

	for(int i=0; i<nocc; ++i){
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
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::K_large_memory(const vecfuncT& vket, const double mul_tol) const {    // Larger memory algorithm ... use i-j sym if psi==f

    auto poisson=set_poisson(world,lo);
    vecfuncT result=compute_K_tile(world,mo_bra,mo_ket,vket,poisson,is_symmetric(),mul_tol);
    truncate(world,result);
    return result;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::compute_K_tile(World &world, const vecfuncT &mo_bra, const vecfuncT &mo_ket,
                                  const vecfuncT &vket, std::shared_ptr<real_convolution_3d> poisson,
                                  const bool symmetric, const double mul_tol) {

    double cpu0=cpu_time();
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
    double cpu1=cpu_time();
    mul1_timer+=long((cpu1-cpu0)*1000l);

    cpu0=cpu_time();
    psif = apply(world, *poisson.get(), psif);
    truncate(world, psif);
    cpu1=cpu_time();
    apply_timer+=long((cpu1-cpu0)*1000l);

    cpu0=cpu_time();
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
    cpu1=cpu_time();
    mul2_timer+=long((cpu1-cpu0)*1000l);
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


template<typename T, std::size_t NDIM>
void Exchange<T,NDIM>::MacroTaskExchange::run(World& subworld, Cloud& cloud,
                                              std::vector<std::shared_ptr<MacroTaskBase> >& taskq) {

    if (not poisson) poisson = Exchange<T,NDIM>::set_poisson(subworld,lo);
    // the argument of the exchange operator is the ket vector

    // load bra and ket if not already loaded
    if (not mo_bra.get()) {
        mo_bra.reset(new vecfuncT(nocc));
        for (int i = 0; i < nocc; ++i) (*mo_bra)[i] = cloud.load<functionT>(subworld, i);
    }
    if (not mo_ket.get()) {
        mo_ket.reset(new vecfuncT(nocc));
        for (int i = 0; i < nocc; ++i) (*mo_ket)[i] = cloud.load<functionT>(subworld, i + nocc);
    }

    if (not vf.get()) {
        vf.reset(new vecfuncT(nf));
        for (int i = 0; i < nf; ++i) (*vf)[i] = cloud.load<functionT>(subworld, i + 2*nocc);
    }

    // make universe-living Kf accessible here in the subworld for result accumulation
    vecfuncT Kf(nf);
    for (int i=0; i<nf; ++i) {
        std::shared_ptr<FunctionImpl<T,NDIM> > rimpl;
        cloud.load(subworld,rimpl,i+2*nocc+nf);
        Kf[i].set_impl(rimpl);
    }

    // compute the tile [column_range,row_range], corresponding to bra[nrow], ket[ncolumn]
    // partition the exchange matrix into tiles
    // result[i]        =                      sum_{k}                ket[k] \int bra[k] vf[i]
    // result[rowbatch] = \sum_{columnbatches} sum_{k in columnbatch} ket[k] \int bra[k] vf[rowbatch]
    vecfuncT ket_batch(mo_ket->begin() + row_range.first, mo_ket->begin() + row_range.second);
    vecfuncT bra_batch(mo_bra->begin() + row_range.first, mo_bra->begin() + row_range.second);
    vecfuncT vf_batch (vf->begin()     + column_range.first,    vf->begin()     + column_range.second);

    if (symmetric and (row_range == column_range)) {
        vecfuncT resultcolumn = compute_diagonal_batch_in_symmetric_matrix(subworld, bra_batch, ket_batch, vf_batch);

        for (int i = row_range.first; i < row_range.second; ++i)
            Kf[i]+=resultcolumn[i - row_range.first];

    } else if (symmetric and !(row_range == column_range)) {
        auto [resultcolumn,resultrow]=compute_offdiagonal_batch_in_symmetric_matrix(subworld, bra_batch, ket_batch, vf_batch);

        for (int i = column_range.first; i < column_range.second; ++i)
            Kf[i]+=resultrow[i - column_range.first];
        for (int i = row_range.first; i < row_range.second; ++i)
            Kf[i]+=resultcolumn[i - row_range.first];
    } else {
        vecfuncT resultcolumn = compute_batch_in_asymmetric_matrix(subworld, bra_batch, ket_batch, vf_batch);

        for (int i = column_range.first; i < column_range.second; ++i)
            Kf[i]+=resultcolumn[i - column_range.first];

    }
    subworld.gop.fence();
}

    /// compute a batch of the exchange matrix, with non-identical ranges

    /// \param subworld     the world we're computing in
    /// \param cloud        where to store the results
    /// \param bra_batch    the bra batch of orbitals (including the nuclear correlation factor square)
    /// \param ket_batch    the ket batch of orbitals, also the orbitals to premultiply with
    /// \param vf_batch     the argument of the exchange operator
    template<typename T, std::size_t NDIM>
    std::pair<std::vector<Function<T,NDIM>>, std::vector<Function<T,NDIM>> >
    Exchange<T,NDIM>::MacroTaskExchange::compute_offdiagonal_batch_in_symmetric_matrix(World& subworld,
                           const vecfuncT& bra_batch, const vecfuncT& ket_batch, const vecfuncT& vf_batch) const {
        // orbital_product is a vector of vectors
        double cpu0 = cpu_time();
        std::vector<vecfuncT> orbital_product=matrix_mul_sparse<T,T,NDIM>(subworld,bra_batch,vf_batch,mul_tol);
        vecfuncT orbital_product_flat=flatten(orbital_product); // convert into a flattened vector
        truncate(subworld, orbital_product_flat);
        double cpu1=cpu_time();
        mul1_timer+=long((cpu1-cpu0)*1000l);

        cpu0=cpu_time();
        vecfuncT Nij = apply(subworld, *poisson.get(), orbital_product_flat);
        truncate(subworld, Nij);
        cpu1=cpu_time();
        apply_timer+=long((cpu1-cpu0)*1000l);

        // accumulate columns:      resultrow(i)=\sum_j j N_ij
        // accumulate rows:      resultcolumn(j)=\sum_i i N_ij
        cpu0=cpu_time();

        // some helper functions
        std::size_t nrow=bra_batch.size();
        std::size_t ncolumn=vf_batch.size();
        auto ij = [&nrow, &ncolumn](const int i, const int j) {return i*ncolumn+j;};

        auto Nslice = [&Nij, &ij, &ncolumn] (const long irow, const Slice s) {
            vecfuncT result;
            MADNESS_CHECK(s.start==0 && s.end==-1 && s.step==1);
            for (std::size_t i=s.start; i<=s.end+ncolumn; ++i) {
                result.push_back(Nij[ij(irow,i)]);
            }
            return result;
        };
        auto Nslice1 = [&Nij, &ij, &nrow] (const Slice s, const long jcolumn) {
            vecfuncT result;
            MADNESS_CHECK(s.start==0 && s.end==-1 && s.step==1);
            for (std::size_t i=s.start; i<=s.end+nrow; ++i) {
                result.push_back(Nij[ij(i,jcolumn)]);
            }
            return result;
        };

        // corresponds to bra_batch and ket_batch, but without the ncf R^2
        // result[i]        =                      sum_{k}                ket[k] \int bra[k] vf[i]
        // result[rowbatch] = \sum_{columnbatches} sum_{k in columnbatch} ket[k] \int bra[k] vf[rowbatch]
        vecfuncT to_dot_with_bra(mo_ket->begin()+row_range.first,mo_ket->begin()+row_range.second);      // regular
        vecfuncT to_dot_with_vf (mo_ket->begin()+column_range.first,mo_ket->begin()+column_range.second);            // transposed

        vecfuncT resultcolumn(nrow);
        for (std::size_t irow=0; irow<nrow; ++irow) {
            resultcolumn[irow]=dot(subworld,to_dot_with_vf,Nslice(irow,_));  // sum over columns result=sum_j ket[j] N[j,i]
        }
        vecfuncT resultrow(ncolumn);
        for (std::size_t icolumn=0; icolumn<ncolumn; ++icolumn) {
            resultrow[icolumn]=dot(subworld,to_dot_with_bra,Nslice1(_,icolumn));  // sum over rows result=sum_i ket[i] N[j,i]
        }

        // !! NO TRUNCATION AT THIS POINT !!
        subworld.gop.fence();
        cpu1=cpu_time();
        mul2_timer+=long((cpu1-cpu0)*1000l);

        return std::make_pair(resultcolumn,resultrow);
    }


    template class Exchange<double_complex,3>;
    template class Exchange<double,3>;

    template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
    template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

//template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
//template <> Spinlock WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);


} /* namespace madness */
