//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <chem/exchangeoperator.h>

#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/orbital_partitioner.h>

using namespace madness;

namespace madness {


template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const SCF* calc, const int ispin)
        : world(world), small_memory_(true), same_(false), lo(calc->param.lo()) {
    if (ispin==0) { // alpha spin
        mo_ket=convert<double,T,NDIM>(world,calc->amo);		// deep copy necessary if T==double_complex
        occ=calc->aocc;
    } else if (ispin==1) {  // beta spin
        mo_ket=convert<double,T,NDIM>(world,calc->bmo);
        occ=calc->bocc;
    }
    mo_bra=conj(world,mo_ket);
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, calc->param.lo(), calc->param.econv()));
}

template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const Nemo* nemo, const int ispin) // @suppress("Class members should be properly initialized")
    : Exchange<T,NDIM>(world,nemo->get_calc().get(),ispin) {

    if (ispin==0) { // alpha spin
        mo_ket=convert<double,T,NDIM>(world,nemo->get_calc()->amo);        // deep copy necessary if T==double_complex
        occ=nemo->get_calc()->aocc;
    } else if (ispin==1) {  // beta spin
        mo_ket=convert<double,T,NDIM>(world,nemo->get_calc()->bmo);        // deep copy necessary if T==double_complex
        occ=nemo->get_calc()->bocc;
    }

    mo_bra=mul(world,nemo->ncf->square(),mo_ket);
    truncate(world,mo_bra);
    poisson = std::shared_ptr<real_convolution_3d>(
            CoulombOperatorPtr(world, nemo->get_calc()->param.lo(),
                    nemo->get_calc()->param.econv()));

}

template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::operator()(
		const std::vector<Function<T,NDIM> >& vket, const double mul_tol) const {
    const bool same = this->same();
    double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb

    reconstruct(world, mo_bra,false);
    reconstruct(world, mo_ket,false);
    world.gop.fence();
    norm_tree(world, mo_bra,false);
    norm_tree(world, mo_ket,false);
    world.gop.fence();

//    if (!same) {
        reconstruct(world, vket);
        norm_tree(world, vket);
//    }

    vecfuncT Kf;
    if (multiworld_) {
        if (efficient_) {
            Kf=K_macrotask_efficient(vket,mul_tol);
        } else {
            Kf=K_macrotask(vket,mul_tol);
        }
    } else if (small_memory_) {
    	Kf=K_small_memory(vket,mul_tol);     // Smaller memory algorithm ... possible 2x saving using i-j sym
    } else {
    	Kf=K_large_memory(vket,mul_tol);
    }
	return Kf;

}

template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::K_macrotask(const vecfuncT& vket, const double mul_tol) const {

	if (world.rank()==0) print("\nentering macrotask version");

    int nf = vket.size();
	const long nocc=mo_ket.size();

	MacroTaskQ taskq(world,world.size());

	double cpu0=cpu_time();
	for (int i=0; i<nocc; ++i) {
		taskq.cloud.store(world,mo_bra[i],i);
		taskq.cloud.store(world,mo_ket[i],i+nocc);
	}

//	taskq.cloud.store(world,mo_bra,0);
//	taskq.cloud.store(world,mo_ket,1);
	taskq.cloud.print_timings(world);
	double cpu1=cpu_time();
	print("cpu time for storing ",cpu1-cpu0);

	double econv=FunctionDefaults<3>::get_thresh();

	long batchsize=std::max(1l,nf/(ntask_per_subworld*world.size()));
	if (world.rank()==0) print("batchsize",batchsize);
	MacroTaskBase::taskqT vtask;
	for (int i=0; i<nf; i+=batchsize) {
		long inputrecord=nocc+nocc+i;
		long outputrecord=nocc+nocc+nf+i;
		auto it1=vket.begin()+i;
		auto it2=std::min(it1+batchsize,vket.end());

		taskq.cloud.store(world,vecfuncT(it1,it2),inputrecord);

		MacroTaskExchange task(inputrecord,outputrecord,lo,nocc, econv, mul_tol);
		vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskExchange(task)));
	}
	world.gop.fence();

	auto current_pmap=FunctionDefaults<3>::get_pmap();
	FunctionDefaults<3>::set_default_pmap(taskq.get_subworld());
	taskq.run_all(vtask);
	FunctionDefaults<3>::set_pmap(current_pmap);
	vecfuncT Kf;
	double cpu2=cpu_time();
	for (int i=0; i<nf; i+=batchsize) {
		vecfuncT dummy;
		taskq.cloud.load<vecfuncT>(world,dummy,nocc+nocc+nf+i);		// spending nearly all i/o time in cloud.h:197
		Kf=append(Kf,dummy);
	}
	double cpu3=cpu_time();
	if (world.rank()==0) printf("loading wall time in collection step %4.1fs\n",cpu3-cpu2);

	taskq.cloud.print_timings(world);
	taskq.get_subworld().gop.fence();
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

        if (world.rank()==0) print("\nentering macrotask_efficient version");

        const long nocc=mo_ket.size();
        const long nsubworld=world.size();

        MacroTaskQ taskq(world,nsubworld);

        double wall0=wall_time();
        for (int i=0; i<nocc; ++i) {
            taskq.cloud.store(world,mo_bra[i],i);
            taskq.cloud.store(world,mo_ket[i],i+nocc);
            taskq.cloud.store(world,vket[i],i+2*nocc);
        }
        double wall1=wall_time();
        if (world.rank()==0) printf("wall time for storing %4.1fs\n",wall1-wall0);

        // partition the exchange matrix into tiles, the tile structure must be symmetric
        std::vector<std::pair<long,long> > ranges=OrbitalPartitioner::partition_for_exchange(5,nsubworld,nocc);
        if (world.rank()==0) {
            print("spliiting nocc into",ranges.size(),"ranges");
            for (auto& r : ranges) printf(" %2ld to %2ld",r.first,r.second);
            print("");
        }


        double econv=FunctionDefaults<3>::get_thresh();

        MacroTaskBase::taskqT vtask;
        for (auto& row_range : ranges) {
            for (auto& column_range : ranges) {
                // compute only the the upper triangular matrix
                if (row_range.first <= column_range.first) {
                    MacroTaskExchangeEfficient task(row_range, column_range, nocc, lo, econv, mul_tol);
                    vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskExchangeEfficient(task)));
                }
            }
        }
        world.gop.fence();
        // sort the tasks according to their tile size, assuming that's a proxy for how much time they take
        std::sort(vtask.begin(),vtask.end(), [](auto a, auto b) {return a->get_priority() > b->get_priority(); });


        auto current_pmap=FunctionDefaults<3>::get_pmap();
        FunctionDefaults<3>::set_default_pmap(taskq.get_subworld());
        taskq.run_all(vtask);
        FunctionDefaults<3>::set_pmap(current_pmap);

        // collect batches of result data
        vecfuncT Kf;
        double cpu2=cpu_time();
        for (auto row_r : ranges) {
            vecfuncT dummy1=zero_functions_compressed<T,NDIM>(world,row_r.second-row_r.first);
        for (auto column_r : ranges) {
                vecfuncT dummy;
                long record = MacroTaskExchangeEfficient::outputrecord(row_r, column_r);
                taskq.cloud.load<vecfuncT>(world, dummy,record);
                dummy1+=dummy;
            }
            Kf = append(Kf, dummy1);
        }
        double cpu3=cpu_time();
        if (world.rank()==0) printf("loading wall time in collection step %4.1fs\n",cpu3-cpu2);

        taskq.cloud.print_timings(world);
        taskq.get_subworld().gop.fence();
        return Kf;
    }


    template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::K_small_memory(const vecfuncT& vket, const double mul_tol) const {

	const long nocc=mo_ket.size();
	vecfuncT Kf=zero_functions_compressed<T,NDIM>(world,nocc);

	for(int i=0; i<nocc; ++i){
		if(occ[i] > 0.0){
		    double cpu0=cpu_time();
			vecfuncT psif = mul_sparse(world, mo_bra[i], vket, mul_tol); /// was vtol
			truncate(world, psif);
            double cpu1=cpu_time();
            double mul1=cpu1-cpu0;
            cpu0=cpu_time();
			psif = apply(world, *poisson.get(), psif);
			truncate(world, psif);
            cpu1=cpu_time();
            double apply1=cpu1-cpu0;
            cpu0=cpu_time();
			psif = mul_sparse(world, mo_ket[i], psif, mul_tol); /// was vtol
			gaxpy(world, 1.0, Kf, occ[i], psif);
            cpu1=cpu_time();
            double dot=cpu1-cpu0;
            printf("timings for mul1, apply, dot: %8.2fs %8.2fs %8.2fs\n",mul1,apply1,dot);

		}
	}
	return Kf;
}

template<typename T, std::size_t NDIM>
std::vector<Function<T,NDIM> > Exchange<T,NDIM>::K_large_memory(const vecfuncT& vket, const double mul_tol) const {    // Larger memory algorithm ... use i-j sym if psi==f

    return compute_K_tile(world,mo_bra,mo_ket,vket,poisson,this->same(),occ,mul_tol);
}

template<typename T, std::size_t NDIM>
std::vector<Function<T, NDIM> >
Exchange<T, NDIM>::compute_K_tile(World &world, const vecfuncT &mo_bra, const vecfuncT &mo_ket,
                                  const vecfuncT &vket, std::shared_ptr<real_convolution_3d> poisson,
                                  const bool same, const Tensor<double> &occ, const double mul_tol) {

    const long nf = vket.size();
    const long nocc = mo_ket.size();
    vecfuncT Kf = zero_functions_compressed<T, NDIM>(world, nocc);
    double tol = FunctionDefaults<3>::get_thresh(); /// Important this is consistent with Coulomb

    vecfuncT psif;
    for (int i = 0; i < nocc; ++i) {
        int jtop = nf;
        if (same)
            jtop = i + 1;
        for (int j = 0; j < jtop; ++j) {
            psif.push_back(mul_sparse(mo_bra[i], vket[j], mul_tol, false));
        }
    }

    world.gop.fence();
    truncate(world, psif, tol);
    psif = apply(world, *poisson.get(), psif);
    truncate(world, psif, tol);
    reconstruct(world, psif);
    norm_tree(world, psif);
    vecfuncT psipsif = zero_functions<T, NDIM>(world, nf * nocc);
    int ij = 0;
    for (int i = 0; i < nocc; ++i) {
        int jtop = nf;
        if (same)
            jtop = i + 1;
        for (int j = 0; j < jtop; ++j, ++ij) {
            psipsif[i * nf + j] = mul_sparse(psif[ij], mo_ket[i], mul_tol, false);
            if (same && i != j) {
                psipsif[j * nf + i] = mul_sparse(psif[ij], mo_ket[j], mul_tol, false);
            }
        }
    }

    world.gop.fence();
    print(occ);
    psif.clear();
    world.gop.fence();
    compress(world, psipsif);
    for (int i = 0; i < nocc; ++i) {
        for (int j = 0; j < nf; ++j) {
            Kf[j].gaxpy(1.0, psipsif[i * nf + j], occ[i], false);
        }
    }
    world.gop.fence();
    psipsif.clear();
    world.gop.fence();
    truncate(world, Kf, tol);
    return Kf;

};

template class Exchange<double_complex,3>;
template class Exchange<double,3>;

template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);


} /* namespace madness */
