//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <chem/exchangeoperator.h>

#include <chem/SCF.h>
#include <chem/nemo.h>
#include <chem/correlationfactor.h>

using namespace madness;

namespace madness {


template<typename T, std::size_t NDIM>
Exchange<T,NDIM>::Exchange(World& world, const SCF* calc, const int ispin)
        : world(world), small_memory_(true), same_(false) {
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
    int nocc = mo_bra.size();
    int nf = vket.size();
    double tol = FunctionDefaults < 3 > ::get_thresh(); /// Important this is consistent with Coulomb
    vecfuncT Kf = zero_functions_compressed<T,NDIM>(world, nf);
    reconstruct(world, mo_bra);
    norm_tree(world, mo_bra);
    reconstruct(world, mo_ket);
    norm_tree(world, mo_ket);
    if (!same) {
        reconstruct(world, vket);
        norm_tree(world, vket);
    }

    // macro task version
    if (multiworld_) {

    	if (world.rank()==0) print("entering macrotask version");

    	MacroTaskQ taskq(world,world.size());

    	double cpu0=cpu_time();
    	const long nocc=mo_ket.size();
		for (int i=0; i<nocc; ++i) {
			taskq.cloud.store(world,mo_bra[i],i);
			taskq.cloud.store(world,mo_ket[i],i+nocc);
		}
		double cpu1=cpu_time();
		print("cpu time for storing ",cpu1-cpu0);

		double lo=1.e-4;
		double econv=FunctionDefaults<3>::get_thresh();

		long batchsize=std::max(1l,nf/(ntask_per_subworld*world.size()));
		if (world.rank()==0) print("batchsize",batchsize);
		MacroTaskBase::taskqT vtask;
		for (int i=0; i<nf; i+=batchsize) {
			long inputrecord=nocc+nocc+i;
			long outputrecord=nocc+nocc+nf+i;
			auto it1=vket.begin()+i;
			auto it2=std::min(it1+batchsize,vket.end());
//			if (world.rank()==0) {
//				print("storing vket: batch(begin,end), inputrecord",it1-vket.begin(),it2-vket.begin(),inputrecord);
//			}
			taskq.cloud.store(world,vecfuncT(it1,it2),inputrecord);

			MacroTaskExchange task(inputrecord,outputrecord,nocc,
							lo, econv, mul_tol);
			vtask.push_back(std::shared_ptr<MacroTaskBase>(new MacroTaskExchange(task)));
		}
		world.gop.fence();

		auto current_pmap=FunctionDefaults<3>::get_pmap();
	    FunctionDefaults<3>::set_default_pmap(taskq.get_subworld());
	    taskq.run_all(vtask);
		FunctionDefaults<3>::set_pmap(current_pmap);
		Kf.clear();
		for (int i=0; i<nf; i+=batchsize) {
//			print("loading results from outputrecord",nocc+nocc+nf+i);
//			taskq.cloud.load<functionT>(world,Kf[i],nocc+nocc+nf+i);
			vecfuncT dummy;
			taskq.cloud.load<vecfuncT>(world,dummy,nocc+nocc+nf+i);
			Kf=append(Kf,dummy);
		}
		taskq.cloud.print_timings(world);


    } else if (small_memory_) {     // Smaller memory algorithm ... possible 2x saving using i-j sym
    	int ii=0;
    	for(int i=0; i<nocc; ++i){
            if(occ[i] > 0.0){
                vecfuncT psif = mul_sparse(world, mo_bra[i], vket, mul_tol); /// was vtol
                truncate(world, psif);
                psif = apply(world, *poisson.get(), psif);
                truncate(world, psif);
                psif = mul_sparse(world, mo_ket[i], psif, mul_tol); /// was vtol
                gaxpy(world, 1.0, Kf, occ[i], psif);
            }
            ii++;
        }
    	print("N (small mem) ii",ii);
    } else {    // Larger memory algorithm ... use i-j sym if psi==f
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
        truncate(world, psif,tol);
        psif = apply(world, *poisson.get(), psif);
        truncate(world, psif, tol);
        reconstruct(world, psif);
        norm_tree(world, psif);
        vecfuncT psipsif = zero_functions<T,NDIM>(world, nf * nocc);
        int ij = 0;
        for (int i = 0; i < nocc; ++i) {
            int jtop = nf;
            if (same)
                jtop = i + 1;
            for (int j = 0; j < jtop; ++j, ++ij) {
                psipsif[i * nf + j] = mul_sparse(psif[ij], mo_ket[i],mul_tol ,false);
                if (same && i != j) {
                    psipsif[j * nf + i] = mul_sparse(psif[ij], mo_ket[j],mul_tol , false);
                }
            }
        }
        world.gop.fence();
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
    }
    truncate(world, Kf, tol);
//    print_size(world,Kf,"Kf");
    return Kf;

}


template class Exchange<double_complex,3>;
template class Exchange<double,3>;

template <> volatile std::list<detail::PendingMsg> WorldObject<MacroTaskQ>::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<MacroTaskQ>::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, Exchange<double,3>::MacroTaskExchange, madness::Hash<long> > >::pending_mutex(0);

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);


} /* namespace madness */
