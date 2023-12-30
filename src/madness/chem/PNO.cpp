/*
 * PNO.cpp
 *
 *  Created on: Oct 22, 2018
 *      Author: kottmanj
 */

#include<madness/chem/PNO.h>

namespace madness {

madness::PairEnergies PNO::compute_projected_mp2_energies(PNOPairs& pairs) const {
	double total = 0.0;
	PAIRLOOP(it)
	{
		double es = 0.0;
		double et = 0.0;
		if(pairs.pno_ij[it.ij()].size()>0){
			const auto& t=pairs.t_ij;
			// make gijab tensors
			const vector_real_function_3d& pno = pairs.pno_ij[it.ij()];
			const real_function_3d& moi = f12.acmos[it.i()];
			const real_function_3d& moj = f12.acmos[it.j()];
			const Tensor<double> tijab = t[it.ij()];
			const Tensor<double> ts = tijab + transpose(tijab);
			const Tensor<double> tt = tijab - transpose(tijab);
			Tensor<double> gijab = matrix_inner(world, moi * pno, apply(world, *poisson, moj * pno));

			double sfactor = 1.0;
			if (it.diagonal())
				sfactor = 0.5;
			const double tfactor = 3.0;

			for (size_t a = 0; a < pno.size(); ++a) {
				for (size_t b = 0; b < pno.size(); ++b) {
					es += sfactor * (ts(a, b) * gijab(a, b));
					et += tfactor * (tt(a, b) * gijab(a, b));
				}
			}
		}
		pairs.energies.eijs[it.ij()] = es;
		pairs.energies.eijt[it.ij()] = et;
		total += (es + et);

	}
	pairs.energies.energy = total;
	pairs.energies.update();
	return pairs.energies;
}

void PNO::solve(std::vector<PNOPairs>& all_pairs) const {
	// set target thresh for correlation computation
	// in general it can be (much?) lower than that of the HF
	// the importance of each PNO is roughly proportional to its norm times the square root of the density
	// (hard to be more precise ... the energy contribution of PNO ~ Lowdin's formula for the increments
	//  in He (n+1/2)^{-6} ... will change when F12 added)
	FunctionDefaults < 3 > ::set_thresh(param.thresh());
	param.print("PNO Parameters\npno","end");
	if (param.f12()){
		f12.param.print("F12 Parameters\nf12", "end");
	}
	if (param.debug()) {
		print("ElectronPairIterator will iterate the pairs: pair_ij, i , j , ij");
		for (ElectronPairIterator it = pit(); it; ++it) {
		  //const size_t ij = it.ij();
		  //const size_t i = it.i();
		  //const size_t j = it.j();
			if (world.rank() == 0)
				std::cout << it.name() << ", " << it.i() << ", " << it.j() << ", " << it.ij() << "\n";
		}
		print("OrbitalIterator iterates");
		for (OrbitalIterator kit = oit(); kit; ++kit) {
			if (world.rank() == 0)
				std::cout << kit.i() << "\n";
		}
	}

	if (param.cispd().empty()){
		if(all_pairs.empty()){
			PNOPairs mp2(MP2_PAIRTYPE,f12.acmos.size());
			solve_mp2(mp2);
			all_pairs.push_back(mp2);
		}
		else{
			solve_mp2(all_pairs);
		}
	}else solve_cispd(all_pairs);
}

vector_real_function_3d PNO::guess_virtuals(const vector_real_function_3d& f, const GuessType& inpgt) const {
	std::string test = std::to_string(inpgt);
	GuessType guesstype = inpgt;
	if (guesstype == UNKNOWN_GUESSTYPE)
		guesstype = param.guesstype();

	if (guesstype == FROM_FILE_GUESSTYPE)
		return (basis.guess_contracted_virtuals_from_file());
	else if (guesstype == PARTIAL_WAVE_GUESSTYPE)
		return (basis.guess_virtuals_internal(param.partial_wave()));
	else if (guesstype == PREDEFINED_GUESSTYPE)
		return basis.predefined_guess(param.predefined_guess());
	else if (guesstype == SCF_GUESSTYPE)
		return Q(nemo.get_calc()->ao);
	else if (guesstype == EXOP_GUESSTYPE)
		return Q(basis.guess_with_exop(f, param.exop(),param.exop_trigo()));
	else if (guesstype == EMPTY_GUESSTYPE)
		return vector_real_function_3d();
	else if (guesstype == PSI4_GUESSTYPE){
		return Q(basis.guess_with_psi4(nemo.get_calc()->amo));
	}
	else
		MADNESS_EXCEPTION(("unknown guesstype" + std::to_string(guesstype)).c_str(), 1);

	MADNESS_EXCEPTION(("unknown guesstype" + std::to_string(guesstype)).c_str(), 1);
	return vector_real_function_3d();
}

PairEnergies PNO::compute_cispd_correction_es(const vector_real_function_3d& xcis, PNOPairs& pairs) const {
	TIMER(timer);
	MADNESS_ASSERT(pairs.type==CISPD_PAIRTYPE);
	double s2b = 0.0;
	double s2c = 0.0;
	for (ElectronPairIterator it = pit(); it; ++it) {
		const Tensor<double> t = (2.0 * pairs.t_ij[it.ij()] - transpose(pairs.t_ij[it.ij()]));
		const real_function_3d moi = f12.acmos[it.i()];
		const real_function_3d moj = f12.acmos[it.j()];
		const real_function_3d xi = xcis[it.i()];
		const real_function_3d xj = xcis[it.j()];
		const auto& pno=pairs.pno_ij[it.ij()];
		double factor = 1.0;
		if (it.diagonal())
			factor = 0.5;

		double s2b_ij = 0.0;
		double s2c_ij = 0.0;
		if(pno.empty()){
			msg << pairs.name(it) << " is empty: no contribution\n";
		}else{
			{
				Tensor<double> gijab;
				if (it.diagonal())
					gijab = 2.0 * matrix_inner(world, apply(world, *poisson, xi * pno), moj * pno); // <xij|g|ab> + <ixj|g|ab> tensor
				else
					gijab = matrix_inner(world, apply(world, *poisson, xi * pno), moj * pno) + matrix_inner(world, moi * pno, apply(world, *poisson, xj * pno));

				for (size_t a = 0; a < pno.size(); ++a) {
					for (size_t b = 0; b < pno.size(); ++b) {
						s2b_ij += (t(a, b) * gijab(a, b));
					}
				}
			}

			{
				// transform pnos according to transposed(Ox)|a>, Ox=|xk><k| -> transposed(Ox)=|k><x_k|
				Projector<double, 3> Oxt(xcis, f12.acmos);
				const vector_real_function_3d tpno = Oxt(pno);
				Tensor<double> gijab;
				if (it.diagonal())
					gijab = 2.0 * matrix_inner(world, apply(world, *poisson, moi * tpno), moj * pno);
				else
					gijab = matrix_inner(world, apply(world, *poisson, moi * tpno), moj * pno) + matrix_inner(world, moi * pno, apply(world, *poisson, moj * tpno));

				for (size_t a = 0; a < pno.size(); ++a) {
					for (size_t b = 0; b < pno.size(); ++b) {
						s2c_ij -= (t(a, b) * gijab(a, b));
					}
				}
			}

		}
		s2b += factor * (s2b_ij);
		s2c += factor * (s2c_ij);
		pairs.energies.eijt[it.ij()]=factor*(s2b_ij+s2c_ij);
	}
	msg << "---------------------\n";
	msg << "CIS(D) ES Correction:\n";
	msg << "S2b=" << s2b << "\n";
	msg << "S2c=" << s2c << "\n";
	msg << "tot=" << s2b+s2c << "\n";
	msg << "---------------------\n";
	pairs.energies.update();
	timer.stop().print("CIS(D) Correction ES");
	return pairs.energies;
}

PairEnergies PNO::compute_cispd_correction_gs(const vector_real_function_3d& xcis,const PNOPairs& pairs) const {
	MADNESS_ASSERT(pairs.type==MP2_PAIRTYPE);
	// compute GS part of CIS(D) energy
	MyTimer timer = MyTimer(world).start();
	double s4a = 0.0;
	double s4b = 0.0;
	double s4c = 0.0;
	PairEnergies energies(pairs.npairs);
	for (ElectronPairIterator it = pit(); it; ++it) {
		const size_t i = it.i();
		const size_t j = it.j();
		const size_t ij = it.ij();
		const real_function_3d& moi = f12.acmos[i];
		const real_function_3d& moj = f12.acmos[j];
		const vector_real_function_3d& acmos = f12.acmos;
		const vector_real_function_3d& pno = pairs.pno_ij[ij];
		const Tensor<double>& t2=pairs.t_ij[ij];
		const real_function_3d& xi = xcis[i];
		const real_function_3d& xj = xcis[j];

		double s4b_ij = 0.0;
		double s4a_ij = 0.0;
		double s4c_ij = 0.0;
		if(pno.empty()){
			msg << pairs.name(it) << " empty: no contribution\n";
		}else{
			// j*<i|g|xk> intermediate
			const vector_real_function_3d igx_j = apply(world, *poisson, moi * xcis) * moj;
			// Integrals <i,j|g|xk,b> : kxb matrix
			const Tensor<double> ijgxb = matrix_inner(world, igx_j, pno);
			// intermediates needed for all parts
			vector_real_function_3d iga = apply(world, (*poisson), (moi * pno));
			vector_real_function_3d jgb = apply(world, (*poisson), (moj * pno)); // simplify later
			const Tensor<double> Sia = inner(world, xi, pno);
			const Tensor<double> Sjb = inner(world, xj, pno); // simplify later
			// 2tijab-tijba part
			const Tensor<double> t = (2.0 * t2 - transpose(t2));

			{
				// transformed a
				vector_real_function_3d tpno;
				for (size_t a = 0; a < pno.size(); ++a) {
					const Tensor<double> Ska = inner(world, pno[a], xcis);
					real_function_3d ta = real_factory_3d(world);
					for (size_t k = 0; k < size_t(Ska.size()); ++k)
						ta += Ska(k) * xcis[k];
					tpno.push_back(ta);
				}
				MADNESS_ASSERT(tpno.size() == pno.size());
				Tensor<double> gijab;
				{
					const vector_real_function_3d ita = moi * tpno;
					const Tensor<double> tmp1 = matrix_inner(world, ita, jgb);
					const vector_real_function_3d jtb = moj * tpno;
					const Tensor<double> tmp2 = matrix_inner(world, iga, jtb);
					gijab = (tmp1 + tmp2);
				}
				for (size_t a = 0; a < pno.size(); ++a) {
					for (size_t b = 0; b < pno.size(); ++b) {
						s4a_ij -= (t(a, b) * gijab(a, b));
					}
				}
			}
			//				std::cout << it.name()<< "s4a_ij=" << s4a_ij << "\ndb1=" << db1 << "\ndb2=" << db2 <<"\n";

			{
				const Tensor<double> Sik = inner(world, xi, xcis);
				real_function_3d ti = real_factory_3d(world);
				for (size_t k = 0; k < size_t(Sik.size()); ++k)
					ti += Sik(k) * acmos[k];
				real_function_3d tj = real_factory_3d(world);
				if (it.diagonal())
					tj = ti;
				else {
					const Tensor<double> Sjk = inner(world, xj, xcis);
					for (size_t k = 0; k < size_t(Sjk.size()); ++k)
						tj += Sjk(k) * acmos[k];
				}
				// transformed gijab with i=Sik_k or/and j=Sjk_k
				Tensor<double> gijab;
				{
					vector_real_function_3d tia = ti * pno;
					const Tensor<double> tmp1 = matrix_inner(world, tia, jgb);
					vector_real_function_3d tjb = tj * pno;
					const Tensor<double> tmp2 = matrix_inner(world, iga, tjb);
					gijab = (tmp1 + tmp2);
				}
				MADNESS_ASSERT(size_t(gijab.size()) == pno.size() * pno.size());
				for (size_t a = 0; a < pno.size(); ++a) {
					for (size_t b = 0; b < pno.size(); ++b) {
						s4b_ij -= (t(a, b) * gijab(a, b));
					}
				}
			}

			{
				Tensor<double> gjb(pno.size()); // sum_k 2<kj|g|xk,b>-<jk|g|xk,b>
				Tensor<double> gia(pno.size()); // sum_k 2<ki|g|xk,a>-<ik|g|xk,a>
				real_function_3d im_kxk = real_factory_3d(world);
				for (size_t k = 0; k < acmos.size(); ++k) {
					im_kxk += acmos[k] * xcis[k];
					gjb -= inner(world, (*poisson)(moj * xcis[k]), acmos[k] * pno);
					gia -= inner(world, (*poisson)(moi * xcis[k]), acmos[k] * pno);
				}
				gjb += 2 * inner(world, im_kxk, jgb);
				gia += 2 * inner(world, im_kxk, iga);
				for (size_t a = 0; a < pno.size(); ++a) {
					for (size_t b = 0; b < pno.size(); ++b) {
						s4c_ij += (t(a, b) * (Sia(a) * gjb(b) + Sjb(b) * gia(a)));
					}
				}
			}

		}
		double factor = 1.0;
		if (it.diagonal()) factor = 0.5;
		s4a += s4a_ij * factor;
		s4b += s4b_ij * factor;
		s4c += s4c_ij * factor;
		energies.eijs[ij]=factor*(s4a_ij+s4b_ij+s4c_ij);
	}
	msg << "---------------------\n";
	msg << "CIS(D) GS Correction:\n";
	msg << "S4a=" << s4a << "\n";
	msg << "S4b=" << s4b << "\n";
	msg << "S4c=" << s4c << "\n";
	msg << "tot=" << s4a+s4b+s4c << "\n";
	msg << "---------------------\n";
	energies.update();
	timer.stop().print("CIS(D) Correction GS");
	return energies;
}

PairEnergies PNO::compute_cispd_f12_correction_gs(const vector_real_function_3d& xcis, PairEnergies& energies) const {
	if (param.f12() == false) return energies;

	MyTimer timer = MyTimer(world).start();
	double s4a_f12 = 0.0;
	double s4b_f12 = 0.0;
	double s4c_f12 = 0.0;
	for (ElectronPairIterator it = pit(); it; ++it) {
		const size_t i = it.i();
		const size_t j = it.j();
		const size_t ij = it.ij();
		const real_function_3d& moi = f12.acmos[i];
		const real_function_3d& moj = f12.acmos[j];
		const vector_real_function_3d& acmos = f12.acmos;
		const real_function_3d& xi = xcis[i];
		const real_function_3d& xj = xcis[j];

		double s4a_f12_ij = 0.0;
		double s4b_f12_ij = 0.0;
		double s4c_f12_ij = 0.0;

		// intermediates
		const Tensor<double> Sik = inner(world, xi, xcis);
		real_function_3d ti = real_factory_3d(world);
		for (size_t k = 0; k < size_t(Sik.size()); ++k)
			ti += Sik(k) * acmos[k];
		real_function_3d tj = real_factory_3d(world);
		if (it.diagonal())
			tj = ti;
		else {
			const Tensor<double> Sjk = inner(world, xj, xcis);
			for (size_t k = 0; k < size_t(Sjk.size()); ++k)
				tj += Sjk(k) * acmos[k];
		}
		real_function_3d im_kxk = real_factory_3d(world);
		for (size_t k = 0; k < acmos.size(); ++k)
			im_kxk += acmos[k] * xcis[k];

		{
			// 2.0<xk,yk|f|ij> - <yk,xk|f|ij> , yk = Q(igxk*j) and similar for ij -> ji
			if (it.diagonal()) {
				const vector_real_function_3d yk = Q(apply(world, *poisson, moi * xcis) * moj);
				s4a_f12_ij = -2.0 * (madness::inner(world, apply(world, *f12.fop, xcis * moi), yk * moj).sum());
			} else {
				{
					const vector_real_function_3d yk = Q(apply(world, *poisson, moi * xcis) * moj);
					s4a_f12_ij = -1.0 * (2.0 * madness::inner(world, apply(world, *f12.fop, xcis * moi), yk * moj).sum() - madness::inner(world, apply(world, *f12.fop, xcis * moj), yk * moi).sum());
				}
				{
					const vector_real_function_3d ykx = Q(apply(world, *poisson, moj * xcis) * moi);
					s4a_f12_ij += -1.0
							* (2.0 * madness::inner(world, apply(world, *f12.fop, xcis * moj), ykx * moi).sum() - madness::inner(world, apply(world, *f12.fop, xcis * moi), ykx * moj).sum());
				}
			}
		}

		{
			// here we have 2.0*(<ti,j|gQf|ij> + <tj,i|gQf|ji>) - (<ti,j|gQf|ji> + <j,ti|gQf|ij>)
			if (it.diagonal()) {
				s4b_f12_ij = -2.0 * (f12.compute_fQg_integral(moi, moj, ti, moj));
			} else
				s4b_f12_ij = -1.0
				* (2.0 * (f12.compute_fQg_integral(moi, moj, ti, moj) + f12.compute_fQg_integral(moj, moi, tj, moi))
						- (f12.compute_fQg_integral(moj, moi, ti, moj) + f12.compute_fQg_integral(moi, moj, tj, moi)));
		}

		{
			// similar than s4a
			double part1 = 0.0;
			{
				const real_function_3d yj = Q((*poisson)(im_kxk) * moj);
				if (it.diagonal())
					part1 = 2.0 * madness::inner((*f12.fop)(xi * moi), yj * moj);
				else {
					const real_function_3d yi = Q((*poisson)(im_kxk) * moi);
					part1 = 2.0 * (madness::inner((*f12.fop)(xi * moi), yj * moj) + madness::inner((*f12.fop)(xj * moj), yi * moi))
																									- (madness::inner((*f12.fop)(xi * moj), yj * moi) + madness::inner((*f12.fop)(xj * moi), yi * moj));
				}
			}
			double part2 = 0.0;
			{
				real_function_3d jgxk_k = real_factory_3d(world);
				for (size_t k = 0; k < acmos.size(); ++k)
					jgxk_k += (*poisson)(moj * xcis[k]) * acmos[k];
				const real_function_3d yj = Q(jgxk_k);
				if (it.diagonal())
					part2 = 2.0 * madness::inner((*f12.fop)(xi * moi), yj * moj);
				else {
					real_function_3d igxk_k = real_factory_3d(world);
					for (size_t k = 0; k < acmos.size(); ++k)
						igxk_k += (*poisson)(moi * xcis[k]) * acmos[k];
					const real_function_3d yi = Q(igxk_k);
					part2 = 2.0 * (madness::inner((*f12.fop)(xi * moi), yj * moj) + madness::inner((*f12.fop)(xj * moj), yi * moi))
																									- (madness::inner((*f12.fop)(xi * moj), yj * moi) + madness::inner((*f12.fop)(xj * moi), yi * moj));
				}
			}
			//std::cout << "part1 = " << part1 << "\npart2=" << part2 << "\n";
			s4c_f12_ij = 2.0 * part1 - part2;
		}
		double factor = 1.0;
		if (it.diagonal())
			factor = 0.5;

		s4a_f12 += s4a_f12_ij * factor;
		s4b_f12 += s4b_f12_ij * factor;
		s4c_f12 += s4c_f12_ij * factor;
		energies.eijs_f12[ij]=factor*(s4a_f12_ij+s4b_f12_ij+s4c_f12_ij);
	}

	msg << "-------------------------\n";
	msg << "CIS(D) GS Correction F12:\n";
	msg << "S4a_f12=" << s4a_f12 << "\n";
	msg << "S4b_f12=" << s4b_f12 << "\n";
	msg << "S4c_f12=" << s4c_f12 << "\n";
	msg << "tot_f12=" << s4a_f12+s4b_f12+s4c_f12 << "\n";
	msg << "-------------------------\n";
	timer.stop().print("CIS(D) F12-Correction GS");

	return energies;
}

PNOPairs PNO::solve_mp2(PNOPairs& mp2) const {
	msg.section("Solve PNO-MP2");
	const size_t nocc=pit().nocc();
	//const size_t npairs=pit().npairs();
	if(mp2.npairs==0) mp2 = PNOPairs(MP2_PAIRTYPE,nocc);
	// compute f12 energies before
	if(param.f12() and f12.param.energytype()==PROJECTED_ENERGYTYPE){
		PairEnergies f12_energies=f12.compute_projected_f12_energies();
		const double old=mp2.energies.energy_f12;
		mp2.energies.eijs_f12=f12_energies.eijs_f12;
		mp2.energies.eijt_f12=f12_energies.eijt_f12;
		mp2.energies.update();
		f12.print_pair_energies(f12_energies.eijs_f12,f12_energies.eijt_f12,"MP2 F12 Energies",mp2.type);
		msg << "difference to old f12 energy = " << mp2.energies.energy_f12-old << "\n";
	}
	mp2 = iterate_pairs(mp2);
	const std::valarray<vector_real_function_3d>& pno_ij = mp2.pno_ij;
	if (param.f12() and f12.param.energytype() == HYLLERAAS_ENERGYTYPE) {
	        //const PairEnergies& energies_from_solve = mp2.energies;
		//const std::valarray<Tensor<double> >& t2_ij = mp2.t_ij;
		// compute energies (or just print them)
		PairEnergies pair_energies = mp2.energies;
		MyTimer timef12 = MyTimer(world).start();
		PairEnergies ef12(f12.npairs());
		ef12 = f12.compute_f12_energies(pno_ij);
		f12.print_pair_energies(pair_energies.eijs, pair_energies.eijt, "Regularized Energies");
		f12.print_f12_energies(ef12, "F12 Energies");
		std::valarray<double> est(f12.npairs());
		std::valarray<double> ett(f12.npairs());
		for (ElectronPairIterator it = pit(); it; ++it) {
			est[it.ij()] = pair_energies.eijs[it.ij()] + ef12.eijs_f12[it.ij()];
			ett[it.ij()] = pair_energies.eijt[it.ij()] + ef12.eijt_f12[it.ij()];
		}
		f12.print_pair_energies(est, ett, "Total Energies");
		timef12.stop().print("Compute F12 Energies (total time)");
		if (world.rank() == 0) {
			std::cout << "\n\n";
			std::cout << "------------------------------------\n";
			std::cout << "MP2 Energy=" << pair_energies.energy << "\n";
			std::cout << "F12 Energy=" << ef12.energy_f12 << "\n";
			std::cout << "Sum Energy=" << pair_energies.energy + ef12.energy_f12 << "\n";
			std::cout << "------------------------------------\n";
			std::cout << "\n\n";
		}
		pair_energies.eijs_f12 = ef12.eijs_f12;
		pair_energies.eijt_f12 = ef12.eijt_f12;
		pair_energies.energy_f12 = ef12.energy_f12;
		mp2.energies=pair_energies;
	} else {
		f12.print_pair_energies(mp2, "Final MP2 Energies");
	}
	// rank information
	std::pair<size_t, size_t> av_rank = get_average_rank(pno_ij);
	if (world.rank() == 0)
		std::cout << "Average | Max Rank: " << av_rank.first << " | " << av_rank.second << "\n";

	mp2.energies.update();
	return mp2;
}

std::vector<PNOPairs> PNO::solve_cispd(std::vector<PNOPairs>& result) const {
	if(result.empty()){
		PNOPairs mp2(MP2_PAIRTYPE,f12.acmos.size());
		result.push_back(mp2);
	}else MADNESS_ASSERT(result.size()>1);
	MyTimer timer_solve_cispd = MyTimer(world).start();
	const size_t nocc=pit().nocc();
	//const size_t npairs=pit().npairs();
	msg.section("Solve PNO-CIS(D)");

	// 1. get CIS vectors
	MyTimer timer_read_cis = MyTimer(world).start();
	std::vector<CISData> cis_result;
	for (const auto& cispd : param.cispd()) {
		// contains excitation number and cis excitation energy
		const size_t ex = cispd.first;
		const double omega = cispd.second;
		if (world.rank() == 0)
			std::cout << "Looking for CIS vectors with omega=" << omega << ", (ex=" << ex << ")\n";

		vector_real_function_3d xvec;
		for (size_t i = 0; i < nocc; ++i) {
			const size_t ii = i + param.freeze(); // this is how the CIS module stores them (easier for me so I dont have to rename)
			const std::string name = std::to_string(ex) + "_" + "x" + std::to_string(ii);
			real_function_3d tmp = real_factory_3d(world);
			load(tmp, name);
			xvec.push_back(tmp);
		}
		MADNESS_ASSERT(xvec.size() == nocc);
		CISData tmpcis(ex,omega,xvec);
		cis_result.push_back(tmpcis);
	}
	timer_read_cis.stop().print("read cis vectors");
	// 2. solve MP2
	MyTimer time_mp2 = MyTimer(world).start();
	if(result.empty()){
		PNOPairs mp2(MP2_PAIRTYPE,f12.acmos.size());
		solve_mp2(mp2);
		result.push_back(mp2);
	}else{
		solve_mp2(result[0]);
	}
	PNOPairs& mp2=result[0];
	double mp2_energy = mp2.energies.total_energy();

	//const std::valarray<Tensor<double> >& t2_mp2_ij = mp2.t_ij;
	const std::valarray<vector_real_function_3d>& pno_mp2_ij = mp2.pno_ij;
	time_mp2.stop().print("MRA-PNO-MP2 Calculation");
	// 4. solve CIS(D) doubles
	size_t count = 1;
	std::vector<std::pair<double, double> > exenergies; // store CIS and CIS(D) energies
	std::vector<std::pair<size_t, size_t> > av_ranks_cispd; // store average and max ranks

	for (const auto& x : cis_result) {
	        //const double omega_cis = x.omega;
		const vector_real_function_3d xcis = x.x;
		//double omega_gs = 0.0;
		PairEnergies cispd_energies(mp2.npairs);
		if (param.no_compute() == MP2_PAIRTYPE || param.no_compute() == ALL_PAIRTYPE) {
			if (world.rank() == 0)
				std::cout << "no compute keyword for MP2: skipping CIS(D)-GS part\n";
		} else {
			if (world.rank() == 0)
				std::cout << "\n\nComputung GS part of CIS(D) energy:\n\n";

			cispd_energies = compute_cispd_correction_gs(xcis, mp2);
			if (param.f12()){
				cispd_energies=compute_cispd_f12_correction_gs(xcis,cispd_energies);
			}
			cispd_energies.update();
			//omega_gs = cispd_energies.total_energy();
		}

		if (param.f12()){
			cispd_energies=compute_cispd_f12_correction_es(xcis,cispd_energies);
		}
		cispd_energies.update();

		MyTimer timer_cispd = MyTimer(world).start();
		PNOPairs cispd(CISPD_PAIRTYPE,nocc,x);
		if(result.size()>count) cispd=result[count];
		else cispd=PNOPairs(CISPD_PAIRTYPE,nocc,x);

		cispd.energies=cispd_energies;
		f12.print_pair_energies(cispd,"CIS(D) with full GS Correction and F12 Parts of ES " + std::to_string(x.number));
		cispd = iterate_pairs(cispd);
		timer_cispd.stop().print("compute CIS(D) doubles");

		//const double cispd_energy_ex = cispd.energies.energy;
		//const std::valarray<Tensor<double> >& t2_cispd_ij = cispd.t_ij;
		const std::valarray<vector_real_function_3d>& pno_cispd_ij = cispd.pno_ij;

		// (re)-compute CIS(D) part of the energy (s2b and s2c)
		cispd.energies=compute_cispd_correction_es(xcis, cispd);

		// print result
		f12.print_pair_energies(cispd,"CIS(D) Correction");

		if (world.rank() == 0) {
			std::cout << "\n\n\n";
			std::cout << "\n================CIS(D) Excitation "+std::to_string(x.number)+"================\n";
			std::cout << "MP2 correlation energy  = " << mp2_energy << "\n";
			std::cout << std::setw(12) << "CIS"  << " | " << std::setw(12) << "CIS(D)" << " | " << std::setw(12) << "delta" << " | " <<  "\n";
			std::cout << std::scientific << std::setprecision(5);
			std::cout << std::setw(12) << cispd.cis.omega  << " | " << std::setw(12) << cispd.cis.omega+cispd.energies.total_energy() << " | " << std::setw(12) << cispd.energies.total_energy() << " | \n";
			std::cout << "\n=================================================\n";
			std::cout << "\n\n\n";
		}


		av_ranks_cispd.push_back(get_average_rank(pno_cispd_ij));
		exenergies.push_back(std::make_pair(cispd.cis.omega,cispd.cis.omega+cispd.energies.total_energy()));
		result.push_back(cispd);
		++count;
	}
	if (world.rank() == 0 and exenergies.size()>1) {
		std::cout << "\n\n\n";
		std::cout << "\n================CIS(D) finished==================\n";
		std::cout << "MP2 correlation energy  = " << mp2_energy << "\n";
		std::cout << std::setw(12) << "CIS"  << " | " << std::setw(12) << "CIS(D)" << " | " << std::setw(12) << "delta" << " | " <<  "\n";
		std::cout << std::scientific << std::setprecision(5);
		for(const auto& o:exenergies)
			std::cout << std::setw(12)<< o.first  << " | " << std::setw(12) << o.second << " | " << std::setw(12) << o.second-o.first  << " | \n";
		std::cout << "\n=================================================\n";
		std::cout << "\n\n\n";
	}
	const std::pair<size_t, size_t> av_rank_mp2 = get_average_rank(pno_mp2_ij);
	if (world.rank() == 0) {
		std::cout << std::setw(12)<< "Average | Max Rank MP2    : " << std::setw(3) << av_rank_mp2.first << " | " << std::setw(3)<< av_rank_mp2.second << "\n";
		for (const auto& rank : av_ranks_cispd)
			std::cout <<             "Average | Max Rank CIS(D) : " << std::setw(3) << rank.first << " | " << std::setw(3) << rank.second << "\n";
	}
	timer_solve_cispd.stop().print("solve_cispd");
	return result;
}

std::pair<size_t, size_t> PNO::get_average_rank(const std::valarray<vector_real_function_3d>& va) const {
	size_t sum = 0;
	size_t vectors = va.size();
	size_t max = 0;
	if (vectors == 0)
		return std::make_pair(0, 0); // failsafe

	for (const auto& v : va) {
		sum += v.size();
		max = std::max(max, v.size());
	}
	return std::make_pair(sum / vectors, max);
}
PNOPairs PNO::truncate_pairs(PNOPairs& pairs) const {
	TIMER(timer);
	const double thresh=param.thresh();
	// truncating with less fences due to shallow copy
	vector_real_function_3d all=pairs.extract(pairs.pno_ij);
	all=append(all,pairs.extract(pairs.Kpno_ij));
	all=append(all,pairs.extract(pairs.W_ij_i));
	//all=append(all,pairs.extract(pairs.W_ij_j)); -> will lead to race condition for diagonal pairs
	PAIRLOOP(it){
		if(not it.diagonal()) all=append(all,pairs.W_ij_j[it.ij()]);
	}

	if(all.empty()) return pairs;
	//const double size_0 = get_size(world,all);
	compress(world,all);
	truncate(world,all,thresh);

	pairs.update_meminfo();
	msg << pairs.meminfo << "\n";

	timer.stop().print("truncate");
	return pairs;
}

/// freeze pairs which do not contribute significantly
/// frozen pairs will not be further optimized
PNOPairs PNO::freeze_insignificant_pairs(PNOPairs& pairs)const{
	pairs.energies.update();
	PAIRLOOP(it){
		const size_t i=it.i()+param.freeze();
		const size_t j=it.j()+param.freeze();
		bool freeze=false;
		const double threshold = std::min(param.thresh()*0.1,param.econv_pairs());
		if(fabs(pairs.energies.eij[it.ij()])<threshold) freeze=true;

		if(freeze){
			msg << "freeze pair " << pairs.name(it) << " (not further optimized: Energy below " << threshold << ")\n";
			pairs.frozen_ij[it.ij()]=true;
		}

		// now check if the pair should be frozen anyway (demanded by parameter)
		// freeze pairs when demanded (but not when not yet initialized and pre-optimized
		// i.e. freeze if demanded and the pair was restarted
		if(not freeze){
			for(const auto orb: param.freeze_pairs_of_orbital()){
			  if(it.i()==size_t(orb) or it.j()==size_t(orb)) freeze=true;
			}
			for(const auto& pair: param.freeze_pairs()){
			  if(it.i()==size_t(pair.first) and it.j()==size_t(pair.second)) freeze=true;
			  if(it.j()==size_t(pair.first) and it.i()==size_t(pair.second)) freeze=true;
			}

			// if a whitelist exists this means that all pairs which are not on that list should be frozen no matter what the other parameters say
			if(param.active_pairs_of_orbital().size()>0){
				freeze=true;
				for(const auto orb: param.active_pairs_of_orbital()){
				  if(i==size_t(orb) or j==size_t(orb)) freeze=false;
				}
			}

			if(freeze){
				msg << "freeze pair " << pairs.name(it) << " demanded by parameters " << "\n";
				pairs.frozen_ij[it.ij()]=true;
			}

		}



	}
	return pairs;
}

PNOPairs PNO::truncate_pair_ranks(PNOPairs& pairs) const {
	// take care this does not truncate intermediates
	std::valarray<vector_real_function_3d>& result = pairs.pno_ij;
	if (pairs.maxranks_ij.max() > 0) {
		PAIRLOOP(it){
			const int maxrank=pairs.maxranks_ij[it.ij()];
			const vector_real_function_3d& tmp = pairs.pno_ij[it.ij()];
			if (tmp.size() < size_t(maxrank) or maxrank<0) continue; // fast cont.
			vector_real_function_3d tr_tmp;
			for (size_t i = 0; i < size_t(maxrank); ++i) tr_tmp.push_back(tmp[tmp.size() - 1 - i]);
			result[it.ij()] = tr_tmp;
			msg << it.name() << " truncated to rank " << maxrank << "\n";
			if(tr_tmp.size()!=tmp.size()) pairs.clear_intermediates(it);
		}
	} else msg << "all maxranks below 0 ... no truncation\n";
	pairs.pno_ij=result;
	return pairs;
}

PNOPairs PNO::initialize_pairs(PNOPairs& pairs, const GuessType& inpgt) const {
	TIMER(timer);
	// convenience
	const PairType& pt=pairs.type;
	const size_t npairs=pit().npairs();
	std::valarray<vector_real_function_3d>& pno_ij=pairs.pno_ij;
	bool do_cispd = (pt==CISPD_PAIRTYPE);

	// determine if the function was speficifally called or if the parameters should be used
	GuessType guesstype = inpgt;
	if (guesstype == UNKNOWN_GUESSTYPE) guesstype = param.guesstype();

	// determine of restart or other options for the speficfic pair type where demanded
	const bool restart = (param.restart() == pt || param.restart() == ALL_PAIRTYPE);
	//const bool no_opt = (param.no_opt() == pt || param.no_opt() == ALL_PAIRTYPE);
	const bool no_guess = (param.no_guess() == pt || param.no_guess() == ALL_PAIRTYPE || guesstype == EMPTY_GUESSTYPE);
	const bool pair_specific_guess = (guesstype == EXOP_GUESSTYPE);

	// 0. look for restart options
	if (restart) {
		msg << "Restart demanded for " << pt << " pairs\n";
		pairs=load_pnos(pairs);
	}

	// 1. compute initial guess for the virtuals
	if (no_guess == false and pair_specific_guess == false) {
		vector_real_function_3d virtuals;
		// use all MOs to seed the guess virtuals
		vector_real_function_3d init_mos = nemo.get_calc()->amo; // only needed for GT_EXOP
		if (not pairs.cis.x.empty())
			init_mos = append(init_mos, pairs.cis.x);
		virtuals = guess_virtuals(init_mos);

		//for(size_t i=0;i<virtuals.size();++i)plot_plane(world,virtuals[i],"guess_virtuals_"+std::to_string(i));
		if (not virtuals.empty()) {
			virtuals = Q(virtuals);
			virtuals = orthonormalize_rrcd(virtuals,0.1*param.thresh());
			TIMER(timer);
			PAIRLOOP(it){
				if(param.diagonal() and not it.diagonal()) continue;
				vector_real_function_3d& pno = pno_ij[it.ij()];
				if (not pno.empty()) {
					msg << it.name() << ": pnos not empty ... project out and assemble\n";
					QProjector<double, 3> Qpno(world, pno, pno);
					pno = append(pno, Qpno(virtuals));
				} else
					pno = append(pno, virtuals);
			}
			timer.stop().print("project restarted pnos out of guess");
		}else msg << "virtuals are empty\n";


	} else if (no_guess == false) {
		// use just the MOs corresponding to the specific pairs to seed the virtual guess
		msg<< "Demanded the pair-specific exop guess\n";
		TIMER(timer);
		PAIRLOOP(it)
		{
			if(param.diagonal() and not it.diagonal()) continue;
			vector_real_function_3d& pno = pno_ij[it.ij()];
			vector_real_function_3d pair_mo;
			pair_mo.push_back(nemo.get_calc()->amo[it.i() + param.freeze()]);
			if (not it.diagonal())
				pair_mo.push_back(nemo.get_calc()->amo[it.j() + param.freeze()]);
			if (do_cispd) {
				const double normi = pairs.cis.x[it.i()].norm2(); // cis vector contains only active orbitals
				const real_function_3d xi = (1.0 / normi) * pairs.cis.x[it.i()];
				pair_mo.push_back(xi);
				if (not it.diagonal()) {
				  //const double normj = pairs.cis.x[it.j()].norm2(); // cis vector contains only active orbitals
					const real_function_3d xj = (1.0 / normi) * pairs.cis.x[it.j()];
					pair_mo.push_back(xj);
				}
			}
			vector_real_function_3d virtij = guess_virtuals(pair_mo, guesstype);
			if (not pno.empty()) {
				QProjector<double, 3> Qpno(world, pno, pno);
				virtij = Qpno(virtij);
			}

			msg << "Created " << virtij.size() << " virtuals for " << it.name() << "\n";
			pno = append(pno, virtij);
		}
		timer.stop().print("Pair specific exop guess");
	} else msg<< "no guess was explicitly demanded for " << pt << "\n\n";

	TIMER(timerQ);
	PAIRLOOP(it)
	{
		if (not pno_ij[it.ij()].empty()) pno_ij[it.ij()] = Q(pno_ij[it.ij()]);
	}
	pairs=orthonormalize_cholesky(pairs);
	timerQ.stop().print("apply Q and orthonormalize PNOs");

	msg << "Created Guess:\n";

	for (ElectronPairIterator it = pit(); it; ++it) {
		msg << it.name() << " with " << pno_ij[it.ij()].size() << " guess functions\n";
	}
	TIMER(trt)

	for (ElectronPairIterator it = pit(); it; ++it) {
		if(param.diagonal() and not it.diagonal()) continue;
		const double size1 = get_size(world, pno_ij[it.ij()]);
		truncate(world, pno_ij[it.ij()], param.thresh());
		const double size2 = get_size(world, pno_ij[it.ij()]);
		if (world.rank() == 0)
			std::cout << "guess " << it.name() << " size: " << size1 << " | " << size2 << " Gbyte (after truncation)\n";
	}
	trt.stop().print("truncate guess");

	pairs.pno_ij=pno_ij;

	// initialize maxranks and frozen pairs
	pairs.frozen_ij=std::valarray<bool>(false,npairs);
	pairs.maxranks_ij=std::valarray<int>(param.maxrank(),npairs);
	// init amplitudes as empty tensors
	pairs.t_ij=std::valarray<Tensor<double> >(npairs);

	// erase off-diagonal pairs and adapt maxranks
	// if diagonal approximation is demanded
	if(param.diagonal()){
		for (ElectronPairIterator it = pit(); it; ++it) {
			if(not it.diagonal()){
				pairs.maxranks_ij[it.ij()]=0;
				pairs.frozen_ij[it.ij()] = true;
				pairs.pno_ij[it.ij()]=vecfuncT();
			}
		}
	}


	timer.stop().print("Initialize Pairs");
	pairs.verify();

	return pairs;
}



PNOPairs PNO::compute_fluctuation_potential(const ElectronPairIterator& it, PNOPairs& pairs) const {
	const size_t i = it.i();
	const size_t j = it.j();
	const auto pno=pairs.pno_ij[it.ij()];
	const real_function_3d& moi = f12.acmos[i];
	const real_function_3d& moj = f12.acmos[j];
	vector_real_function_3d& W_ij_i=pairs.W_ij_i[it.ij()];
	vector_real_function_3d& W_ij_j=pairs.W_ij_j[it.ij()];
	if(pno.empty()) return pairs;

	if (param.f12()) {
		// for f12 reuse the Kpno intermediate for both W potentials
		const vector_real_function_3d& Kpno=pairs.Kpno_ij[it.ij()];
		W_ij_i = compute_Vreg_aj_i(i, j, pno, Kpno);
		if (it.diagonal())
			W_ij_j = W_ij_i;
		else
			W_ij_j = compute_Vreg_aj_i(j, i, pno, Kpno);
	} else {
		W_ij_i = compute_V_aj_i(moi, moj, pno, Q);
		if (it.diagonal())
			W_ij_j = W_ij_i;
		else
			W_ij_j = compute_V_aj_i(moj, moi, pno, Q);
	}
	return pairs;
}

PNOPairs PNO::compute_cispd_fluctuation_potential(const ElectronPairIterator& it, PNOPairs& pairs) const {
  //const size_t ij = it.ij();
	const size_t i = it.i();
	const size_t j = it.j();
	const auto& pno=pairs.pno_ij[it.ij()];
	const auto& cis=pairs.cis.x;
	const auto& Vx=pairs.cis.Vx;
	const real_function_3d& moi = f12.acmos[i];
	const real_function_3d& moj = f12.acmos[j];
	vector_real_function_3d& W_ij_i=pairs.W_ij_i[it.ij()];
	vector_real_function_3d& W_ij_j=pairs.W_ij_j[it.ij()];
	if(pno.empty()) return pairs;

	Projector<double, 3> Ox(f12.acmos, cis); // adjoint Ox projector
	Projector<double, 3> Oxt(cis, f12.acmos); // adjoint Ox projector
	const real_function_3d& Vxi = Vx[i];
	const real_function_3d& Vxj = Vx[j];
	if (param.f12()) {
		const real_function_3d& xi = cis[i];
		const real_function_3d& xj = cis[j];
		const auto& Kpno=pairs.Kpno_ij[it.ij()];

		// avoid double computation in off-diag pairs
		TIMER(ktimer);
		const vector_real_function_3d Ox_pno = Oxt(pno);
		const auto KOx_pno = K(Ox_pno);
		ktimer.stop().print("K(Ox(pno))");

		W_ij_i = compute_Vreg_aj_i(xi, moj, pno, Q, Kpno)
										+ compute_Vreg_aj_i(moi, xj, pno, Q, Kpno)
										- compute_Vreg_aj_i(moi, moj, Ox_pno, Q, KOx_pno) // can not use Kpno here (would need K(Ox_pno)
										- compute_Vreg_aj_i(moi, moj, pno, Ox, Kpno)
										- compute_Vreg_aj_i_fock_residue(Vxi, moj, pno)
										- compute_Vreg_aj_i_fock_residue(moi, Vxj, pno)
										+ compute_Vreg_aj_i_commutator_response(moi, moj, pno, Vx);
		//truncate(world,W_ij_i[ij],param.thresh);
		if (it.diagonal())
			W_ij_j = W_ij_i;
		else {
			W_ij_j = compute_Vreg_aj_i(xj, moi, pno, Q, Kpno)
											+ compute_Vreg_aj_i(moj, xi, pno, Q, Kpno)
											- compute_Vreg_aj_i(moj, moi, Ox_pno, Q, KOx_pno) // can not use Kpno here (would need K(Ox_pno)
											- compute_Vreg_aj_i(moj, moi, pno, Ox, Kpno)
											- compute_Vreg_aj_i_fock_residue(Vxj, moi, pno)
											- compute_Vreg_aj_i_fock_residue(moj, Vxi, pno)
											+ compute_Vreg_aj_i_commutator_response(moj, moi, pno, Vx);
			//truncate(world,W_ij_j[ij],param.thresh);
		}
	} else {
		const vector_real_function_3d Ox_pno = Oxt(pno);
		const real_function_3d& xi = cis[i];
		const real_function_3d& xj = cis[j];
		W_ij_i = compute_V_aj_i(xi, moj, pno, Q) + compute_V_aj_i(moi, xj, pno, Q) - compute_V_aj_i(moi, moj, Ox_pno, Q) - compute_V_aj_i(moi, moj, pno, Ox);
		if (it.diagonal())
			W_ij_j = W_ij_i;
		else {
			W_ij_j = compute_V_aj_i(xj, moi, pno, Q) + compute_V_aj_i(moj, xi, pno, Q) - compute_V_aj_i(moj, moi, Ox_pno, Q) - compute_V_aj_i(moj, moi, pno, Ox);
		}
	}
	return pairs;
}

PNOPairs PNO::iterate_pairs(PNOPairs & pairs) const {
	msg.section("Iterate Pairs");
	if (param.adaptive_solver() == pairs.type or param.adaptive_solver() == ALL_PAIRTYPE){
		if(param.guesstype() != EXOP_GUESSTYPE){
			msg << "demanded adaptive_solver=" << param.adaptive_solver() << " and guesstype=" << param.guesstype() << " change one of them!\n";
			MADNESS_EXCEPTION("adaptive solver only works for guesstype = exop, deactivate adaptive_solver for the demanded guesstype",1)
		}
		return adaptive_solver(pairs);
	}

	msg.section("Standard Solver");
	if(pairs.empty() or pairs.npairs==0){
		msg.output("Initializing Pairs:");
		pairs = initialize_pairs(pairs,param.guesstype());
		// set the maxranks to maxrank (or to pno_ij size if there was a restart)
		PAIRLOOP(it) pairs.maxranks_ij[it.ij()] = param.maxrank();
		print_ranks(pairs);
	}
	else{
		msg.warning("Standard solver expects empty pairs!");
		MADNESS_EXCEPTION("Standard solver expects empty pairs",1);
	}

	pairs = iterate_pairs_internal(pairs, param.maxiter_micro(), param.econv_micro());

	return pairs;
}


PNOPairs PNO::adaptive_solver(PNOPairs& pairs)const{
	msg.section("Adaptive Solver");
	if(pairs.empty() or pairs.npairs==0){
		msg.output("Initializing empty Pairs:");
		pairs = initialize_pairs(pairs,EMPTY_GUESSTYPE);
	}
	else{
		msg.output("Set back intermediates of pairs:");
		PAIRLOOP(it) pairs.frozen_ij[it.ij()]=false;
		PAIRLOOP(it) pairs.clear_intermediates(it);
	}

	double energy = pairs.energies.total_energy();

	// set the maxranks to zero (or to pno_ij size if there was a restart)
	PAIRLOOP(it) pairs.maxranks_ij[it.ij()] = pairs.pno_ij[it.ij()].size();

	// find out how many virtuals are created for diagonal pairs (this is used to get the same growing for all pairs)
	int rank_increase = (madness::guessfactory::make_predefined_exop_strings(param.exop())).size();
	if(param.rank_increase()>0) rank_increase=param.rank_increase();

	// setup the adaptive protocol
	std::vector<std::string> multipoles(param.maxiter_macro(),param.exop());
	if(param.exop()=="multipole"){
		// multipole protocol needs at least 3 macroiterions therefore the +2
		multipoles=std::vector<std::string>(param.maxiter_macro()+2,"octopole");
		multipoles[0]="dipole+";
		multipoles[1]="quadrupole";
		multipoles = std::vector<std::string>(multipoles.begin(), multipoles.begin()+param.maxiter_macro());
	}else if(param.exop()=="octopole"){
		// fix that stupid name in madness TDHF at one point
		std::vector<std::string> multipoles(param.maxiter_macro(),"octopole");
	}

	msg << "Guess protocol for adaptively growing ranks is: " << multipoles << "\n";

	for (size_t i = 0; i < size_t(param.maxiter_macro()); ++i) {

		// fast return if possible
		if(param.exop()=="none" and i>0){
			msg << "exop=none ... no further macroiterations needed\n";
			break;
		}

		if (i>0){
			bool maxrank_reached=true;
			PAIRLOOP(it)
			{
			  const bool frozen = pairs.frozen_ij[it.ij()];
			  if (not frozen and pairs.pno_ij[it.ij()].size() < size_t(param.maxrank())) maxrank_reached = false;
			}
			if (maxrank_reached){
				msg << "Maximum Rank of" << param.maxrank() << " reached for all non-frozen pairs, no need to iterate further";
				break;
			}
		}

		if (world.rank() == 0)
			std::cout << "\n\n--------Adaptive Solver Macroiteration " << i << " from max " << param.maxiter_macro() << " ----------\n\n";

		if(param.exop()!="none"){
			pairs = grow_rank(pairs,multipoles[i]);
		}
		else {
			msg << "no growing of ranks demadned\n";
		}

		PAIRLOOP(it) pairs.maxranks_ij[it.ij()] += rank_increase;
		// never exceed the global maxrank parameter
		PAIRLOOP(it) pairs.maxranks_ij[it.ij()] = std::min(pairs.maxranks_ij[it.ij()], param.maxrank());

		const PairEnergies oldE=pairs.energies;
		pairs = iterate_pairs_internal(pairs, param.maxiter_micro(), param.econv_micro());
		const double deltaE = pairs.energies.total_energy()-energy;
		energy = pairs.energies.total_energy();

		// test convergence of individual pairs
		bool all_pairs_converged=true;
		double maxdE=0.0;
		size_t maxij=0;
		PAIRLOOP(it){
			const double dE=pairs.energies.eij[it.ij()]-oldE.eij[it.ij()];
			if(fabs(dE)>maxdE){
				maxdE=dE;
				maxij=it.ij();
			}
			if(fabs(dE)>param.econv_pairs()) all_pairs_converged=false;
			else{
				msg << pairs.name(it) << " converged!";
				pairs.frozen_ij[it.ij()]=true;
			}
		}
		msg << "\n\n\nAdaptive Sovler:" <<
				"\n Current Energy=" << energy <<
				"\n DeltaE        =" << deltaE <<
				"\n max dE (pair) =" << maxdE << " (" << maxij << ")" <<
				"\n all_pairs_conv=" << all_pairs_converged << "\n\n\n";
		if (fabs(deltaE) < param.econv_macro()) {
			msg<< "\n\ndeltaE below threshold Macroiterations converged!\n\n";
			break;
		}else if(all_pairs_converged){
			msg << "\n\nall pairs converged but deltaE above given threshold -> check your parameters\n\n";
			break;
		}else msg << "\n\nnot converged yet\n\n";
		if (i == size_t(param.maxiter_macro()))
			msg << "Maximum number of Macroiterations reached\n";

		// do not further optimize and grow ranks of pairs which do not really contribute
		pairs=freeze_insignificant_pairs(pairs);
	}
	return pairs;
}

PNOPairs PNO::iterate_pairs_internal(PNOPairs& pairs, const int maxiter, const double econv) const {
	MyTimer timer_solve = MyTimer(world).start();

	pairs.verify();
	// convenience
	const PairType& pt=pairs.type;
	const CISData& cis=pairs.cis;
	const std::valarray<int>& maxrank=pairs.maxranks_ij;
	std::valarray<vector_real_function_3d>& pno_ij=pairs.pno_ij;
	std::valarray<Tensor<double> >& t2_ij=pairs.t_ij; // ij -> amplitudes in PNO basis
	std::valarray<vector_real_function_3d>& W_ij_i=pairs.W_ij_i; // fluctuation potential
	std::valarray<vector_real_function_3d>& W_ij_j=pairs.W_ij_j; // fluctuation potential
	std::valarray<Tensor<double> >& W_ij=pairs.W_ij; // ij -> <ij|W|ab>
	//std::valarray<Tensor<double> >& F_ij=pairs.F_ij; // ij -> Fock matrix in PNO basis
	PNOTensors::Tensor_IJ_IK<double>& S_ij_ik=pairs.S_ij_ik; // f12.pit().tridx(ij,ik) -> <a_ij|b_ik>
	PNOTensors::Tensor_IJ_KJ<double>& S_ij_kj=pairs.S_ij_kj; // f12.pit().tridx(ij,kj) -> <a_ij|b_kj>
	// see if this is the CIS(D) solver

	bool do_cispd = false;
	if (pt == CISPD_PAIRTYPE) do_cispd = true;
	// assign parameters
	bool no_compute = (param.no_compute() == pt || param.no_compute() == ALL_PAIRTYPE);
	bool restart = (param.restart() == pt || param.restart() == ALL_PAIRTYPE);
	bool no_opt = (param.no_opt() == pt || param.no_opt() == ALL_PAIRTYPE);
	bool no_opt_in_first_iteration = param.no_opt_in_first_iteration();
	if (no_opt_in_first_iteration) msg << "no optimization in first iteration!\n";

	double omega = cis.omega;
	if (do_cispd) {
		TIMER(timerX);
		pairs.cis.Vx = compute_CIS_potentials(cis.x);
		pairs.cis.Kx = K(cis.x);
		MADNESS_ASSERT(cis.Vx.size() == cis.x.size() && cis.x.size() == f12.acmos.size());
		if (omega == 0.0)
			MADNESS_EXCEPTION(("could not find omega for cis excitation number " + std::to_string(cis.number)).c_str(), 1);
		msg<< "CIS excitation energy is " << omega << "\n";
		timerX.stop().print("compute CIS intermediates\n");
	}

	if(econv<0){
		no_compute=false;
		restart=true;
		no_opt=true;
		no_opt_in_first_iteration=true;
	}

	// information about the current iteration cycle
	msg << "SolveInternal: Parameters are assigned as\n";
	msg << "no_compute        = " << no_compute << "\n";
	msg << "no_opt            = " << no_opt << "\n";
	msg << "restart           = " << restart << "\n";
	msg << "do_cispd          = " << do_cispd << "\n";
	msg << "econv             ="  << econv <<"\n";
	msg << "maxrank (min,max) = (" << maxrank.min() << "," << maxrank.max() << ")\n";

	if (do_cispd && world.rank() == 0)
		std::cout << "Solving CIS(D) doubles\n";
	else if (world.rank() == 0)
		std::cout << "Solving MP2 doubles\n";

	// occupied orbitals
	const vector_real_function_3d& amo = nemo.get_calc()->amo;
	// compute the occ-occ block of the Fock matrix
	vector_real_function_3d actmo;
	for (size_t i = param.freeze(); i < amo.size(); ++i)
		actmo.push_back(amo[i]);
	Tensor<double> fmat = F(actmo, actmo);

	auto econverged = false;
	auto dconverged = false;
	auto energy = pairs.energies.total_energy();

	//const auto nocc_act = amo.size() - param.freeze();
	//auto npairs = nocc_act * (nocc_act + 1) / 2; // # of {i>=j} pairs

	size_t iter = 0;
	bool converged = false;
	if (no_compute) {
		print("Found no_compute keyword!\n");
	}
	std::valarray<Tensor<double> > rdm_evals_ij(f12.npairs());
	while (((iter < size_t(maxiter))) and (no_compute == false)) {
		TIMER(timer_iter);
		print_ranks(pairs);

		auto last_iter_energy = energy;
		S_ij_ik.reset();
		S_ij_kj.reset();

		// determines if the W_ij matrix should be computed directly
		bool no_opt_in_this_iteration=false;
		if(no_opt) no_opt_in_this_iteration=true;
		if(iter==0 and no_opt_in_first_iteration) no_opt_in_this_iteration=true;
		if(converged) no_opt_in_this_iteration=true;
		if(iter==size_t(maxiter-1)) no_opt_in_this_iteration=true;

		// if the key word is not set everything else does not matter at all
		if(no_opt_in_first_iteration==false) no_opt_in_this_iteration=false;

		if(no_opt_in_this_iteration) msg << "iter "<<std::setfill(' ')<<std::setw(2) << iter <<  " no opt in this iteration\n";


		truncate_pairs(pairs);

		// 1 Compute K intermediate and Fock matrices
		if(no_opt_in_this_iteration){
			msg << "no opt in this iteration: K intermediate is neither computed nor stored\n";
			// do not store K intermediate (not needed and large in first iteration
			TIMER(timer2);
			PAIRLOOP(it){
				const size_t ij=it.ij();
				if(pairs.pno_ij[ij].empty()) pairs.F_ij[ij]=Tensor<double>(std::vector<long>(2,0));
				else if(pairs.frozen_ij[it.ij()]){
					MADNESS_ASSERT(pairs.F_ij[it.ij()].ndim()==2);
					MADNESS_ASSERT(size_t(pairs.F_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
					msg << "frozen pair " + pairs.name(it)<< ", no recomputation of F_ij\n";
				}
				else{
					pairs.F_ij[ij]=F(pno_ij[ij],pno_ij[ij]);
				}

			}
			timer2.stop().print("compute Fock matrices");
		}else{
			TIMER(timer1);
			vector_real_function_3d allP=pairs.extract(pno_ij);
			// apply K to max 500 functions at a time
			int chunk=param.chunk();
			vector_real_function_3d allK;
			auto allPit =allP.begin();
			while(allK.size()<allP.size()){
				MADNESS_ASSERT(allPit!=allP.end());
				const auto dist=(std::distance(allPit,allP.end()));
				if(dist<chunk) chunk=dist;
				const auto tmp=vector_real_function_3d(allPit,allPit+chunk);
				vector_real_function_3d Ktmp=K(tmp);
				allK.insert(allK.end(),Ktmp.begin(),Ktmp.end());
				allPit+=chunk;
				const double allKs=get_size(world,allK);
				msg << "wall time=" << wall_time()  << " chunk=" << chunk << " allK.size()=" << allK.size() << " memory=" << allKs << " GB\n";
			}
			pairs.Kpno_ij=pairs.reassemble(allK);
			timer1.stop().print("compute KPNO");

			TIMER(timerj);
			vector_real_function_3d allJ;
			allPit =allP.begin();
			while(allJ.size()<allP.size()){
				MADNESS_ASSERT(allPit!=allP.end());
				const auto dist=(std::distance(allPit,allP.end()));
				if(dist<chunk) chunk=dist;
				const auto tmp=vector_real_function_3d(allPit,allPit+chunk);
				vector_real_function_3d Jtmp=J(tmp);
				allJ.insert(allJ.end(),Jtmp.begin(),Jtmp.end());
				allPit+=chunk;
				const double allJs=get_size(world,allJ);
				msg << "wall time=" << wall_time()  << " chunk=" << chunk << " allJ.size()=" << allJ.size() << " memory=" << allJs << " GB\n";
			}
			const auto Jpno_ij=pairs.reassemble(allJ);
			timerj.stop().print("compute JPNO");


			TIMER(timer2);
			PAIRLOOP(it){
				const size_t ij=it.ij();
				if(pairs.pno_ij[it.ij()].empty()){
					pairs.F_ij[it.ij()]=Tensor<double>(std::vector<long>(2,0));
				}
				else if(pairs.frozen_ij[it.ij()]){
					if(pairs.F_ij[it.ij()].size()==0) pairs.F_ij[it.ij()]=F(pairs.pno_ij[it.ij()],pairs.pno_ij[it.ij()]);
					else{
					  MADNESS_ASSERT(pairs.F_ij[it.ij()].dim(0)==pairs.F_ij[it.ij()].dim(1));
					  MADNESS_ASSERT(size_t(pairs.F_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
						continue;
					}
				}
				else if(pno_ij[ij].size()>0){
					Tensor<double> tmat=T(pno_ij[ij],pno_ij[ij]);
					Tensor<double> vmat=V(pno_ij[ij],pno_ij[ij]);
					Tensor<double> jmat=matrix_inner(world,pno_ij[ij],Jpno_ij[ij]);
					Tensor<double> kmat=matrix_inner(world,pno_ij[ij],pairs.Kpno_ij[ij]);
					pairs.F_ij[ij]=tmat+vmat+jmat-kmat;
				}
				else pairs.F_ij[ij]=Tensor<double>(std::vector<long>(2,0));
			}
			timer2.stop().print("compute Fock matrices");
		}

		// 2. (optional) canonicalize PNOs

		if (param.canonicalize_pno()) {
			TIMER(timer2);
			canonicalize(pairs);
			timer2.stop().print("Canonicalize");
		}

		truncate_pairs(pairs);


		// 3.a compute fluctuation potential for each pair
		// in the first iteration: No optimization just pno truncation
		// so no explicit fluctuation potential needed
		if (no_opt_in_this_iteration==false) {

			update_fluctuation_potentials(pairs);

			// 3.b compute fluctuation integrals for each pair
			TIMER(timer3b);
			for (ElectronPairIterator it = pit(); it; ++it) {
				const size_t ij = it.ij();
				//const size_t i = it.i();
				//const size_t j = it.j();
				if(pno_ij[ij].empty()){
					W_ij[ij]=Tensor<double>(std::vector<long>(2,0));
					continue;
				}
				// make sure the matrix has been calculated in the first iteration
				if (pairs.frozen_ij[it.ij()]) {
					if(pairs.W_ij[ij].size()==0 and do_cispd) pairs.W_ij[ij]=compute_cispd_fluctuation_matrix(it, pairs);
					else if(pairs.W_ij[ij].size()==0) pairs.W_ij[ij]=compute_fluctuation_matrix(it, pno_ij[it.ij()], pairs.Kpno_ij[it.ij()]);
					continue;
				}
				// evaluate as average of <a|i W_jb> and <j W_ia| b>
				W_ij[ij] = 0.5 * (matrix_inner(world, pno_ij[ij], W_ij_i[ij]) + matrix_inner(world, W_ij_j[ij], pno_ij[ij]));

			}
			timer3b.stop().print("Fluctuation Matrix");

		} else {

			// for cases without optimization it is cheaper to calculate the fluctuation matrix directly (since then no potentials are needed)
			// in the first iteration this is done for various reasons like
			// 1. The guess is usually the largest and some functions will be gone in the next iterations
			// 2. If some pairs are frozen the matrix is computed here in the first iteration and is then kept
			TIMER(timer3c);
			// no matrix routine yet, we have do to the expensive thing
			for (ElectronPairIterator it = pit(); it; ++it) {
				if(pno_ij[it.ij()].size()>0){
					if(pairs.frozen_ij[it.ij()]){
						MADNESS_ASSERT(pairs.W_ij[it.ij()].ndim()==2);
						MADNESS_ASSERT(size_t(pairs.W_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
						msg << pairs.name(it) << " frozen pair: no recomputation of fluctuation matrix\n";
					}
					else if(pairs.W_ij[it.ij()].size()>0){
						MADNESS_ASSERT(pairs.W_ij[it.ij()].ndim()==2);
						MADNESS_ASSERT(size_t(pairs.W_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
						msg <<  pairs.name(it) << " fluctuation matrix already computed\n";
					}
					else if(do_cispd) W_ij[it.ij()] = compute_cispd_fluctuation_matrix(it, pairs);
					else W_ij[it.ij()]=compute_fluctuation_matrix(it, pno_ij[it.ij()], pairs.Kpno_ij[it.ij()]);
				}
				else W_ij[it.ij()] = Tensor<double>(std::vector<long>(2,0));
			}

			timer3c.stop().print("Whole Direct Fluctuation Matrix");
		}

		pairs.update_meminfo();
		msg << pairs.meminfo << "\n";
		truncate_pairs(pairs);

		// 3.d compute PNO overlaps
		TIMER(timer3d);
		for (ElectronPairIterator it = pit(); it; ++it) {
			const size_t ij = it.ij();
			const size_t i = it.i();
			const size_t j = it.j();
			for (OrbitalIterator kit = oit(); kit; ++kit) {
				const size_t k = kit.i();
				if (not S_ij_kj.is_initialized(i, j, k)) {
					const auto kj = f12.pit().tridx(k, j);
					Tensor<double> s_ij_kj(std::vector<long>(2,0));
					if(not(pno_ij[ij].empty() or pno_ij[kj].empty())){
						s_ij_kj = matrix_inner(world, pno_ij[ij], pno_ij[kj], false);
					}
					S_ij_kj.set(i, j, k, s_ij_kj);
				}
				if (not S_ij_ik.is_initialized(i, j, k)) {
					const auto ik = f12.pit().tridx(i, k);
					Tensor<double> s_ij_ik(std::vector<long>(2,0));
					if(not(pno_ij[ij].empty() or pno_ij[ik].empty())){
						s_ij_ik= matrix_inner(world, pno_ij[ij], pno_ij[ik], false);
					}
					S_ij_ik.set(i, j, k, s_ij_ik);
				}
			}

		}
		timer3d.stop().print("PNO Overlaps");

		// 3. solve for the amplitudes
		TIMER(timer3);
		t_solve(pairs,fmat);
		energy=pairs.energies.total_energy();

		timer3.stop().print("Solve Amplitudes");

		// 4. compute the energy -> done in the T solver loop

		// 5. check convergence
		econverged = (std::abs(last_iter_energy - energy) < econv);
		const double delta_energy = std::abs(last_iter_energy - energy);

		bool do_compression = true;

		// 6. compress the amplitudes
		if (do_compression) {
			TIMER(timer6);
			if (iter == 0 or (iter == 1 and no_opt_in_first_iteration)) {
				msg << "No optimization untill now: Compressing pno ranks with tighter threshold\n";
				msg << param.tpno_tight() << " instead of " << param.tpno() << "\n";
				rdm_evals_ij = pno_compress(pairs,param.tpno_tight());
			} else {
				rdm_evals_ij = pno_compress(pairs,param.tpno());
			}
			pairs.rdm_evals_ij = rdm_evals_ij;

			// test if ampltides are transformed correctly
			//PairEnergies edb=compute_projected_mp2_energies(t2_ij, pno_ij);
			//f12.print_pair_energies(edb.eijs,edb.eijt,std::to_string(energytype()) + " MP2 pair energies after compression:");
			timer6.stop().print("PNO Compression");
		} else if (world.rank() == 0)
			std::cout << "Convergence reached or last iteration. No compression\n";

		print_ranks(pairs);
		save_pnos(pairs);
		// see if we can stop
		if (iter == size_t(maxiter - 1)) {
			if (world.rank() == 0)
				std::cout << "exiting: Last iteration\n";
			break;
		}
		if (no_opt) {
			if (world.rank() == 0)
				std::cout << "exiting: Found no_opt keyword\n";
			break;
		}
		if (converged) {
			if (world.rank() == 0)
				std::cout << "exiting: Iterations converged\n";
			break;
		}

		// update PNOs, unless converged
		bool do_optimization = true;
		if (no_opt_in_this_iteration){
			do_optimization = false;
		}
		if (do_optimization) {

			truncate_pairs(pairs);

			// 8. update the PNO subspaces
			dconverged = update_pno(pairs, rdm_evals_ij,fmat);

			// 9. Orthonormalize
			pairs = orthonormalize_cholesky(pairs);

		} else {
			if (world.rank() == 0) {
				std::cout << "PNO optimization not performed ";
				if (iter == 0)
					std::cout << "in the first iteration\n";
				else if (no_opt)
					std::cout << "since no_opt keyword was given\n";
				else if (econverged and dconverged)
					std::cout << "since the PNOs are already converged\n";
				else if (iter == size_t(param.maxiter() - 1))
					std::cout << "since this is the last iteration\n";
				else
					std::cout << "for unknown reasons ... should not happen!!!!!!!\n";
			}
		}

		if ((no_opt_in_first_iteration and iter<2) and (no_opt==false)) {
			msg << "Set back econverged since there was no optimization in the first iteration\n";
			econverged = false; // since there is no optimization in iteration 0
			converged=false;
		}
		save_pnos(pairs);
		++iter;
		if ((econverged == true) and (dconverged == true)) {
			converged = true;
			print("\n\nCONVERGED!");
			if (world.rank() == 0)
				std::cout << "energy change = " << delta_energy << "\n\n";
			if (world.rank() == 0)
				std::cout << "re-entering loop to compute the correct amplitudes and energies\n";
		} else {
			if (world.rank() == 0)
				std::cout << "\n\n" << "No convergence yet. Energy change = " << delta_energy << "\n\n";
		}
		timer_iter.stop().print("Whole Iteration " + std::to_string(iter));
	}
	if (no_compute == true) {
		for (ElectronPairIterator it = pit(); it; ++it)
			t2_ij[it.ij()] = madness::Tensor<double>(pno_ij[it.ij()].size(), pno_ij[it.ij()].size());
	}
	// give information about the rank
	{
		if (world.rank() == 0)
			std::cout << "PNO ranks:\n";

		double average_rank = 0.0;
		for (ElectronPairIterator it = pit(); it; ++it) {
			if (world.rank() == 0)
				std::cout << "PNO Rank for pair " << it.name() << " = " << pno_ij[it.ij()].size() << "\n";

			average_rank += pno_ij[it.ij()].size();
		}
		average_rank = average_rank / double(pno_ij.size());
		if (world.rank() == 0)
			std::cout << "average rank " << average_rank << "\n";
	}

	timer_solve.stop().print("PNO::solve");
	msg.subsection("ITERATE-PAIRS-INTERNAL-ENDED");

	pairs.t_ij=t2_ij;
	pairs.pno_ij=pno_ij;
	return pairs;
}

PNOPairs PNO::grow_rank(PNOPairs& pairs, std::string exop)const{
	msg << "Growing ranks with exop=" << exop << "\n";
	bool no_need=true;
	PAIRLOOP(it)
	{
	  if (pairs.pno_ij[it.ij()].size() < size_t(param.maxrank())) no_need = false;
	}
	if (no_need) {
		msg << "maxrank reached for all pairs, no need to grow ranks\n";
	} else {
		TIMER(timer);
		PAIRLOOP(it)
		{
		  if (pairs.pno_ij[it.ij()].size() >= size_t(param.maxrank())) {
				msg << "maxrank reached for " << it.name() << " no further growing\n";
				continue;
			} else if (pairs.frozen_ij[it.ij()]) {
				msg << it.name() << " is frozen: no growing\n";
				continue;
			} else {
				vector_real_function_3d pair_mo;
				pair_mo.push_back(nemo.get_calc()->amo[it.i() + param.freeze()]);
				if (not it.diagonal())
					pair_mo.push_back(nemo.get_calc()->amo[it.j() + param.freeze()]);
				bool do_cispd = true;
				if (pairs.cis.x.empty())
					do_cispd = false;
				if (do_cispd) {
					const double normi = pairs.cis.x[it.i()].norm2(); // cis vector contains only active orbitals
					const real_function_3d xi = (1.0 / normi) * pairs.cis.x[it.i()];
					pair_mo.push_back(xi);
					if (not it.diagonal()) {
					        //const double normj = pairs.cis.x[it.j()].norm2(); // cis vector contains only active orbitals
						const real_function_3d xj = (1.0 / normi) * pairs.cis.x[it.j()];
						pair_mo.push_back(xj);
					}
				}

				vector_real_function_3d virtij = Q(basis.guess_with_exop(pair_mo, exop,param.exop_trigo()));// guess_virtuals(pair_mo, EXOP_TYPE);
				// project out already existing pno pairs
				if (not pairs.pno_ij[it.ij()].empty()) {
					QProjector<double, 3> Qpno(world, pairs.pno_ij[it.ij()], pairs.pno_ij[it.ij()]);
					virtij = Qpno(virtij);

				}

				msg << "Created " << virtij.size() << " virtuals for " << it.name() << "\n";
				pairs.pno_ij[it.ij()] = append(pairs.pno_ij[it.ij()], virtij);

				// delete all intermediates and potentials
				MADNESS_ASSERT(pairs.frozen_ij[it.ij()]==false);
				pairs.clear_intermediates(it);
			}
		}
		pairs=orthonormalize_cholesky(pairs);
		timer.stop().print("grow PNO subspace");
	}

	return pairs;
}

PNOPairs PNO::orthonormalize_cholesky(PNOPairs& pairs) const {
	TIMER(timer);
	std::valarray<vector_real_function_3d> result(pairs.pno_ij.size());
	std::valarray<Tensor<double> > U_ij(pairs.pno_ij.size());
	PAIRLOOP(it)
	{
		auto& pno=pairs.pno_ij[it.ij()];
		if(pairs.frozen_ij[it.ij()]){
			msg << "frozen pair: no orthonormalization\n";
			Tensor<double> unit(pno.size(),pno.size());
			for(size_t x=0;x<pno.size();++x) unit(x,x)=1.0;
			U_ij[it.ij()]=unit;
			continue;
		}
		else if(pno.empty()){
			pairs.t_ij[it.ij()]=Tensor<double>(std::vector<long>(2,0));
			pairs.W_ij[it.ij()]=Tensor<double>(std::vector<long>(2,0));
			pairs.F_ij[it.ij()]=Tensor<double>(std::vector<long>(2,0));
			U_ij[it.ij()]=Tensor<double>(std::vector<long>(2,0));
		}
		else {
			int rank=0;
			const double tol=FunctionDefaults<3>::get_thresh()*0.1;
			Tensor<integer> piv;
			auto ovlp=matrix_inner(world,pno,pno);
			madness::rr_cholesky(ovlp,tol,piv,rank); // destroys ovlp and gives back Upper ∆ Matrix from CCD

			// rearrange and truncate the functions according to the pivoting of the rr_cholesky
			vector_real_function_3d pp(rank);
			vector_real_function_3d pk(rank);
			vector_real_function_3d pw1(rank);
			vector_real_function_3d pw2(rank);
			for(integer i=0;i<rank;++i){
				pp[i] =pairs.pno_ij[it.ij()][piv[i]];
				if(not pairs.Kpno_ij[it.ij()].empty()) pk[i] =pairs.Kpno_ij[it.ij()][piv[i]];
				if(not pairs.W_ij_i[it.ij()].empty())  pw1[i]=pairs.W_ij_i[it.ij()][piv[i]];
				if(not pairs.W_ij_j[it.ij()].empty())  pw2[i]=pairs.W_ij_j[it.ij()][piv[i]];
			}
			pairs.pno_ij[it.ij()] =pp;
			if(not pairs.Kpno_ij[it.ij()].empty())pairs.Kpno_ij[it.ij()]=pk;
			if(not pairs.W_ij_i[it.ij()].empty())pairs.W_ij_i[it.ij()] =pw1;
			if(not pairs.W_ij_j[it.ij()].empty())pairs.W_ij_j[it.ij()] =pw2;
			ovlp=ovlp(Slice(0,rank-1),Slice(0,rank-1));

			Tensor<double> L = transpose(ovlp);
			Tensor<double> Linv = inverse(L);
			Tensor<double> U = transpose(Linv);
			U_ij[it.ij()]=U;

			if(size_t(rank)<pno.size()) msg << "RRCD: truncate ranks of " << std::setw(25) << pairs.name(it) << " from " << std::setw(3) << pno.size() <<  " to " << std::setw(3) << rank << "\n";
		}
	}
	timer.stop().print("RRCD");
	pairs=transform_pairs(pairs,U_ij);
	return pairs;
}

PairEnergies PNO::compute_cispd_f12_correction_es(const vector_real_function_3d& xcis, PairEnergies& energies) const {
	TIMER(timer);
	double s2b_f12 = 0.0;
	double s2c_f12 = 0.0;
	for (ElectronPairIterator it = pit(); it; ++it) {
		const real_function_3d moi = f12.acmos[it.i()];
		const real_function_3d moj = f12.acmos[it.j()];
		const real_function_3d xi = xcis[it.i()];
		const real_function_3d xj = xcis[it.j()];
		double factor = 1.0;
		if (it.diagonal())
			factor = 0.5;

		double s2b_f12_ij = 0.0;
		{
			// functional part |ph> + |hp>
			if (it.diagonal())
				s2b_f12_ij += 2.0 * (f12.compute_fQg_integral(xi, moj, xi, moj) + f12.compute_fQg_integral(moi, xj, xi, moj));
			else
				s2b_f12_ij += 2.0 * (f12.compute_fQg_integral(xi, moj, xi, moj) + f12.compute_fQg_integral(moi, xj, xi, moj))
				+ 2.0 * (f12.compute_fQg_integral(xj, moi, xj, moi) + f12.compute_fQg_integral(moj, xi, xj, moi))
				- 1.0 * (f12.compute_fQg_integral(xi, moj, xj, moi) + f12.compute_fQg_integral(moi, xj, xj, moi))
				- 1.0 * (f12.compute_fQg_integral(xj, moi, xi, moj) + f12.compute_fQg_integral(moj, xi, xi, moj));

			// projector-response part: OxQ + QOx
			if (it.diagonal()) {
				vector_real_function_3d p1 = Q(apply(world, *f12.fop, f12.acmos * moj) * moi);
				vector_real_function_3d p2 = xcis;
				const double pr_part = 2.0 * (madness::inner(world, xi * p1, apply(world, *poisson, p2 * moj)).sum() + madness::inner(world, xi * p2, apply(world, *poisson, p1 * moj)).sum());
				s2b_f12_ij -= pr_part;
			} else {
				vector_real_function_3d p1i = Q(apply(world, *f12.fop, f12.acmos * moj) * moi);
				vector_real_function_3d p1j = Q(apply(world, *f12.fop, f12.acmos * moi) * moj);
				vector_real_function_3d p2 = xcis;
				s2b_f12_ij -= 2.0
						* (madness::inner(world, xi * p1i, apply(world, *poisson, p2 * moj)).sum() + madness::inner(world, xj * p1j, apply(world, *poisson, p2 * moi)).sum()
								+ madness::inner(world, xi * p2, apply(world, *poisson, p1j * moj)).sum() + madness::inner(world, xj * p2, apply(world, *poisson, p1i * moi)).sum())
								- 1.0
								* (madness::inner(world, moj * p1i, apply(world, *poisson, p2 * xi)).sum() + madness::inner(world, moi * p1j, apply(world, *poisson, p2 * xj)).sum()
										+ madness::inner(world, moj * p2, apply(world, *poisson, p1j * xi)).sum() + madness::inner(world, moi * p2, apply(world, *poisson, p1i * xj)).sum());
			}
		}
		double s2c_f12_ij = 0.0;
		{
			// transform pnos according to transposed(Ox)|a>, Ox=|xk><k| -> transposed(Ox)=|k><x_k|
			Projector<double, 3> Oxt(xcis, f12.acmos);
			//double pr_part = 0.0;
			if (it.diagonal()) {
				const vector_real_function_3d& p1 = xcis;
				const vector_real_function_3d p2 = Q(apply(world, *poisson, moi * f12.acmos) * moj);
				s2c_f12_ij -= 2.0 * (madness::inner(world, p1 * xi, apply(world, *f12.fop, p2 * moj)).sum() + madness::inner(world, apply(world, *f12.fop, p1 * moi), p2 * xj).sum());
				//const double tmp = s2c_f12_ij;
				// projector response
				//Projector<double,3> Oxt(xcis,f12.acmos);
				const vector_real_function_3d Oxtp1 = Oxt(p1);
				const vector_real_function_3d Oxtp2 = Oxt(p2);
				s2c_f12_ij += 2.0 * (madness::inner(world, Oxtp1 * moi, apply(world, *f12.fop, p2 * moj)).sum() + madness::inner(world, apply(world, *f12.fop, p1 * moi), Oxtp2 * moj).sum());
			} else {
				const vector_real_function_3d& p1 = xcis;
				const vector_real_function_3d p2j = Q(apply(world, *poisson, moi * f12.acmos) * moj);
				const vector_real_function_3d p2i = Q(apply(world, *poisson, moj * f12.acmos) * moi);
				s2c_f12_ij -= 2.0 * (madness::inner(world, p1 * xi, apply(world, *f12.fop, p2j * moj)).sum() + madness::inner(world, apply(world, *f12.fop, p1 * moi), p2j * xj).sum())
																								+ 2.0 * (madness::inner(world, p1 * xj, apply(world, *f12.fop, p2i * moi)).sum() + madness::inner(world, apply(world, *f12.fop, p1 * moj), p2i * xi).sum())
																								- 1.0 * (madness::inner(world, p2j * xi, apply(world, *f12.fop, p1 * moj)).sum() + madness::inner(world, apply(world, *f12.fop, p2j * moi), p1 * xj).sum())
																								- 1.0 * (madness::inner(world, p2i * xj, apply(world, *f12.fop, p1 * moi)).sum() + madness::inner(world, apply(world, *f12.fop, p2i * moj), p1 * xi).sum());
				;
				// projector response
				//Projector<double,3> Oxt(xcis,f12.acmos);
				const vector_real_function_3d Op1 = Oxt(p1);
				const vector_real_function_3d Op2j = Oxt(p2j);
				const vector_real_function_3d Op2i = Oxt(p2i);
				s2c_f12_ij += 2.0 * (inner(world, Op1 * moi, apply(world, *f12.fop, p2j * moj)).sum() + inner(world, p1 * moi, apply(world, *f12.fop, Op2j * moj)).sum());
				s2c_f12_ij -= 1.0 * (inner(world, Op1 * moj, apply(world, *f12.fop, p2j * moi)).sum() + inner(world, p1 * moj, apply(world, *f12.fop, Op2j * moi)).sum());
				// conjugated pair
				s2c_f12_ij += 2.0 * (inner(world, Op1 * moj, apply(world, *f12.fop, p2i * moi)).sum() + inner(world, p1 * moj, apply(world, *f12.fop, Op2i * moi)).sum());
				s2c_f12_ij -= 1.0 * (inner(world, Op1 * moi, apply(world, *f12.fop, p2i * moj)).sum() + inner(world, p1 * moi, apply(world, *f12.fop, Op2i * moj)).sum());
			}
		}
		s2b_f12 += factor * (s2b_f12_ij);
		s2c_f12 += factor * (s2c_f12_ij);
		energies.eijt_f12[it.ij()]=factor*(s2b_f12_ij+s2c_f12_ij);
	}
	msg << "---------------------\n";
	msg << "CIS(D) ES Correction:\n";
	msg << "S2b_f12=" << s2b_f12 << "\n";
	msg << "S2c_f12=" << s2c_f12 << "\n";
	msg << "tot_f12=" << s2b_f12+s2c_f12 << "\n";
	msg << "---------------------\n";
	energies.update();
	timer.stop().print("CIS(D) F12-Correction ES");
	return energies;
}

madness::PairEnergies PNO::t_solve(PNOPairs& pairs, const Tensor<double>& F_occ,
		const double R_convergence_threshold, const size_t max_niter) const {

	// convenience assignement
	auto& t2_ij=pairs.t_ij;
	const auto& F_ij=pairs.F_ij;
	const auto& W_ij=pairs.W_ij;
	const auto& S_ij_ik=pairs.S_ij_ik;
	const auto& S_ij_kj=pairs.S_ij_kj;
	//const auto& pno_ij=pairs.pno_ij;
	//const auto& nocc_act=pairs.nocc;
	double omega=0.0;
	if(pairs.type==CISPD_PAIRTYPE) omega=pairs.cis.omega;
	double energy=0.0;

	//	const auto npno_total = std::accumulate(std::begin(W_ij), std::end(W_ij), 0.0, [](size_t partial_sum, const Tensor<double>& elem) {
	//		return partial_sum + elem.size();
	//	});
	//  do consistency check while at it
	auto npno_total=0;
	PAIRLOOP(it){
		npno_total+=pairs.W_ij[it.ij()].size();
		MADNESS_ASSERT(W_ij[it.ij()].ndim()==2);
		MADNESS_ASSERT(W_ij[it.ij()].dim(0)==W_ij[it.ij()].dim(1));
		MADNESS_ASSERT(size_t(W_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
	}

	size_t tdim = 0;
	for (const auto& t : W_ij) tdim += (t.dim(0) * t.dim(0));
	valarray_allocator valloc(tdim);
	XNonlinearSolver<std::valarray<double>, double, valarray_allocator> kain(valloc, true);
	kain.do_print = true;
	kain.set_maxsub(param.kain_subspace());
	const auto npairs = t2_ij.size();
	std::valarray<double> energy_ij_s(npairs); // singlet pair energies
	std::valarray<double> energy_ij_t(npairs); // triplet pair energies
	bool converged = false;
	auto iter = 0;
	while (!converged && size_t(iter) < max_niter) {
		std::valarray<Tensor<double> > R_ij(npairs); // ij -> <ij|R|ab>
		auto flatten = [npno_total](const std::valarray<Tensor<double> >& arr) {
			std::valarray<double> flattened_arr(npno_total);
			auto iter = std::begin(flattened_arr);
			for (const auto& tensor: arr) {
				std::copy(tensor.ptr(), tensor.ptr() + tensor.size(), iter);
				iter += tensor.size();
			}
			return flattened_arr;
		};
		auto unflatten = [&t2_ij,&pairs](const std::valarray<double>& flattened_arr) {
			std::valarray<Tensor<double> > arr(t2_ij.size());
			auto iter = std::begin(flattened_arr);
			for (size_t ij = 0;ij != t2_ij.size();++ij) {
				const auto npno = pairs.pno_ij[ij].size();
				arr[ij]=Tensor<double>(std::vector<long>(2,npno)); // this constructor also works for npno=0
				const auto npno2 = npno * npno;
				std::copy(iter, iter + npno2, arr[ij].ptr());
				iter += npno2;
			}
			return arr;
		};
		auto kain_update_t2 = [&kain, &t2_ij, &R_ij, &flatten, &unflatten]() {
			auto t2_ij_flattened = flatten(t2_ij);
			auto R_ij_flattened = flatten(R_ij);
			auto new_t2_ij_flattened = kain.update(t2_ij_flattened, R_ij_flattened);
			t2_ij = unflatten(new_t2_ij_flattened);
		};
		auto R_norm2 = 0.0;
		// W + F_vv terms deal with 1 pair at a time
		for (ElectronPairIterator it = pit(); it; ++it) {
			if(pairs.frozen_ij[it.ij()]){
				MADNESS_ASSERT(pairs.t_ij[it.ij()].ndim()==2);
				MADNESS_ASSERT(pairs.t_ij[it.ij()].dim(0)==pairs.t_ij[it.ij()].dim(1));
				MADNESS_ASSERT(size_t(pairs.t_ij[it.ij()].dim(0))==pairs.pno_ij[it.ij()].size());
				R_ij[it.ij()]=Tensor<double>(std::vector<long>(2,pairs.pno_ij[it.ij()].size()));
				continue;
			}
			const size_t ij = it.ij();
			//const size_t i = it.i();
			//const size_t j = it.j();
			const auto npno = W_ij[ij].dim(0);
			if(npno==0){
				t2_ij[ij]=Tensor<double>(std::vector<long>(2,0));
				continue;
			}
			// zero out the amplitudes if this is the first iteration
			if (iter == 0)
				t2_ij[ij] = Tensor<double>(npno, npno);
			R_ij[ij] = W_ij[ij] + inner(F_ij[ij], t2_ij[ij]) + inner(t2_ij[ij], F_ij[ij]);
		}
		// F_oo terms couple pairs
		for (ElectronPairIterator it = pit(); it; ++it) {
			if(W_ij[it.ij()].size()==0) continue;
			const size_t ij = it.ij();
			const size_t i = it.i();
			const size_t j = it.j();
			if(pairs.frozen_ij[it.ij()]){
				R_ij[ij]=Tensor<double>(std::vector<long>(2,pairs.pno_ij[it.ij()].size()));
				continue;
			}
			for (OrbitalIterator kit = oit(); kit; ++kit) {
				const size_t k = kit.i();
				const auto ik = f12.pit().tridx(i, k);
				const auto swap_ik = i < k;
				const auto T_ik = swap_ik ? t2_ij[ik].swapdim(0, 1) : t2_ij[ik];
				const auto kj = f12.pit().tridx(k, j);
				const auto swap_kj = k < j;
				const auto T_kj = swap_kj ? t2_ij[kj].swapdim(0, 1) : t2_ij[kj];
				auto s_ij_ik = S_ij_ik.get(i, j, k);
				auto s_ij_kj = S_ij_kj.get(i, j, k);
				if(s_ij_ik.size()>0 or T_ik.size()>0) R_ij[ij]-= F_occ(k, j) * inner(inner(s_ij_ik, T_ik), s_ij_ik.swapdim(0, 1));
				if(s_ij_kj.size()>0 or T_kj.size()>0) R_ij[ij]-= F_occ(i, k) * inner(inner(s_ij_kj, T_kj), s_ij_kj.swapdim(0, 1));
			}
			R_ij[ij] -= omega * t2_ij[it.ij()];
			// keep track of the residual norm
			const auto R_ij_norm = R_ij[ij].normf();
			R_norm2 += R_ij_norm * R_ij_norm;
		}
		converged = (R_norm2 <= R_convergence_threshold * R_convergence_threshold);
		// energy evaluation
		{
			energy = 0.0;
			for (ElectronPairIterator it = pit(); it; ++it) {
				if(pairs.frozen_ij[it.ij()]) continue;
				const size_t ij = it.ij();
				const size_t i = it.i();
				const size_t j = it.j();
				const auto perm_ij_factor = i == j ? 1 : 2;
				auto& tpno = t2_ij[ij];
				const auto npno = tpno.dim(0);
				if(npno==0) continue;
				auto& rpno = R_ij[ij];
				const auto& wpno = W_ij[ij];
				const auto& fpno = F_ij[ij];
				double eps_s = 0.0, eps_t = 0.0; // singlet + triplet pair energies
				for (size_t a = 0; a != size_t(npno); ++a) {
					for (size_t b = 0; b <= a; ++b) {
						const auto perm_ab_factor = a == b ? 1 : 2;
						const auto t_ab_s = (tpno(a, b) + tpno(b, a)) / 2;
						const auto wr_ab_s = ((wpno(a, b) + wpno(b, a)) + (rpno(a, b) + rpno(b, a))) / 2;
						eps_s += perm_ij_factor * perm_ab_factor * t_ab_s * wr_ab_s;
						if (i != j && a != b) {
							const auto t_ab_t = (tpno(a, b) - tpno(b, a)) / 2;
							const auto wr_ab_t = ((wpno(a, b) - wpno(b, a)) + (rpno(a, b) - rpno(b, a))) / 2;
							eps_t += 3 * perm_ij_factor * perm_ab_factor * t_ab_t * wr_ab_t;
						}
					}
				}
				// if we calculated CIS(D) we have ( 1/(1+dij) instead of 2/(1+dij))
				double factor = 1.0;
				if (omega > 0.0)
					factor = 0.5;

				energy += factor * (eps_s + eps_t);
				energy_ij_s[ij] = factor * eps_s;
				energy_ij_t[ij] = factor * eps_t;
				// Jacobi step on first iteration, KAIN afterwards
				for (size_t a = 0; a != size_t(npno); ++a) {
				  for (size_t b = 0; b != size_t(npno); ++b) {
						const auto denom = (F_occ(i, i) + F_occ(j, j) + omega - fpno(a, a) - fpno(b, b));
						rpno(a, b) = rpno(a, b) / denom; // rescale residue for kain update
						if (iter == 0 || !(param.kain())) {
							tpno(a, b) += rpno(a, b);
						} // jacobi update, += is correct
					} // b
				} // a
			} // it
			if (iter != 0 && param.kain()) {
				kain_update_t2();
			}
			if (param.debug()) {
				for (ElectronPairIterator it = pit(); it; ++it) {
					if (world.rank() == 0)
						std::cout << it.name() << " amplitudes after iter=" << iter << "\n" << t2_ij[it.ij()] << "\n";
				}
			}
			if (omega > 0.0)
				energy *= 0.5;

			if (world.rank() == 0)
				print("t_iter = ", iter, " Hylleraas energy = ", energy, " ||R||_F = ", sqrt(R_norm2), " conv = ", (converged ? "yes" : "no"));
		}
		++iter;
		if (!converged && iter == param.maxiter_t() && R_norm2 > 1.0)
			throw "MP1 equations did not converge at all";
	} // T solver loop
	if (!converged) {
		if (false) {
			throw "MP1 equations did not converge";
		} else {
			if (world.rank() == 0) {
				std::cout << "!!!No convergence in ampltiudes!!!: try to proceed\n";
			}
		}
	}

	if(pairs.type==MP2_PAIRTYPE and energytype()==HYLLERAAS_ENERGYTYPE){
		pairs.energies.energy = energy;
		pairs.energies.eijs = energy_ij_s;
		pairs.energies.eijt = energy_ij_t;
		f12.print_pair_energies(pairs,"Hylleraas MP2 Energies");
	}else if(pairs.type==CISPD_PAIRTYPE){
		pairs.energies=compute_cispd_correction_es(pairs.cis.x,pairs);
		f12.print_pair_energies(pairs,"CIS(D) Correction");
	}else if(pairs.type==MP2_PAIRTYPE){
		pairs.energies=compute_projected_mp2_energies(pairs);
		f12.print_pair_energies(pairs,"Projected MP2 Energies");
	}else msg.warning("unknown pairtype in t_solve\n");
	pairs.energies.update();
	return pairs.energies;
}

std::valarray<Tensor<double> > PNO::pno_compress(PNOPairs& pairs, const double tpno) const {

	TIMER(timer);

	// convenience assigmenet
	auto& t2_ij=pairs.t_ij;
	auto& pno_ij=pairs.pno_ij;
	//auto& W_ij_i=pairs.W_ij_i;
	//auto& W_ij_j=pairs.W_ij_j;
	//auto& W_ij=pairs.W_ij;
	//auto& S_ij_ik=pairs.S_ij_ik;
	//auto& S_ij_kj=pairs.S_ij_kj;
	const auto& nocc_act=pairs.nocc;
	const auto& maxranks=pairs.maxranks_ij;

	const auto npairs = PNOTensors::ntri(nocc_act);
	std::valarray<Tensor<double> > evals_ij(npairs);
	std::valarray<Tensor<double> > U_ij(npairs); // "unitary" rotations for each pair
	for (ElectronPairIterator it = pit(); it; ++it) {
		const size_t ij = it.ij();
		const size_t i = it.i();
		const size_t j = it.j();
		if(pno_ij[ij].empty()) continue;
		const int maxrank = maxranks[it.ij()];
		if (pairs.frozen_ij[it.ij()]) {
			msg << it.name() << " is frozen: no pno-compression" << "\n";
			const int pno_size = pno_ij[it.ij()].size();
			U_ij[ij] = Tensor<double>(pno_size, pno_size);
			U_ij[ij] = 0.0;
			for (int x = 0; x < pno_size; ++x)
				U_ij[ij](x, x) = 1.0;
			continue;
		}
		// to be consistent with other codes defining density here as in
		// Eq. 23, JCP 138 034106 (2013)
		const auto delta_ij = i == j ? 1. : 0.;
		const auto& T = t2_ij[ij];
		const auto T_t = T.swapdim(0, 1);
		const auto T_tilde = 4 * T - 2 * T_t;
		const auto T_tilde_t = T_tilde.swapdim(0, 1);
		auto Dij = (inner(T_tilde_t, T) + inner(T_tilde, T_t)) / (1 + delta_ij);
		const auto n = Dij.dim(0);
		// compute PNOs
		Tensor<double> U, evals;
		syev(Dij, U, evals);
		if (world.rank() == 0) print(it.name() + ": density evals:\n", evals);

		// truncate PNOs
		size_t npno = 0;
		for (long k = n - 1; k >= 0; --k) {
			if (evals(k) > tpno)
				++npno;
		}
		if (tpno < 0) {
			msg << "tpno<0: no truncation\n";
			npno = n;
		}
		if (npno > size_t(maxrank) and maxrank>0)
			npno = maxrank;

		if (npno != size_t(n) && world.rank() == 0)
			print(it.name() + ": # of PNOs reduced from ", n, " to ", npno, " maxrank=", maxrank);

		bool emptypair = false;
		if (npno == 0)
			emptypair = true;

		Tensor<double> U_tr;
		if (!emptypair)
			U_tr = U(_, Slice(n - npno, -1, 1));

		Tensor<double> r_tr;
		if (!emptypair)
			r_tr = evals(Slice(n - npno, -1, 1));

		evals_ij[ij] = r_tr;
		U_ij[ij] = U_tr;

	} // it
	// NOW rotate all pnos and intermediates
	timer.stop().print("compress PNOs");
	// has its own timer
	transform_pairs(pairs,U_ij);


	return evals_ij;
}

PNOPairs PNO::transform_pairs(PNOPairs& pairs, const std::valarray<Tensor<double> >& U_ij) const {
	// compress all functions to avoid fences
	TIMER(timer0);
	auto allf = pairs.extract(pairs.pno_ij);
	allf = append(allf, pairs.extract(pairs.Kpno_ij));
	allf = append(allf, pairs.extract(pairs.W_ij_i));
	allf = append(allf, pairs.extract(pairs.W_ij_j));
	compress(world, allf);
	timer0.stop().print("MRA compress");
	TIMER(time);
	PAIRLOOP(it)
	{
		const size_t ij = it.ij();
		const auto& U = U_ij[ij];
		if (U.size() == 0) {
			pairs.W_ij_i[ij] = vector_real_function_3d();
			pairs.W_ij_j[ij] = vector_real_function_3d();
			pairs.pno_ij[ij] = vector_real_function_3d();
			pairs.pno_ij[ij] = vector_real_function_3d();
			pairs.Kpno_ij[ij] = vector_real_function_3d();
			pairs.F_ij[ij] = Tensor<double>(std::vector<long>(2, 0));
			pairs.t_ij[ij] = Tensor<double>(std::vector<long>(2, 0));
			pairs.W_ij[ij] = Tensor<double>(std::vector<long>(2, 0));
		} else {
			if (pairs.F_ij[ij].size() > 0)
				pairs.F_ij[ij] = inner(inner(U.swapdim(0, 1), pairs.F_ij[ij]), U);
			if (pairs.t_ij[ij].size() > 0)
				pairs.t_ij[ij] = inner(inner(U.swapdim(0, 1), pairs.t_ij[ij]), U);
			if (pairs.W_ij[ij].size() > 0)
				pairs.W_ij[ij] = inner(inner(U.swapdim(0, 1), pairs.W_ij[ij]), U);
			if (pairs.W_ij_i[ij].size() > 0)
				pairs.W_ij_i[ij] = transform(world, pairs.W_ij_i[ij], U);
			if (pairs.W_ij_j[ij].size() > 0)
				pairs.W_ij_j[ij] = transform(world, pairs.W_ij_j[ij], U);
			if (pairs.pno_ij[ij].size() > 0)
				pairs.pno_ij[ij] = transform(world, pairs.pno_ij[ij], U);
			if (pairs.Kpno_ij[ij].size() > 0)
				pairs.Kpno_ij[ij] = transform(world, pairs.Kpno_ij[ij], U);
		}
		const double thresh = FunctionDefaults<3>::get_thresh();
		truncate(pairs.pno_ij[ij], thresh);
		truncate(pairs.W_ij_i[ij], thresh);
		truncate(pairs.W_ij_j[ij], thresh);
		truncate(pairs.Kpno_ij[ij], thresh);
	}
	// rotate overlaps
	for (ElectronPairIterator it = pit(); it; ++it) {
		const size_t ij = it.ij();
		const size_t i = it.i();
		const size_t j = it.j();
		auto& S_ij_kj = pairs.S_ij_kj;
		auto& S_ij_ik = pairs.S_ij_ik;
		for (OrbitalIterator kit = oit(); kit; ++kit) {
			const size_t k = kit.i();
			{
				const auto kj = f12.pit().tridx(k, j);
				//const auto ijk = f12.pit().tridx(ij, kj);
				//const auto swap_ij_kj = ij < kj;
				if (U_ij[kj].size() == 0 || U_ij[ij].size() == 0) {
					S_ij_kj.set(i, j, k, Tensor<double>(std::vector<long>(2, 0)));
				} else if (S_ij_kj.is_unique(i, j, k)) {
					if (!S_ij_kj.is_initialized(i, j, k))
						continue;

					S_ij_kj.set(i, j, k, inner(inner(U_ij[ij].swapdim(0, 1), S_ij_kj.get(i, j, k)), U_ij[kj]));
				}
			}
			{
				const auto ik = f12.pit().tridx(i, k);
				//const auto ijk = f12.pit().tridx(ij, ik);
				//const auto swap_ij_ik = ij < ik;
				if (U_ij[ik].size() == 0 || U_ij[ij].size() == 0) {
					S_ij_ik.set(i, j, k, Tensor<double>(std::vector<long>(2, 0)));
				} else if (S_ij_ik.is_unique(i, j, k)) {
					if (!S_ij_ik.is_initialized(i, j, k))
						continue;

					S_ij_ik.set(i, j, k, inner(inner(U_ij[ij].swapdim(0, 1), S_ij_ik.get(i, j, k)), U_ij[ik]));
				}
			}
		}
	}
	time.stop().print("Transform Pairs");
	return pairs;
}

/// Print information about the PNO rank of all pairs (i.e. the number of PNO functions)
void PNO::print_ranks(const PNOPairs& pairs)const{
	msg <<"Information About All Pairs:\n";
	size_t sum=0;
	size_t maxi=0;
	size_t s=0;
	PAIRLOOP(it){
		const std::string name=pairs.name(it);
		msg << name << std::setw(std::max(0,int(20-name.size())))  <<  " " <<" number="<<it.ij()<<" ";
		const size_t rank= pairs.pno_ij[it.ij()].size();
		msg<<"rank=" <<std::setw(3) <<rank<<" frozen=" << pairs.frozen_ij[it.ij()] << "\n";
		maxi=std::max(maxi,rank);
		sum+=rank;
		++s;
	}
	//const size_t lr=maxi;
	const double ar=double(sum)/double(s);
	msg << "number of pairs="  << s << "\n";
	msg << "average rank   ="  << sum/s << " ("<<ar<<")" << "\n";
	msg << "largest rank   ="  << maxi << "\n";

}

bool PNO::update_pno(PNOPairs& pairs, const std::valarray<Tensor<double> >& rdm_evals_ij, const Tensor<double>& F_occ) const {
	TIMER(timer);
	// convenience assignement
	const auto nocc_act = pairs.nocc;
	auto& pno_ij=pairs.pno_ij;
	auto& W_ij_i=pairs.W_ij_i;
	auto& W_ij_j=pairs.W_ij_j;
	//auto& W_ij=pairs.W_ij;
	auto& F_ij=pairs.F_ij;
	auto& t2_ij=pairs.t_ij;
	auto& S_ij_ik=pairs.S_ij_ik;
	auto& S_ij_kj=pairs.S_ij_kj;
	const double omega=pairs.cis.omega;
	if(pairs.type!=CISPD_PAIRTYPE) MADNESS_ASSERT(omega==0.0);

	bool converged = true;
	// compute densities between coupled pairs
	PNOTensors::Tensor_IJ_IK<double> D_ij_ik(nocc_act); // f12.pit().tridx(ij,ik) -> D_{a_ij,b_ik}
	PNOTensors::Tensor_IJ_KJ<double> D_ij_kj(nocc_act); // f12.pit().tridx(ij,kj) -> D_{a_ij,b_kj}
	{
		TIMER(time);
		for (ElectronPairIterator it = pit(); it; ++it) {
			const size_t ij = it.ij();
			const size_t i = it.i();
			const size_t j = it.j();
			const auto& T_ij = t2_ij[ij];
			if (T_ij.normf() == 0.0)
				continue;

			const auto Tt_ij = 4 * T_ij - 2 * T_ij.swapdim(0, 1);
			for (OrbitalIterator kit = oit(); kit; ++kit) {
				const size_t k = kit.i();
				// D_ij_kj
				const auto kj = f12.pit().tridx(k, j);
				if (t2_ij[kj].size() == 0.0) {
					D_ij_kj.set(i, j, k, Tensor<double>(std::vector<long>(2,0)));
				} else if (!D_ij_kj.is_initialized(i, j, k)) {
					const auto swap_kj = k < j;
					const auto T_kj = swap_kj ? t2_ij[kj].swapdim(0, 1) : t2_ij[kj];
					const auto s_ij_kj = S_ij_kj.get(i, j, k);
					const auto d_ij_kj = inner(inner(Tt_ij.swapdim(0, 1), s_ij_kj), T_kj) + inner(inner(Tt_ij, s_ij_kj), T_kj.swapdim(0, 1));
					if (param.debug() && world.rank() == 0)
						print("D{i,j}{k,j} (i=", i, ",j=", j, ",k=", k, "):", d_ij_kj);

					D_ij_kj.set(i, j, k, d_ij_kj);
				}

				// D_ij_ik
				const auto ik = f12.pit().tridx(i, k);
				if (t2_ij[ik].size() == 0.0) {
					D_ij_ik.set(i, j, k, Tensor<double>(std::vector<long>(2,0)));
				} else if (!D_ij_ik.is_initialized(i, j, k)) {
					const auto ik = f12.pit().tridx(i, k);
					const auto swap_ik = i < k;
					const auto T_ik = swap_ik ? t2_ij[ik].swapdim(0, 1) : t2_ij[ik];
					const auto s_ij_ik = S_ij_ik.get(i, j, k);
					const auto d_ij_ik = inner(inner(Tt_ij.swapdim(0, 1), s_ij_ik), T_ik) + inner(inner(Tt_ij, s_ij_ik), T_ik.swapdim(0, 1));
					if (param.debug() && world.rank() == 0)
						print("D{i,j}{i,k} (i=", i, ",j=", j, ",k=", k, "):", d_ij_ik);

					D_ij_ik.set(i, j, k, d_ij_ik);
				}
			}
		}
		time.stop().print("compute densities for G application");
	}

	// for each pair assemble the potentials
	vector_real_function_3d allVP;
	std::vector<poperatorT> allG;
	{
		TIMER(timer2);

		for (ElectronPairIterator it = pit(); it; ++it) {
			const size_t ij = it.ij();
			const size_t i = it.i();
			const size_t j = it.j();
			const auto& pno = pno_ij[ij];

			if (pairs.frozen_ij[it.ij()]) {
				msg << it.name() << " is frozen: no optimization" << "\n";
				continue;
			}
			if (pno.empty()){
				msg << it.name() << " is empty: no optimization" << "\n";
				continue;
			}
			const auto npno = pno.size();
			// compute density and energy density
			const auto& T = t2_ij[ij];
			if (T.normf() == 0.0)
				continue;

			const auto FT = inner(F_ij[ij], T);
			const auto T_t = T.swapdim(0, 1);
			// (T F)^\dagger = F^\dagger T^\dagger = F T^\dagger
			const auto TF_t = inner(F_ij[ij], T_t);
			const auto T_tilde = 4 * T - 2 * T_t;
			const auto T_tilde_t = T_tilde.swapdim(0, 1);
			auto Dij = inner(T_tilde_t, T) + inner(T_tilde, T_t);
			auto Eij = inner(T_tilde_t, FT) + inner(T_tilde, TF_t);
			//        print("pair {", i, ",", j, "}: density:\n", Dij);
			//        print("pair {", i, ",", j, "}: energy density:\n", Eij);
			// compute (J-K+V) | a >  (J here = 2J in the notes)
			//vector_real_function_3d Ja = J(pno);
			//vector_real_function_3d Ka = K(pno);
			//vector_real_function_3d Va = V(pno);
			//vector_real_function_3d JKVa = add(world, sub(world, Ja, Ka), Va);
			//truncate(world, JKVa, param.thresh);
			//const auto JKVa=J(pno_ij[ij])+V(pno_ij[ij])-pairs.Kpno_ij[ij];
			// this will contain all contributions to apply Green's function to
			vector_real_function_3d Vphi(npno);
			Tensor<double> Vphi_energy(npno);
			// compute Vphi and matching Green's function energy
			for (size_t e = 0; e != npno; ++e) {
				if (omega == 0.0)
					Vphi_energy(e) = -1 * ((Eij(e, e) / Dij(e, e)) - F_occ(i, i) - F_occ(j, j));
				else {
					Vphi_energy(e) = -1 * (((Eij(e, e)) / Dij(e, e)) - F_occ(i, i) - F_occ(j, j) - omega);
				}
				// QW terms
				for (size_t b = 0; b != npno; ++b) {
					if (b == 0)
						Vphi[e] = T_tilde(e, b) * W_ij_i[ij][b] + T_tilde(b, e) * W_ij_j[ij][b];
					else
						Vphi[e] += T_tilde(e, b) * W_ij_i[ij][b] + T_tilde(b, e) * W_ij_j[ij][b];
				}
				// Fock terms
				if (nocc_act != 1) {
					// need to include pair-pair terms if more than 1 pair
					for (size_t k = 0; k != nocc_act; ++k) {
						// ik
						if (k != j) {
							const auto ik = f12.pit().tridx(i, k);
							const auto& pno_ik = pno_ij[ik];
							const auto d_ij_ik = D_ij_ik.get(i, j, k);
							for (size_t c = 0; c != pno_ik.size(); ++c) {
								Vphi[e] -= pno_ik[c] * (d_ij_ik(e, c) * F_occ(k, j));
							}
						}
						// kj
						if (k != i) {
							const auto kj = f12.pit().tridx(k, j);
							const auto& pno_kj = pno_ij[kj];
							const auto d_ij_kj = D_ij_kj.get(i, j, k);
							for (size_t c = 0; c != pno_kj.size(); ++c) {
								Vphi[e] -= pno_kj[c] * (d_ij_kj(e, c) * F_occ(i, k));
							}
						}
					}
				}
				for (size_t c = 0; c != npno; ++c) {
					if (c != e) {
						Vphi[e] += Eij(e, c) * pno[c];
					}
				}
				Vphi[e] *= 1 / Dij(e, e);
				// diag term
				//Vphi[e] += JKVa[e];
			}

			allVP=append(allVP,Vphi);
			std::vector<poperatorT> bsh3 = make_bsh_operators(world, Vphi_energy);
			allG.insert(allG.end(),bsh3.begin(),bsh3.end());

			// set back all matrices for pairs which will be optimized
			pairs.t_ij[ij]=Tensor<double>(std::vector<long>(2,0));
			pairs.W_ij[ij]=Tensor<double>(std::vector<long>(2,0));
			pairs.F_ij[ij]=Tensor<double>(std::vector<long>(2,0));
		}
		timer2.stop().print("Prepare for G application");
	}
	if(K.type != "neglect"){
		// compute J potetial and add it with the precomputed K potential
		TIMER(timerK);
		const auto allK=pairs.extract(pairs.Kpno_ij);
		allVP -= allK;
		timerK.stop().print("add K");
	}
	TIMER(timer3);
	// clear intermediates (for frozen pairs the matrices will stay to avoid recomputation)
	PAIRLOOP(it) pairs.clear_intermediates(it);
	double tmps=get_size(world,allVP);
	msg << "size of all Vpsi potentials (without Vnuc and J) " << tmps << " GB \n";
	vector_real_function_3d allP=pairs.extract(pairs.pno_ij);
	vector_real_function_3d newP;
	size_t chunk=param.chunk();
	auto itVP=allVP.begin();
	auto itP=allP.begin();
	auto itG=allG.begin();
	std::vector<double> allR;
	while(allR.size()<allP.size()){
		const auto dist=std::distance(itVP,allVP.end());
		if(chunk>size_t(dist)) chunk=dist;
		// include the nuclear potential
		auto Vtmp=vector_real_function_3d(itVP,itVP+chunk);
		Vtmp+=(V(vector_real_function_3d(itP,itP+chunk)));
		Vtmp+=(J(vector_real_function_3d(itP,itP+chunk)));
		{
			double tmps=get_size(world,Vtmp);
			msg << "time=" << wall_time() << " chunk=" << chunk <<  " Vpsi+Vnuc+J=" << tmps << " GB\n";
		}
		truncate(world,Vtmp);
		{
			double tmps=get_size(world,Vtmp);
			msg << "time=" << wall_time() << " chunk=" << chunk <<  " Vpsi+Vnuc+J=" << tmps << " GB (after truncation)\n";
		}
		scale(world,Vtmp,-2.0);
		auto G=std::vector<poperatorT>(itG,itG+chunk);
		auto Gtmp=Q(apply(world,G,Vtmp));
		Vtmp.clear();
		auto Rnorm = norm2s<double,3>(world, Gtmp-vector_real_function_3d(itP,itP+chunk));
		newP.insert(newP.end(), Gtmp.begin(), Gtmp.end());
		// keep the norms
		allR.insert(allR.end(), Rnorm.begin(), Rnorm.end());
		itVP+=chunk;
		itP+=chunk;
		itG+=chunk;
	}

	timer3.stop().print("Apply G to all");


	TIMER(timer4);
	// reassemble without overwriting frozen pairs
	pairs.pno_ij=pairs.reassemble(newP,pairs.pno_ij);
	pairs.S_ij_ik.reset();
	pairs.S_ij_kj.reset();
	timer4.stop().print("apply Q and update");
	// print some detailed information about the errors
	auto iterr=allR.begin();
	PAIRLOOP(it){
		const size_t ij = it.ij();
		const size_t ps=pno_ij[ij].size();
		if (pairs.frozen_ij[it.ij()]) continue;
		if (ps==0) continue;
		MADNESS_ASSERT(iterr!=allR.end());
		const auto errij=std::vector<double>(iterr,iterr+ps);
		msg << pairs.name(it) << " residual norms:\n";
		msg << errij << "\n";
		msg << pairs.name(it) << " scaled residual norms:\n";
		std::vector<double> scerrij(errij.size());
		for(size_t a=0;a<scerrij.size();++a) scerrij[a]=errij[a]*sqrt(rdm_evals_ij[ij](a));
		msg << scerrij << "\n";
		iterr+=ps;
		if(*std::max_element(scerrij.begin(),scerrij.end())>param.dconv()) converged=false;
	}

	timer.stop().print("update all PNOs");
	return converged;
}

std::vector<PNO::poperatorT> PNO::make_bsh_operators(World& world, const tensorT& evals) const
{
	PROFILE_MEMBER_FUNC (SCF);
	//int nmo = evals.dim(0);
	std::vector<poperatorT> ops(evals.size());
	double tol = param.op_thresh();
	for (int i = 0;i < evals.size();++i) {
		double eps = evals(i);
		if (eps > 0) {
			MADNESS_EXCEPTION("Positive Eigenvalue for BSH Operator?????", 1);
		}
		ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps), nemo.get_calc()->param.lo(), tol));
	}
	return ops;
}

vector_real_function_3d PNO::compute_CIS_potentials(const vector_real_function_3d& xcis) const {
	real_function_3d tH = real_factory_3d(world);
	for (size_t i = 0; i < xcis.size(); ++i) {
		tH += f12.acmos[i] * xcis[i];
	}
	const vector_real_function_3d J = (*poisson)(tH) * f12.acmos;
	vector_real_function_3d K;
	for (size_t i = 0; i < xcis.size(); ++i) {
		real_function_3d Ki = real_factory_3d(world);
		for (size_t j = 0; j < xcis.size(); ++j)
			Ki += (*poisson)(f12.acmos[j] * f12.acmos[i]) * xcis[j];
		K.push_back(Ki);
	}
	vector_real_function_3d result = Q(2.0 * J - K);
	truncate(world, result, param.thresh());
	return result;
}

Tensor<double> PNO::compute_fluctuation_matrix(const ElectronPairIterator& it, const vector_real_function_3d& pnos, const vector_real_function_3d& Kpnos_in) const {
	if (param.f12()) {
		vector_real_function_3d Kpnos = Kpnos_in;
		if (Kpnos_in.empty())
			Kpnos = K(pnos);

		const double sizetmp = get_size(world, Kpnos);
		if (world.rank() == 0)
			std::cout << "Kpnos Intermediate:" << sizetmp << " Gbyte\n";

		return (f12.compute_regularized_fluctuation_matrix(it, pnos, Kpnos));
	} else {
		const vector_real_function_3d vi = mul(world, f12.acmos[it.i()], pnos, false);
		vector_real_function_3d vj;
		if (it.diagonal())
			vj = vi;
		else
			vj = mul(world, f12.acmos[it.j()], pnos, false);

		world.gop.fence();
		const vector_real_function_3d gvi = apply(world, *poisson, vi);
		return transpose(matrix_inner(world, vj, gvi));
	}
}

Tensor<double> PNO::compute_cispd_fluctuation_matrix(const ElectronPairIterator& it, PNOPairs& pairs) const {
	TIMER(timer);
	Tensor<double> result;
	// get intermediates
	vector_real_function_3d& pnos = pairs.pno_ij[it.ij()];
	vector_real_function_3d& Kpnos = pairs.Kpno_ij[it.ij()];
	if(Kpnos.empty()){
		TIMER(timerK);
		const auto Kp=K(pnos);
		Kpnos=Kp;
		timerK.stop().print("recompute Kpno");
	}
	std::pair<vector_real_function_3d, vector_real_function_3d> empty = std::make_pair(vector_real_function_3d(), vector_real_function_3d());
	const real_function_3d xi = pairs.cis.x[it.i()];
	const real_function_3d xj = pairs.cis.x[it.j()];
	const real_function_3d Kxi = pairs.cis.Kx[it.i()];
	const real_function_3d Kxj = pairs.cis.Kx[it.j()];
	const real_function_3d Vxi = pairs.cis.Vx[it.i()];
	const real_function_3d Vxj = pairs.cis.Vx[it.j()];
	const real_function_3d moi = f12.acmos[it.i()];
	const real_function_3d moj = f12.acmos[it.j()];
	const real_function_3d Ki = f12.acKmos[it.i()];
	const real_function_3d Kj = f12.acKmos[it.j()];
	const auto& x=pairs.cis.x;
	const auto& Vx=pairs.cis.Vx;
	if (param.f12()) {

		TIMER(time1);
		// xo part
		const auto PNO=std::make_pair(pnos,Kpnos);
		Tensor<double> xo_part = (f12.compute_regularized_fluctuation_matrix(PNO, empty, std::make_pair(xi, Kxi), std::make_pair(moj, Kj)));
		// ox part
		Tensor<double> ox_part;
		if (it.diagonal())
			ox_part = transpose(xo_part);
		else
			ox_part = (f12.compute_regularized_fluctuation_matrix(PNO, empty, std::make_pair(moi, Ki), std::make_pair(xj, Kxj)));
		time1.stop().print("CIS(D)-Wij: xo/ox").start();

		// fock residue
		Tensor<double> fr_xo_part = (f12.compute_xyab_integrals(Vxi, moj, pnos, pnos, f12.fop));
		Tensor<double> fr_ox_part;
		if (it.diagonal())
			fr_ox_part = transpose(fr_xo_part);
		else
			fr_ox_part = (f12.compute_xyab_integrals(moi, Vxj, pnos, pnos, f12.fop));
		time1.stop().print("CIS(D)-Wij: fr").start();

		// projector response
		Tensor<double> pr_part_1(pnos.size(), pnos.size());
		Tensor<double> pr_part_2(pnos.size(), pnos.size());
		{
			const Tensor<double> xpno = matrix_inner(world, x, pnos);
			Tensor<double> im1 = (f12.compute_regularized_fluctuation_matrix(std::make_pair(f12.acmos, f12.acKmos), PNO, std::make_pair(moi, Ki), std::make_pair(moj, Kj)));
			Tensor<double> im2;
			if (it.diagonal())
				im2 = transpose(im1);
			else
				im2 = (f12.compute_regularized_fluctuation_matrix(PNO, std::make_pair(f12.acmos, f12.acKmos), std::make_pair(moi, Ki), std::make_pair(moj, Kj)));
			for (size_t a = 0; a < pnos.size(); ++a) {
				for (size_t b = 0; b < pnos.size(); ++b) {
					for (size_t k = 0; k < f12.acmos.size(); ++k) {
						// projector response
						pr_part_1(a, b) += xpno(k, a) * im1(k, b);
						pr_part_2(a, b) += xpno(k, b) * im2(a, k);
					}
				}
			}
		}
		time1.stop().print("CIS(D)-Wij: pr").start();
		// commutator response
		Tensor<double> cr_part(pnos.size(), pnos.size());
		{
			const Tensor<double> Vxpno = matrix_inner(world, Vx, pnos);
			Tensor<double> im1 = (f12.compute_xyab_integrals(moi, moj, f12.acmos, pnos, f12.fop));
			Tensor<double> im2;
			if (it.diagonal())
				im2 = transpose(im1);
			else
				im2 = (f12.compute_xyab_integrals(moi, moj, pnos, f12.acmos, f12.fop));
			for (size_t a = 0; a < pnos.size(); ++a) {
				for (size_t b = 0; b < pnos.size(); ++b) {
					for (size_t k = 0; k < f12.acmos.size(); ++k) {
						// commutator response
						cr_part(a, b) += Vxpno(k, a) * im1(k, b);
						cr_part(a, b) += Vxpno(k, b) * im2(a, k);
					}
				}
			}
		}
		time1.stop().print("CIS(D)-Wij: cr").start();
		result = xo_part - fr_xo_part + ox_part - fr_ox_part - (pr_part_1 + pr_part_2) + cr_part;
	} else {
		TIMER(time1)
																				// xo part
																				Tensor<double> xo_part = (f12.compute_xyab_integrals(xi, moj, pnos, pnos, poisson));
		// ox part
		Tensor<double> ox_part = (f12.compute_xyab_integrals(moi, xj, pnos, pnos, poisson));
		time1.stop().print("CIS(D)-Wij: xo/ox").start();
		// projector response
		Tensor<double> pr_part(pnos.size(), pnos.size());
		{
			const Tensor<double> xpno = matrix_inner(world, x, pnos);
			Tensor<double> im1 = (f12.compute_xyab_integrals(moi, moj, f12.acmos, pnos, poisson));
			Tensor<double> im2 = (f12.compute_xyab_integrals(moi, moj, pnos, f12.acmos, poisson));
			for (size_t a = 0; a < pnos.size(); ++a) {
				for (size_t b = 0; b < pnos.size(); ++b) {
					for (size_t k = 0; k < f12.acmos.size(); ++k) {
						// projector response
						pr_part(a, b) += xpno(k, a) * im1(k, b);
						pr_part(a, b) += xpno(k, b) * im2(a, k);
					}
				}
			}
		}
		time1.stop().print("CIS(D)-Wij: pr").start();
		result = (xo_part + ox_part - pr_part);
		Kpnos.clear();
	}
	timer.stop().print("CIS(D)-Wij: total");
	pairs.W_ij[it.ij()]=result;
	return result;
}

void PNO::canonicalize(PNOPairs& pairs)const{
	pairs.rdm_evals_ij = std::valarray<Tensor<double> >(); // make clear that they are not valid anymore
	// diagonalize all Fock matrices
	std::valarray<Tensor<double> > U_ij(Tensor<double>(std::vector<long>(2,0)),pairs.npairs);
	PAIRLOOP(it){
		const size_t n=pairs.pno_ij[it.ij()].size();
		if(n==0) continue;

		Tensor<double> U, evals;
		syev(pairs.F_ij[it.ij()], U, evals);
		// F is laster transformed in transform_pairs
		U_ij[it.ij()]=U;
	}
	// transform all pairs and potenials
	transform_pairs(pairs, U_ij);
}

void PNO::save_pnos(const PNOPairs& pairs) const {
	TIMER(timer);
	for (ElectronPairIterator it = pit(); it; ++it) {
		std::string name = pairs.name(it);
		save_function(pairs.pno_ij[it.ij()], name);
	}
	timer.stop().print("save pnos on disc");
}

PNOPairs PNO::load_pnos(PNOPairs& pairs) const {
	MADNESS_ASSERT(pairs.pno_ij.size() == pit().npairs());
	std::valarray<vector_real_function_3d>& pno_ij = pairs.pno_ij;
	for (ElectronPairIterator it = pit(); it; ++it) {
		vector_real_function_3d tmp;
		std::string name = pairs.name(it);
		try {
			load_function(world, tmp, name);
			msg << "Found PNOs for " << name << " with rank " << tmp.size()
					<< "\n";
		} catch (...) {
			msg << "failed to find PNO " << name
					<< " initialize as zero function\nn";
		}
		pno_ij[it.ij()] = tmp;
	}
	return pairs;
}

void PNO::update_fluctuation_potentials(PNOPairs& pairs) const {
	TIMER(timer);
	PAIRLOOP(it)
	{
		if (pairs.frozen_ij[it.ij()]) {
			msg << pairs.name(it) << " is frozen: potential not computed\n";
			continue;
		}
		std::pair<vector_real_function_3d, vector_real_function_3d> W;
		if (pairs.type == CISPD_PAIRTYPE)
			compute_cispd_fluctuation_potential(it, pairs);
		else
			compute_fluctuation_potential(it, pairs);
	}
	timer.stop().print("Fluctuation Potentials");
}

/// the terms are expanded as follows:
/// Q (-J1 +K1) | i(1) >  < a(2) | j(2) >
///  +  Q | i(1) > < a(2) | -J(2) + K(2) | j(2) >
///  +  i(1) * \int \dr2 1/|r12| a(2) j(2)
/// the first line is zero due to orthogonality of a and j, the second line is zero due to action of Q on i
template<typename projector>
vector_real_function_3d PNO::compute_V_aj_i(
		const real_function_3d& moi, const real_function_3d& moj,
		const vector_real_function_3d& virtuals, const projector& Qpr) const {
	MADNESS_ASSERT(not param.f12());
	const vector_real_function_3d aj = mul(world, moj, virtuals); // multiply a \times j
	vector_real_function_3d gaj = apply(world, *poisson, aj); // \int \dr2 aj(2)/r12
	vector_real_function_3d Vaj_i = mul(world, moi, gaj);
	vector_real_function_3d Vaj_i1 = Qpr(Vaj_i);
	truncate(world, Vaj_i1, param.thresh());
	return Vaj_i1;
}

// used by CIS(D) (could also be used by MP2)
template<typename projector>
vector_real_function_3d PNO::compute_Vreg_aj_i(
		const real_function_3d& moi, const real_function_3d& moj,
		const vector_real_function_3d& virtuals, const projector& Qpr,
		const vector_real_function_3d& Kpno) const {
	MyTimer time = MyTimer(world).start();
	MADNESS_ASSERT(param.f12());
	const real_function_3d Ki = K(moi); // can be much faster with intermediates of Ki and Kxi
	const real_function_3d Kj = K(moj);
	vector_real_function_3d tmp = f12.apply_regularized_potential(moj, moi, Kj,
			Ki, virtuals, Kpno);
	vector_real_function_3d Vaj_i = Qpr(tmp);
	//truncate(world,Vaj_i,param.thresh);
	time.stop().print("compute_Vreg_aj_i");
	return Vaj_i;
}

vector_real_function_3d PNO::compute_Vreg_aj_i(const size_t& i,
		const size_t& j, const vector_real_function_3d& virtuals,
		const vector_real_function_3d& Kpno) const {
	MADNESS_ASSERT(param.f12());
	const real_function_3d& moi = f12.acmos[i];
	const real_function_3d& moj = f12.acmos[j];
	const real_function_3d& Ki = f12.acKmos[i];
	const real_function_3d& Kj = f12.acKmos[j];
	vector_real_function_3d tmp = f12.apply_regularized_potential(moj, moi, Kj,
			Ki, virtuals, Kpno);
	vector_real_function_3d Vaj_i = Q(tmp);
	truncate(world, Vaj_i, param.thresh());
	return Vaj_i;
}

vector_real_function_3d PNO::compute_Vreg_aj_i_fock_residue(
		const real_function_3d& ket1, const real_function_3d& ket2,
		const vector_real_function_3d& virtuals) const {
	MyTimer time = MyTimer(world).start();
	const vector_real_function_3d result = Q(
			apply(world, *f12.fop, ket2 * virtuals) * ket1);
	time.stop().print("compute_Vreg_aj_i_fock_residue");
	return result;
}

vector_real_function_3d PNO::compute_Vreg_aj_i_commutator_response(
		const real_function_3d& moi, const real_function_3d& moj,
		const vector_real_function_3d& virtuals,
		const vector_real_function_3d& Vx) const {
	MyTimer time = MyTimer(world).start();
	Projector<double, 3> OVx(f12.acmos, Vx);
	Projector<double, 3> OVxt(Vx, f12.acmos);
	const vector_real_function_3d part1 = OVx(
			apply(world, *f12.fop, virtuals * moj) * moi);
	const vector_real_function_3d part2 = Q(
			apply(world, *f12.fop, OVxt(virtuals) * moj) * moi);
	time.stop().print("compute_Vreg_aj_i_commutator_response");
	return (part1 + part2);
}

void PNO::check_orthonormality(const vector_real_function_3d& v) const {
	Tensor<double> ovlp = matrix_inner(world, v, v);
	for (int i = 0; i < ovlp.dim(0); ++i)
		ovlp(i, i) -= 1.0;
	double error = ovlp.normf() / ovlp.size();
	if (error > 1.e-14 && world.rank() == 0)
		print("orthonormality error: ", error);
}

} /* namespace madness */
