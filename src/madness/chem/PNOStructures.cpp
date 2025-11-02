/*
 * PNOStructures.cpp
 *
 *  Created on: Oct 24, 2018
 *      Author: kottmanj
 */

#include <PNOStructures.h>

namespace madness {

/// inner product for std::valarray
double inner(const std::valarray<double>& bra, const std::valarray<double>& ket) {
	assert(bra.size() == ket.size());
	const auto n = bra.size();
	double result = 0;
	const auto* bra_data = &bra[0];
	const auto* ket_data = &ket[0];
	for (size_t i = 0; i != n; ++i, ++bra_data, ++ket_data)
		result += *bra_data * *ket_data;
	return result;
}

PairEnergies PairEnergies::operator +=(const PairEnergies& right) {
	MADNESS_ASSERT(eijs.size() == right.eijs.size());
	MADNESS_ASSERT(eijt.size() == right.eijt.size());
	MADNESS_ASSERT(eijs_f12.size() == right.eijs_f12.size());
	MADNESS_ASSERT(eijt_f12.size() == right.eijt_f12.size());
	MADNESS_ASSERT(eij.size() == right.eij.size());
	eijs += right.eijs;
	eijt += right.eijs;
	eijs_f12 += right.eijs_f12;
	eijt_f12 += right.eijt_f12;
	eij += right.eij;
	energy += right.energy;
	energy_f12 += right.energy_f12;
	return *this;
}

PairEnergies PairEnergies::operator +(const PairEnergies& right) const {
	PairEnergies result(eij.size());
	MADNESS_ASSERT(result.eijs.size() == right.eijs.size());
	MADNESS_ASSERT(result.eijt.size() == right.eijt.size());
	MADNESS_ASSERT(result.eijs_f12.size() == right.eijs_f12.size());
	MADNESS_ASSERT(result.eijt_f12.size() == right.eijt_f12.size());
	MADNESS_ASSERT(result.eij.size() == right.eij.size());
	result.eijs = eijs + right.eijs;
	result.eijt = eijt + right.eijt;
	result.eijs_f12 = eijs_f12 + right.eijs_f12;
	result.eijt_f12 = eijt_f12 + right.eijt_f12;
	result.eij = eij + right.eij;
	result.energy = energy + right.energy;
	result.energy_f12 = energy_f12 + right.energy_f12;
	return result;
}

ElectronPairIterator& ElectronPairIterator::operator ++() {
	if (finished_)
		return *this;
	else {
		if (j_ < i_)
			++j_;
		else if (i_ < (stop_ - 1)) {
			++i_;
			j_ = start_;
		} else
			finished_ = true;
	}
	++ij_;
	// check consistency with pno-mp2.cc code
	if (!finished_)
		MADNESS_ASSERT((ij_ = tridx(i_, j_))); //-offset_);

	return *this;
}

bool PNOPairs::is_consistent(std::string& errm) const {
	if (npairs!=nocc*(nocc+1)/2) {
		errm = "PNOPairs: size inconsistency between nocc and npairs";
		return false;
	}
	if (pno_ij.size() != npairs) {
		errm = "PNOPairs: size inconsistency between pno_ij and npairs";
		return false;
	}
	if (pno_ij.size() != t_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and t_ij";
		return false;
	}
	if (pno_ij.size() != frozen_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and frozen_ij";
		return false;
	}
	if (pno_ij.size() != maxranks_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and maxranks_ij";
		return false;
	}
	if (type == UNKNOWN_PAIRTYPE) {
		errm = "PNOPair: pairtype is not determined";
		return false;
	}
	if (type == CISPD_PAIRTYPE && !cis.initialized()) {
		errm = "PNOPair: CIS(D) pairtype but CIS was not initialized";
		return false;
	}
	if (type != CISPD_PAIRTYPE && cis.initialized()) {
		errm = "PNOPair: CIS initialized but pairtype is not CIS(D)";
		return false;
	}
	if (pno_ij.size() != Kpno_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and Kpno_ij";
		return false;
	}
	if (pno_ij.size() != W_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and W_ij";
		return false;
	}
	if (pno_ij.size() != F_ij.size()) {
		errm = "PNOPairs: size inconsistency between pno_ij and F_ij";
		return false;
	}
	return true;
}

vecfuncT ParametrizedExchange::operator ()(const vecfuncT& vket,
		const double& mul_tol) const {
	if (type == "neglect") {
		return zero_functions<double, 3>(world, vket.size());
	} else if (type == "full") {
//		return K(vket, mul_tol);
        print("set mul_tol =0");
		return K(vket);

	} else {
		// Occupation numbers
		const Tensor<double> occ = nemo.get_calc()->get_aocc();
		// Closed shell full density of the nemo orbitals (without any nuclear cusps)
		const real_function_3d nemo_density = 2.0
				* nemo.make_density(occ, nemo.get_calc()->amo);
		// Real Alpha density (with nuclear cusps)
		const real_function_3d alpha_density = 0.5 * nemo.R_square
				* nemo_density;
		std::string xc_data = type;
		xc_data = xc_data.erase(0, xc_data.find_first_not_of(" "));
		xc_data = xc_data.erase(xc_data.find_last_not_of(" ") + 1);
		const XCOperator<double,3> xc(world, xc_data,
				!nemo.get_calc()->param.spin_restricted(), alpha_density,
				alpha_density);
		real_function_3d xc_pot = xc.make_xc_potential();
		return xc_pot * vket;
	}
}

void PNOPairs::initialize(const std::size_t& nocc) {
	const size_t n = nocc * (nocc + 1) / 2;
	pno_ij = std::valarray<vector_real_function_3d>(n);
	t_ij = std::valarray<Tensor<double> >(
			Tensor<double>(std::vector<long>(2, 0)), n);
	W_ij = std::valarray<Tensor<double> >(
			Tensor<double>(std::vector<long>(2, 0)), n);
	F_ij = std::valarray<Tensor<double> >(
			Tensor<double>(std::vector<long>(2, 0)), n);
	maxranks_ij = std::valarray<int>(-1, n);
	frozen_ij = std::valarray<bool>(false, n);
	energies = PairEnergies(n);
	Kpno_ij = std::valarray<vector_real_function_3d>(n);
	W_ij_i = std::valarray<vector_real_function_3d>(n);
	W_ij_j = std::valarray<vector_real_function_3d>(n);
}

PNOPairs PNOPairs::operator =(const PNOPairs& other) {
	MADNESS_ASSERT(type == other.type);
	pno_ij = other.pno_ij;
	t_ij = other.t_ij;
	frozen_ij = other.frozen_ij;
	maxranks_ij = other.maxranks_ij;
	energies = other.energies;
	return *this;
}

std::string PNOPairs::name(
		const ElectronPairIterator& it) const {
	std::string pre;
	std::stringstream ss;
	ss << type;
	ss >> pre;
	return pre + "_" + it.name();
}

vector_real_function_3d PNOPairs::extract(const vfT& vf) const {
	vector_real_function_3d result;
	for (size_t ij = 0; ij < pno_ij.size(); ++ij) {
		if (frozen_ij[ij])
			continue;
		else
			result = append(result, vf[ij]);
	}
	return result;
}

PNOPairs::vfT PNOPairs::reassemble(
		const vector_real_function_3d& v, vfT& result) const {
	MADNESS_ASSERT(result.size() == npairs);
	if (v.empty())
		return result;

	auto itv = v.begin();
	for (size_t ij = 0; ij < pno_ij.size(); ++ij) {
		if (frozen_ij[ij])
			continue;

		const auto& x = pno_ij[ij];
		MADNESS_ASSERT(itv < v.end());
		result[ij] = vector_real_function_3d(itv, itv + x.size());
		itv += x.size();
	}
	MADNESS_ASSERT(result.size() == npairs);
	return result;
}

void PNOPairs::clear_intermediates(const ElectronPairIterator& it) {
	S_ij_ik.reset();
	S_ij_kj.reset();
	Kpno_ij[it.ij()].clear();
	W_ij_i[it.ij()].clear();
	W_ij_j[it.ij()].clear();
	if (frozen_ij[it.ij()] == false)
		W_ij[it.ij()] = Tensor<double>(std::vector<long>(2, 0));

	if (frozen_ij[it.ij()] == false)
		F_ij[it.ij()] = Tensor<double>(std::vector<long>(2, 0));

	if (frozen_ij[it.ij()] == false)
		t_ij[it.ij()] = Tensor<double>(std::vector<long>(2, 0));

	update_meminfo();
}

PNOPairs::MemInfo PNOPairs::update_meminfo() const {
	meminfo = MemInfo();
	if (pno_ij.size() > 0 && pno_ij[0].size() > 0
			&& pno_ij[0][0].is_initialized()) {
		World& world = pno_ij[0][0].world();
		meminfo.pno = get_size(world, extract(pno_ij));
		meminfo.Kpno = get_size(world, extract(Kpno_ij));
		meminfo.W = get_size(world, extract(W_ij_i))
				+ get_size(world, extract(W_ij_j));
	}
	return meminfo;
}


} /* namespace madness */
