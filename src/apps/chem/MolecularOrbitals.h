/*
 * MolecularOrbitals.h
 *
 *  Created on: 11 Jul 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_MOLECULARORBITALS_H_
#define SRC_APPS_CHEM_MOLECULARORBITALS_H_

namespace madness {

template<typename T, std::size_t NDIM>
class MolecularOrbitals {
public:

	std::vector<Function<T,NDIM> > get_mos() const {
		return mos;
	}

	Tensor<double> get_eps() const {
		return eps;
	}

	std::vector<std::string> get_irreps() const {
		return irreps;
	}

	std::vector<int> get_localize_sets() const {
		return localize_sets;
	}

	std::vector<double> get_occ() const {
		return occ;
	}

	/// setters will always invalidate all other member variables
	void set_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		invalidate_all();
		mo=mo_new;
	}

	/// updates will keep other member variables
	void update_mos(const std::vector<Function<T,NDIM> >& mo_new) {
		mo=mo_new;
	}

	void update_occ(const std::vector<int>& occ_new) {
		occ=occ_new;
	}

	/// updates will keep other member variables
	void update_mos_and_eps(const std::vector<Function<T,NDIM> >& mo_new,
			const Tensor<double>& eps_new) {
		mo=mo_new;
		eps=eps_new;
	}

	void recompute_irreps() {
	}

	void recompute_localize_sets() {
	}

	void invalidate_all() {
		invalidate_mos();
		invalidate_eps();
		invalidate_irreps();
		invalidate_occ();
		invalidate_localize_sets();
	}

	void invalidate_mos() {
		mo.clear();
	}

	void invalidate_eps() {
		eps.clear();
	}

	void invalidate_irreps() {
		irreps.clear();
	}

	void invalidate_occ() {
		occ.clear();
	}

	void invalidate_localize_sets() {
		localize_sets.clear();
	}

	template <typename Archive>
	void serialize (Archive& ar) {
		ar & mo & eps & irreps & occ & localize_sets;
	}

private:
	std::vector<Function<T,NDIM> > mo;
	std::vector<double> eps;
	std::vector<std::string> irreps;
	std::vector<double> occ;
	std::vector<int> localize_sets;

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_MOLECULARORBITALS_H_ */
