/*
 * CISOperators.cc
 *
 *  Created on: Aug 28, 2015
 *      Author: kottmanj
 */

#include "CISOperators.h"

namespace madness {


CIS_Operators::CIS_Operators(World&world, const Nemo &nemo, const vecfuncT &mos) :
		world(world),
		xcoperator(world,&nemo),
		use_nuclear_correlation_factor_(true),
		mo_ket_(mos),
		orbital_energies_(nemo.get_calc()->aeps),
		R2(init_R2(nemo))
				{
	if (nemo.nuclear_correlation->type()
			== NuclearCorrelationFactor::None) {
		std::cout << "No nuclear correlation factor used" << std::endl;
		use_nuclear_correlation_factor_ = false;
	}
	poisson = std::shared_ptr < real_convolution_3d
			> (CoulombOperatorPtr(world, nemo.get_calc()->param.lo,
					1.e-6));
	if (use_nuclear_correlation_factor_)
		mo_bra_ = mul(world, R2, mo_ket_);
	else
		mo_bra_ = mo_ket_;
	dR2 = get_gradient(R2);
	plot_plane(world, dR2[0], "dxR2");
	set_thresh(world, mo_bra_, FunctionDefaults<3>::get_thresh());
	set_thresh(world, mo_ket_, FunctionDefaults<3>::get_thresh());
	truncate(world, mo_ket_);
	truncate(world, mo_bra_);
	exchange_intermediate_ = make_exchange_intermediate();
	sanitycheck();
	//		START_TIMER();
	//		CC_Operators tmp(world,nemo);
	//		END_TIMER("CC TEST INIT OPERATORS");
}

/// Make shure that R2 gets the right thresh and is constant
real_function_3d CIS_Operators::init_R2(const Nemo &nemo) const {
	if (nemo.nuclear_correlation) {
		real_function_3d tmp = copy(nemo.nuclear_correlation->square());
		tmp.set_thresh(FunctionDefaults<3>::get_thresh());
		tmp.truncate();
		tmp.verify();
		return tmp;
	}
	real_function_3d constant = real_factory_3d(world);
	return (constant + 1.0);
}

// Make the derivative of R2
vecfuncT CIS_Operators::get_gradient(const real_function_3d f) const {
	std::vector < std::shared_ptr<real_derivative_3d> > gradop =
			gradient_operator<double, 3>(world);
	f.verify();
	vecfuncT gradf;
	for (size_t i = 0; i < 3; i++) {
		real_function_3d dfi = (*gradop[i])(f);
		gradf.push_back(dfi);
	}
	set_thresh(world, gradf, FunctionDefaults<3>::get_thresh());
	truncate(world, gradf);
	return gradf;
}

void CIS_Operators::sanitycheck() const {
	if (mo_ket_.empty())
		error("mo_ket_ is empty");
	if (mo_bra_.empty())
		error("mo_bra_ is empty");
	if (exchange_intermediate_.empty())
		error("exchange intermediate is empty");
	for (auto x : mo_ket_) {
		if (x.thresh() != FunctionDefaults<3>::get_thresh())
			error("Wrong thresh in mo_ket_ functions");
	}
	for (auto x : mo_bra_) {
		if (x.thresh() != FunctionDefaults<3>::get_thresh())
			error("Wrong thresh in mo_bra_ functions");
	}
	for (auto tmp : exchange_intermediate_) {
		if (tmp.empty())
			error("Exchange Intermediate contains empty vectors");
		for (auto x : tmp) {
			if (x.thresh() != FunctionDefaults<3>::get_thresh())
				error("Wrong thresh in Exchange Intermediate");
		}
	}
	R2.verify();
}

void CIS_Operators::memory_information(const vecfuncT &v, const std::string &msg) const {
	const double x = get_size(world, v);
	if (world.rank() == 0)
		std::cout << msg << "(" << x << " GB)" << " for " << v.size()
				<< " functions\n";
}

std::vector<vecfuncT> CIS_Operators::make_exchange_intermediate() const {
	std::vector<vecfuncT> intermediate;
	double memory = 0.0;
	for (size_t i = 0; i < mo_bra_.size(); i++) {
		const vecfuncT integrant = mul(world, mo_bra_[i], mo_ket_);
		const vecfuncT intermediate_i = apply(world, *poisson, integrant);
		intermediate.push_back(intermediate_i);
		memory += get_size(world, intermediate_i);
	}
	if (world.rank() == 0)
		std::cout << "Created exchange intermediate of dimension "
				<< intermediate.size() << "x" << intermediate.front().size()
				<< " and size (" << memory << " GB)\n";
	return intermediate;

}

// The Fock operator is partitioned into F = T + Vn + R
// the fock residue R= 2J-K for closed shell is computed here
// J_j = \sum_i <i|r12|i> |tau>
// K_j = \sum_i <i|r12|tau_j> |i>
vecfuncT CIS_Operators::fock_residue_closed_shell(const vecfuncT &tau) const {
	START_TIMER();
	vecfuncT J = mul(world, (*poisson)(make_density()), tau);
	truncate(world, J);
	scale(world, J, 2.0);
	END_TIMER("J");
	START_TIMER();
	vecfuncT K;
	for (size_t j = 0; j < tau.size(); j++) {
		real_function_3d Kj = real_factory_3d(world);
		for (size_t i = 0; i < mo_bra_.size(); i++) {
			Kj += (*poisson)(mo_bra_[i] * tau[j]) * mo_ket_[i];
		}
		K.push_back(Kj);
	}
	truncate(world, K);
	scale(world, K, -1);
	END_TIMER("K");
	return add(world, J, K);
}


// The same residue for the case that the Fock operator is the Kohn-Sham Operator
vecfuncT CIS_Operators::KS_residue_closed_shell(const vecfuncT &tau) const {
	START_TIMER();
	vecfuncT J = mul(world, (*poisson)(make_density()), tau);
	truncate(world, J);
	scale(world, J, 2.0);
	END_TIMER("J");
	START_TIMER();
	real_function_3d vxc = xcoperator.make_xc_potential();
	vxc.truncate();
	vecfuncT applied_vxc = mul(world, vxc, tau);
	truncate(world, applied_vxc);
	plot_plane(world,vxc,"vxc");
	plot_plane(world,applied_vxc.back(),"vxc_x_xhomo");
	plot_plane(world,add(world,J,applied_vxc).back(),"KS_residue_homo");
	END_TIMER("Vxc");
	return add(world, J, applied_vxc);
}

// Kinetik energy
// -1/2 <x|R2Nabla2|x> = +1/2 <Nabla R2 x | Nabla x> = grad(R2x)*grad(x)
// grad(R2x) = GradR2*x + R2*gradx
// grad(R2x)*grad(y) = GradR2*x*Grady + R2*Gradx*Grady
double CIS_Operators::get_matrix_element_kinetic_energy(const vecfuncT &ket,
		const vecfuncT &bra) const {
	//std::cout << " Making Kintic Energy Matrix Element 1:\n";
	double value = 0.0;
	vecfuncT R2bra = mul(world, R2, bra);
	truncate(world, R2bra);
	std::vector<vecfuncT> dket, dbra;
	std::vector < std::shared_ptr<real_derivative_3d> > gradop;
	gradop = gradient_operator<double, 3>(world);
	for (size_t axis = 0; axis < 3; axis++) {
		START_TIMER();
		const vecfuncT gradbra = apply(world, *gradop[axis], R2bra);
		END_TIMER("Gradient of R2Bra");
		START_TIMER();
		const vecfuncT gradket = apply(world, *gradop[axis], ket);
		END_TIMER("Gradient of Ket");
		START_TIMER();
		value += 0.5 * inner(world, gradbra, gradket).sum();
		END_TIMER("Inner Product");
	}
	return value;
}
double CIS_Operators::get_matrix_element_kinetic_2(const vecfuncT &bra,
		const vecfuncT &ket) const {
	std::cout << "Making Kinetic 2 Element\n";
	double value = 0.0;
	std::vector < std::shared_ptr<real_derivative_3d> > gradop =
			gradient_operator<double, 3>(world);
	for (size_t axis = 0; axis < 3; axis++) {
		START_TIMER();
		vecfuncT gradbra = apply(world, *gradop[axis], bra);
		truncate(world, gradbra);
		END_TIMER("make gradbra");
		START_TIMER();
		vecfuncT gradket = apply(world, *gradop[axis], ket);
		truncate(world, gradket);
		END_TIMER("make gradket");
		START_TIMER();
		vecfuncT gradR2bra = mul(world, dR2[axis], bra);
		truncate(world, gradR2bra);
		END_TIMER("multiply dR2 and bra");
		START_TIMER();
		value += (inner(world, gradR2bra, gradket).sum());
		END_TIMER("Inner product1");
		START_TIMER();
		value += make_inner_product(gradbra, gradket);
		END_TIMER("Inner product2");

	}
	return 0.5 * value;
}

// Kinetic part of the CIS perturbed fock matrix
Tensor<double> CIS_Operators::get_matrix_kinetic(const std::vector<vecfuncT> &x) const {
	Tensor<double> result(x.size(), x.size());
	// make the x,y and z parts of the matrix and add them
	std::vector < std::shared_ptr<real_derivative_3d> > gradop =
			gradient_operator<double, 3>(world);
	for (size_t axis = 0; axis < 3; axis++) {
		// make all gradients
		std::vector<vecfuncT> dx, dR2x;
		for (auto xi : x) {
			const vecfuncT dxi = apply(world, *(gradop[axis]), xi);
			dx.push_back(dxi);
			if (use_nuclear_correlation_factor_) {
				const vecfuncT dR2xi = apply(world, *(gradop[axis]),
						mul(world, R2, xi));
				dR2x.push_back(dR2xi);
			}
		}
		for (size_t i = 0; i < x.size(); i++) {
			for (size_t j = 0; j < x.size(); j++) {
				if (use_nuclear_correlation_factor_)
					result(i, j) += 0.5
							* inner(world, dR2x[j], dx[i]).sum();
				else
					result(i, j) += 0.5 * inner(world, dx[j], dx[i]).sum();
			}
		}

	}
	return result;
}

// Diagrammatic Potentials:

// The coulomb Term of the S3C diagram: Positive sign
// \     /
//  \---/  = Q\sum_j(<j|g12|tau_j>)|i>
//  _\_/_
vecfuncT CIS_Operators::S3C_C(const vecfuncT &tau) const {
	START_TIMER();
	vecfuncT result = mul(world, (*poisson)(make_density(mo_ket_, tau)),
			mo_ket_);
	Q(result);
	truncate(world, result);
	END_TIMER("S3C_C");
	return result;
}

// The Exchange Term of the S3C diagram: Negative sign
// \  /
//  \/...   = Q\sum_j(<j|g12|i>|tau_j>)
//     / \
//    _\_/_
vecfuncT CIS_Operators::S3C_X(const vecfuncT &tau) const {
	START_TIMER();
	vecfuncT result;
	for (size_t i = 0; i < tau.size(); i++) {
		real_function_3d tmp = real_factory_3d(world);
		vecfuncT vectmp = mul(world, exchange_intermediate_[i], tau);
		for (size_t j = 0; j < tau.size(); j++)
			tmp += vectmp[j];
		tmp.truncate();
		result.push_back(tmp);
	}
	Q(result);
	truncate(world, result);
	scale(world, result, -1.0);
	END_TIMER("S3C_X");
	return result;
}

// Project out the occupied space
void CIS_Operators::Q(vecfuncT &f) const {
	for (size_t i = 0; i < f.size(); i++)
		Q(f[i]);
}
void CIS_Operators::Q(real_function_3d &f) const {
	for (size_t i = 0; i < mo_ket_.size(); i++) {
		f -= mo_bra_[i].inner(f) * mo_ket_[i];
	}
}

// Make an inner product between vecfunctions
double CIS_Operators::make_inner_product(const vecfuncT &bra, const vecfuncT &ket) const {
	if (use_nuclear_correlation_factor_)
		return inner(world, mul(world, R2, bra), ket).sum();
	else
		return inner(world, bra, ket).sum();
}
// inner product between functions
double CIS_Operators::make_inner_product(const real_function_3d &bra,
		const real_function_3d &ket) const {
	if (use_nuclear_correlation_factor_)
		return (bra * R2).inner(ket);
	else
		return bra.inner(ket);
}
// inner product between function and vecfunction
double CIS_Operators::make_inner_product(const real_function_3d &bra,
		const vecfuncT &ket) const {
	if (use_nuclear_correlation_factor_)
		return inner(world, bra * R2, ket).sum();
	else
		return inner(world, bra, ket).sum();
}


} /* namespace madness */
