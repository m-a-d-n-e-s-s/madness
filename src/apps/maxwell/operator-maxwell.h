#ifndef __operator_maxwell__
#define __operator_maxwell__

#include <string>
#include <constants.h>
#include <mra/mra.h>
#include <linalg/gmres.h>

using namespace std;
using namespace madness;

typedef double_complex complexd;
typedef Function<complexd, 3> function;
typedef Vector<function, 3> vecfunc;

/// A class for storing all the pertinent data for an enveloped incident
/// pulse.
///
/// NOTE that all units are in atomic units, unless specified otherwise.
class EFieldOperator : public Operator<vecfunc> {
protected:
	const function &epshat;
	const complexd prefact;
	const vecfunc &gradlnepshat;
	const Function<double, 3> &box_mask; // mask to make the boundary 0
	const Function<double, 3> &grad_mask; // mask to damp out derivatives
	const SeparatedConvolution<double, 3> &G;

	void action(const vecfunc &in, vecfunc &out) const {
		function dotp = (gradlnepshat[0]*in[0] + gradlnepshat[1]*in[1]
			+ gradlnepshat[2]*in[2]) * box_mask;
		vecfunc preop;
		int i;

		for(i = 0; i < 3; ++i)
			preop[i] = diff(dotp, i) * grad_mask;
		dotp.clear();

		for(i = 0; i < 3; ++i) {
			dotp = epshat*in[i];
			dotp.compress();
			preop[i].compress();
			preop[i].gaxpy(complexd(1.0, 0.0), dotp, prefact);
			out[i] = in[i] - apply(G, preop[i]);
			out[i].truncate();
			dotp.clear();
			preop[i].clear();
		}
	}

public:
	/// Needs the epsilon-hat complex dielectric, the frequency omega,
	/// mu0, eps0, the gradient of ln(epshat), the box_mask function to make
	/// the function 0 at the boundary, the grad_mask function to smooth
	/// out derivatives near the boundaries, and the Poisson Green's function
	EFieldOperator(const function &_epshat, const double _omega,
		const double _mu0, const double _eps0, const vecfunc &_gradlnepshat,
		const Function<double, 3> &_box_mask,
		const Function<double, 3> &_grad_mask,
		const SeparatedConvolution<double, 3> &_G) : epshat(_epshat),
		prefact(complexd(0.0, -_omega*_eps0*_mu0)), gradlnepshat(_gradlnepshat),
		box_mask(_box_mask), grad_mask(_grad_mask), G(_G) {}
};

#endif
