#ifndef __frequency_handler__
#define __frequency_handler__

#include <string>
#include <constants.h>
#include <mra/mra.h>
#include "envelopedpulse.h"

using namespace std;
using namespace madness;

typedef double_complex complexd;
typedef Function<complexd, 3> function;
typedef Vector<function, 3> vecfunc;

/// A class for storing all the pertinent data for each frequency-specific
/// field component.
///
/// NOTE that all units are in atomic units, unless specified otherwise.
///
/// This class also extends FunctionFunctorInterface to convert the 1-D
/// incident pulse information into the 3-D incident pulse components.
/// The polarization and propagation directions are stored in the
/// EnvelopedPulse object.
class FrequencyHandler {
public:
	const int freq_index;
	Vector<Function<double, 3>, 3> inc; //< the 3-D incident pulse
	vecfunc scat; //< the 3-D scattered field
	const EnvelopedPulse &pulse; //< The incident pulse / simulation parameters
	World &world; //< The World to use when making functions

	/// Constructor: needs the frequency index, the simulation parameters,
	/// and the MPI communicator (World)
	FrequencyHandler(const int freq, const EnvelopedPulse &_pulse,
		World &_world) : freq_index(freq), pulse(_pulse), world(_world) {}

	/// Read the *.freq.n file to get the parameters, etc.
	/// The specific file name is obtained from pulse.
	/// Returns true on successful read, false otherwise
	bool read_file();
};

/// A class for turning the frequency-dependent incident field function into
/// a 3-D madness function (the dot product is used to convert the 3-D point
/// into the 1-D argument)
class FrequencyIncident : public FunctionFunctorInterface<double, 3> {
protected:
	const EnvelopedPulse &pulse; //< The incident pulse / simulation parameters
	const CubicInterpolationTable<double> &interp; //< The interpolation data

public:
	/// Constructor
	FrequencyIncident(const EnvelopedPulse &_pulse,
		const CubicInterpolationTable<double> &_interp) : pulse(_pulse),
		interp(_interp) {}

	/// The interface for turning the 1-D incident function into the 3-D
	/// vector.  This really works best if the E-field polarization is
	/// only in one of the vector directions.
	double operator() (const Vector<double, 3> &pt) const;
};

#endif
