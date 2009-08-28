/// Implementation of functions in frequencyhandler.h

#include "frequencyhandler.h"
#include <misc/interpolation_1d.h>

bool FrequencyHandler::read_file() {
	char filename[80];
	FILE *file;
	int i;

	// open the file
	sprintf(filename, "%s.freq.%d", pulse.basefile.c_str(), freq_index);
	file = fopen(filename, "r");
	if(file == NULL) {
		fprintf(stderr, "Unable to open file %s\n", filename);
		return false;
	}

	// read the parameters:
	// sanity check: the frequency
	double read1;
	if(fscanf(file, "%le", &read1) < 1) {
		fprintf(stderr, "Error reading file\n");
		fclose(file);
		return false;
	}

	if(fabs(read1 - constants::pi*freq_index / pulse.T) > 1.0e-8) {
		fprintf(stderr, "Error: corrupted file\n");
		fclose(file);
		return false;
	}

	// read in the dot-product function
	int npts;
	if(fscanf(file, "%d", &npts) < 1) {
		fprintf(stderr, "Error reading file\n");
		fclose(file);
		return false;
	}
	vector<double> pts(npts);
	for(i = 0; i < npts; ++i) {
		if(fscanf(file, "%le", &pts[i]) < 1) {
			fprintf(stderr, "Error reading file\n");
			fclose(file);
			return false;
		}
	}

	fclose(file);

	CubicInterpolationTable<double> interp(pulse.dotp_min, pulse.dotp_max,
		npts, pts);

	// construct the incident fields
	inc[0] = FunctionFactory<double, 3>(world).functor(
		SharedPtr<FunctionFunctorInterface<double, 3> >(
		new FrequencyIncident(pulse, interp)));
	inc[1] = copy(inc[0]);
	inc[2] = copy(inc[0]);

	for(i = 0; i < 3; ++i) {
		inc[i].scale(pulse.e_polar[i]);
		inc[i].truncate();
	}

	return true;
}

double FrequencyIncident::operator() (const Vector<double, 3> &pt) const {
	double dotp = 0.0;
	for(int i = 0; i < 3; ++i)
		dotp += pulse.prop_dir[i]*pt[i];

	return interp(dotp);
}
