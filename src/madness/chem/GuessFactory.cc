/*
 * guessfactory.cc
 *
 *  Created on: Sep 27, 2018
 *      Author: kottmanj
 */

#include "GuessFactory.h"

namespace madness {
namespace guessfactory{

/// compute the centroid of a function i.e. c[xi]=<f|xi|f>/<f|f> i.e. position expectation value
coord_3d compute_centroid(const real_function_3d& f) {
	coord_3d result(0.0);
	real_function_3d density = f * f;
	const double integral = density.trace();
	for (size_t x = 0; x < 3; ++x) {
		const auto mf = PolynomialFunctor(x);
		real_function_3d m = real_factory_3d(f.world()).functor(mf);
		result[x] = (m * density).trace() / integral;
	}
	return result;
}

template<typename T, std::size_t NDIM>
std::vector<coord_3d> compute_centroids(const std::vector<Function<T,NDIM> > & vf){
	std::vector<coord_3d> results(vf.size());
	if(vf.empty()) return results;

	World& world=vf.front().world();
	const std::vector<double> denom=norm2s(world,vf);
	std::vector<Tensor<double> >nums(3);
	for (size_t x = 0; x < 3; ++x) {
		const auto mf = PolynomialFunctor(x);
		real_function_3d m = real_factory_3d(world).functor(mf);
		nums[x]=real(inner(world,vf,m*vf));
	}

	// create correct data structure
	for(size_t i=0;i<vf.size();++i){
		for(size_t x=0;x<3;++x){
			results[i][x]=(nums[x][i]/(denom[i]*denom[i]));
		}
	}

	return results;
}



/// excite a vector of functions with a specific excitation operator
/// @param[in/out] vf the function which gets excited, exop*f on return
/// @param[in] exop_input , the excitation operator defined by a string (see the polynomial_functor class for details)
/// @return exop*vf i.e. result[i]=exop*vf[i]
vector_real_function_3d apply_polynomial_exop(vector_real_function_3d& vf, const std::string& exop_input, std::vector<coord_3d> centers, const bool& fence) {
	if (vf.empty())
		return vf;

	//recompute centers if necessary
	if (centers.empty())
		centers = compute_centroids(vf);

	ExopUnaryOpStructure exop(std::make_shared<PolynomialFunctor>(PolynomialFunctor(exop_input)));
	for (auto& f : vf) {
		f.unaryop(exop, false);
	}
	if (fence)
		vf.front().world().gop.fence();

	return vf;
}

/// convenience wrapper
real_function_3d apply_polynomial_exop(real_function_3d& f, const std::string& exop_input, coord_3d center, const bool& fence) {
	vector_real_function_3d vf(1, f);
	std::vector<coord_3d> centers(1, center);
	return apply_polynomial_exop(vf, exop_input, centers, fence).front();
}

///// excite a vector of functions with a specific excitation operator
///// @param[in/out] vf the function which gets excited, exop*f on return
///// @param[in] exop_input, the excitation operator defined by a string (see the polynomial_functor class for details)
///// @param[in] the centers of the vf functions, if none were given they are recomputed
///// @return exop*vf i.e. result[i]=exop*vf[i]
//template<typename T, std::size_t NDIM>
//std::vector<Function<T,NDIM> > apply_trigonometric_exop(std::vector<Function<T,NDIM> >& vf,
//		const std::string& exop_input, std::vector<coord_3d> centers, const bool& fence) {
//	if (vf.empty())
//		return vf;
//
//	//recompute centers if necessary
//	if (centers.empty())
//		centers = compute_centroids(vf);
//
//	ExopUnaryOpStructure exop(std::make_shared<PolynomialTrigonometricsFunctor>(PolynomialTrigonometricsFunctor(exop_input)));
//	for (auto& f : vf) {
//		f.unaryop(exop, false);
//	}
//	if (fence)
//		vf.front().world().gop.fence();
//
//	return vf;
//}
//
///// convenience wrapper
//template<typename T, std::size_t NDIM>
//Function<T,NDIM> apply_trigonometric_exop(Function<T,NDIM>& f, const std::string& exop_input, coord_3d center, const bool& fence) {
//	std::vector<Function<T,NDIM> > vf(1, f);
//	std::vector<coord_3d> centers(1, center);
//	return apply_trigonometric_exop(vf, exop_input, centers, fence).front();
//}

void PolynomialFunctor::test() {
	std::cout << "Test polynomial functor " << "\n input string is " << input_string_ << std::endl;
	for(const auto& mono : data_){
		for(const auto entry : mono){
			std::cout << entry << ",";
		}
		std::cout << "\n";
	}
}

void ExopUnaryOpStructure::operator ()(const Key<3>& key, Tensor<double>& t) const {
	Tensor<double> exop(t.ndim(), t.dims());
	const Tensor<double>& qp = cdata.quad_x;
	fcube(key, (*exfunc), qp, exop);
	t.emul(exop);
}

void ExopUnaryOpStructure::operator ()(const Key<3>& key, Tensor<double_complex>& t) const {
	Tensor<double> exop(t.ndim(), t.dims());
	const Tensor<double>& qp = cdata.quad_x;
	fcube(key, (*exfunc), qp, exop);
	t.emul(exop);
}

double PolynomialFunctor::operator ()(const coord_3d& rr) const {
	coord_3d r;
	r[0] = rr[0] - center[0];
	r[1] = rr[1] - center[1];
	r[2] = rr[2] - center[2];
	return dampf(r) * compute_value(r);
}



double PolynomialFunctor::compute_value(const coord_3d& r) const {
	double result = 0.0;
	for (size_t i = 0; i < data_.size(); i++) {
		if (data_[i].size() != 4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, data is faulty", 1);
		result += (data_[i][3] * pow(r[0], data_[i][0]) * pow(r[1], data_[i][1]) * pow(r[2], data_[i][2]));
	}

	return result;
}

std::vector<std::vector<double> > PolynomialFunctor::read_string(const std::string string) const
{
	std::stringstream line(string);
	std::string name;
	size_t counter = 0;
	std::vector<double> current_data = vector_factory(0.0, 0.0, 0.0, 1.0);
	std::vector<std::vector<double> > read_data;
	while (line >> name) {
		if (name == "c")
		line >> current_data[3];
		else
		if (name == "x")
		line >> current_data[0];
		else
		if (name == "y")
		line >> current_data[1];
		else
		if (name == "z")
		line >> current_data[2];
		else
		if (name == ",") {
			counter++;
			read_data.push_back(current_data);
			current_data = vector_factory(0.0, 0.0, 0.0, 1.0);
		}
	}
	// dont forget the last read polynomial
	read_data.push_back(current_data);
	return read_data;
}

double PolynomialTrigonometricsFunctor::compute_value(const coord_3d& r) const {
       double result = 0.0;
       for (size_t i = 0; i < data_.size(); i++) {
       	if (data_[i].size() != 4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, data is faulty", 1);
       	result += (data_[i][3] * pow(sin(r[0]), data_[i][0]) * pow(sin(r[1]), data_[i][1]) * pow(sin(r[2]), data_[i][2]));
       }
       throw "CONTROL REACHES END OF NON-VOID FUNCTION JUST GUESSING THAT RESULT IS SUPPOSED TO BE RETURNED!";
       return result; 
}

double GaussFunctor::operator ()(const coord_3d& rr) const {
	coord_3d r;
	r[0] = rr[0] - center[0];
	r[1] = rr[1] - center[1];
	r[2] = rr[2] - center[2];
	if (width_ <= 0.0)
		return 1.0;
	else {
		const double prefactor = 0.06349363593424097 / (width_ * width_ * width_);
		const double exponent = 0.5 / (width_ * width_);
		return prefactor * exp(-exponent * (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
	}
	return 1.0;
}

double PlaneWaveFunctor::operator ()(const coord_3d& r) const {
    //SPLITCOORD(x, y, z, r);
	double result = 1.0;
	for (int i = 0; i < 3; ++i) {
		result = result * (*this)(r[i], i);
	}
	return result;
}

double PlaneWaveFunctor::operator ()(const double& x, const int& dim) const {
	const double argument = (2.0 * n[dim]) * M_PI / L[dim];
	if (cosinus[dim]) {
		return cos(argument * x);
	} else {
		return sin(argument * x);
	}
}

std::string PlaneWaveFunctor::name(const bool& compact) const {
	std::string result = "";
	for (size_t i = 0; i < 3; ++i) {
		if (compact && n[i] > 0) {
			if (cosinus[i])
				result += "c_" + std::to_string(n[i]) + "_";
			else
				result += "s_" + std::to_string(n[i]) + "_";
		} else if (n[i] > 0) {
			if (cosinus[i])
				result += "cos(n=" + std::to_string(n[i]) + " l=" + std::to_string(L[i]) + ")";
			else
				result += "sin(n=" + std::to_string(n[i]) + " l=" + std::to_string(L[i]) + ")";
		} else
			result += "0";
	}
	return result;
}

/// Makes an excitation operator string based on predefined keywords
std::vector<std::string>make_predefined_exop_strings(const std::string what){
	std::vector<std::string> exop_strings;
	if (what == "dipole") {
		exop_strings.resize(3);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
	} else if (what == "x") {
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0";
	} else if (what == "y") {
		exop_strings.resize(1);
		exop_strings[0] = "y 1.0";
	} else if (what == "z") {
		exop_strings.resize(1);
		exop_strings[0] = "z 1.0";
	} else if (what == "r2") {
		exop_strings.resize(1);
		exop_strings[0] = "x 2.0 , y 2.0 , z 2.0";
	} else if (what == "quadrupole") {
		exop_strings.resize(9);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0";
		exop_strings[4] = "y 2.0";
		exop_strings[5] = "z 2.0";
		exop_strings[6] = "x 1.0 y 1.0";
		exop_strings[7] = "x 1.0 z 1.0";
		exop_strings[8] = "y 1.0 z 1.0";
	} else if (what == "dipole+") {
		exop_strings.resize(4);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
	} else if (what == "dipole+diffuse") {
		exop_strings.resize(4);
		exop_strings[0] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[1] = "x 2.0 y 1.0 , y 3.0 , y 1.0 z 2.0";
		exop_strings[2] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[3] = "x 4.0 , y 4.0 , z 4.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	} else if (what == "dipole+diffuse_big") {
		exop_strings.resize(8);
		exop_strings[0] = "x 1.0";
		exop_strings[1] = "y 1.0";
		exop_strings[2] = "z 1.0";
		exop_strings[3] = "x 2.0 , y 2.0 , z 2.0";
		exop_strings[4] = "x 3.0 , x 1.0 y 2.0 , x 1.0 z 2.0";
		exop_strings[5] = "x 2.0 y 1.0 , y 3.0 ,y 1.0 z 2.0";
		exop_strings[6] = "x 2.0 z 1.0 , y 2.0 z 1.0 , z 3.0";
		exop_strings[7] = "x 4.0 , 4 2.0 , 4 2.0, x 2.0 y 2.0, x 2.0 z 2.0, y 2.0 z 2.0";
	} else if (what == "c2v") {
		exop_strings.resize(4);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0 , x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[1] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0 , x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[2] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0 , x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[3] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0, x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0 ";
	} else if (what == "water_first") {
		exop_strings.resize(1);
		exop_strings[0] = "x 1.0 y 1.0, x 1.0 y 1.0 z 1.0";
	} else if (what == "c2v_big") {
		exop_strings.resize(8);
		exop_strings[0] = "z 1.0 , z 3.0 , x 2.0 z 1.0 , y 2.0 z 1.0";
		exop_strings[1] = "x 2.0 , y 2.0 , z 2.0 , x 4.0 , y 4.0 , z 4.0 , x 2.0 y 2.0 , x 2.0 z 2.0 , y 2.0 z 2.0";
		exop_strings[2] = "x 1.0 y 1.0 , x 3.0 y 1.0 , x 1.0 y 3.0";
		exop_strings[3] = "x 1.0 y 1.0 z 1.0 , x 1.0 y 1.0 z 2.0";
		exop_strings[4] = "x 1.0 , x 1.0 z 1.0 , x 1.0 z 2.0 , x 3.0 , x 3.0 z 1.0 , x 1.0 z 3.0";
		exop_strings[5] = "x 1.0 y 2.0 , x 1.0 y 2.0 z 1.0";
		exop_strings[6] = "y 1.0 , y 1.0 z 1.0 , y 1.0 z 2.0 , y 3.0 z 1.0 , y 1.0 z 3.0 , y 3.0";
		exop_strings[7] = "x 2.0 y 1.0 , x 2.0 y 1.0 z 1.0";
	} else if (what == "auto-1") {
		exop_strings = guessfactory::make_auto_polynom_strings(1);
	} else if (what == "auto-2") {
		exop_strings = guessfactory::make_auto_polynom_strings(2);
	} else if (what == "auto-3" or what == "octopole") {
		exop_strings = guessfactory::make_auto_polynom_strings(3);
	} else if (what == "auto-4") {
		exop_strings = guessfactory::make_auto_polynom_strings(4);
	} else if (what == "auto-5") {
		exop_strings = guessfactory::make_auto_polynom_strings(5);
	} else if (what == "auto-6") {
		exop_strings = guessfactory::make_auto_polynom_strings(6);
	} else if (what == "auto-7") {
		exop_strings = guessfactory::make_auto_polynom_strings(7);
	} else if (what == "auto-8") {
		exop_strings = guessfactory::make_auto_polynom_strings(8);
	} else if (what == "auto-9") {
		exop_strings = guessfactory::make_auto_polynom_strings(9);
	} else
		std::cout << "Keyword " << what << " is not known" << std::endl;

	return exop_strings;
}

/// Makes an automated excitation operator string for the excitation operators needed to create virtuals from the reference orbitals
std::vector<std::string> make_auto_polynom_strings(const size_t order){
	std::vector<std::string> exop_strings;
	for (size_t i = 0; i < order + 1; i++) {
		for (size_t j = 0; j < order + 1; j++) {
			for (size_t k = 0; k < order + 1; k++) {
				if (i + j + k > order)
					MADNESS_ASSERT(i + j + k > order); // do nothing
				else if (i == 0 && j == 0 && k == 0)
					MADNESS_ASSERT(i == 0); // do nothing
				else {
					if (i == 0 && j != 0 && k != 0)
						exop_strings.push_back(" y " + madness::stringify(j) + " z " + madness::stringify(k));
					else if (j == 0 && i != 0 && k != 0)
						exop_strings.push_back("x " + madness::stringify(i) + " z " + madness::stringify(k));
					else if (k == 0 && i != 0 && j != 0)
						exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j));
					else if (i == 0 && j == 0)
						exop_strings.push_back(" z " + madness::stringify(k));
					else if (i == 0 && k == 0)
						exop_strings.push_back(" y " + madness::stringify(j));
					else if (j == 0 && k == 0)
						exop_strings.push_back("x " + madness::stringify(i));
					else
						exop_strings.push_back("x " + madness::stringify(i) + " y " + madness::stringify(j) + " z " + madness::stringify(k));
				}
			}
		}
	}
	return exop_strings;
}


// explicit instantiation
template std::vector<coord_3d> compute_centroids(const std::vector<Function<double,3> > & vf);
template std::vector<coord_3d> compute_centroids(const std::vector<Function<double_complex,3> > & vf);



} /* namespace guessfactory */

} /* namespace madness */
