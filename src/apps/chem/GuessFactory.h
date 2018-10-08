/*
 * GuessFactory.h
 *
 *  Created on: Sep 27, 2018
 *      Author: kottmanj
 *
 *  Guess Factory for TDHF which can be used standalone for other purposes
 */

#ifndef SRC_APPS_CHEM_GUESSFACTORY_H_
#define SRC_APPS_CHEM_GUESSFACTORY_H_

#include <madness.h>

namespace madness {

// convenience macros
#define SPLITCOORD(x,y,z,r) double x=r[0]-origin[0]; double y=r[1]-origin[1];double z=r[2]-origin[2];

/// compute the centroid of a function i.e. c[xi]=<f|xi|f>/<f|f> i.e. position expectation value
coord_3d compute_centroid(const real_function_3d& f);
std::vector<coord_3d> compute_centroids(const vector_real_function_3d & vf);
/// little helper for coord (Vector<3>) and Tensor data formats
template<typename T, size_t NDIM>
Vector<T,NDIM> tensor_to_coord(const Tensor<T>& t) {
	Vector<T,NDIM> result;
	MADNESS_ASSERT(t.size() >= NDIM);
	for (size_t i = 0; i < NDIM; ++i)
		result[i] = t[i];
	return result;
}



/// create excitation operators with unaryop (faster as explicit construction and multiplication)
/// Guess function do not need to be perfectly refined
struct ExopUnaryOpStructure {

	ExopUnaryOpStructure(const std::shared_ptr<FunctionFunctorInterface<double, 3> >& f) :
			exfunc(f), cdata(FunctionCommonData<double, 3>::get(FunctionDefaults<3>::get_k())) {
	}

	/// multiplies the target function by the excitation operator defined by exfunc
	void operator ()(const Key<3>& key, Tensor<double>& t) const;
    /// shared pointer to object of excitation operator
    std::shared_ptr<FunctionFunctorInterface<double,3> >  exfunc;
    FunctionCommonData<double,3> cdata;
    template <typename Archive> void serialize(Archive& ar) {}
};

/// creates a plane-wave: sin (or cos) with argument (npi/L*x)
struct PlaneWaveFunctor : public FunctionFunctorInterface<double,3> {

	PlaneWaveFunctor(std::vector<double> vn,std::vector<bool> vc, const coord_3d& c) : L(FunctionDefaults<3>::get_cell_width()), n(vn), cosinus(vc), origin(c) {}

	typedef double resultT;

    /// for explicit construction of this plane-wave function
	double operator ()(const coord_3d& r) const;
    /// operator for the 1D plane waves
	double operator ()(const double& x, const int& dim) const;
    /// in case this is needed at some point
    double operator()(const coord_1d & x, const int& dim)const{
    	return (*this)(x[0],dim);
    }

	std::string name(const bool& compact = false) const;

    const Tensor<double> L;
    const std::vector<double> n;
    const std::vector<bool> cosinus;
    const coord_3d origin;


};

/// GaussFunctor to let the exciation operators go to zero at the boundaries
/// totally symmetric
struct gauss_functor : public FunctionFunctorInterface<double,3> {
public:
	gauss_functor();
	gauss_functor(const double& width): width_(width){
		MADNESS_ASSERT(not(width<0.0));
	}
	gauss_functor(const double& width, const coord_3d c): width_(width), center(c){
		MADNESS_ASSERT(not(width<0.0));
	}
	gauss_functor(const double& width, const Tensor<double> c): width_(width), center(tensor_to_coord<double,3>(c)){
		MADNESS_ASSERT(not(width<0.0));
	}
	const double width_;
	const coord_3d center=coord_3d();

	/// explicit construction
	double operator ()(const coord_3d& rr) const;


};

/// Project a general 3D polynomial to the MRA Grid
struct polynomial_functor : public FunctionFunctorInterface<double,3> {
public :
	polynomial_functor(const std::string input, const double& damp_width=0.0, const coord_3d& c=coord_3d()) : input_string_(input), data_(read_string(input)), dampf(damp_width), center(c) {}
	polynomial_functor(const std::string input,const double& damp_width, const Tensor<double>& c) : input_string_(input), data_(read_string(input)), dampf(damp_width), center(tensor_to_coord<double,3>(c)) {}

	/// construction by coordinates
	double operator ()(const coord_3d& rr) const;

	/// create the value of the polynomial according to the data in the data_ structure
	double compute_value(const coord_3d& r) const;

protected:
	const std::string input_string_;
	/// The data for the construction of the polynomial chain
	/// every entry of data_ is vector containing the threee exponents and the coefficient of a monomial dx^ay^bz^c , data_[i] = (a,b,c,d)
	const std::vector<std::vector<double>> data_;
	/// damping function
	gauss_functor dampf;
	coord_3d center=coord_3d();
public:
	std::vector<std::vector<double> > read_string(const std::string string) const;
	void test();
	std::vector<std::vector<double> > give_data(){return data_;}
};

/// helper struct for computing moments
struct xyz {
	int direction;
	xyz(int direction) : direction(direction) {}
	double operator()(const coord_3d& r) const {
		return r[direction];
	}
};

/// instead of x,y,z use sin(x), sin(y), sin(z)
struct polynomial_trigonometrics_functor : public polynomial_functor {
	/// c++11 constructor inheritance
	using polynomial_functor::polynomial_functor;
	// overload
	/// create the value of the polynomial according to the data in the data_ structure
	/// instead of x,y,z use sin(x), sin(y), sin(z)
	double compute_value(const coord_3d& r) const;
};


/// excite a vector of functions with a specific excitation operator
/// @param[in/out] vf the function which gets excited, exop*f on return
/// @param[in] exop_input , the excitation operator defined by a string (see the polynomial_functor class for details)
/// @return exop*vf i.e. result[i]=exop*vf[i]
vector_real_function_3d apply_polynomial_exop(vector_real_function_3d& vf, const std::string& exop_input, std::vector<coord_3d> centers = std::vector<coord_3d>(), const bool& fence = false);
/// convenience wrapper
real_function_3d apply_polynomial_exop(real_function_3d& f, const std::string& exop_input, coord_3d center = coord_3d(), const bool& fence = false);
/// excite a vector of functions with a specific excitation operator
/// @param[in/out] vf the function which gets excited, exop*f on return
/// @param[in] exop_input, the excitation operator defined by a string (see the polynomial_functor class for details)
/// @param[in] the centers of the vf functions, if none were given they are recomputed
/// @return exop*vf i.e. result[i]=exop*vf[i]
vector_real_function_3d apply_trigonometric_exop(vector_real_function_3d& vf, const std::string& exop_input, std::vector<coord_3d> centers = std::vector<coord_3d>(), const bool& fence = false);
/// convenience wrapper
real_function_3d apply_trigonometric_exop(real_function_3d& f, const std::string& exop_input, coord_3d center = coord_3d(), const bool& fence = false);

/// Makes an excitation operator string based on predefined keywords
std::vector<std::string> make_predefined_exop_strings(const std::string what);
/// Makes an automated excitation operator string for the excitation operators needed to create virtuals from the reference orbitals
std::vector<std::string> make_auto_polynom_strings(const size_t order);



} /* namespace madness */

#endif /* SRC_APPS_CHEM_GUESSFACTORY_H_ */
