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

// typedefs
typedef std::vector<Function<double, 3> > vecfuncT;

// convenience macros
#define SPLITCOORD(x,y,z,r) double x=r[0]-origin[0]; double y=r[1]-origin[1];double z=r[2]-origin[2];

/// compute the centroid of a function i.e. c[xi]=<f|xi|f>/<f|f> i.e. position expectation value
coord_3d compute_centroid(const real_function_3d& f);
std::vector<coord_3d> compute_centroids(const vecfuncT & vf);
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
    double operator()(const coord_3d& r)const{
    	SPLITCOORD(x,y,z,r);

    	double result=1.0;
    	for(int i=0;i<3;++i){
    		result=result*(*this)(r[i],i);
    	}
    	return result;

    }
    /// operator for the 1D plane waves
    double operator()(const double & x, const int& dim)const{
    	const double argument = (2.0*n[dim])*M_PI/L[dim];
    	if(cosinus[dim]){
    		return cos(argument*x);
    	}else{
    		return sin(argument*x);
    	}
    }
    /// in case this is needed at some point
    double operator()(const coord_1d & x, const int& dim)const{
    	return (*this)(x[0],dim);
    }

    std::string name(const bool& compact=false)const{
    	std::string result="";

    	for(size_t i=0;i<3;++i){
    		if(compact and n[i]>0){
    			if(cosinus[i])	result+="c_" + std::to_string(n[i])+"_";
    			else 			result+="s_" + std::to_string(n[i])+"_";
    		}
    		else if(n[i]>0){
    			if(cosinus[i])	result+="cos(n=" + std::to_string(n[i]) +" l="+ std::to_string(L[i]) + ")";
    			else 			result+="sin(n=" + std::to_string(n[i]) +" l="+ std::to_string(L[i]) + ")";
    		}else result +="0";

    	}

    	return result;
    }

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
static vecfuncT apply_polynomial_exop(vecfuncT& vf, const std::string& exop_input, std::vector<coord_3d> centers=std::vector<coord_3d>(), const bool& fence=false){
	if(vf.empty()) return vf;
	//recompute centers if necessary
	if(centers.empty()) centers=compute_centroids(vf);

	ExopUnaryOpStructure exop(std::make_shared<polynomial_functor>(polynomial_functor(exop_input)));

	for(auto& f:vf){
		f.unaryop(exop,false);
	}
	if(fence) vf.front().world().gop.fence();
	return vf;
}
/// convenience wrapper
static real_function_3d apply_polynomial_exop(real_function_3d & f, const std::string& exop_input, coord_3d center=coord_3d(), const bool& fence=false){
	vecfuncT vf(1,f);
	std::vector<coord_3d> centers(1,center);
	return apply_polynomial_exop(vf,exop_input,centers,fence).front();
}
/// excite a vector of functions with a specific excitation operator
/// @param[in/out] vf the function which gets excited, exop*f on return
/// @param[in] exop_input, the excitation operator defined by a string (see the polynomial_functor class for details)
/// @param[in] the centers of the vf functions, if none were given they are recomputed
/// @return exop*vf i.e. result[i]=exop*vf[i]
static vecfuncT apply_trigonometric_exop(vecfuncT& vf, const std::string& exop_input, std::vector<coord_3d> centers=std::vector<coord_3d>()  , const bool& fence=false){
	if(vf.empty()) return vf;
	//recompute centers if necessary
	if(centers.empty()) centers=compute_centroids(vf);

	ExopUnaryOpStructure exop(std::make_shared<polynomial_trigonometrics_functor>(polynomial_trigonometrics_functor(exop_input)));

	for(auto& f:vf){
		f.unaryop(exop,false);
	}
	if(fence) vf.front().world().gop.fence();
	return vf;
}
/// convenience wrapper
static real_function_3d apply_trigonometric_exop(real_function_3d & f, const std::string& exop_input, coord_3d center=coord_3d(), const bool& fence=false){
	vecfuncT vf(1,f);
	std::vector<coord_3d> centers(1,center);
	return apply_trigonometric_exop(vf,exop_input,centers,fence).front();
}



class GuessFactory {

};

inline std::vector<std::vector<double> > polynomial_functor::read_string(const std::string string) const
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

inline double polynomial_functor::operator ()(const coord_3d& rr) const {
	coord_3d r;
	r[0] = rr[0] - center[0];
	r[1] = rr[1] - center[1];
	r[2] = rr[2] - center[2];
	return dampf(r) * compute_value(r);
}

inline void polynomial_functor::test() {
	std::cout << "Test polynomial functor " << "\n input string is " << input_string_ << std::endl;
	std::cout << "\n read data is \n" << data_ << std::endl;
}

inline double polynomial_functor::compute_value(const coord_3d& r) const {
	double result = 0.0;
	for (size_t i = 0; i < data_.size(); i++) {
		if (data_[i].size() != 4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, data is faulty", 1);
		result += (data_[i][3] * pow(r[0], data_[i][0]) * pow(r[1], data_[i][1]) * pow(r[2], data_[i][2]));
	}

	return result;
}

inline double polynomial_trigonometrics_functor::compute_value(const coord_3d& r) const {
	double result = 0.0;
	for (size_t i = 0; i < data_.size(); i++) {
		if (data_[i].size() != 4) MADNESS_EXCEPTION("ERROR in polynomial exop functor, data is faulty", 1);
		result += (data_[i][3] * pow(sin(r[0]), data_[i][0]) * pow(sin(r[1]), data_[i][1]) * pow(sin(r[2]), data_[i][2]));
	}
}

inline double gauss_functor::operator ()(const coord_3d& rr) const {
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

} /* namespace madness */

#endif /* SRC_APPS_CHEM_GUESSFACTORY_H_ */
