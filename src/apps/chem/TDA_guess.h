/*
 * TDA_guess.h
 *
 *  Created on: Jan 16, 2015
 *      Author: kottmanj
 */



#ifndef TDA_GUESS_H_
#define TDA_GUESS_H_

#include <math.h>

namespace madness{

/// Gaussian analoges of atomic orbitals
struct general_ao : public FunctionFunctorInterface<double,3> {
private:
	const double exponent_;
	const std::string type_;
	double rad2(const coord_3d &r)const{
		return r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
	}
	double s_orbital(const coord_3d &r)const{
		return exp(-rad2(r)*exponent_);
	}
	double px_orbital(const coord_3d &r)const{
		return r[0]*s_orbital(r);
	}
	double py_orbital(const coord_3d &r)const{
		return r[1]*s_orbital(r);
	}
	double pz_orbital(const coord_3d &r)const{
		return r[2]*s_orbital(r);
	}
	double dxx_orbital(const coord_3d &r) const{
		return r[0]*r[0]*s_orbital(r);
	}
	double dyy_orbital(const coord_3d &r) const{
		return r[1]*r[1]*s_orbital(r);
	}
	double dzz_orbital(const coord_3d &r) const{
		return r[2]*r[2]*s_orbital(r);
	}
	double dxy_orbital(const coord_3d &r) const{
		return r[0]*r[1]*s_orbital(r);
	}
	double dxz_orbital(const coord_3d &r) const{
		return r[0]*r[2]*s_orbital(r);
	}
	double dyz_orbital(const coord_3d &r) const{
		return r[1]*r[2]*s_orbital(r);
	}
public:
	general_ao(const double ceta, const std::string typ) : exponent_(ceta),type_(typ) {}

	double operator()(const coord_3d &r)const{
		if(type_=="s") return s_orbital(r);
		else if(type_=="px") return px_orbital(r);
		else if(type_=="py") return py_orbital(r);
		else if(type_=="pz") return pz_orbital(r);
		else if(type_=="dxx") return dxx_orbital(r);
		else if(type_=="dyy") return dyy_orbital(r);
		else if(type_=="dzz") return dzz_orbital(r);
		else if(type_=="dxy") return dxy_orbital(r);
		else if(type_=="dxz") return dxz_orbital(r);
		else if(type_=="dyz") return dyz_orbital(r);
		else MADNESS_EXCEPTION("ERROR in creating atomic orbital guess: Type unknown, use s, px, py or pz" ,1);
	}
};

/// Creates a general polynomial C*x^aY^bz^c with origin at xpos,ypos,zpos which is 1 at the borders (discontinous, but at the dyadic points)
struct general_polynomial : public FunctionFunctorInterface<double,3> {
public:
	general_polynomial(const double &c,const double &x, const double &y, const double &z, const coord_3d &pos, const double & box) :
		coefficient(c), exponent_x(x), exponent_y(y), exponent_z(z), position(pos), box_size_(box) {}
	double operator()(const coord_3d &r)const{
		if(fabs(r[0])>box_size_) return 1.0;
		else if(fabs(r[0])>box_size_) return 1.0;
		else if(fabs(r[0])>box_size_) return 1.0;
		else {
			const coord_3d r2 = r-position;
			return coefficient*pow(r2[0],exponent_x)*pow(r2[1],exponent_y)*pow(r2[2],exponent_z);
		}
	}
private:
	const double coefficient;
	const double exponent_x;
	const double exponent_y;
	const double exponent_z;
	const coord_3d position;
	const double box_size_;
};
class guess{
private :
// size of the simulation box
const double L;
public :
/// Constructor
guess(const double box_size) : L(box_size) {}
/// Make a guess excitation around a given point with a generalized polynomial
vecfuncT make_guess_xfunction(World &world, const std::string input_info, const vecfuncT &mos)const{
	double c=1.0,cx=0.0,cy=0.0,cz=0.0;
	coord_3d position(0.0);
	double box_size = L; // default value -> no "smoothing"
	// Read the input line given in the string input_info
			std::stringstream line(input_info);
			std::string name;
			while(line>>name){
				std::transform(name.begin(), name.end(), name.begin(), ::tolower);
				if(name=="c") line >> c;
				else if(name=="cx") line >> cx;
				else if(name=="cy") line >> cy;
				else if(name=="cz") line >> cz;
				else if(name=="xpos") line >> position[0];
				else if(name=="ypos") line >> position[1];
				else if(name=="zpos") line >> position[2];
				else if(name=="range") line >> box_size;
			}

	// Determine the box size to be at a dyadic point
	double range = L;
	double smallest_step = L/16.0;
	double current = smallest_step;
	for(size_t i=0;i<16;i++){
		current+=smallest_step;
		if(current > box_size){
			range = current;
			break;
		}
	}
	// Generate the general polynomial
	if(world.rank()==0){
		std::cout << "Making excitation operator " << std::fixed << std::setprecision(1) << c << "x^" << cx << "y^" << cy << "z^" << cz << " at position " << position << " with range " << range << std::endl;
	}
	if(cx==0.0 and cy==0.0 and cz==0.0) MADNESS_EXCEPTION("ERROR IN make_guess_xfunction of TDA_guess.h ... no input given",1);
	std::shared_ptr<FunctionFunctorInterface<double, 3> > polynom_functor(new general_polynomial(c,cx,cy,cz,position,range));
	real_function_3d poly_exop = real_factory_3d(world).functor(polynom_functor);
	poly_exop.truncate();

	// Excite the mos
	vecfuncT guess_x;
	for(size_t i=0;i<mos.size();i++){
		real_function_3d tmp = poly_exop * mos[i];
		tmp.truncate();
		guess_x.push_back(tmp);
	}

	return guess_x;
}

/// Creates guess functions like big atomic orbitals (for anions and extreme rydberg states)
vecfuncT make_atomic_guess(World & world,const std::string input_info, const size_t number_of_active_mos_){
	// Read the input line given in the string input_info
			std::stringstream line(input_info);
			std::string name;
			double exponent=1.0;
			std::string type = "s";
			while(line>>name){
				std::transform(name.begin(), name.end(), name.begin(), ::tolower);
				if(name=="exponent") line >> exponent;
				else if(name=="type") line >> type;
			}
			// Generate atomic orbital
			if(world.rank()==0){
				std::cout << "Making " << type<< " type atomic orbital as guess (exponent is " << exponent << ")" << std::endl;
			}
			std::shared_ptr<FunctionFunctorInterface<double, 3> > ao_functor(new general_ao(exponent,type));
			real_function_3d ao_guess = real_factory_3d(world).functor(ao_functor);
			ao_guess.truncate();

			vecfuncT guess_xfunction_x;
			for(size_t i=0;i<number_of_active_mos_;i++){
				guess_xfunction_x.push_back(ao_guess);
			}
			return guess_xfunction_x;
}

};
}// namespace madness



#endif /* TDA_GUESS_H_ */
