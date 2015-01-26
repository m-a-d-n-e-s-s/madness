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
	/// Molecular orbitals, the MRA mos or the MRA mos projected on the STO basis
	const std::string mo_mode_;
	bool projected_;
	vecfuncT mos_;
	/// Maximum order of polynomials for excitation operators
	size_t order_;
public :
	/// Constructor
	guess(World & world,const double box_size, const std::string mo_mode,const vecfuncT &active_mo, const vecfuncT & ao_basis ) :
		L(box_size),
		mo_mode_(mo_mode),
		projected_(false)
	{
		if(world.rank()==0) std::cout << "Making guess excitation operators, guess_mode is " << mo_mode_ << std::endl;
		if(mo_mode_=="projected" or mo_mode_=="projected_mo" or mo_mode_=="projected_mos" or mo_mode_=="ao"){
			if(world.rank()==0) std::cout << "Molecular orbitals are getting projected to AO-Basis ..." << std::endl;
			projected_ = true;
			mos_ = project_to_ao_basis(world,active_mo,ao_basis);
		}else mos_ = active_mo;
		order_ =6;
		make_polynomial_exops_key(world);
		if(world.rank()==0) std::cout << "\n\n" << std::endl;
	}

	/// The key and the exponents for custom guess exops
	/// assignment goes as follows : if(input_line == key[i]) exponents for exop are in exop_exopnents[i]
	/// construct with make_polynomial_exops_key
	std::vector<std::string> exop_key_;
	std::vector<std::vector<double> > exop_exponents_;




	/// Initialize a polynomial
	real_function_3d make_polynomial(World & world,const std::vector<double> exponents, const coord_3d &position, const double range)const{
		std::shared_ptr<FunctionFunctorInterface<double, 3> > polynom_functor(new general_polynomial(1.0,exponents[0],exponents[1],exponents[2],position,range));
		real_function_3d poly_exop = real_factory_3d(world).functor(polynom_functor);
		poly_exop.truncate();
		return poly_exop;
	}

	/// Construct the polynomial exops key
	void make_polynomial_exops_key(World & world){
		std::vector<std::string> key;
		vecfuncT polynomials;
		std::vector<std::vector<double>> exponents;
		for(size_t x=0;x<order_;x++){
			for(size_t y=0;y<order_;y++){
				for(size_t z=0;z<order_;z++){
					if(x+y+z <order_){
					double exx = (double) x; double exy = (double) y ; double exz = (double) z;
					std::string key_tmp =  "x" + stringify(x) + "y" + stringify(y) + "z" + stringify(z);
					std::vector<double> ex_tmp; ex_tmp.push_back(exx); ex_tmp.push_back(exy); ex_tmp.push_back(exz);
					exponents.push_back(ex_tmp);
					key.push_back(key_tmp);
					}
				}
			}
		}
		exop_key_ = key;
		exop_exponents_ = exponents;
	}

	/// Make custom exoperator from a given input line
	/// input format must be xaybzc C , with a,b,c,C free doubles
	/// x as excitationoperator would be x1y0z0 1.0
	real_function_3d make_custom_excitation_operator(World & world, const std::string input_info,bool silent = false)const{
		// Read the input lines and assign the key
		coord_3d position(0.0);
		double range = L;
		std::vector<size_t> key_numbers_;
		std::vector<double> coefficients;
		std::stringstream line(input_info);
		std::string name;
		while(line>>name){
			std::transform(name.begin(), name.end(), name.begin(), ::tolower);
			for(size_t i=0;i<exop_key_.size();i++){
				if(name==exop_key_[i]){
					key_numbers_.push_back(i);
					double c;
					line>>c;
					coefficients.push_back(c);
				}
				else if(name == "pos" or name == "position"){
					line >> position[0];
					line >> position[1];
					line >> position[2];
				}
				else if(name=="range" or name=="box_size"){
					line>>range;
				}
			}
		}
		// failsafe
		if(coefficients.size()!=key_numbers_.size()) MADNESS_EXCEPTION("ERROR in make_custom_excitation_operator, coefficients and polynomials not the same size",1);
		if(coefficients.size()==0) MADNESS_EXCEPTION("ERROR in make_custom_excitation_operator, no input or wrong format: Format is xaybzc d with doubles a,b,c,d",1);
		// Control output
		if(world.rank()==0 and not silent) std::cout << "Make Custom Operator at position "
				<< position << " with range " << range << " and polynomials: \n"<< std::endl;
		// Make the exop
		real_function_3d exop = real_factory_3d(world);
		for(size_t i=0;i<key_numbers_.size();i++){
			real_function_3d tmp = make_polynomial(world,exop_exponents_[key_numbers_[i]],position,range);
			exop += coefficients[i]*tmp;
			if(world.rank()==0 and not silent) std::cout << coefficients[i] << " "<< exop_key_[key_numbers_[i]] << " ";
		}
		if(world.rank()==0 and not silent) std::cout << "\n-----\n" << std::endl;
		exop.truncate();
		return exop;
	}

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
		double norm = ao_guess.norm2();
		ao_guess.scale(1.0/norm);
		double xnorm = (double) number_of_active_mos_;
		ao_guess.scale(1.0/xnorm);
		vecfuncT guess_xfunction_x;
		for(size_t i=0;i<number_of_active_mos_;i++){
			guess_xfunction_x.push_back(ao_guess);
		}
		return guess_xfunction_x;
	}

	/// Project to AO basis
	vecfuncT project_to_ao_basis(World & world, const vecfuncT & active_mo, const vecfuncT & ao_basis)const{
		//if(world.rank()==0)std::cout << "MO size " << active_mo.size() << " AO size " << ao_basis.size() << std::endl;
		// Make AOAO overlap Matrix
		Tensor<double> Saa = matrix_inner(world,ao_basis,ao_basis,true);
		// Make AOMO overlap Matrix
		Tensor<double> Sam = matrix_inner(world,ao_basis,active_mo,false);
		// Get the coefficients : solve Sam = C*Sao
		Tensor<double> C;
		gesv(Saa, Sam, C);
		C = transpose(C);
		//if(world.rank()==0)std::cout << "AOAO overlap Matrix:\n" << std::setprecision(2) << Saa << std::endl;
		//if(world.rank()==0)std::cout << "AOMO overlap Matrix:\n" << std::setprecision(2) << Sam << std::endl;
		//if(world.rank()==0)std::cout << "Transformation Matrix:\n" << std::setprecision(2) << C << std::endl;

		// make projected mos
		vecfuncT projected_mos;
		for(size_t i=0;i<active_mo.size();i++){
			//if(world.rank()==0)std::cout << "Projecting MO " << i << " to AO Basis ...\n Coefficients: ";
			real_function_3d tmp = real_factory_3d(world);
			for(size_t j=0;j<ao_basis.size();j++){
				tmp += C(i,j)*ao_basis[j];
				if(world.rank()==0 and fabs(C(i,j))>FunctionDefaults<3>::get_thresh()*100.0)std::cout << C(i,j) <<  " ";
			}
			if(world.rank()==0)std::cout << std::endl;
			tmp.truncate();
			projected_mos.push_back(tmp);
			plot_plane(world,tmp,"projected_mo_"+stringify(i));
		}
		return projected_mos;
	}

	/// Make a custom guess excitation vector
	vecfuncT make_custom_guess(World & world, const std::string input_info, bool silent=false)const{
		if(mos_.empty())MADNESS_EXCEPTION("Failure in projection of MOs, no MOs given, maybe too much frozen ?",1);
		// Get the excitation operators
		real_function_3d exop = make_custom_excitation_operator(world,input_info,silent);
		// excite the projected mos
		vecfuncT guess_x;
		for(size_t i=0;i<mos_.size();i++){
			real_function_3d xtmp = exop*mos_[i];
			xtmp.truncate();
			double norm = xtmp.norm2();
			xtmp.scale(1.0/norm);
			guess_x.push_back(xtmp);
		}
		// normalize
		Tensor<double> x_self_overlap = inner(world,guess_x,guess_x);
		double xnorm2 = x_self_overlap.sum();
		scale(world,guess_x,1.0/sqrt(xnorm2));
		return guess_x;
	}

	/// Project the solution vectors back to the polynomial basis
	void project_to_guess_basis(World &world, const vecfuncT &solution,const double tol)const{
		if(world.rank()==0) std::cout << "Projecting back the solution vectors to the polynomial guess basis up to order " << order_-1 << std::endl;
		if(world.rank()==0) std::cout << "Projection of MOs to AO basis is " << projected_ << std::endl;
		std::vector<std::string> dominant_c;
		for(size_t i=0;i<exop_exponents_.size();i++){
			// make the corresponding guess_excitation_vector to every key element
			std::string input = exop_key_[i] + " 1.0";
			vecfuncT guess_x = make_custom_guess(world,input,true);
			// calculate overlap with the solution vector
			Tensor<double> ctensor = inner(world,solution,guess_x);
			double c = ctensor.sum();
			std::string tmp = exop_key_[i] + " " + stringify(c);
			dominant_c.push_back(tmp);
			if(world.rank()==0 and fabs(c)>tol) std::cout << std::setprecision(2) << exop_key_[i] << " " << c << " ";
		}

	}
	/// Project chosen components the solution vectors back to the polynomial basis
	void project_to_guess_basis(World &world, const real_function_3d &solution,const double tol,const size_t mo_number)const{
		if(world.rank()==0) std::cout << " Projection corresponds to MO " << mo_number << std::endl;
		std::vector<std::string> dominant_c;
		for(size_t i=0;i<exop_exponents_.size();i++){
			// make the corresponding guess_excitation_vector to every key element
			std::string input = exop_key_[i] + " 1.0";
			vecfuncT guess_x = make_custom_guess(world,input,true);
			// calculate overlap with the solution vector
			double c = inner(solution,guess_x[mo_number]);
			std::string tmp = exop_key_[i] + " " + stringify(c);
			dominant_c.push_back(tmp);
			if(world.rank()==0 and fabs(c)>tol) std::cout << std::setprecision(2) << exop_key_[i] << " " << c << " ";
		}

	}
};
}// namespace madness



#endif /* TDA_GUESS_H_ */
