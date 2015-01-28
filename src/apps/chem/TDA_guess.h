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
	// The world
	World &world;
	// size of the simulation box
	const double L;
	/// Molecular orbitals, the MRA mos or the MRA mos projected on the STO basis
	const std::string mo_mode_;
	bool projected_;
	vecfuncT mos_;
	/// Maximum order of polynomials for excitation operators
	size_t order_;
	const coord_3d zero_position_;
public :
	/// Constructor
	guess(World & world,const double box_size, const std::string mo_mode,const vecfuncT &active_mo, const vecfuncT & ao_basis ) :
		world(world),
		L(box_size),
		mo_mode_(mo_mode),
		projected_(false),
		zero_position_(initv(0.0,0.0,0.0))
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
				if(world.rank()==0 and fabs(C(i,j))>FunctionDefaults<3>::get_thresh()*100.0 and active_mo.size()<6)std::cout << C(i,j) <<  " ";
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

	/// Make a symmetric guess
	std::vector<vecfuncT> make_symmetry_guess(World & world, const std::string point_group)const{
		if(world.rank()==0) std::cout << "Making symmetric guess excitation vectors with point group: " << point_group << std::endl;
		coord_3d position; position[0]=0.0; position[1]=0.0; position[2]=0.0;
		std::vector<vecfuncT> guess_vectors;
		if(point_group=="c2v"){
			guess_vectors = make_c2v_guess(position);
		}else if(point_group=="d2h"){
			guess_vectors = make_d2h_guess(position);
		}else{
			MADNESS_EXCEPTION("ERROR in guess calculation: Unknown point group",1);
		}

		MADNESS_EXCEPTION("ERROR in guess calculation: Should not end up here ... unknown point group ?",1);
	}

	std::vector<vecfuncT> make_c2v_guess(const coord_3d &position, const bool big=false)const{
		std::vector<vecfuncT> guess_vectors;
		{
			// A1 (odd and even)
			real_function_3d a1_odd_exop = real_factory_3d(world);
			real_function_3d a1_eve_exop = real_factory_3d(world);
			std::vector<std::vector<double> > a1_odd_exponents;
			std::vector<std::vector<double> > a1_eve_exponents;
			a1_odd_exponents.push_back(initv(0,0,1));
			a1_odd_exponents.push_back(initv(2,0,1));
			a1_odd_exponents.push_back(initv(0,2,1));
			a1_odd_exponents.push_back(initv(0,0,3));
			a1_odd_exponents.push_back(initv(0,0,1));

			a1_eve_exponents.push_back(initv(2,0,0));
			a1_eve_exponents.push_back(initv(0,2,0));
			a1_eve_exponents.push_back(initv(0,0,2));
			a1_eve_exponents.push_back(initv(4,0,0));
			a1_eve_exponents.push_back(initv(0,4,0));
			a1_eve_exponents.push_back(initv(0,0,4));
			a1_eve_exponents.push_back(initv(2,2,0));
			a1_eve_exponents.push_back(initv(2,0,2));
			a1_eve_exponents.push_back(initv(0,2,2));

			for(size_t i=0;i<a1_odd_exponents.size();i++){
				a1_odd_exop += make_polynomial(world,a1_odd_exponents[i],position,L);
			}
			for(size_t i=0;i<a1_eve_exponents.size();i++){
				a1_eve_exop += make_polynomial(world,a1_eve_exponents[i],position,L);
			}
			vecfuncT a1_odd_guess = mul(world,a1_odd_exop,mos_);
			vecfuncT a1_eve_guess = mul(world,a1_eve_exop,mos_);
			if(big){
				guess_vectors.push_back(a1_odd_guess);
				guess_vectors.push_back(a1_eve_guess);
			}else{
				vecfuncT small_guess = add(world,a1_odd_guess,a1_eve_guess);
				guess_vectors.push_back(small_guess);
			}
		}
		{
			// A2 (with and without z component)
			real_function_3d a2_1_exop = real_factory_3d(world);
			real_function_3d a2_2_exop = real_factory_3d(world);
			std::vector<std::vector<double> > a2_1_exponents;
			std::vector<std::vector<double> > a2_2_exponents;
			a2_1_exponents.push_back(initv(1,1,0));
			a2_1_exponents.push_back(initv(3,1,0));
			a2_1_exponents.push_back(initv(1,3,0));

			a2_2_exponents.push_back(initv(1,1,1));

			for(size_t i=0;i<a2_1_exponents.size();i++){
				a2_1_exop += make_polynomial(world,a2_1_exponents[i],position,L);
			}
			for(size_t i=0;i<a2_2_exponents.size();i++){
				a2_2_exop += make_polynomial(world,a2_2_exponents[i],position,L);
			}
			vecfuncT a2_1_guess = mul(world,a2_1_exop,mos_);
			vecfuncT a2_2_guess = mul(world,a2_1_exop,mos_);
			if(big){
				guess_vectors.push_back(a2_1_guess);
				guess_vectors.push_back(a2_2_guess);
			}else{
				guess_vectors.push_back(add(world,a2_1_guess,a2_2_guess));
			}
		}
		{
			// B1 (with and without y component)
			real_function_3d b1_1_exop = real_factory_3d(world);
			real_function_3d b1_2_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b1_1_exponents;
			std::vector<std::vector<double> > b1_2_exponents;
			b1_1_exponents.push_back(initv(1,0,0));
			b1_1_exponents.push_back(initv(1,0,1));
			b1_1_exponents.push_back(initv(3,0,0));
			b1_2_exponents.push_back(initv(1,2,0));
			b1_1_exponents.push_back(initv(1,0,2));
			b1_1_exponents.push_back(initv(3,0,1));
			b1_2_exponents.push_back(initv(1,2,1));
			b1_1_exponents.push_back(initv(1,0,3));


			for(size_t i=0;i<b1_1_exponents.size();i++){
				b1_1_exop += make_polynomial(world,b1_1_exponents[i],position,L);
			}
			for(size_t i=0;i<b1_2_exponents.size();i++){
				b1_1_exop += make_polynomial(world,b1_2_exponents[i],position,L);
			}
			vecfuncT b1_1_guess = mul(world,b1_1_exop,mos_);
			vecfuncT b1_2_guess = mul(world,b1_1_exop,mos_);
			if(big){
				guess_vectors.push_back(b1_1_guess);
				guess_vectors.push_back(b1_2_guess);
			}else{
				guess_vectors.push_back(add(world,b1_1_guess,b1_2_guess));
			}
		}
		{
			// B2 (with and without x component)
			real_function_3d b2_1_exop = real_factory_3d(world);
			real_function_3d b2_2_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b2_1_exponents;
			std::vector<std::vector<double> > b2_2_exponents;
			b2_1_exponents.push_back(initv(0,1,0));
			b2_1_exponents.push_back(initv(0,1,1));
			b2_1_exponents.push_back(initv(0,1,2));
			b2_2_exponents.push_back(initv(2,1,0));
			b2_1_exponents.push_back(initv(0,3,0));
			b2_2_exponents.push_back(initv(2,1,1));
			b2_1_exponents.push_back(initv(0,3,1));
			b2_1_exponents.push_back(initv(0,1,3));


			for(size_t i=0;i<b2_1_exponents.size();i++){
				b2_1_exop += make_polynomial(world,b2_1_exponents[i],position,L);
			}
			for(size_t i=0;i<b2_2_exponents.size();i++){
				b2_2_exop += make_polynomial(world,b2_2_exponents[i],position,L);
			}
			vecfuncT b2_1_guess = mul(world,b2_1_exop,mos_);
			vecfuncT b2_2_guess = mul(world,b2_1_exop,mos_);
			if(big){
				guess_vectors.push_back(b2_1_guess);
				guess_vectors.push_back(b2_2_guess);
			}else{
				guess_vectors.push_back(add(world,b2_1_guess,b2_2_guess));
			}
		}
		// check
		if(big and not guess_vectors.size()==8) MADNESS_EXEPTION("ERROR in c2v guess, big guess demanded but not 8 guess vectors created ... ???",1);
		if(not big and not guess_vectors.size()==4) MADNESS_EXEPTION("ERROR in c2v guess, normal guess demanded but not 4 guess vectors created ... ???",1);
		return guess_vectors;

	}
	std::vector<vecfuncT> make_d2h_guess(const coord_3d &position, const bool big=false)const{
		std::vector<vecfuncT> guess_vectors;
		{
			// Ag
			real_function_3d ag_exop = real_factory_3d(world);
			std::vector<std::vector<double> > ag_exponents;
			ag_exponents.push_back(initv(2,0,0));
			ag_exponents.push_back(initv(0,2,0));
			ag_exponents.push_back(initv(0,0,2));
			ag_exponents.push_back(initv(4,0,0));
			ag_exponents.push_back(initv(0,4,0));
			ag_exponents.push_back(initv(0,0,4));
			ag_exponents.push_back(initv(2,2,0));
			ag_exponents.push_back(initv(0,2,2));
			ag_exponents.push_back(initv(2,0,2));

			for(size_t i=0;i<ag_exponents.size();i++){
				ag_exop += make_polynomial(world,ag_exponents[i],position,L);
			}
			vecfuncT ag_guess = mul(world,ag_exop,mos_);
			guess_vectors.push_back(ag_guess);

		}
		{
			// B1g
			real_function_3d b1g_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b1g_exponents;
			b1g_exponents.push_back(initv(1,1,0));
			b1g_exponents.push_back(initv(3,1,0));
			b1g_exponents.push_back(initv(1,3,0));
			b1g_exponents.push_back(initv(1,1,2));

			for(size_t i=0;i<b1g_exponents.size();i++){
				b1g_exop += make_polynomial(world,b1g_exponents[i],position,L);
			}
			vecfuncT b1g_guess = mul(world,b1g_exop,mos_);
			guess_vectors.push_back(b1g_guess);
		}
		{
			// b2g
			real_function_3d b2g_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b2g_exponents;
			b2g_exponents.push_back(initv(1,0,1));
			b2g_exponents.push_back(initv(3,0,1));
			b2g_exponents.push_back(initv(1,2,1));
			b2g_exponents.push_back(initv(1,0,3));

			for(size_t i=0;i<b2g_exponents.size();i++){
				b2g_exop += make_polynomial(world,b2g_exponents[i],position,L);
			}
			vecfuncT b2g_guess = mul(world,b2g_exop,mos_);
			guess_vectors.push_back(b2g_guess);
		}
		{
			// b3g
			real_function_3d b3g_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b3g_exponents;
			b3g_exponents.push_back(initv(0,1,1));
			b3g_exponents.push_back(initv(2,1,1));
			b3g_exponents.push_back(initv(0,3,1));
			b3g_exponents.push_back(initv(0,1,3));


			for(size_t i=0;i<b3g_exponents.size();i++){
				b3g_exop += make_polynomial(world,b3g_exponents[i],position,L);
			}
			vecfuncT b3g_guess = mul(world,b3g_exop,mos_);
			guess_vectors.push_back(b3g_guess);
		}
		{
			// au
			real_function_3d au_exop = real_factory_3d(world);
			std::vector<std::vector<double> > au_exponents;
			au_exponents.push_back(initv(1,1,1));

			for(size_t i=0;i<au_exponents.size();i++){
				au_exop += make_polynomial(world,au_exponents[i],position,L);
			}
			vecfuncT au_guess = mul(world,au_exop,mos_);
			guess_vectors.push_back(au_guess);
		}
		{
			// b1u
			real_function_3d b1u_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b1u_exponents;
			b1u_exponents.push_back(initv(0,0,1));
			b1u_exponents.push_back(initv(0,0,3));
			b1u_exponents.push_back(initv(0,2,1));
			b1u_exponents.push_back(initv(2,0,1));

			for(size_t i=0;i<b1u_exponents.size();i++){
				b1u_exop += make_polynomial(world,b1u_exponents[i],position,L);
			}
			vecfuncT b1u_guess = mul(world,b1u_exop,mos_);
			guess_vectors.push_back(b1u_guess);
		}
		{
			// b2u
			real_function_3d b2u_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b2u_exponents;
			b2u_exponents.push_back(initv(0,1,2));
			b2u_exponents.push_back(initv(2,1,0));
			b2u_exponents.push_back(initv(0,3,0));
			b2u_exponents.push_back(initv(0,1,0));

			for(size_t i=0;i<b2u_exponents.size();i++){
				b2u_exop += make_polynomial(world,b2u_exponents[i],position,L);
			}
			vecfuncT b2u_guess = mul(world,b2u_exop,mos_);
			guess_vectors.push_back(b2u_guess);
		}
		{
			// b3u
			real_function_3d b3u_exop = real_factory_3d(world);
			std::vector<std::vector<double> > b3u_exponents;
			b3u_exponents.push_back(initv(1,0,0));
			b3u_exponents.push_back(initv(1,0,2));
			b3u_exponents.push_back(initv(1,2,0));
			b3u_exponents.push_back(initv(3,0,0));

			for(size_t i=0;i<b3u_exponents.size();i++){
				b3u_exop += make_polynomial(world,b3u_exponents[i],position,L);
			}
			vecfuncT b3u_guess = mul(world,b3u_exop,mos_);
			guess_vectors.push_back(b3u_guess);
		}
		//check
		if(not guess_vectors.size()==8) MADNESS_EXCEPTION("ERROR in creating d2h guess, not 8 guess vectors created",1);
		return guess_vectors;
	}
	std::vector<vecfuncT> make_c2h_guess(const coord_3d &position, const bool big=false)const{
		std::vector<vecfuncT> guess_vectors;
		{
			//Ag
			std::vector<std::vector<double> > exponents1;
			std::vector<std::vector<double> > exponents2;
			real_function_3d exop1 = real_factory_3d(world);
			real_function_3d exop2 = real_factory_3d(world);
			exponents1.push_back(initv(2.0,0.0,0.0));
			exponents1.push_back(initv(0.0,2.0,0.0));
			exponents1.push_back(initv(0.0,0.0,2.0));
			exponents1.push_back(initv(4.0,0.0,0.0));
			exponents1.push_back(initv(0.0,4.0,0.0));
			exponents1.push_back(initv(0.0,0.0,4.0));
			exponents1.push_back(initv(2.0,2.0,0.0));
			exponents1.push_back(initv(0.0,2.0,2.0));
			exponents1.push_back(initv(2.0,0.0,2.0));
			exponents2.push_back(initv(1.0,1.0,0.0));
			exponents2.push_back(initv(3.0,1.0,0.0));
			exponents2.push_back(initv(1.0,3.0,0.0));
			exponents2.push_back(initv(1.0,1.0,2.0));

			for(size_t i=0;i<exponents1.size();i++){
				exop1 += make_polynomial(world,exponents1[i],position,L);
			}
			for(size_t i=0;i<exponents2.size();i++){
				exop2 += make_polynomial(world,exponents2[i],position,L);
			}
			vecfuncT guess1 = mul(world,exop1,mos_);
			vecfuncT guess2 = mul(world,exop2,mos_);
			if(big){
				guess_vectors.push_back(guess1);
				guess_vectors.push_back(guess2);
			}else{
				guess_vectors.push_back(add(world,guess1,guess2));
			}
		}
		{
			//Au
			std::vector<std::vector<double> > exponents1;
			std::vector<std::vector<double> > exponents2;
			real_function_3d exop1 = real_factory_3d(world);
			real_function_3d exop2 = real_factory_3d(world);
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));

			for(size_t i=0;i<exponents1.size();i++){
				exop1 += make_polynomial(world,exponents1[i],position,L);
			}
			for(size_t i=0;i<exponents2.size();i++){
				exop2 += make_polynomial(world,exponents2[i],position,L);
			}
			vecfuncT guess1 = mul(world,exop1,mos_);
			vecfuncT guess2 = mul(world,exop2,mos_);
			if(big){
				guess_vectors.push_back(guess1);
				guess_vectors.push_back(guess2);
			}else{
				guess_vectors.push_back(add(world,guess1,guess2));
			}
		}
		{
			//B1
			std::vector<std::vector<double> > exponents1;
			std::vector<std::vector<double> > exponents2;
			real_function_3d exop1 = real_factory_3d(world);
			real_function_3d exop2 = real_factory_3d(world);
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents1.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));
			exponents2.push_back(initv(0.0,0.0,0.0));

			for(size_t i=0;i<exponents1.size();i++){
				exop1 += make_polynomial(world,exponents1[i],position,L);
			}
			for(size_t i=0;i<exponents2.size();i++){
				exop2 += make_polynomial(world,exponents2[i],position,L);
			}
			vecfuncT guess1 = mul(world,exop1,mos_);
			vecfuncT guess2 = mul(world,exop2,mos_);
			if(big){
				guess_vectors.push_back(guess1);
				guess_vectors.push_back(guess2);
			}else{
				guess_vectors.push_back(add(world,guess1,guess2));
			}
		}
		//check
		if(big and not guess_vectors.size()==8) MADNESS_EXEPTION("ERROR in c2h guess, big guess demanded but not 8 guess vectors created ... ???",1);
		if(not big and not guess_vectors.size()==4) MADNESS_EXEPTION("ERROR in c2h guess, normal guess demanded but not 4 guess vectors created ... ???",1);
		return guess_vectors;
	}
	/// Needed by symmetry guess function
	std::vector<double> initv(const double &x,const double &y, const double &z)const{
		std::vector<double> tmp;
		tmp.push_back(x); tmp.push_back(y); tmp.push_back(z);
		return tmp;
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
