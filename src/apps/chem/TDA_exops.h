/*
 * TDA_exops.h
 *
 *  Created on: Sep 2, 2014
 *      Author: kottmanj
 */

#ifndef TDA_EXOPS_H_
#define TDA_EXOPS_H_

namespace madness{
struct exoperators{
private:
	/// The point of excitation for the custom guess
	coord_3d excitation_point_;
	vecfuncT polynom_basis_;
public:
	exoperators(World &world,const coord_3d excitation_point) : excitation_point_(excitation_point){}
	std::vector<std::string> key_;
	// Creates the polynom basis
	void get_polynom_basis(World &world){
		if(polynom_basis_.empty()){
		if(world.rank()==0) std::cout << "Creating the polynomial guess basis ... " << std::endl;
		vecfuncT creators;
		std::vector<std::string> key_creator;
		std::vector<std::string> key;
		 real_function_3d x   = real_factory_3d(world).f(x_function); x.add_scalar(-excitation_point_[0]); creators.push_back(x);
		 real_function_3d y   = real_factory_3d(world).f(y_function); y.add_scalar(-excitation_point_[1]); creators.push_back(y);
		 real_function_3d z   = real_factory_3d(world).f(z_function); z.add_scalar(-excitation_point_[2]); creators.push_back(z);
		 std::string xstring = "x"; key_creator.push_back(xstring);
		 std::string ystring = "y"; key_creator.push_back(ystring);
		 std::string zstring = "z"; key_creator.push_back(zstring);

		vecfuncT basis;
		//dipole
		for(size_t i=0;i<creators.size();i++){
			basis.push_back(creators[i]);
			key.push_back(key_creator[i]);
		}
		std::cout << "Dipole Basis size " << key.size();
		//quadrupole
		for(size_t i=0;i<creators.size();i++){
			for(size_t j=0;j<=i;j++){
				real_function_3d tmp = creators[i]*creators[j];
				basis.push_back(tmp);
				std::string stmp = key_creator[i]+key_creator[j];
				key.push_back(stmp);
			}
		}
		std::cout << "Quadrupole Basis size " << key.size();
		// cubics
		//quadrupole
		for(size_t i=0;i<creators.size();i++){
			for(size_t j=0;j<=i;j++){
				for(size_t k=0;k<=j;k++){
					real_function_3d tmp = creators[i]*creators[j]*creators[k];
					basis.push_back(tmp);
					std::string stmp = key_creator[i]+key_creator[j]+key_creator[k];
					key.push_back(stmp);
				}
			}
		}
		std::cout << "cubics Basis size " << key.size();
		// quartics
		//quadrupole
		for(size_t i=0;i<creators.size();i++){
			for(size_t j=0;j<=i;j++){
				for(size_t k=0;k<=j;k++){
					for(size_t l=0;l<=k;l++){
						real_function_3d tmp = creators[i]*creators[j]*creators[k]*creators[l];
						basis.push_back(tmp);
						std::string stmp = key_creator[i]+key_creator[j]+key_creator[k]+key_creator[l];
						key.push_back(stmp);
					}
				}
			}
		}
		std::cout << "quartics Basis size " << key.size();
		//quints
		// quartics
		//quadrupole
		for(size_t i=0;i<creators.size();i++){
			for(size_t j=0;j<=i;j++){
				for(size_t k=0;k<=j;k++){
					for(size_t l=0;l<=k;l++){
						for(size_t m=0;m<=l;m++){
							real_function_3d tmp = creators[i]*creators[j]*creators[k]*creators[l]*creators[m];
							basis.push_back(tmp);
							std::string stmp = key_creator[i]+key_creator[j]+key_creator[k]+key_creator[l]+key_creator[m];
							key.push_back(stmp);
					}
				}
			}
		}
	}
	std::cout << "quints Basis size " << key.size();
		for(size_t i=0;i<creators.size();i++){
			for(size_t j=0;j<=i;j++){
				for(size_t k=0;k<=j;k++){
					for(size_t l=0;l<=k;l++){
						for(size_t m=0;m<=l;m++){
							for(size_t n=0;n<=m;n++){
							real_function_3d tmp = creators[i]*creators[j]*creators[k]*creators[l]*creators[m]*creators[n];
							basis.push_back(tmp);
							std::string stmp = key_creator[i]+key_creator[j]+key_creator[k]+key_creator[l]+key_creator[m]+key_creator[n];
							key.push_back(stmp);
					}
				}
			}
		}
	}}
//	std::cout << "6th order Basis size " << key.size();
//		for(size_t i=0;i<creators.size();i++){
//			for(size_t j=0;j<=i;j++){
//				for(size_t k=0;k<=j;k++){
//					for(size_t l=0;l<=k;l++){
//						for(size_t m=0;m<=l;m++){
//							for(size_t n=0;n<=m;n++){
//								for(size_t u=0;u<=n;u++){
//							real_function_3d tmp = creators[j]*creators[i]*creators[k]*creators[l]*creators[m]*creators[n]*creators[u];
//							basis.push_back(tmp);
//							std::string stmp = key_creator[j]+key_creator[i]+key_creator[k]+key_creator[l]+key_creator[m]+key_creator[n]+key_creator[u];
//							key.push_back(stmp);}}}}}}}
//		for(size_t i=0;i<creators.size();i++){
//			for(size_t j=0;j<=i;j++){
//				for(size_t k=0;k<=j;k++){
//					for(size_t l=0;l<=k;l++){
//						for(size_t m=0;m<=l;m++){
//							for(size_t n=0;n<=m;n++){
//								for(size_t u=0;u<=n;u++){
//									for(size_t q=0;q<=u;q++){
//							real_function_3d tmp = creators[j]*creators[i]*creators[k]*creators[l]*creators[m]*creators[n]*creators[u]*creators[q];
//							basis.push_back(tmp);
//							std::string stmp = key_creator[j]+key_creator[i]+key_creator[k]+key_creator[l]+key_creator[m]+key_creator[n]+key_creator[u]+key_creator[q];
//							key.push_back(stmp);}}}}}}}}
//		std::cout << "7th order Basis size " << key.size();
//		for(size_t i=0;i<creators.size();i++){
//			for(size_t j=0;j<=i;j++){
//				for(size_t k=0;k<=j;k++){
//					for(size_t l=0;l<=k;l++){
//						for(size_t m=0;m<=l;m++){
//							for(size_t n=0;n<=m;n++){
//								for(size_t u=0;u<=n;u++){
//									for(size_t q=0;q<=u;q++){
//										for(size_t w=0;w<=q;w++){
//							real_function_3d tmp = creators[j]*creators[i]*creators[k]
//							        *creators[l]*creators[m]*creators[n]*creators[u]*creators[q]*creators[w];
//							basis.push_back(tmp);
//							std::string stmp = key_creator[j]+key_creator[i]+key_creator[k]+key_creator[l]
//							               +key_creator[m]+key_creator[n]+key_creator[u]+key_creator[q]+key_creator[w];
//							key.push_back(stmp);}}}}}}}}}
	polynom_basis_ = basis;
	key_=key;
		}
	}
	void plot_polynom_basis(World &world){
		// check if polynom_basis_ is created, if not do it
		get_polynom_basis(world);

		for(size_t i=0;i<polynom_basis_.size();i++){
			plot_plane(world,polynom_basis_[i],"polynom_basis_"+stringify(i));
		}
	}

vecfuncT atomic_excitation_dipole(World& world){
	vecfuncT exops;
	 real_function_3d x   = real_factory_3d(world).f(x_function); x.add_scalar(-excitation_point_[0]); exops.push_back(x);
	 real_function_3d y   = real_factory_3d(world).f(y_function); y.add_scalar(-excitation_point_[1]); exops.push_back(y);
	 real_function_3d z   = real_factory_3d(world).f(z_function); z.add_scalar(-excitation_point_[2]); exops.push_back(z);
	 real_function_3d rr  = x*x + y*y +z*z; exops.push_back(rr);
	 truncate(world,exops);
	 return exops;
}
	std::vector<double> get_overlaps_with_guess(World &world,vecfuncT &excitation,const vecfuncT &mo,const real_function_3d smoothing){
		if(excitation.size()!=mo.size()) MADNESS_EXCEPTION("Error in calculating overlap with the guess operators, excitation vector and mo vector have different sizes",1);
		get_polynom_basis(world);
		vecfuncT quadbas = polynom_basis_;
		 std::vector<double> overlaps;
		 for(size_t i=0;i< quadbas.size();i++){
			 vecfuncT guess_ex = mul(world,quadbas[i],mo);
			 guess_ex = mul(world,smoothing,guess_ex);
			 double norm2= inner(world,guess_ex,guess_ex).sum();
			 scale(world,guess_ex,1.0/sqrt(norm2));
			 double norm2_2 = inner(world,excitation,excitation).sum();
			 scale(world,excitation,1.0/sqrt(norm2_2));
			 double overlap_tmp = inner(world, guess_ex,excitation).sum();
			 if(overlap_tmp < 1.e-5) overlap_tmp=0.0;
			 overlaps.push_back(overlap_tmp);
		 }
		 return overlaps;
	}

	std::vector<std::vector<double> > get_coefficients(World &world, const std::vector< std::string > input_lines){
		if(world.rank()==0)std::cout << "Creating custom guess operators ... " << input_lines.size() << " demanded" << std::endl;
		if(input_lines.size()==0) MADNESS_EXCEPTION("Custom guess demanded but no input was given",1);
		//std::cout << input_lines << std::endl;

		// get the coefficients from the strings
		std::vector<std::vector<double> > coefficients;
		for(size_t i=0;i<input_lines.size();i++){
		std::stringstream line(input_lines[i]);
		std::string name;
		if(key_.empty()) get_polynom_basis(world);
		std::vector<double> icoeff(key_.size(),0);

		while(line>>name){
			std::transform(name.begin(), name.end(), name.begin(), ::tolower);
			for(size_t i=0;i<key_.size();i++){
				if(name==key_[i]) line>>icoeff[i];
			}
		}
		coefficients.push_back(icoeff);
		}

		// controll output
		if(world.rank()==0){
			std::cout << "\n\n---------CUSTOM OPERATORS ARE----------\n\n" << std::endl;
			for(size_t i=0;i<coefficients.size();i++){
				std::cout << "\n exop" << i << " (x,y,z,xx,yy,zz,xy,xz,yz)" << std::endl;
				for(size_t j=0;j<coefficients[i].size();j++){
					if(coefficients[i][j]!=0.0){
						std::cout << std::fixed <<std::setprecision(3)
						<<" " <<  key_[j]<< " " <<coefficients[i][j];
					}
				}
				std::cout << std::endl;
			}
			std::cout << "\n\n---------------------------------------\n\n" << std::endl;
		}
		return coefficients;
	}
	vecfuncT get_custom_exops(World &world,const std::vector< std::string > input_lines){
		std::vector<std::vector<double> > coefficients = get_coefficients(world, input_lines);
		return make_custom_exops(world,coefficients);
	}
	vecfuncT make_custom_exops(World &world,std::vector<std::vector<double> > coefficients){
		get_polynom_basis(world);
		vecfuncT exops;
		vecfuncT quadbas=polynom_basis_;

		 for(size_t i=0;i<coefficients.size();i++){
			 real_function_3d tmp = real_factory_3d(world);
			 for(size_t j=0;j<coefficients[i].size();j++){
				  tmp += coefficients[i][j]*quadbas[j];
			 }
			 exops.push_back(tmp);
		 }
		 return exops;
	}

	std::vector<vecfuncT> make_custom_guess(World &world,const std::vector< std::string > input_lines,const vecfuncT &mos, const real_function_3d &smoothing){
		get_polynom_basis(world);
		std::vector<std::vector<double> > coefficients = get_coefficients(world, input_lines);
		std::vector<vecfuncT> guess;
		vecfuncT quadbas = polynom_basis_;

		// make guess
		for(size_t i=0;i<coefficients.size();i++){
			int n = mos.size();
			vecfuncT xtmp = zero_functions<double,3>(world,n);
			for(size_t j=0;j<coefficients[i].size();j++){
				if(coefficients[i][j]!=0){
					real_function_3d smoothed_exop = quadbas[j]*smoothing;
					vecfuncT scaled_basis_function = mul(world,smoothed_exop,mos);
					double norm = inner(world,scaled_basis_function,scaled_basis_function).sum();
					norm = sqrt(norm);
					scale(world,scaled_basis_function,coefficients[i][j]/norm);
					xtmp = add(world,xtmp,scaled_basis_function);
				}
			}
			guess.push_back(xtmp);
		}
		return guess;
	}

	vecfuncT get_exops(World &world,const std::string exop){

		vecfuncT xoperators;
		if(exop == "dipole"){
			real_function_3d fx = real_factory_3d(world).f(x_function);
			xoperators.push_back(fx);
			real_function_3d fy = real_factory_3d(world).f(y_function);
			xoperators.push_back(fy);
			real_function_3d fz = real_factory_3d(world).f(z_function);
			xoperators.push_back(fz);
		}
		else if(exop == "dipole+"){
			real_function_3d fr = real_factory_3d(world).f(r_function);
			xoperators.push_back(fr);
			real_function_3d fx = real_factory_3d(world).f(x_function);
			xoperators.push_back(fx);
			real_function_3d fy = real_factory_3d(world).f(y_function);
			xoperators.push_back(fy);
			real_function_3d fz = real_factory_3d(world).f(z_function);
			xoperators.push_back(fz);
		}
		else if(exop == "quadrupole"){
			real_function_3d fx = real_factory_3d(world).f(x_function);
			xoperators.push_back(fx);
			real_function_3d fy = real_factory_3d(world).f(y_function);
			xoperators.push_back(fy);
			real_function_3d fz = real_factory_3d(world).f(z_function);
			xoperators.push_back(fz);
			real_function_3d fxy = real_factory_3d(world).f(xy_function);
			xoperators.push_back(fxy);
			real_function_3d fyz = real_factory_3d(world).f(yz_function);
			xoperators.push_back(fyz);
			real_function_3d fxz = real_factory_3d(world).f(xz_function);
			xoperators.push_back(fxz);
		}
		else if(exop == "quadrupole+"){
			real_function_3d fx = real_factory_3d(world).f(x_function);
			xoperators.push_back(fx);
			real_function_3d fy = real_factory_3d(world).f(y_function);
			xoperators.push_back(fy);
			real_function_3d fz = real_factory_3d(world).f(z_function);
			xoperators.push_back(fz);
			real_function_3d fxx = real_factory_3d(world).f(xx_function);
			xoperators.push_back(fxx);
			real_function_3d fyy = real_factory_3d(world).f(yy_function);
			xoperators.push_back(fyy);
			real_function_3d fzz = real_factory_3d(world).f(zz_function);
			xoperators.push_back(fzz);
			real_function_3d fxy = real_factory_3d(world).f(xy_function);
			xoperators.push_back(fxy);
			real_function_3d fyz = real_factory_3d(world).f(yz_function);
			xoperators.push_back(fyz);
			real_function_3d fxz = real_factory_3d(world).f(xz_function);
			xoperators.push_back(fxz);
		}
		else if(exop == "hybrids"){
			real_function_3d tmp_sp30=real_factory_3d(world).f(sp30);xoperators.push_back(tmp_sp30);
			real_function_3d tmp_sp31=real_factory_3d(world).f(sp31);xoperators.push_back(tmp_sp31);
			real_function_3d tmp_sp32=real_factory_3d(world).f(sp32);xoperators.push_back(tmp_sp32);
			real_function_3d tmp_sp33=real_factory_3d(world).f(sp33);xoperators.push_back(tmp_sp33);
		}
		else if(exop == "C2v"){
			real_function_3d A1=real_factory_3d(world).f(c2v_A1);xoperators.push_back(A1);
			real_function_3d A2=real_factory_3d(world).f(c2v_A2);xoperators.push_back(A2);
			real_function_3d B1=real_factory_3d(world).f(c2v_B1);xoperators.push_back(B1);
			real_function_3d B2=real_factory_3d(world).f(c2v_B2);xoperators.push_back(B2);
		}
		else if(exop == "C2h"){
			real_function_3d tmp0 =real_factory_3d(world).f(c2v_Ag); xoperators.push_back(tmp0);
			real_function_3d tmp1 =real_factory_3d(world).f(c2v_Bg); xoperators.push_back(tmp1);
			real_function_3d tmp2 =real_factory_3d(world).f(c2v_Au); xoperators.push_back(tmp2);
			real_function_3d tmp3 =real_factory_3d(world).f(c2v_Bu); xoperators.push_back(tmp3);

		}
		else if(exop == "D2h"){
			real_function_3d tmp0 = real_factory_3d(world).f(d2h_Ag );xoperators.push_back(tmp0);
			real_function_3d tmp1 = real_factory_3d(world).f(d2h_B1g);xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(d2h_B2g);xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(d2h_B3g);xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(d2h_Au );xoperators.push_back(tmp4);
			real_function_3d tmp5 = real_factory_3d(world).f(d2h_B1u);xoperators.push_back(tmp5);
			real_function_3d tmp6 = real_factory_3d(world).f(d2h_B2u);xoperators.push_back(tmp6);
			real_function_3d tmp7 = real_factory_3d(world).f(d2h_B3u);xoperators.push_back(tmp7);

		}
		else if(exop == "C2"){
			real_function_3d tmp0 = real_factory_3d(world).f(c2_A );xoperators.push_back(tmp0);
			real_function_3d tmp1 = real_factory_3d(world).f(c2_B) ;xoperators.push_back(tmp1);
		}
		else if (exop == "Cs"){
			 real_function_3d tmp0 = real_factory_3d(world).f(cs_A);xoperators.push_back(tmp0);
			 real_function_3d tmp1 = real_factory_3d(world).f(cs_B);xoperators.push_back(tmp1);
		}
		else if(exop == "D3"){
			real_function_3d tmp0  = real_factory_3d(world).f(d3_A1);xoperators.push_back(tmp0);
			real_function_3d tmp1  = real_factory_3d(world).f(d3_A2);xoperators.push_back(tmp1);
			real_function_3d tmp2  = real_factory_3d(world).f(d3_E1);xoperators.push_back(tmp2);
			real_function_3d tmp3  = real_factory_3d(world).f(d3_E2);xoperators.push_back(tmp3);
		}
		else if(exop=="C3v"){
			real_function_3d tmp0  = real_factory_3d(world).f(c3v_A1);xoperators.push_back(tmp0);
			real_function_3d tmp1  = real_factory_3d(world).f(c3v_A2);xoperators.push_back(tmp1);
			real_function_3d tmp2  = real_factory_3d(world).f(c3v_E11);xoperators.push_back(tmp2);
			real_function_3d tmp3  = real_factory_3d(world).f(c3v_E12);xoperators.push_back(tmp3);
		}
		else if(exop=="benzene"){
			real_function_3d tmp0  = real_factory_3d(world).f(benzene_0);xoperators.push_back(tmp0);
			real_function_3d tmp1  = real_factory_3d(world).f(benzene_1);xoperators.push_back(tmp1);
			real_function_3d tmp2  = real_factory_3d(world).f(benzene_2);xoperators.push_back(tmp2);
			real_function_3d tmp3  = real_factory_3d(world).f(benzene_3);xoperators.push_back(tmp3);
			real_function_3d tmp4  = real_factory_3d(world).f(benzene_4);xoperators.push_back(tmp4);
			real_function_3d tmp5  = real_factory_3d(world).f(benzene_5);xoperators.push_back(tmp5);
			real_function_3d tmp6  = real_factory_3d(world).f(benzene_6);xoperators.push_back(tmp6);
		}
		else if(exop=="uracil"){
			 real_function_3d tmp1 = real_factory_3d(world).f(ura_1);xoperators.push_back(tmp1);
			 real_function_3d tmp2 = real_factory_3d(world).f(ura_2);xoperators.push_back(tmp2);
			 real_function_3d tmp3 = real_factory_3d(world).f(ura_3);xoperators.push_back(tmp3);
			 real_function_3d tmp4 = real_factory_3d(world).f(ura_4);xoperators.push_back(tmp4);
			 real_function_3d tmp5 = real_factory_3d(world).f(ura_5);xoperators.push_back(tmp5);
		}
		else if(exop=="D2d"){
			real_function_3d tmp0 = real_factory_3d(world).f(d2d_A1); xoperators.push_back(tmp0);
			real_function_3d tmp1 = real_factory_3d(world).f(d2d_A2); xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(d2d_B1); xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(d2d_B2); xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(d2d_E1); xoperators.push_back(tmp4);
			real_function_3d tmp5 = real_factory_3d(world).f(d2d_E2); xoperators.push_back(tmp5);
		}
		else if(exop=="reduced-D6h"){
			real_function_3d tmp0 = real_factory_3d(world).f(d6h_A1g  );xoperators.push_back(tmp0);
			real_function_3d tmp1 = real_factory_3d(world).f(d6h_E1g1 );xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(d6h_E1g2 );xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(d6h_E2g1 );xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(d6h_E2g2 );xoperators.push_back(tmp4);
			real_function_3d tmp5 = real_factory_3d(world).f(d6h_A2u  );xoperators.push_back(tmp5);
			real_function_3d tmp6 = real_factory_3d(world).f(d6h_B1u  );xoperators.push_back(tmp6);
			real_function_3d tmp7 = real_factory_3d(world).f(d6h_B2u  );xoperators.push_back(tmp7);
			real_function_3d tmp8 = real_factory_3d(world).f(d6h_E1u1 );xoperators.push_back(tmp8);
			real_function_3d tmp9 = real_factory_3d(world).f(d6h_E1u2 );xoperators.push_back(tmp9);
		}
		else if(exop=="D6h"){
			real_function_3d tmp0 = real_factory_3d(world).f(d6h_A1g  );xoperators.push_back(tmp0);
			real_function_3d tmp1 = real_factory_3d(world).f(d6h_E1g1 );xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(d6h_E1g2 );xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(d6h_E2g1 );xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(d6h_E2g2 );xoperators.push_back(tmp4);
			real_function_3d tmp5 = real_factory_3d(world).f(d6h_A2u  );xoperators.push_back(tmp5);
			real_function_3d tmp6 = real_factory_3d(world).f(d6h_B1u  );xoperators.push_back(tmp6);
			real_function_3d tmp7 = real_factory_3d(world).f(d6h_B2u  );xoperators.push_back(tmp7);
			real_function_3d tmp8 = real_factory_3d(world).f(d6h_E1u1 );xoperators.push_back(tmp8);
			real_function_3d tmp9 = real_factory_3d(world).f(d6h_E1u2 );xoperators.push_back(tmp9);
			std::cout << "Random guess exops for irreps A2g,B1g,B2g,A2u " << std::endl;
			real_function_3d tmp = real_factory_3d(world);
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			for(size_t k=0;k<polynom_basis_.size();k++){
				// dont take linear, quadratic and cubic functions
				double range = 0.0;
				if(k>19) range = 0.1;
				double c = random_number()*range;
				tmp += c*polynom_basis_[k];
				std::cout <<std::fixed << std::setprecision(2) << " "<< c << " " <<  key_[k];
			}
			xoperators.push_back(tmp);
			std::cout << std::endl;

		}
		else if(exop=="sequential_full"){
			real_function_3d tmp = real_factory_3d(world).f(seq_full); xoperators.push_back(tmp);
		}
		else if(exop=="sequential_dipole"){
			real_function_3d tmp = real_factory_3d(world).f(seq_dipole); xoperators.push_back(tmp);
		}
		else if(exop=="sequential_dipole+"){
			real_function_3d tmp = real_factory_3d(world).f(seq_dipolep); xoperators.push_back(tmp);
		}
		else if (exop=="big_3"){
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			xoperators = polynom_basis_;
			if(world.rank()==0){
				std::cout << "Guess exop is big_3 operators are:" << std::endl;
				for(size_t i=0;i<3+6+10;i++) std::cout << " " <<key_[i] << " ";
				std::cout << std::endl;
			}
			xoperators.erase(xoperators.begin()+3+6+10+1,xoperators.end());
		}
		else if (exop=="big_4"){
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			xoperators = polynom_basis_;
			if(world.rank()==0){
				std::cout << "Guess exop is big_4 excitation operators are:" << std::endl;
				for(size_t i=0;i<3+6+10+15;i++) std::cout << " " <<key_[i] << " ";
				std::cout << std::endl;
			}
			xoperators.erase(xoperators.begin()+3+6+10+15+1,xoperators.end());
		}
		else if (exop=="random_4"){
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			for(size_t i=0;i<4;i++){
				std::cout << "Random guess exops " << i << std::endl;
				real_function_3d tmp = real_factory_3d(world);
				for(size_t k=0;k<polynom_basis_.size();k++){
					double range = 0.1;
					if(k>3) range =0.1;
					if(k>10) range = 0.1;
					if(k>20) range =0.1;
					if(k>34) range =0.1;
					double c = random_number()*range;
					tmp += c*polynom_basis_[k];
					std::cout <<std::fixed << std::setprecision(2) << " "<< c << " " <<  key_[k];
				}
				xoperators.push_back(tmp);
				std::cout << std::endl;
			}
			//std::cout << "testing random numbers" << std::endl;
			//for(size_t i=0;i<100;i++) std::cout<< std::fixed << std::setprecision(1) << random_number()<< std::endl;
		}
		else if (exop=="random_8"){
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			for(size_t i=0;i<8;i++){
				std::cout << "Random guess exops " << i << std::endl;
				real_function_3d tmp = real_factory_3d(world);
				for(size_t k=0;k<polynom_basis_.size();k++){
					double range = 0.1;
					if(k>3) range =0.1;
					if(k>10) range = 0.1;
					if(k>20) range =0.1;
					if(k>34) range =0.1;
					double c = random_number()*range;
					tmp += c*polynom_basis_[k];
					std::cout <<std::fixed << std::setprecision(2) << " "<< c << " " <<  key_[k];
				}
				xoperators.push_back(tmp);
				std::cout << std::endl;
			}
			//std::cout << "testing random numbers" << std::endl;
			//for(size_t i=0;i<100;i++) std::cout<< std::fixed << std::setprecision(1) << random_number()<< std::endl;
		}
		else if (exop=="random"){
			// check if polynom_basis is created, if not, do so
			get_polynom_basis(world);
			real_function_3d tmp = real_factory_3d(world);
			for(size_t k=0;k<polynom_basis_.size();k++){
				double range = 0.1;
				double c = random_number()*range;
				tmp += c*polynom_basis_[k];
			}
			xoperators.push_back(tmp);
		}
		else if (exop=="oscillating"){
			real_function_3d tmp1 = real_factory_3d(world).f(sinus_x);  xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(sinus_y);  xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(sinus_z);  xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(cosinus_x);  xoperators.push_back(tmp4);
			real_function_3d tmp5 = real_factory_3d(world).f(cosinus_y);  xoperators.push_back(tmp5);
			real_function_3d tmp6 = real_factory_3d(world).f(cosinus_z);  xoperators.push_back(tmp6);
		}
		else if (exop=="oscillating_2"){
			real_function_3d tmp1 = real_factory_3d(world).f(sinus_x);  xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(sinus_y);  xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(sinus_z);  xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(cosinus_r); xoperators.push_back(tmp4);
		}
		else if (exop=="oscillating_small"){
			real_function_3d tmp1 = real_factory_3d(world).f(sinus_r); xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(cosinus_r); xoperators.push_back(tmp2);
		}
		else if (exop=="oscillating_small_2"){
			real_function_3d tmp1 = real_factory_3d(world).f(sinus_x); xoperators.push_back(tmp1);
			real_function_3d tmp2 = real_factory_3d(world).f(sinus_y); xoperators.push_back(tmp2);
			real_function_3d tmp3 = real_factory_3d(world).f(sinus_z); xoperators.push_back(tmp3);
			real_function_3d tmp4 = real_factory_3d(world).f(cosinus_2r); xoperators.push_back(tmp4);
		}
		else {
			std::cout << "exop keyword " << exop << "is not known" << std::endl;
			MADNESS_EXCEPTION("Unknown keyword in exop struct",1);
		}

		return xoperators;

	}


private:

	// Random Number generator
	double random_number(){
		int signa = rand()% 10 + 1;
		int signb = rand()% 10 +1;
		double tmp = (double) signa - (double) signb;
		double  sign = 1.0;
		if(tmp<0.0) sign = -1.0;
		int a = rand() % 10 + 1;
		double aa = (double) a * sign;
		return aa;
	}

	// oscillating functions
	static double sinus_x  (const coord_3d &r){return sin(r[0]);}
	static double sinus_y  (const coord_3d &r){return sin(r[1]);}
	static double sinus_z  (const coord_3d &r){return sin(r[2]);}
	static double cosinus_x(const coord_3d &r){return cos(r[0]);}
	static double cosinus_y(const coord_3d &r){return cos(r[1]);}
	static double cosinus_z(const coord_3d &r){return cos(r[2]);}

	static double sinus_r (const coord_3d &r){return sin(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));}
	static double cosinus_r (const coord_3d &r){return cos(sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));}

	static double sinus_2r (const coord_3d &r){return sin(2*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));}
	static double cosinus_2r (const coord_3d &r){return cos(2*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]));}

	static double sinus_2x  (const coord_3d &r){return sin(2*r[0]);}
	static double sinus_2y  (const coord_3d &r){return sin(2*r[1]);}
	static double sinus_2z  (const coord_3d &r){return sin(2*r[2]);}
	static double cosinus_2x(const coord_3d &r){return cos(2*r[0]);}
	static double cosinus_2y(const coord_3d &r){return cos(2*r[1]);}
	static double cosinus_2z(const coord_3d &r){return cos(2*r[2]);}

	// Excitation operator creators for custom guess
	double xc_function(const coord_3d &r){return r[0];}
	double yc_function(const coord_3d &r){return r[1];}
    double zc_function(const coord_3d &r){return r[2];}

	// Dipole operators
	static double x_function(const coord_3d &r){return r[0];}
	static double y_function(const coord_3d &r){return r[1];}
	static double z_function(const coord_3d &r){return r[2];}
	// quadrupole operators
	static double xy_function(const coord_3d &r){return r[0]*r[1];}
	static double xz_function(const coord_3d &r){return r[0]*r[2];}
	static double yz_function(const coord_3d &r){return r[1]*r[2];}
	static double xx_yy_function(const coord_3d &r){return r[0]*r[0]-r[1]*r[1];}
	static double xx_yy_2zz_function(const coord_3d &r){return r[0]*r[0]+r[1]*r[1]-2.0*r[3]*r[3];}
	static double zz_function(const coord_3d &r){return r[2]*r[2];}
	static double xx_function(const coord_3d &r){return r[0]*r[0];}
	static double yy_function(const coord_3d &r){return r[1]*r[1];}
	// dark states
	static double r_function(const coord_3d &r){return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);}
	static double rr_function(const coord_3d &r){return r[0]*r[0]+r[1]*r[1]+r[2]*r[2];}
	// Hybrids (just for testing)
	static double sp30(const coord_3d &r){return r_function(r)+r[0]+r[1]+r[2];}
	static double sp31(const coord_3d &r){return r_function(r)+r[0]-r[1]-r[2];}
	static double sp32(const coord_3d &r){return r_function(r)-r[0]+r[1]-r[2];}
	static double sp33(const coord_3d &r){return r_function(r)-r[0]-r[1]+r[2];}
	// Excitation operator for sequential guess as linear combination of x y z r
	static double f_seqop(const coord_3d &r) {return (x_function(r)+y_function(r)+z_function(r)+r_function(r));}
	// cubic excitations
	static double xxx_function(const coord_3d &r){ return r[0]*r[0]*r[0];   }
	static double yyy_function(const coord_3d &r){ return r[1]*r[1]*r[1];   }
	static double zzz_function(const coord_3d &r){ return r[2]*r[2]*r[2];   }
	static double xyy_function(const coord_3d &r){ return r[0]*r[1]*r[1];   }
	static double xxy_function(const coord_3d &r){ return r[0]*r[0]*r[1];   }
	static double xzz_function(const coord_3d &r){ return r[0]*r[2]*r[2];   }
	static double xxz_function(const coord_3d &r){ return r[0]*r[0]*r[2];   }
	static double yzz_function(const coord_3d &r){ return r[1]*r[2]*r[2];   }
	static double yyz_function(const coord_3d &r){ return r[1]*r[1]*r[2];   }
	static double xyz_function(const coord_3d &r){ return r[0]*r[1]*r[2];   }
// quartic excitations
	static double xxxx_function(const coord_3d &r){return r[0]*r[0]*r[0]*r[0];}
	static double xxxy_function(const coord_3d &r){return r[0]*r[0]*r[0]*r[1];}
	static double xxxz_function(const coord_3d &r){return r[0]*r[0]*r[0]*r[2];}
	static double xxyy_function(const coord_3d &r){return r[0]*r[0]*r[1]*r[1];}
	static double xxyz_function(const coord_3d &r){return r[0]*r[0]*r[1]*r[2];}
	static double xxzz_function(const coord_3d &r){return r[0]*r[0]*r[2]*r[2];}
	static double xyyy_function(const coord_3d &r){return r[0]*r[1]*r[1]*r[1];}
	static double xyyz_function(const coord_3d &r){return r[0]*r[1]*r[1]*r[2];}
	static double xyzz_function(const coord_3d &r){return r[0]*r[1]*r[2]*r[2];}
	static double xzzz_function(const coord_3d &r){return r[0]*r[2]*r[2]*r[2];}
	static double yyyy_function(const coord_3d &r){return r[1]*r[1]*r[1]*r[1];}
	static double yyyz_function(const coord_3d &r){return r[1]*r[1]*r[1]*r[2];}
	static double yyzz_function(const coord_3d &r){return r[1]*r[1]*r[2]*r[2];}
	static double yzzz_function(const coord_3d &r){return r[1]*r[2]*r[2]*r[2];}
	static double zzzz_function(const coord_3d &r){return r[2]*r[2]*r[2]*r[2];}


	/// Ci

	/// C2h
	static double c2v_Ag(const coord_3d &r){return rr_function(r)+xy_function(r);}
	static double c2v_Bg(const coord_3d &r){return xy_function(r)+yz_function(r);}
	static double c2v_Au(const coord_3d &r){return z_function(r)+zzz_function(r)+xyz_function(r)+xxz_function(r)+yyz_function(r);}
	static double c2v_Bu(const coord_3d &r){return x_function(r)+y_function(r)+xzz_function(r)+yzz_function(r)+xxy_function(r)+xyy_function(r)+xxx_function(r)+yyy_function(r);}


	/// D2h
	static double d2h_Ag   (const coord_3d &r){return rr_function(r) ;}
	static double d2h_B1g  (const coord_3d &r){return xy_function(r) ;}
	static double d2h_B2g  (const coord_3d &r){return xz_function(r) ;}
	static double d2h_B3g  (const coord_3d &r){return yz_function(r) ;}
	static double d2h_Au   (const coord_3d &r){return xyz_function(r);}
	static double d2h_B1u  (const coord_3d &r){return z_function(r)+zzz_function(r)+yyz_function(r)+xzz_function(r) ;}
	static double d2h_B2u  (const coord_3d &r){return y_function(r)+yzz_function(r)+xxy_function(r)+yyy_function(r) ;}
	static double d2h_B3u  (const coord_3d &r){return x_function(r)+xzz_function(r)+xyy_function(r)+xxx_function(r) ;}
	/// C1

	/// C2
	static double c2_A(const coord_3d &r){return z_function(r)+rr_function(r)+xy_function(r)+zzz_function(r)+xyz_function(r)+yyz_function(r)+xxz_function(r);}
	static double c2_B(const coord_3d &r){return x_function(r)+y_function(r)+yz_function(r)+xz_function(r)+xzz_function(r)+yzz_function(r)+xxy_function(r)+xyy_function(r)+xxx_function(r)+yyy_function(r);}

	/// Cs
	static double cs_A(const coord_3d &r){return x_function(r)+y_function(r)+xx_function(r)+yy_function(r)+zz_function(r)+xy_function(r)
			+xzz_function(r)+yzz_function(r)+xxy_function(r)+xyy_function(r)+xxx_function(r)+yyy_function(r);}
	static double cs_B(const coord_3d &r){return z_function(r)+yz_function(r)+xz_function(r)+zzz_function(r)+xyz_function(r)+yyz_function(r)+xxz_function(r)  ;}

	/// C2v
	static double c2v_A1(const coord_3d &r){return z_function(r)+xx_function(r)+yy_function(r)+zz_function(r)+xxz_function(r)+yyz_function(r)+zzz_function(r) ;}
	static double c2v_A2(const coord_3d &r){return xy_function(r)+xyz_function(r) ;}
	static double c2v_B1(const coord_3d &r){return x_function(r)+xz_function(r)+xzz_function(r)+xxx_function(r)+xyy_function(r);}
	static double c2v_B2(const coord_3d &r){return y_function(r)+yz_function(r)+yzz_function(r)+yyy_function(r)+xxy_function(r);}

	/// D2

	/// D3
	static double d3_A1(const coord_3d &r){return rr_function(r)+xxx_function(r)-3.0*xyy_function(r);}
	static double d3_A2(const coord_3d &r){return z_function(r)+zzz_function(r)+3.0*xxy_function(r)-yyy_function(r)+xxz_function(r)+yyz_function(r);}
	static double d3_E1(const coord_3d &r){return x_function(r)+xx_function(r)-yy_function(r)+xz_function(r)+xzz_function(r)+xyz_function(r)+xxx_function(r)+xyy_function(r);}
	static double d3_E2(const coord_3d &r){return y_function(r)+xy_function(r)+yz_function(r)+yzz_function(r)+xxz_function(r)-yyz_function(r)+xxy_function(r)+yyy_function(r);}

	/// C3v
	static double c3v_A1 (const coord_3d &r){return z_function(r)+rr_function(r)+zzz_function(r)+xxx_function(r)-3.0*xyy_function(r)+xxz_function(r)+yyz_function(r) ;}
	static double c3v_A2 (const coord_3d &r){return 3.0*xxy_function(r)-yyy_function(r) ;}
	static double c3v_E11(const coord_3d &r){return x_function(r)+xx_function(r)-yy_function(r)+xz_function(r)+xzz_function(r)+xyz_function(r)+xxx_function(r)+xyy_function(r);}
	static double c3v_E12(const coord_3d &r){return y_function(r)+xy_function(r)+yz_function(r)+yzz_function(r)+xxz_function(r)-yyz_function(r)+xxy_function(r)+yyy_function(r);}

	/// Benzene
	 static double benzene_0(const coord_3d &r){return xxx_function(r) -3.0*xyy_function(r) ;}
	 static double benzene_1(const coord_3d &r){return 3.0*xxy_function(r) -yyy_function(r);}
	 static double benzene_2(const coord_3d &r){return xz_function(r);}
	 static double benzene_3(const coord_3d &r){return yz_function(r);}
	 static double benzene_4(const coord_3d &r){return z_function(r);}
	 static double benzene_5(const coord_3d &r){return xyz_function(r);}
	 static double benzene_6(const coord_3d &r){return xxz_function(r) - yyz_function(r);}

	 /// D6h without irreps A2g, B1g, B2g, A2u because no linear, quadratic or cubic polynomial transforms like them
	 // Suggestion to include them: Make random guess functions with polynomials of order 4 and higher and orthonormalize
	 static double d6h_A1g (const coord_3d &r){return xx_function(r)+yy_function(r)+zz_function(r) ;}
	 static double d6h_E1g1(const coord_3d &r){return xz_function(r);}
	 static double d6h_E1g2(const coord_3d &r){return yz_function(r);}
	 static double d6h_E2g1(const coord_3d &r){return xx_function(r)-yy_function(r);}
	 static double d6h_E2g2(const coord_3d &r){return xy_function(r);}
	 static double d6h_A2u (const coord_3d &r){return z_function(r) + zzz_function(r) + xxz_function(r) + yyz_function(r);}
	 static double d6h_B1u (const coord_3d &r){return xxx_function(r)-3.0*xyy_function(r);}
	 static double d6h_B2u (const coord_3d &r){return 3.0*xxy_function(r)-yyy_function(r);}
	 static double d6h_E1u1(const coord_3d &r){return x_function(r) + xzz_function(r) + xxx_function(r) +xyy_function(r);}
	 static double d6h_E1u2(const coord_3d &r){return y_function(r) + yzz_function(r) + xxy_function(r) +yyy_function(r);}
	 //static double rest(const coord_3d &r){return ;}

	 /// D2d
	 static double d2d_A1(const coord_3d &r){return rr_function(r) +xyz_function(r)   ;}
	 static double d2d_A2(const coord_3d &r){return xxz_function(r) - yyz_function(r)   ;}
	 static double d2d_B1(const coord_3d &r){return xx_function(r) -yy_function(r)   ;}
	 static double d2d_B2(const coord_3d &r){return z_function(r) +xy_function(r) +zzz_function(r) + xxz_function(r) +yyz_function(r)   ;}
	 static double d2d_E1(const coord_3d &r){return x_function(r)+xz_function(r)+xzz_function(r)+xyy_function(r)+xxx_function(r);}
	 static double d2d_E2(const coord_3d &r){return y_function(r)+yz_function(r)+yzz_function(r)+xxy_function(r)+yyy_function(r);}

	 /// Sequential guesses:
	 /// full
	 static double seq_full(const coord_3d&r){return x_function(r)+y_function(r)+z_function(r)+rr_function(r)+xy_function(r)+xz_function(r)+yz_function(r)
	 + xxx_function(r)+yyy_function(r)+zzz_function(r)+xxz_function(r)+xxy_function(r)+xyy_function(r)+xyz_function(r)+xzz_function(r)+yyz_function(r)
	 + yzz_function(r);}
	 // dipole
	 static double seq_dipole(const coord_3d&r){return x_function(r)+y_function(r)+z_function(r);}
	 // dipole+
	 static double seq_dipolep(const coord_3d&r){return seq_dipole(r)+rr_function(r);}

	 /// Uracil custom guess for testing purposes
	 static double ura_1(const coord_3d &r){return 0.5*x_function(r)-1.5*y_function(r)+2*xx_function(r)-2.5*yy_function(r)-0.05*zz_function(r)+xy_function(r) ;}
	 static double ura_2(const coord_3d &r){return -2.3*x_function(r) -2*y_function(r)-0.44*xx_function(r)+yy_function(r)-1.3*zz_function(r)-xy_function(r);}
	 static double ura_3(const coord_3d &r){return 0.004*z_function(r)+0.12*xz_function(r)+0.09*yz_function(r);}
	 static double ura_4(const coord_3d &r){return 0.14*z_function(r)+xz_function(r)+yz_function(r);}
	 static double ura_5(const coord_3d &r){return -0.003*z_function(r)+0.1*xz_function(r)-0.01*yz_function(r);}



};
}



#endif /* TDA_EXOPS_H_ */
