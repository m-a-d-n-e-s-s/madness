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

public:
	exoperators(){}

	vecfuncT get_polynom_basis(World &world){
		vecfuncT basis;
		 //dipole
		 real_function_3d x   = real_factory_3d(world).f(x_function);  basis.push_back(x);
		 real_function_3d y   = real_factory_3d(world).f(y_function);  basis.push_back(y);
		 real_function_3d z   = real_factory_3d(world).f(z_function);  basis.push_back(z);
		 //quadrupole
		 real_function_3d xx  = real_factory_3d(world).f(xx_function); basis.push_back(xx);
		 real_function_3d yy  = real_factory_3d(world).f(yy_function); basis.push_back(yy);
		 real_function_3d zz  = real_factory_3d(world).f(zz_function); basis.push_back(zz);
		 real_function_3d xy  = real_factory_3d(world).f(xy_function); basis.push_back(xy);
		 real_function_3d xz  = real_factory_3d(world).f(xz_function); basis.push_back(xz);
		 real_function_3d yz  = real_factory_3d(world).f(yz_function); basis.push_back(yz);

		 real_function_3d xxx  = real_factory_3d(world).f(xxx_function); basis.push_back(xxx);
		 real_function_3d yyy  = real_factory_3d(world).f(yyy_function); basis.push_back(yyy);
		 real_function_3d zzz  = real_factory_3d(world).f(zzz_function); basis.push_back(zzz);
		 real_function_3d xxy  = real_factory_3d(world).f(xxy_function); basis.push_back(xxy);
		 real_function_3d xxz  = real_factory_3d(world).f(xxz_function); basis.push_back(xxz);
		 real_function_3d xyy  = real_factory_3d(world).f(xyy_function); basis.push_back(xyy);
		 real_function_3d xyz  = real_factory_3d(world).f(xyz_function); basis.push_back(xyz);
		 real_function_3d xzz  = real_factory_3d(world).f(xzz_function); basis.push_back(xzz);
		 real_function_3d yyz  = real_factory_3d(world).f(yyz_function); basis.push_back(yyz);
		 real_function_3d yzz  = real_factory_3d(world).f(yzz_function); basis.push_back(yzz);

		 real_function_3d xxxx  = real_factory_3d(world).f(xxxx_function); basis.push_back(xxxx);
		 real_function_3d xxxy  = real_factory_3d(world).f(xxxy_function); basis.push_back(xxxy);
		 real_function_3d xxxz  = real_factory_3d(world).f(xxxz_function); basis.push_back(xxxz);
		 real_function_3d xxyy  = real_factory_3d(world).f(xxyy_function); basis.push_back(xxyy);
		 real_function_3d xxyz  = real_factory_3d(world).f(xxyz_function); basis.push_back(xxyz);
		 real_function_3d xxzz  = real_factory_3d(world).f(xxzz_function); basis.push_back(xxzz);
		 real_function_3d xyyy  = real_factory_3d(world).f(xyyy_function); basis.push_back(xyyy);
		 real_function_3d xyyz  = real_factory_3d(world).f(xyyz_function); basis.push_back(xyyz);
		 real_function_3d xyzz  = real_factory_3d(world).f(xyzz_function); basis.push_back(xyzz);
		 real_function_3d xzzz  = real_factory_3d(world).f(xzzz_function); basis.push_back(xzzz);
		 real_function_3d yyyy  = real_factory_3d(world).f(yyyy_function); basis.push_back(yyyy);
		 real_function_3d yyyz  = real_factory_3d(world).f(yyyz_function); basis.push_back(yyyz);
		 real_function_3d yyzz  = real_factory_3d(world).f(yyzz_function); basis.push_back(yyzz);
		 real_function_3d yzzz  = real_factory_3d(world).f(yzzz_function); basis.push_back(yzzz);
		 real_function_3d zzzz  = real_factory_3d(world).f(zzzz_function); basis.push_back(zzzz);
		 return basis;
	}

	std::vector<double> get_overlaps_with_guess(World &world,const vecfuncT &excitation,const vecfuncT &mo){
		if(excitation.size()!=mo.size()) MADNESS_EXCEPTION("Error in calculating overlap with the guess operators, excitation vector and mo vector have different sizes",1);
		vecfuncT quadbas = get_polynom_basis(world);
		 std::vector<double> overlaps;
		 for(size_t i=0;i< quadbas.size();i++){
			 vecfuncT guess_ex = mul(world,quadbas[i],mo);
			 double norm2= inner(world,guess_ex,guess_ex).sum();
			 scale(world,guess_ex,1.0/sqrt(norm2));
			 double overlap_tmp = inner(world, guess_ex,excitation).sum();
			 if(overlap_tmp < 1.e-5) overlap_tmp=0.0;
			 overlaps.push_back(overlap_tmp);
		 }
		 return overlaps;
	}

	vecfuncT get_custom_exops(World &world,const std::vector< std::string > input_lines){

		if(world.rank()==0)std::cout << "Creating custom guess operators ... " << input_lines.size() << " demanded" << std::endl;
		if(input_lines.size()==0) MADNESS_EXCEPTION("Custom guess demanded but no input was given",1);
		//std::cout << input_lines << std::endl;

		// get the coefficients from the strings
		std::vector<std::vector<double> > coefficients;
		for(size_t i=0;i<input_lines.size();i++){
		std::stringstream line(input_lines[i]);
		std::string name;
		std::vector<double> icoeff(19,0);

		while(line>>name){
			std::transform(name.begin(), name.end(), name.begin(), ::toupper);
			if(name == "X") line >> icoeff[0];
			if(name == "Y") line >> icoeff[1];
			if(name == "Z") line >> icoeff[2];
			if(name == "XX") line >> icoeff[3];
			if(name == "YY") line >> icoeff[4];
			if(name == "ZZ") line >> icoeff[5];
			if(name == "XY") line >> icoeff[6];
			if(name == "XZ") line >> icoeff[7];
			if(name == "YZ") line >> icoeff[8];
			if(name == "XXX") line >> icoeff[9];
			if(name == "YYY") line >> icoeff[10];
			if(name == "ZZZ") line >> icoeff[11];
			if(name == "XXY") line >> icoeff[12];
			if(name == "XXZ") line >> icoeff[13];
			if(name == "XYY") line >> icoeff[14];
			if(name == "XYZ") line >> icoeff[15];
			if(name == "XZZ") line >> icoeff[16];
			if(name == "YYZ") line >> icoeff[17];
			if(name == "YZZ") line >> icoeff[18];

		}
		coefficients.push_back(icoeff);
		}

		// controll output
		if(world.rank()==0){
			std::cout << "\n\n---------CUSTOM OPERATORS ARE----------\n\n" << std::endl;
			for(size_t i=0;i<coefficients.size();i++){
				std::cout << "\n exop" << i << " (x,y,z,xx,yy,zz,xy,xz,yz)" << std::endl;
					std::cout << std::fixed <<std::setprecision(2) << coefficients[i] << std::endl;
			}
			std::cout << "\n\n---------------------------------------\n\n" << std::endl;
		}


		vecfuncT exops;
		vecfuncT quadbas=get_polynom_basis(world);

		 for(size_t i=0;i<coefficients.size();i++){
			 real_function_3d tmp = real_factory_3d(world);
			 for(size_t j=0;j<coefficients[i].size();j++){
				  tmp += coefficients[i][j]*quadbas[j];
			 }
			 exops.push_back(tmp);
		 }
		 return exops;
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
		else if(exop=="sequential_full"){
			real_function_3d tmp = real_factory_3d(world).f(seq_full); xoperators.push_back(tmp);
		}
		else if(exop=="sequential_dipole"){
			real_function_3d tmp = real_factory_3d(world).f(seq_dipole); xoperators.push_back(tmp);
		}
		else if(exop=="sequential_dipole+"){
			real_function_3d tmp = real_factory_3d(world).f(seq_dipolep); xoperators.push_back(tmp);
		}
		else {
			std::cout << "exop keyword " << exop << "is not known" << std::endl;
			MADNESS_EXCEPTION("Unknown keyword in exop struct",1);
		}

		return xoperators;

	}


private:
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
	static double xxxx_function(const coord_3d &r){ r[0]*r[0]*r[0]*r[0];}
	static double xxxy_function(const coord_3d &r){ r[0]*r[0]*r[0]*r[1];}
	static double xxxz_function(const coord_3d &r){ r[0]*r[0]*r[0]*r[2];}
	static double xxyy_function(const coord_3d &r){ r[0]*r[0]*r[1]*r[1];}
	static double xxyz_function(const coord_3d &r){ r[0]*r[0]*r[1]*r[2];}
	static double xxzz_function(const coord_3d &r){ r[0]*r[0]*r[2]*r[2];}
	static double xyyy_function(const coord_3d &r){ r[0]*r[1]*r[1]*r[1];}
	static double xyyz_function(const coord_3d &r){ r[0]*r[1]*r[1]*r[2];}
	static double xyzz_function(const coord_3d &r){ r[0]*r[1]*r[2]*r[2];}
	static double xzzz_function(const coord_3d &r){ r[0]*r[2]*r[2]*r[2];}
	static double yyyy_function(const coord_3d &r){ r[1]*r[1]*r[1]*r[1];}
	static double yyyz_function(const coord_3d &r){ r[1]*r[1]*r[1]*r[2];}
	static double yyzz_function(const coord_3d &r){ r[1]*r[1]*r[2]*r[2];}
	static double yzzz_function(const coord_3d &r){ r[1]*r[2]*r[2]*r[2];}
	static double zzzz_function(const coord_3d &r){ r[2]*r[2]*r[2]*r[2];}


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
	static double d2h_B1u  (const coord_3d &r){return x_function(r)+zzz_function(r)+yyz_function(r)+xzz_function(r) ;}
	static double d2h_B2u  (const coord_3d &r){return y_function(r)+yzz_function(r)+xxy_function(r)+yyy_function(r) ;}
	static double d2h_B3u  (const coord_3d &r){return z_function(r)+xzz_function(r)+xyy_function(r)+xxx_function(r) ;}
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
