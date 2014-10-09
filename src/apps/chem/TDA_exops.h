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

	vecfuncT get_custom_exops(World &world,const std::vector< std::string > input_lines){

		if(world.rank()==0)std::cout << "Creating custom guess operators ... " << input_lines.size() << " demanded" << std::endl;
		if(input_lines.size()==0) MADNESS_EXCEPTION("Custom guess demanded but no input was given",1);
		//std::cout << input_lines << std::endl;

		// get the coefficients from the strings
		std::vector<std::vector<double> > coefficients;
		for(size_t i=0;i<input_lines.size();i++){
		std::stringstream line(input_lines[i]);
		std::string name;
		std::vector<double> icoeff(9,0);

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
		vecfuncT qudadbas;
		 real_function_3d x   = real_factory_3d(world).f(x_function);  qudadbas.push_back(x);
		 real_function_3d y   = real_factory_3d(world).f(y_function);  qudadbas.push_back(y);
		 real_function_3d z   = real_factory_3d(world).f(z_function);  qudadbas.push_back(z);
		 real_function_3d xx  = real_factory_3d(world).f(xx_function); qudadbas.push_back(xx);
		 real_function_3d yy  = real_factory_3d(world).f(yy_function); qudadbas.push_back(yy);
		 real_function_3d zz  = real_factory_3d(world).f(zz_function); qudadbas.push_back(zz);
		 real_function_3d xy  = real_factory_3d(world).f(xy_function); qudadbas.push_back(xy);
		 real_function_3d xz  = real_factory_3d(world).f(xz_function); qudadbas.push_back(xz);
		 real_function_3d yz  = real_factory_3d(world).f(yz_function); qudadbas.push_back(yz);

		 for(size_t i=0;i<coefficients.size();i++){
			 real_function_3d tmp = real_factory_3d(world);
			 for(size_t j=0;j<coefficients[i].size();j++){
				  tmp += coefficients[i][j]*qudadbas[j];
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
