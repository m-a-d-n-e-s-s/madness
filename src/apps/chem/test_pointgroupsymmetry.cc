/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#include <madness/mra/mra.h>
#include <madness/mra/functypedefs.h>
#include <chem/pointgroupoperator.h>
#include <chem/pointgroupsymmetry.h>

using namespace madness;

namespace madness {
extern std::vector<std::string> cubefile_header(std::string filename="input",
		const bool& no_orient=false);
}


double dgaussian(const coord_3d& r) {
    double x=r[0]-1, y=r[1]-2, z=r[2]-3;
    return x*y*exp(-(x*x + y*y + z*z));
}

double xyyygaussian(const coord_2d& r) {
    double x=r[0], y=r[1];
    return x*y*y*y*exp(-(x*x + y*y ));
}


double gaussian_shift_2d(const coord_2d& r) {
    double x=r[0]-2., y=r[1]-1.0;
    return exp(-(x*x + y*y ));
}

double gaussian_shift_3d(const coord_3d& r) {
    double x=r[0]-2., y=r[1]-1.0, z=r[2]-3.0;
    return exp(-(x*x + y*y + z*z ));
}

void test_symmetry_operators(World& world) {

	const real_function_2d f=real_factory_2d(world).f(gaussian_shift_2d);
	pg_operator sx=pg_sigma_x();
    pg_operator sy=pg_sigma_y();
    pg_operator c4=pg_c4();
    print("applying",sx);
	real_function_2d opf1=sx(f);
    real_function_2d opf2=sy(f);
    real_function_2d opf3=c4(f);
    std::vector<real_function_2d> ff=vector_factory<real_function_2d>(f,opf1,opf2,opf3);
	plot_plane(world,ff,"f");


	real_function_3d p=real_factory_3d(world).f(dgaussian);
	pg_operator sxy=pg_sigma_xy();
	real_function_3d opp1=sxy(p);
	pg_operator syz=pg_sigma_yz();
	real_function_3d opp2=syz(p);

	// plot the Gaussian cube file
	std::vector<std::string> molecular_info=cubefile_header("input",false);

	std::string filename="dgaussian.cube";
	plot_cubefile<3>(world,p,filename,molecular_info);

	std::string filename1="opp1.cube";
	plot_cubefile<3>(world,opp1,filename1,molecular_info);

	std::string filename2="opp2.cube";
	plot_cubefile<3>(world,opp2,filename2,molecular_info);

	pg_operator c4x=pg_c4x();
	pg_operator c4y=pg_c4y();
	pg_operator c4z=pg_c4z();
	real_function_3d c4x_p=c4x(p);
	real_function_3d c4y_p=c4y(p);
	real_function_3d c4z_p=c4z(p);
	plot_cubefile<3>(world,c4x_p,"c4x_p.cube",molecular_info);
	plot_cubefile<3>(world,c4y_p,"c4y_p.cube",molecular_info);
	plot_cubefile<3>(world,c4z_p,"c4z_p.cube",molecular_info);

	pg_operator c2x=pg_c2x();
	pg_operator c2y=pg_c2y();
	pg_operator c2z=pg_c2z();

	print("test1");
	real_function_3d gonex=c2x(p) - c4x(c4x(p));
	real_function_3d goney=c2y(p) - c4y(c4y(p));
	real_function_3d gonez=c2z(p) - c4z(c4z(p));
	double nx=gonex.norm2();
	double ny=goney.norm2();
	double nz=gonez.norm2();
	print("norm(gone)",nx,ny,nz);

}

void test_projector(World& world) {
	const real_function_2d ff=real_factory_2d(world).f(gaussian_shift_2d);

	const vector_real_function_2d f(1,ff);
	projector_irrep proj_1("cs","a'");
	const std::vector<vector_real_function_2d> result1=proj_1.apply_symmetry_operators(f);

	projector_irrep proj_2("cs","a''");
	const std::vector<vector_real_function_2d> result2=proj_2.apply_symmetry_operators(f);

	proj_2.set_irrep("all");
	const std::vector<vector_real_function_2d> result3=proj_2.apply_symmetry_operators(f);

	plot_plane(world,ff,result1[0][0],result1[1][0],"result1");
	plot_plane(world,ff,result2[0][0],result2[1][0],"result2");
	plot_plane(world,ff,result3[0][0],result3[1][0],"result3");

}

void test_irreps(World& world) {
	const real_function_2d f=real_factory_2d(world).f(gaussian_shift_2d);

	// 3d
	pg_operator sxy=pg_sigma_xy();
	real_function_3d p1=real_factory_3d(world).f(gaussian_shift_3d);
	real_function_3d p2=real_factory_3d(world).f(gaussian_shift_3d);
	real_function_3d p3=sxy(p1);
	std::vector<real_function_3d> pp=vector_factory<real_function_3d>(p1,p2,p3);

//	projector_irrep proj_c2v_a1("c2v","a1");
//	std::vector<real_function_3d> a1=proj_c2v_a1(pp);
//	projector_irrep proj_c2v_a2("c2v","a2");
//	std::vector<real_function_3d> a2=proj_c2v_a2(pp);
//	projector_irrep proj_c2v_b1("c2v","b1");
//	std::vector<real_function_3d> b1=proj_c2v_b1(pp);
//	projector_irrep proj_c2v_b2("c2v","b2");
//	std::vector<real_function_3d> b2=proj_c2v_b2(pp);

//	proj_c2v_a1.print_character_table();
//	std::vector<std::string> reduced=proj_c2v_a1.reduce("a1","b2","b1");


//	std::vector<std::string> molecular_info=cubefile_header("input",false);
//
//	plot_cubefile<3>(world,a1,"a1.cube",molecular_info);
//	plot_cubefile<3>(world,a2,"a2.cube",molecular_info);
//	plot_cubefile<3>(world,b1,"b1.cube",molecular_info);
//	plot_cubefile<3>(world,b2,"b2.cube",molecular_info);

}

void test_nemo(World& world) {

	// load nemos from file
	vector_real_function_3d nemos;
	for (int i=0; i<5; ++i) {
		real_function_3d tmp=real_factory_3d(world);
		std::string name="nemo"+stringify(i);
		load(tmp,name);
		nemos.push_back(tmp);
	}
	int k=nemos[0].get_impl()->get_k();
	FunctionDefaults<3>::set_k(k);

	// project nemos on irreps
	projector_irrep proj_c2v("c2v");

//	for (const std::string& irrep : proj_c2v.get_all_irreps()) {
//		proj_c2v.set_irrep(irrep);
//		vector_real_function_3d a1new=proj_c2v(nemos);
//
//		if (a1new.size()==0) continue;
//
//		Tensor<double> ovlp2=matrix_inner(world,a1new,a1new);
//		print("overlap of", irrep,"orbitals");
//		print(ovlp2);
//
//
//		// plot a1 orbitals
//		std::vector<std::string> molecular_info=cubefile_header("input",false);
//		for (int i=0; i<a1new.size(); ++i) {
//			std::string name="irrep"+irrep+stringify(i)+".cube";
//			plot_cubefile<3>(world,a1new[i],name,molecular_info);
//		}
//	}

}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    FunctionDefaults<2>::set_cubic_cell(-6.0,6.0);
    FunctionDefaults<3>::set_cubic_cell(-6.0,6.0);


    test_projector(world);
//    test_symmetry_operators(world);
//    test_irreps(world);
//    test_nemo(world);

    madness::finalize();
    return 0;
}
