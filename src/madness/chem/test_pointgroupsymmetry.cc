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
#include<madness/chem/pointgroupoperator.h>
#include<madness/chem/pointgroupsymmetry.h>

using namespace madness;

namespace madness {
extern std::vector<std::string> cubefile_header(std::string filename="input",
		const bool& no_orient=false);
}


double dgaussian(const coord_3d& r) {
    double x=r[0]-1, y=r[1]-2, z=r[2]-3;
    return x*y*exp(-(x*x + y*y + z*z));
}

// same as dgaussian, but reflected at the xz plane (y -> -y)
double s_xz_dgaussian(const coord_3d& r) {
    double x=r[0]-1, y=-r[1]-2, z=r[2]-3;
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

// same as gaussian_shift_2d, but reflected at the x plane (y -> -y)
double s_x_gaussian_shift_2d(const coord_2d& r) {
    double x=r[0]-2., y=-r[1]-1.0;
    return exp(-(x*x + y*y ));
}


double gaussian_shift_3d(const coord_3d& r) {
    double x=r[0]-2., y=r[1]-1.0, z=r[2]-3.0;
    return exp(-(x*x + y*y + z*z ));
}

double one1(const coord_3d& r) {
    return 1.0;
}

void plot_symmetry_operators(World& world) {

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


	// test function
	real_function_3d p=real_factory_3d(world).f(dgaussian);

	// plot the Gaussian cube file
	std::vector<std::string> molecular_info(1,"0 0 0.0 0.0 0.0\n");

	// mirroring
	pg_operator e=pg_identity();
	pg_operator i=pg_inversion();
	pg_operator sxy=pg_sigma_xy();
	pg_operator sxz=pg_sigma_xz();
	pg_operator syz=pg_sigma_yz();
	pg_operator c4x=pg_c4x();
	pg_operator c4y=pg_c4y();
	pg_operator c4z=pg_c4z();
	pg_operator c2x=pg_c2x();
	pg_operator c2y=pg_c2y();
	pg_operator c2z=pg_c2z();

	real_function_3d opp1=sxy(p);
	real_function_3d opp2=syz(p);

	plot_cubefile<3>(world,p,"dgaussian.cube",molecular_info);
	plot_cubefile<3>(world,opp1,"opp1.cube",molecular_info);
	plot_cubefile<3>(world,opp2,"opp2.cube",molecular_info);

	real_function_3d c4x_p=c4x(p);
	real_function_3d c4y_p=c4y(p);
	real_function_3d c4z_p=c4z(p);
	plot_cubefile<3>(world,c4x_p,"c4x_p.cube",molecular_info);
	plot_cubefile<3>(world,c4y_p,"c4y_p.cube",molecular_info);
	plot_cubefile<3>(world,c4z_p,"c4z_p.cube",molecular_info);
}


/// test a small number of operator multiplications, including a hard-wire test
int check_operator_multiplications_2d(World& world) {

	print("test operator multiplications 2D");

	pg_operator e=pg_identity();
	pg_operator i=pg_inversion();
	pg_operator sx=pg_sigma_x();
	pg_operator sy=pg_sigma_y();
	pg_operator c4=pg_c4();
	pg_operator c2=pg_c2();

	real_function_2d p=real_factory_2d(world).f(gaussian_shift_2d);
	real_function_2d sp=real_factory_2d(world).f(s_x_gaussian_shift_2d);

	// hard-wired test
	double h1=(sx(p)-sp).norm2();
//	print("norm(hardwire) ",h1);

	double n1=(c2(p) - c4(c4(p))).norm2();
	double n2=(sx(sy(p))-c2(p)).norm2();

	int result=0;
	if (h1+n1+n2>1.e-13) {
		result=1;
		print("large error norm in the operator multiplications");
		print(h1);
		print(n1,n2);
	} else {
		print("  .. all good");
	}
	return result;

}

/// test a small number of operator multiplications, including a hard-wire test
int check_operator_multiplications_3d(World& world) {

	print("test operator multiplications 3D");

	pg_operator e=pg_identity();
	pg_operator i=pg_inversion();
	pg_operator sxy=pg_sigma_xy();
	pg_operator sxz=pg_sigma_xz();
	pg_operator syz=pg_sigma_yz();
	pg_operator c4x=pg_c4x();
	pg_operator c4y=pg_c4y();
	pg_operator c4z=pg_c4z();
	pg_operator c2x=pg_c2x();
	pg_operator c2y=pg_c2y();
	pg_operator c2z=pg_c2z();

	real_function_3d p=real_factory_3d(world).f(dgaussian);
	real_function_3d sp=real_factory_3d(world).f(s_xz_dgaussian);

	// hard-wired test
	double h1=(sxz(p)-sp).norm2();
//	print("norm(hardwire) ",h1);

//	print("testing the norm of the differences c2-c4*c4 for x,y,z");
	real_function_3d gonex=c2x(p) - c4x(c4x(p));
	real_function_3d goney=c2y(p) - c4y(c4y(p));
	real_function_3d gonez=c2z(p) - c4z(c4z(p));
	double nx=gonex.norm2();
	double ny=goney.norm2();
	double nz=gonez.norm2();
//	print("norm(gone)",nx,ny,nz);


	double n1=(sxy(sxz(syz(p)))-i(p)).norm2();
//	print("norm s(s(s(f))) - i(f):",n1);

	double n2=(sxz(sxy(p))-c2x(p)).norm2();
//	print("norm sxz(sxy(p))-c2x(p):",n2);

	int result=0;
	if (h1+nx+ny+nz+n1+n2>1.e-13) {
		result=1;
		print("large error norm in the operator multiplications");
		print(h1);
		print(nx,ny,nz);
		print(n1,n2);
	} else {
		print("  .. all good");
	}
	return result;

}

/// perform all operator multiplication of the c2v group, subtract the result
int check_multiplication_table_c2v(World& world) {

	print("test multiplication table c2v");

	real_function_3d p=real_factory_3d(world).f(dgaussian);

	pg_operator e=pg_identity();
	pg_operator i=pg_inversion();
	pg_operator sxy=pg_sigma_xy();
	pg_operator sxz=pg_sigma_xz();
	pg_operator syz=pg_sigma_yz();
	pg_operator c4x=pg_c4x();
	pg_operator c4y=pg_c4y();
	pg_operator c4z=pg_c4z();
	pg_operator c2x=pg_c2x();
	pg_operator c2y=pg_c2y();
	pg_operator c2z=pg_c2z();

	// check multiplication table of c2v
	double a1=(e(e(p))-e(p)).norm2();
	double a2=(e(c2z(p))-c2z(p)).norm2();
	double a3=(e(sxz(p))-sxz(p)).norm2();
	double a4=(e(syz(p))-syz(p)).norm2();

	double b1=(c2z(e(p))-c2z(p)).norm2();
	double b2=(c2z(c2z(p))-e(p)).norm2();
	double b3=(c2z(sxz(p))-syz(p)).norm2();
	double b4=(c2z(syz(p))-sxz(p)).norm2();

	double c1=(sxz(e(p))-sxz(p)).norm2();
	double c2=(sxz(c2z(p))-syz(p)).norm2();
	double c3=(sxz(sxz(p))-e(p)).norm2();
	double c4=(sxz(syz(p))-c2z(p)).norm2();

	double d1=(syz(e(p))-syz(p)).norm2();
	double d2=(syz(c2z(p))-sxz(p)).norm2();
	double d3=(syz(sxz(p))-c2z(p)).norm2();
	double d4=(syz(syz(p))-e(p)).norm2();


	int result=0;
	if (a1+a2+a3+a4+b1+b2+b3+b4+c1+c2+c3+c4+d1+d2+d3+d4 > 1.e-13) {
		print("large error norm in the c2v multiplication table");
		print(a1,a2,a3,a4);
		print(b1,b2,b3,b4);
		print(c1,c2,c3,c4);
		print(d1,d2,d3,d4);
		result=1;
	} else {
		print("  .. all good");
	}
	return result;
}


/// symmetrize a trial function, check its behavior according to the group table
int test_projector(World& world) {

	// the trial function
	const real_function_3d f=real_factory_3d(world).f(gaussian_shift_3d);

	double error=0.0;
	std::string all_pg[]={"c1","cs","c2","ci","c2v","c2h","d2","d2h"};
	for (const std::string& pg : all_pg) {
		print("point group",pg);
		projector_irrep proj(pg);

		for (const std::string& irrep : proj.get_all_irreps()) {  // loop over all irreps
			print(" irrep",irrep);
			proj.set_irrep(irrep);
			real_function_3d f1=proj(f)[0];	// result is the first element of result vector

			charactertable table=proj.get_table();
			for (int i=0; i<table.order_; ++i) {	// loop over all symmetry operations
				const pg_operator syop = table.operators_[i];
				const int character=table.irreps_[irrep][i];
				double n1=(f1 - character * syop(f1)).norm2();
				print("  operator, character, norm ",syop.name(),character, n1);
				error+=n1;
			}
		}
	}

	int result=0;
	if (error > 1.e-13) {
		print("large error norm test_projector:", error);
		result=1;
	} else {
		print("  .. all good");
	}
	return result;
}

int test_orthogonalization(World& world) {

	int result=0;

	pg_operator e=pg_identity();
	pg_operator i=pg_inversion();
	pg_operator sxy=pg_sigma_xy();
	pg_operator sxz=pg_sigma_xz();
	pg_operator syz=pg_sigma_yz();
	pg_operator c4x=pg_c4x();
	pg_operator c4y=pg_c4y();
	pg_operator c4z=pg_c4z();
	pg_operator c2x=pg_c2x();
	pg_operator c2y=pg_c2y();
	pg_operator c2z=pg_c2z();
	std::vector<std::string> molecular_info(1,"0 0 0.0 0.0 0.0\n");

	// the trial function
	real_function_3d f=real_factory_3d(world).f(gaussian_shift_3d);
	double fnorm=f.norm2();
	f.scale(1.0/fnorm);
	const real_function_3d one=real_factory_3d(world).f(one1);

	// create a set of symmetry-equivalent functions
	std::string all_pg[]={"c1","cs","c2","ci","c2v","c2h","d2","d2h"};
	for (const std::string& pg : all_pg) {
		print("point group",pg);
		projector_irrep proj(pg);
		proj.set_verbosity(1);
		charactertable table=proj.get_table();

		// create a set of symmetry-adapted functions from one function
		std::vector<std::string> sirrep1;
		vector_real_function_3d sym_adapt1=proj.create_symmetry_adapted_basis(f,sirrep1);
		                if (sym_adapt1.size()!= (size_t) table.order_) result++;
		// create a set of symmetry-equivalent functions
		vector_real_function_3d sym_equiv;
		for (auto& syop : table.operators_) sym_equiv.push_back(syop(f));
		print("sym_equiv.size()",sym_equiv.size());

		// create a set of symmetry-adapted functions from the set of symm-equiv functions
		std::vector<std::string> sirrep;
		vector_real_function_3d sym_adapt=proj(sym_equiv,one,sirrep);
		print("sym_adapt.size()",sym_adapt.size());
		if (sym_adapt.size()!=sym_equiv.size()) result++;
		print("sirrep",sirrep);

		// double check the two sets of functions
		Tensor<double> ovlp=matrix_inner(world,sym_adapt,sym_adapt);

		for (int i=0; i<ovlp.dim(0); ++i) ovlp(i,i)-=1.0;
		double norm=ovlp.normf()/ovlp.size();
		if (norm>1.e-5) {
			result++;
			print(ovlp);
			print("norm of the overlap matrix (should be zero)",norm);
		}
	}

	double error=0;
	if (error > 1.e-13) {
		print("large error norm test_orthogonalization:", error);
		result=1;
	} else {
		print("  .. all good");
	}
	return result;
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    FunctionDefaults<2>::set_cubic_cell(-6.0,6.0);
    FunctionDefaults<3>::set_cubic_cell(-6.0,6.0);


    int result=0;
    result+=check_operator_multiplications_2d(world);
    result+=check_operator_multiplications_3d(world);
    result+=check_multiplication_table_c2v(world);
    result+=test_projector(world);
    result+=test_orthogonalization(world);
//    plot_symmetry_operators(world);

    print("result",result);
    madness::finalize();
    return result;
}
