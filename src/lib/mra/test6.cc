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


  $Id$
*/

/// \file test6.cc
/// \brief test various functionality for 6d functions

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include "mra/mra.h"

using namespace madness;

bool is_small(const double& val, const double& eps) {
	return (val<eps);
}

std::string ok(const bool b) {if (b) return "ok   "; return "fail ";};

bool is_large(const double& val, const double& eps) {
	return (val>eps);
}

static double gauss_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= sqrt(x*x + y*y + z*z);
    const double norm=0.712705695388313;
    return norm*exp(-r2);
}

static double gauss_6d(const coord_6d& r) {
    coord_3d r1, r2;
    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
    return gauss_3d(r1)*gauss_3d(r2);
}

static double r2r(const coord_6d& r) {
    coord_3d r1, r2;
    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
    double g1=gauss_3d(r1);
    return g1*g1*gauss_3d(r2);
}

static double add_test(const coord_6d& r) {
    coord_3d r1, r2;
    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
    double g1=gauss_3d(r1);
    double g2=gauss_3d(r2);

    return g1*g2 + g1*g1*g2;
}

/// test f(1,2) = g(1) h(2)
int test_hartree_product(World& world, const long& k, const double thresh) {

	print("entering hartree_product");
	int nerror=0;
	bool good;

    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    real_function_3d phisq=phi*phi;

    {
        real_function_6d ij=hartree_product(phi,phi);
        ij.truncate();

        double norm=ij.norm2();
        print("norm(ij)",norm);

        double err=ij.err(gauss_6d);
        good=is_small(err,thresh);
        print(ok(good), "hartree_product(phi,phi) error:",err);
        if (not good) nerror++;

    }

    {
        real_function_6d ij=hartree_product(phisq,phi);
        double err=ij.err(r2r);
        good=is_small(err,thresh);
        print(ok(good), "hartree_product(phi^2,phi) error:",err);
        if (not good) nerror++;
    }

	print("all done\n");
	return nerror;
}

/// test f(1,2)*g(1)
int test_multiply(World& world, const long& k, const double thresh) {

    print("entering multiply f(1,2)*g(1)");
    int nerror=0;
    bool good;

    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    real_function_3d phisq=phi*phi;

    real_function_6d ij=hartree_product(phi,phi);
    real_function_6d iij2=multiply(ij,phi,1);

    double err=iij2.err(r2r);
    good=is_small(err,thresh);
    print(ok(good), "multiply f(1,2)*g(1) error:",err);
    if (not good) nerror++;

    print("all done\n");
    return nerror;
}

/// test f(1,2) + g(1,2) for both f and g reconstructed
int test_add(World& world, const long& k, const double thresh) {

    print("entering add f(1,2)*g(1)");
    int nerror=0;
    bool good;

    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    real_function_3d phisq=phi*phi;

    real_function_6d f=hartree_product(phi,phi);
    real_function_6d g=hartree_product(phisq,phi);

    // this is in reconstructed form
    real_function_6d h1=f+g;
    double err1=h1.err(add_test);
    good=is_small(err1,thresh);
    print(ok(good), "add error1:",err1);
    if (not good) nerror++;

    // more explicitly;
    real_function_6d h2=gaxpy_oop_reconstructed(1.0,f,1.0,g,true);
    double err2=h2.err(add_test);
    good=is_small(err2,thresh);
    print(ok(good), "add error2:",err2);
    if (not good) nerror++;


    print("all done\n");
    return nerror;
}


/// test f(1,2) + g(1,2) for both f and g reconstructed
int test_exchange(World& world, const long& k, const double thresh) {

    print("entering exchange f(1,2)*g(1)");
    int nerror=0;
    bool good;

    double norm;
    real_function_3d phi=real_factory_3d(world).f(gauss_3d);


    real_convolution_3d poisson = CoulombOperator(world,0.0001,thresh);

    real_function_3d rho = 2.0*phi*phi;
    real_function_3d coulombpot=poisson(rho);


    real_function_6d f=2.0*hartree_product(phi,phi);
    real_function_6d f2=multiply(f,phi,1);
    f2.print_size("f2 after apply");
    norm=f2.norm2();
    if (world.rank()==0) print("f2 norm",norm);
    real_function_6d x=poisson(f2);

    x.print_size("x after apply");
    norm=x.norm2();
    if (world.rank()==0) print("x norm",norm);
    x=multiply(x,phi,1);
    x.print_size("x after multiply");
    norm=x.norm2();
    if (world.rank()==0) print("x norm",norm);


    real_function_6d tmp=0.5*multiply(f,coulombpot,1);
    tmp.print_size("tmp after multiply");
    norm=tmp.norm2();
    if (world.rank()==0) print("tmp norm",norm);

    real_function_6d diff=tmp-x;
    diff.print_size("diff");
    norm=diff.norm2();
    if (world.rank()==0) print("diff norm",norm);
    good=is_small(norm,thresh);
    print(ok(good), "exchange error:",norm);
    if (not good) nerror++;

    // do only orbital
    real_function_3d tmp2=phi*coulombpot;
    tmp2.print_size("J phi after multiply");
    norm=tmp2.norm2()*phi.norm2();
    if (world.rank()==0) print("J phi norm",norm);



    print("all done\n");
    return nerror;
}

/// test inner product using redundant wave functions
int test_inner(World& world, const long& k, const double thresh) {

    print("entering inner");
    int nerror=0;
    bool good;

    double norm;
    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    real_function_3d phi2=phi*phi;

    real_function_6d f1=hartree_product(phi,phi);
    real_function_6d f2=hartree_product(phi,phi2);

    double a1=inner(f1,f2);
    double a2=inner(phi,phi) * inner(phi,phi2);
    norm=a1-a2;

    if (world.rank()==0) print("diff norm",norm);
    good=is_small(norm,thresh);
    print(ok(good), "inner error:",norm);
    if (not good) nerror++;

    print("all done\n");
    return nerror;
}



int test(World& world, const long& k, const double thresh) {

    print("entering test");
    int nerror=0;


    double norm;
    real_function_3d phi=real_factory_3d(world).f(gauss_3d);


    real_convolution_3d poisson = CoulombOperator(world,0.0001,thresh);
    poisson.modified()=false;

    real_function_3d rho = 2.0*phi*phi;
    real_function_3d coulombpot=poisson(rho);
    norm=coulombpot.norm2();
    print("coulombpot",norm);
    coulombpot.print_size("coulombpot");


    real_function_6d f=2.0*hartree_product(phi,phi);
    real_function_6d f2=multiply(f,phi,1);
    f2.print_size("f2 after apply");
    norm=f2.norm2();
    if (world.rank()==0) print("f2 norm",norm);
    real_function_6d x=poisson(f2);

    x.print_size("x after apply");
    norm=x.norm2();
    if (world.rank()==0) print("x norm",norm);
    x=multiply(x,phi,1);
    x.print_size("x after multiply");
    norm=x.norm2();
    if (world.rank()==0) print("x norm",norm);


    real_function_6d tmp=0.5*multiply(f,coulombpot,1);
    tmp.print_size("tmp after multiply");
    norm=tmp.norm2();
    if (world.rank()==0) print("tmp norm",norm);

    real_function_6d diff=tmp-x;
    diff.print_size("diff");
    norm=diff.norm2();
    if (world.rank()==0) print("diff norm",norm);

    // do only orbital
    real_function_3d tmp2=phi*coulombpot;
    tmp2.print_size("J phi after multiply");
    norm=tmp2.norm2()*phi.norm2();
    if (world.rank()==0) print("J phi norm",norm);



    print("all done\n");
    return nerror;
}


int main(int argc, char**argv) {

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    srand(time(NULL));
    startup(world,argc,argv);

    // the parameters
    long k=6;
    double thresh=1.e-3;
    double L=16;
    TensorType tt=TT_2D;

    // override the default parameters
    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="size") L=atof(val.c_str());               // usage: size=10
        if (key=="k") k=atoi(val.c_str());                  // usage: k=5
        if (key=="thresh") thresh=atof(val.c_str());        // usage: thresh=1.e-3
        if (key=="TT") {
            if (val=="TT_2D") tt=TT_2D;
            else if (val=="TT_3D") tt=TT_3D;
            else if (val=="TT_FULL") tt=TT_FULL;
            else {
                print("arg",arg, "key",key,"val",val);
                MADNESS_EXCEPTION("confused tensor type",0);
            }

        }
    }

    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_tensor_type(tt);

    print("entering testsuite for 6-dimensional functions\n");
    print("k            ",k);
    print("thresh       ",thresh);
    print("boxsize      ",L);
    print("tensor type: ", FunctionDefaults<6>::get_tensor_type());
    print("");


    int error=0;

    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    double norm=phi.norm2();
    if (world.rank()==0) printf("phi.norm2()   %12.8f\n",norm);

    real_function_3d phi2=2.0*phi*phi;
    norm=phi2.norm2();
    if (world.rank()==0) printf("phi2.norm2()  %12.8f\n",norm);

//    test(world,k,thresh);
    error+=test_hartree_product(world,k,thresh);
//    error+=test_multiply(world,k,thresh);
//    error+=test_add(world,k,thresh);
//    error+=test_exchange(world,k,thresh);
    error+=test_inner(world,k,thresh);


    print(ok(error==0),error,"finished test suite\n");

    world.gop.fence();
    finalize();

    return 0;
}


