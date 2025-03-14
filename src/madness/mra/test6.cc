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

#include <madness/mra/mra.h>
#include <madness/world/test_utilities.h>

using namespace madness;

std::string ok(const bool b) {if (b) return "ok   "; return "fail ";};

int check_small(const double val, const double eps, const std::string message) {
	bool is_small=(fabs(val)<eps);
	print(ok(is_small),val,message);
	return (is_small) ? 0 : 1;
}

int check(bool b, const std::string message) {
	print(ok(b),message);
	return (b) ? 0 : 1;
}

bool is_small(const double& val, const double& eps) {
	return (val<eps);
}


bool is_large(const double& val, const double& eps) {
	return (val>eps);
}

template<size_t NDIM>
void load_function(World& world, Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0)  print("loading function ", name);

    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, name.c_str());
    ar & pair;

    FunctionDefaults<3>::set_k(pair.k());
    FunctionDefaults<6>::set_k(pair.k());

    FunctionDefaults<3>::set_thresh(pair.thresh());
    FunctionDefaults<6>::set_thresh(pair.thresh());

    std::string line="loaded function "+name;
    pair.print_size(line);

}

template<size_t NDIM>
void save_function(World& world, Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0)  print("saving function ", name);

    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str());
    ar & pair;

    std::string line="saved function "+name;
    pair.print_size(line);

}


static double one_3d(const coord_3d& r) {
	return 1.0;
}

static double zero_3d(const coord_3d& r) {
	return 0.0;
}

static double one_6d(const coord_6d& r) {
	return 1.0;
}

static double gauss_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= sqrt(x*x + y*y + z*z);
    const double norm=0.712705695388313;
    return norm*exp(-r2);
}

static double tightgauss_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= sqrt(x*x + y*y + z*z);
    const double norm=0.712705695388313;
    return norm*exp(-2.0*r2);
}

static double gauss_plus_one_3d(const coord_3d& r) {
    return gauss_3d(r)+one_3d(r);
}
static double gauss_plus_tight_3d(const coord_3d& r) {
    return gauss_3d(r)+tightgauss_3d(r);
}


// static double gauss_6d(const coord_6d& r) {
//     coord_3d r1, r2;
//     r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
//     r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
//     return gauss_3d(r1)*gauss_3d(r2);
// }
//
//
// static double r2r(const coord_6d& r) {
//     coord_3d r1, r2;
//     r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
//     r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
//     double g1=gauss_3d(r1);
//     return g1*g1*gauss_3d(r2);
// }

//static double add_test(const coord_6d& r) {
//    coord_3d r1, r2;
//    r1[0]=r[0],    r1[1]=r[1],    r1[2]=r[2];
//    r2[0]=r[3],    r2[1]=r[4],    r2[2]=r[5];
//    double g1=gauss_3d(r1);
//    double g2=gauss_3d(r2);
//
//    return g1*g2 + g1*g1*g2;
//}

static double V(const Vector<double,3>& r) {
  return -1.0/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-12);
}

/// test f(1,2) = g(1) h(2)
template<typename T, std::size_t NDIM>
int test_hartree_product(World& world, const long& k, const double thresh) {
    test_output t1(std::string("testing hartree_product for dimension "+std::to_string(NDIM)));
//    t1.set_cout_to_terminal();
    constexpr std::size_t LDIM=NDIM/2;
    static_assert(LDIM*2==NDIM);
    double loose=1.1*std::pow(2.0,1.5*LDIM);
    print("loosening factor",loose);
    FunctionDefaults<NDIM>::set_tensor_type(TT_2D);
    FunctionDefaults<LDIM>::set_tensor_type(TT_FULL);

    {
        Function<T,LDIM> phi=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r) {return exp(-inner(r,r));});
        auto reference=[](const Vector<double,NDIM>& r) {return exp(-inner(r,r));};

        Function<T,NDIM> ij=hartree_product(phi,phi);
        double err1=ij.err(reference);
        ij.truncate();
        double err2=ij.err(reference);
        print("err1, err2",err1,err2);
        t1.checkpoint(is_small(err2,thresh*loose),"hartree_product(gauss,gauss) error:");
    }

    {
        Function<T,LDIM> phi=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r) {return exp(-sqrt(inner(r,r)+1.e-1));});
        Function<T,NDIM> ij=hartree_product(phi,phi);
        auto reference=[](const Vector<double,NDIM>& r) {
            double r1=0.0, r2=0.0;
           	for (std::size_t i=0; i<LDIM; ++i) {
                 r1+=r[i]*r[i];
                 r2+=r[i+LDIM]*r[i+LDIM];
            }
            return exp(-sqrt(r1+1.e-1)-sqrt(r2+1.e-1));
        };
        double err=ij.err(reference);
        print("err",err);
        t1.checkpoint(is_small(err,thresh*loose*3.0),"hartree_product(slater,slater) error:");
    }

	print("all done\n");
    return t1.end();
}

/// test f(1,2)*g(1)
template<typename T, std::size_t NDIM>
int test_multiply(World& world, const long& k, double thresh) {
	constexpr std::size_t LDIM=NDIM/2;
	FunctionDefaults<NDIM>::set_tensor_type(TT_2D);
	FunctionDefaults<LDIM>::set_tensor_type(TT_FULL);
	thresh*=50;

    test_output t1(std::string("testing multiply f(1,2)*g(1) for dimension "+std::to_string(NDIM)));
    t1.set_cout_to_terminal();

	/// Slater factor
    Function<T,NDIM> f12=FunctionFactory<T,NDIM>(world).is_on_demand().f([](const Vector<double,NDIM>& r) {
        Vector<double,LDIM> r1, r2;
        for (std::size_t i=0; i<LDIM; ++i) {
            r1[i]=r[i];
            r2[i]=r[i+LDIM];
        }
        return exp(- (r1-r2).normf());
    });


    const Function<T,LDIM> phi=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r){return exp(-r.normf());});
    const Function<T,LDIM> phisq=phi*phi;

    Function<T,NDIM> fii=CompositeFactory<T,NDIM,LDIM>(world)
    	    	.particle1(copy(phi))
    	    	.particle2(copy(phi))
    	    	.g12(f12);
    Function<T,NDIM> fi2i=CompositeFactory<T,NDIM,LDIM>(world)
    	    	.particle1(copy(phisq))
    	    	.particle2(copy(phi))
    	    	.g12(f12);

    Function<T,NDIM> fii2=CompositeFactory<T,NDIM,LDIM>(world)
                .particle1(copy(phi))
                .particle2(copy(phisq))
                .g12(f12);


    fii.fill_tree();
    fi2i.fill_tree();
    fii2.fill_tree();

    fii.print_size("f12 |phi phi>");
    fi2i.print_size("f12 |phi^2 phi>");
    fii2.print_size("f12 |phi phi^2>");


 	for (int particle=1; particle<3; ++particle) {
 		Function<T,NDIM> iij1=multiply(copy(fii),phi,particle);
 		std::string sp=std::to_string(particle);
     	iij1.print_size("multiply");
     	iij1.truncate().reduce_rank();
     	iij1.print_size("multiply truncated");
     	double err=(fi2i-iij1).norm2();
 		t1.checkpoint(err,thresh,"multiply f(1,2)*g("+sp+"1) error:");
 	}

 	for (int particle=1; particle<3; ++particle) {
 		std::string sp=std::to_string(particle);
 		Function<T,NDIM> iij3;
 		if (particle==1) iij3=CompositeFactory<T,NDIM,LDIM>(world)
 					.ket(copy(fii)).V_for_particle1(copy(phi));
 		if (particle==2) iij3=CompositeFactory<T,NDIM,LDIM>(world)
 					.ket(copy(fii)).V_for_particle2(copy(phi));

     	iij3.fill_tree();
     	iij3.print_size("CompositeFactory");
     	iij3.truncate();
     	iij3.print_size("CompositeFactory truncated");
 		double err=(fi2i-iij3).norm2();

 		t1.checkpoint(err,thresh,"CompositeFactory (1,2)*g("+sp+") error:");
 	}

	return t1.end();
}

/// test f(1,2) + g(1,2) for both f and g reconstructed
int test_add(World& world, const long& k, const double thresh) {

    print("entering add");
    int nerror=0;

    // simple test for 3d functions and different adding schemes
    real_function_3d one3=real_factory_3d(world).f(one_3d);
    real_function_3d gauss3=real_factory_3d(world).f(gauss_3d);
    real_function_3d gauss_plus_one3=real_factory_3d(world).f(gauss_plus_one_3d);

    {
    	real_function_3d r1=one3+gauss3;
    	double error1=r1.err(gauss_plus_one_3d);
    	nerror+=check_small(error1,thresh,"operator+");
    }
    {
    	real_function_3d r1=one3+gauss3;	// this has been checked before
    	real_function_3d r2=r1-gauss_plus_one3;
    	double error2=r2.err(zero_3d);
    	nerror+=check_small(error2,thresh,"operator-");
    }
    {
    	one3.compress(); gauss3.compress();
    	real_function_3d r3=gaxpy_oop(1.0,one3,1.0,gauss3);
    	nerror+=check(r3.is_compressed(),"is compressed");
    	double error3=r3.err(gauss_plus_one_3d);
    	nerror+=check_small(error3,thresh,"gaxpy_oop add");
    }
    {
    	one3.reconstruct(); gauss3.reconstruct();
    	real_function_3d r4=gaxpy_oop_reconstructed(1.0,one3,1.0,gauss3);
    	nerror+=check(!r4.is_compressed(),"is reconstructed");
    	double error4=r4.err(gauss_plus_one_3d);
    	nerror+=check_small(error4,thresh,"gaxpy_oop_reconstructed");
    }
    {
    	real_function_3d r=copy(one3);
    	r.compress(); gauss3.compress();
    	r+=gauss3;
    	nerror+=check(r.is_compressed(),"is reconstructed");
    	double error1=r.err(gauss_plus_one_3d);
    	nerror+=check_small(error1,thresh,"operator+=, compressed");
    }
    {
    	real_function_3d r=copy(one3);
    	r.reconstruct(); gauss3.reconstruct();
    	r+=gauss3;
    	nerror+=check(r.is_compressed(),"is reconstructed");
    	double error1=r.err(gauss_plus_one_3d);
    	nerror+=check_small(error1,thresh,"operator+=, reconstructed");
    }
    {
    	one3.reconstruct(); gauss3.reconstruct();
    	real_function_3d r=gaxpy_oop_reconstructed(1.0,one3,-1.0,gauss3);
    	nerror+=check(!r.is_compressed(),"is reconstructed");
    	real_function_3d r2=one3-gauss3;
    	double error=(r-r2).norm2();
    	nerror+=check_small(error,thresh,"gaxpy_oop_reconstructed subtract");
    }
    {
    	real_function_3d r=copy(gauss3);
    	r.reconstruct(); gauss3.reconstruct();
    	r.add_scalar(1.0);
    	nerror+=check(!r.is_compressed(),"is reconstructed");
    	double error1=r.err(gauss_plus_one_3d);
    	nerror+=check_small(error1,thresh,"add_scalar");
    }


    // 6d tests; mainly consistency checks since err() function is very expensive in 6d
    real_function_6d f=hartree_product(gauss3,gauss3);
    real_function_6d one6=real_factory_6d(world).f(one_6d);

    {
    	real_function_6d r=copy(f);
    	r.add_scalar(1.0);
    	real_function_6d r2=one6+f;
    	double error=(r2-r).norm2();
    	nerror+=check_small(error,thresh,"6d add_scalar/operator+/-");
    }
    {
        real_function_3d tightgauss3=real_factory_3d(world).f(tightgauss_3d);
        real_function_3d gauss_plus_tight3=real_factory_3d(world).f(gauss_plus_tight_3d);

    	real_function_6d r=hartree_product(gauss_plus_tight3,gauss_plus_tight3);
    	real_function_6d r1=hartree_product(gauss3,gauss3);
    	real_function_6d r2=hartree_product(gauss3,tightgauss3);
    	real_function_6d r22=hartree_product(tightgauss3,gauss3);
    	r22+=r2;
    	r22.scale(0.5);
    	real_function_6d r3=hartree_product(tightgauss3,tightgauss3);
    	real_function_6d r4=gaxpy_oop_reconstructed(1.0,r,-2.0,r22)-r1-r3;
    	double error=r4.norm2();
    	nerror+=check_small(error,1.5*thresh,"6d gaxpy_oop_reconstructed/operator+=/operator- note loosened threshold");
    }

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


/// test 6D convolution
int test_convolution(World& world, const long& k, const double thresh) {

    print("entering convolution");
    int nerror=0;
    bool good;

    double eps=-0.5;

    // solve the 3D H-atom
    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    const real_function_3d v=real_factory_3d(world).f(V);
	real_function_3d vphi=v*phi;

    v.print_size("v");
//    real_convolution_3d poisson = CoulombOperator(world,0.0001,thresh);
    real_convolution_3d green = BSHOperator<3>(world, sqrt(-2.0*eps), 1.e-8, 1.e-6);

    for (int i=0; i<10; ++i) {

    	vphi=v*phi;
    	double PE=inner(vphi,phi);
    	print("<phi | V | phi>: ",PE);

    	real_derivative_3d Dx = free_space_derivative<double,3>(world,0);
    	Function<double,3> du = Dx(phi);
    	double KE = 3*0.5*(du.inner(du));
    	print("<phi | T | phi>: ",KE);
    	print("<phi | H | phi>: ",KE + PE);


    	phi=green(-2.0*vphi).truncate();
    	double norm=phi.norm2();
    	phi.scale(1.0/norm);
    	print("phi.norm2()",norm);
    }

    vphi=v*phi;
    double PE=inner(vphi,phi);
   	print("<phi | V | phi>: ",PE);

    // solve the H-atom in 6D
    real_convolution_6d green6 = BSHOperator<6>(world, sqrt(-2.0*(eps+eps)), 1.e-8, 1.e-6);

	real_function_6d result1=-2.0*green6(copy(vphi),copy(phi)).truncate().reduce_rank();
	result1=result1-2.0*green6(copy(phi),copy(vphi)).truncate().reduce_rank();
	world.gop.fence();

	double a=result1.norm2();
	if (world.rank()==0) print("<GVphi | GVphi> ",a);

	real_function_6d diff=result1-hartree_product(phi,phi);
	double norm=diff.norm2();

    if (world.rank()==0) print("diff norm",norm);
    good=is_small(norm,thresh);
    print(ok(good), "inner error:",norm);
    if (not good) nerror++;

    print("all done\n");
    return nerror;
}


int test_replicate(World& world, const long& k, const double thresh) {
    real_function_3d phi=real_factory_3d(world).f(gauss_3d);
    auto map=phi.get_pmap();
    map->print_data_sizes(world,"before replication");
    phi.replicate();
    phi.get_pmap()->print_data_sizes(world,"replicated");
    map->print_data_sizes(world,"after replication");
    phi.distribute(map);
    map->print_data_sizes(world,"after distribution");
    return 0;
}


int test(World& world, const long& k, const double thresh) {

    print("entering test");
    int nerror=0;

    typedef Key<6> keyT;

    real_function_3d phi=real_factory_3d(world).f(gauss_3d);

    real_function_6d ij=hartree_product(phi,phi);
//    real_function_3d ij=phi;
    ij.compress();

    // get the root NS coeffs
    keyT key0=ij.get_impl()->get_cdata().key0;
    GenTensor<double> NScoeff=(ij.get_impl()->get_coeffs().find(key0)).get()->second.coeff();

    {
		// convert NS coeffs to values directly
		GenTensor<double> val1=ij.get_impl()->NScoeffs2values(key0,NScoeff,false);

		// convert NS coeffs to S coeffs, and then to values
		Tensor<double> Scoeff=ij.get_impl()->unfilter(NScoeff).full_tensor_copy();
		Tensor<double> val2(ij.get_impl()->get_cdata().v2k);

		for (KeyChildIterator<6> kit(key0); kit; ++kit) {
			const keyT& child = kit.key();
			std::vector<Slice> cp = ij.get_impl()->child_patch(child);
			Tensor<double> child_s_coeff=Scoeff(cp);
			val2(cp)=ij.get_impl()->coeffs2values(child,child_s_coeff);
		}

		Tensor<double> diff=val2-val1.full_tensor_copy();
		double error=diff.normf();
		print("error in NScoeff2values",error);
    }

    {
		// convert S coeffs to values directly
		const std::vector<Slice>& s0=ij.get_impl()->get_cdata().s0;
		GenTensor<double> val1=ij.get_impl()->NScoeffs2values(key0,GenTensor<double>(NScoeff(s0)),true);

		// convert NS coeffs to S coeffs, and then to values
		Tensor<double> Scoeff(ij.get_impl()->get_cdata().v2k);
		Scoeff(s0)=(NScoeff.full_tensor_copy()(s0));
		Scoeff=ij.get_impl()->unfilter(Scoeff);
		Tensor<double> val2(ij.get_impl()->get_cdata().v2k);

		for (KeyChildIterator<6> kit(key0); kit; ++kit) {
			const keyT& child = kit.key();
			std::vector<Slice> cp = ij.get_impl()->child_patch(child);
			Tensor<double> child_s_coeff=Scoeff(cp);
			val2(cp)=ij.get_impl()->coeffs2values(child,child_s_coeff);
		}

		Tensor<double> diff=val2-val1.full_tensor_copy();
		double error=diff.normf();
		print("error in Scoeff2values",error);
    }

    {
		// convert S coeffs to values directly
		const std::vector<Slice>& s0=ij.get_impl()->get_cdata().s0;
		GenTensor<double> val1=ij.get_impl()->NS_fcube_for_mul(key0,key0,GenTensor<double>(NScoeff(s0)),true);

		// convert NS coeffs to S coeffs, and then to values
		Tensor<double> Scoeff(ij.get_impl()->get_cdata().v2k);
		Scoeff(s0)=(NScoeff.full_tensor_copy()(s0));
		Scoeff=ij.get_impl()->unfilter(Scoeff);
		Tensor<double> val2(ij.get_impl()->get_cdata().v2k);

		for (KeyChildIterator<6> kit(key0); kit; ++kit) {
			const keyT& child = kit.key();
			std::vector<Slice> cp = ij.get_impl()->child_patch(child);
			Tensor<double> child_s_coeff=Scoeff(cp);
			val2(cp)=ij.get_impl()->coeffs2values(child,child_s_coeff);
		}

		Tensor<double> diff=val2-val1.full_tensor_copy();
		double error=diff.normf();
		print("error in NS_fcube_for_mul",error);
    }

    {
		// convert NS coeffs to values directly
		GenTensor<double> val1=ij.get_impl()->NScoeffs2values(key0,NScoeff,false);
		GenTensor<double> coeff1=ij.get_impl()->values2NScoeffs(key0,val1);
		GenTensor<double> val2=ij.get_impl()->NScoeffs2values(key0,coeff1,false);

		Tensor<double> diff=val2.full_tensor_copy()-val1.full_tensor_copy();
		double error=diff.normf();
		print("error in values2NScoeffs(NScoeff2values)",error);
    }

    print("all done\n");
    return nerror;
}

template<typename T, std::size_t NDIM>
int test_vector_composite(World& world, const long& k, const double thresh) {
	constexpr std::size_t LDIM=NDIM/2;
	MADNESS_CHECK(NDIM==2*LDIM);

	FunctionDefaults<NDIM>::set_tensor_type(TT_2D);

	test_output t1("test_vector_composite for dimension "+std::to_string(NDIM));
//	t1.set_cout_to_terminal();

	Function<T,LDIM> l1_func=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r) {return exp(-inner(r,r));});
	Function<T,LDIM> l2_func=FunctionFactory<T,LDIM>(world).f([](const Vector<double,LDIM>& r) {return exp(-2.0*inner(r,r));});
	Function<T,NDIM> h1_func=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r) {return exp(-inner(r,r));});
	Function<T,NDIM> h2_func=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r) {return exp(-2.0*inner(r,r));});
	Function<T,NDIM> h12_func=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r) {return exp(-inner(r,r)) + exp(-2.0*inner(r,r));});
    Function<T,NDIM> h1_eri_func=FunctionFactory<T,NDIM>(world).f([](const Vector<double,NDIM>& r) {
        Vector<double,LDIM> r1, r2;
        for (std::size_t i=0; i<LDIM; ++i) {
             r1[i]=r[i];
             r2[i]=r[i+LDIM];
        }
        return exp(-inner(r1,r1) - inner(r2,r2)- (r1-r2).normf());
    }).thresh(thresh*0.01); // high accuracy for the reference number

    Function<T,NDIM> eri=FunctionFactory<T,NDIM>(world).is_on_demand().f([](const Vector<double,NDIM>& r) {
        Vector<double,LDIM> r1, r2;
        for (std::size_t i=0; i<LDIM; ++i) {
            r1[i]=r[i];
            r2[i]=r[i+LDIM];
        }
        return exp(- (r1-r2).normf());
    });

    // the analytical norm of a 1D Gauss functions with f(x) = exp(-alpha x^2) : ||f(x)||_2= (pi/(2 alpha))^1/4
    double nh_analytical=std::pow(M_PI/2.0,0.25*NDIM);
    double nl=l1_func.norm2();
    double nh=h1_func.norm2();
    print("l1/h1_func.norm2(),n_analytical",nl,nh,nh_analytical);
	/// test construction with a single pair of functions
	{
		Function<T,NDIM> h1_1=CompositeFactory<T,NDIM,LDIM>(world) .ket(h1_func);
		Function<T,NDIM> h1_2=CompositeFactory<T,NDIM,LDIM>(world).particle1(l1_func).particle2(l1_func);
		h1_1.fill_tree();
		h1_2.fill_tree();
		double n1=h1_1.norm2();
		double n2=h1_2.norm2();
        print("n1, n2, nh", n1, n2, nh);
		t1.checkpoint(std::abs(n1-nh)<thresh,"construction with a single pair -- direct");
		t1.checkpoint(std::abs(n2-nh)<thresh,"construction with a single pair -- direct");
		t1.checkpoint(std::abs(n1-n2)<thresh,"construction with a single pair -- consistent");
	}

    /// test construction with a single pair of functions and eri
    {
        Function<T,NDIM> h1_1=CompositeFactory<T,NDIM,LDIM>(world).ket(h1_func).g12(eri);
        Function<T,NDIM> h1_2=CompositeFactory<T,NDIM,LDIM>(world).g12(eri).particle1(l1_func).particle2(l1_func);
        h1_1.fill_cuspy_tree();
        h1_2.fill_cuspy_tree();
        double n1=h1_1.norm2();
        double n2=h1_2.norm2();
        double n3=h1_eri_func.norm2();
        print("n1, n2, nh", n1, n2, n3);
        t1.checkpoint(std::abs(n1-n3)<thresh*30,"construction with a single pair and eri -- direct REDUCED ACC");
        t1.checkpoint(std::abs(n2-n3)<thresh*30,"construction with a single pair and eri -- direct REDUCED ACC");
        t1.checkpoint(std::abs(n1-n2)<thresh*30,"construction with a single pair and eri -- consistent REDUCED ACC");
    }

	/// test construction with a vector of functions:
	{
		Function<T,NDIM> h1_1=CompositeFactory<T,NDIM,LDIM>(world).ket(h12_func);   // reference
		Function<T,NDIM> h1_2=CompositeFactory<T,NDIM,LDIM>(world).particle1({l1_func,l2_func}).particle2({l1_func,l2_func});
		Function<T,NDIM> h1_3=CompositeFactory<T,NDIM,LDIM>(world).ket({h1_func,h2_func});
		h1_1.fill_tree();
		h1_2.fill_tree();
		h1_3.fill_tree();
		double n1=h1_1.norm2();
		double n2=(h1_2).norm2();
		double n3=(h1_3).norm2();
		h1_2+=h1_3;
		double n12=(h1_2).norm2();
		print("n1, n2, n3, n12", n1, n2, n3, n12);
		t1.checkpoint(std::abs(n1-n2)<thresh,"construction with a vector of functions -- particles");
        t1.checkpoint(std::abs(n1-n3)<thresh,"construction with a vector of functions -- ket");
        t1.checkpoint(std::abs(n2-n3)<thresh,"construction with a vector of functions -- consistent");
	}

    Function<T,NDIM> h1_1=CompositeFactory<T,NDIM,LDIM>(world) .ket(copy(h1_func));




	return t1.end();
}


int main(int argc, char**argv) {

//    initialize(argc,argv);
//    World world(SafeMPI::COMM_WORLD);
    World& world= initialize(argc,argv);
    srand(time(nullptr));
    startup(world,argc,argv);

    // the parameters
    long k=7;
    double thresh=1.e-5;
    double L=40;
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
            else if (val=="TT_FULL") tt=TT_FULL;
            else if (val=="TT_TENSORTRAIN") tt=TT_TENSORTRAIN;
            else {
                print("arg",arg, "key",key,"val",val);
                MADNESS_EXCEPTION("confused tensor type",0);
            }

        }
    }

    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<1>::set_truncate_mode(1);
    FunctionDefaults<2>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<4>::set_truncate_mode(1);

    FunctionDefaults<1>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<2>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<4>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<5>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);

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

	error+=test_vector_composite<double,2>(world,k,thresh);
//    test(world,k,thresh);
    error+=test_hartree_product<double,2>(world,k,thresh);
    error+=test_hartree_product<double,4>(world,k,thresh);
    error+=test_convolution(world,k,thresh);
    error+=test_multiply<double,4>(world,k,thresh);
    error+=test_add(world,k,thresh);
    error+=test_exchange(world,k,thresh);
    error+=test_inner(world,k,thresh);
    error+=test_replicate(world,k,thresh);

    print(ok(error==0),error,"finished test suite\n");
    world.gop.fence();
    finalize();

    return error;
}


