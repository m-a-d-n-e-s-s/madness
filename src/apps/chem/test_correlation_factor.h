/*
 * test_correlation_factor.h
 *
 *  Created on: Aug 28, 2015
 *      Author: kottmanj
 */
#ifndef TEST_CORRELATION_FACTOR
#define TEST_CORRELATION_FACTOR

#include<chem/electronic_correlation_factor.h>


struct test_correlation_factor{
	test_correlation_factor(World &world, const real_function_3d &f, const double eps) : world(world), corrfac(world), f(f), eps(eps){
		corrfac = CorrelationFactor(world, 1.0, 1.e-7, 1.e-6);
	}
	test_correlation_factor(World &world,const CorrelationFactor & corrfactor, const real_function_3d &f, const double eps) : world(world), corrfac(corrfactor), f(f), eps(eps){

	}

	World &world;
	CorrelationFactor corrfac;
	real_function_3d f;
	const double eps;


	void make_test(){
		output("Testing Correlation Factor");
		real_function_6d ftest = corrfac.apply_U(f,f,eps);
		output("Ended Testing for the given function f\n the norm is: " + stringify(ftest.norm2()));
		plot_plane(world,ftest,"Uef");
	}

	void output(const std::string &msg)const{
		if(world.rank()==0){
			std::cout << msg << std::endl;
		}
	}

};

#endif
