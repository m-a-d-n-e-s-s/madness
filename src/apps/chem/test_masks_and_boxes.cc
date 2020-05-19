/*
 * test_masks_and_boxes.cc
 *
 *  Created on: 16 May 2019
 *      Author: fbischoff
 */


#include <chem/masks_and_boxes.h>
#include <madness/mra/mra.h>

using namespace madness;

bool test_smooth_maxout(World& world) {

	bool success=true;
	for (int itight=3; itight<4; ++itight) { // was 10 ... reduced to speed up the test
		double radius=10.0;
		double deviation=std::pow(10,-double(itight));
		double tightness=max_of_x_1_smooth::compute_tightness(deviation,radius);
		// plot
		std::string filename="line_flatten_out_tightness_"+std::to_string(deviation);
		std::function<double(const coord_1d&)> op = [&tightness, &radius] (const coord_1d& r)
						{return max_of_x_1_smooth::compute_factor(r[0],tightness,radius);};
		plot_line<std::function<double(const coord_1d&)>,1 >(world, filename.c_str(), 100, {0.0}, {2.0*radius}, op);

		// compare values
		double r=0.8*radius;
		double actual_deviation=fabs(r-r*op({r}));
		print("deviation, tightness, r, maxed(r)",deviation,tightness,r,r*op({r}), actual_deviation);
		if (actual_deviation>deviation) {
			print("test failed for demanded deviation ",deviation);
			success=false;
		}
	}
	return success;
}


bool test_spherical_box(World& world) {

	print("entering test_spherical_box");
	bool success=true;
	for (int ideviation=2; ideviation<3; ++ideviation) { // was 4 ... reduced to speed up test
		double radius=10.0;
		double deviation=std::pow(10,-double(ideviation));

		Vector<double,3> offset={0.0,0.0,0.0};
		Vector<double,3> direction={0.0,0.0,1.0};

		spherical_box<3> sbox1(radius,deviation,offset,direction);

		// project
		real_function_3d sbox=real_factory_3d(world).functor(sbox1);
		sbox.print_size("sbox with deviation 1e.-"+stringify(ideviation));

		// test accuracy
		double actual=sbox(0.0,0.8*radius,0.0);
		double expect=1.0-std::pow(10.0,-double(ideviation));
		double error=std::abs(actual-expect);

		// compare values
		print("deviation, actual, expect,error",deviation,actual,expect,error);
		if (error>deviation*1.2) {
			print("test failed for demanded deviation ",deviation);
			success=false;
		}
	}
	return success;
}



int main(int argc, char** argv) {

	World& world=initialize(argc, argv);
	if (world.rank() == 0) {
		print("\n  test masks_and_boxes \n");
		printf("starting at time %.1f\n", wall_time());

	}
	startup(world,argc,argv);
	std::cout.precision(6);


	FunctionDefaults<3>::set_k(8); // was 8 but lower order can be faster for deeper refinement
	FunctionDefaults<3>::set_thresh(1.e-4); // was 1e-5 but increaesd to speed up test
	FunctionDefaults<3>::set_refine(true);
	FunctionDefaults<3>::set_initial_level(5);
	FunctionDefaults<3>::set_truncate_mode(1);
	FunctionDefaults<3>::set_cubic_cell(-20, 20);


	try {

		test_smooth_maxout(world);
		test_spherical_box(world);


	} catch (const SafeMPI::Exception& e) {
		print(e);
		error("caught an MPI exception");
	} catch (const madness::MadnessException& e) {
		print(e);
		error("caught a MADNESS exception");
	} catch (const madness::TensorException& e) {
		print(e);
		error("caught a Tensor exception");
	} catch (const char* s) {
		print(s);
		error("caught a string exception");
	} catch (const std::string& s) {
		print(s);
		error("caught a string (class) exception");
	} catch (const std::exception& e) {
		print(e.what());
		error("caught an STL exception");
	} catch (...) {
		error("caught unhandled exception");
	}

	finalize();
	return 0;
}
