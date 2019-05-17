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
	for (int itight=3; itight<10; ++itight) {
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


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  test masks_and_boxes \n");
    	printf("starting at time %.1f\n", wall_time());

    }
    startup(world,argc,argv);
    std::cout.precision(6);

    try {

    	test_smooth_maxout(world);


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
