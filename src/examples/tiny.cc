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
/*!
  \file tiny.cc
  \brief Solves the Hartree-Fock and MP2 equations for the helium atom
  \defgroup examplehehf Hartree-Fock and MP2 for the helium atom
  \ingroup examples

  The source is
  <a href=https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/examples/helium_mp2.cc>here</a>.


*/


#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/lbdeux.h>

#include <iostream>

using namespace madness;

namespace madness{
extern std::vector<std::string> cubefile_header(std::string filename="input", const bool& no_orient=false);
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
    if (world.rank()==0)  print("loading function ", name);

    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, name.c_str());
    ar & pair;

    std::string line="saved function "+name;
    pair.print_size(line);

}


template<size_t NDIM>
void draw_line(World& world, Function<double,NDIM>& pair, const std::string restart_name) {

    Vector<double,NDIM> lo(0.0), hi(0.0);
    lo[2]=-8.0;
    hi[2]=8.0;

    {
        std::string filename="line_"+restart_name;
        trajectory<NDIM> line=trajectory<NDIM>::line2(lo,hi,601);
        plot_along<NDIM>(world,line,pair,filename);
    }

}

template<size_t NDIM>
void draw_circle(World& world, Function<double,NDIM>& pair, const std::string restart_name) {

	std::string filename="circle_"+restart_name;
	coord_3d el2(0.0);
	el2[1]=0.5;
	trajectory<NDIM> circ(0.5,el2,601);
	plot_along<NDIM>(world,circ,pair,filename);

}


void dostuff(World& world) {
    real_function_3d rho=real_factory_3d(world),rhonemo=real_factory_3d(world),
            vsigaa,vsigaanemo;
    real_function_3d sigmanemo=real_factory_3d(world);
    real_function_3d dens_pt0=real_factory_3d(world);

    load(rho,"rho");
    load(rhonemo,"rhonemo");
    load(sigmanemo,"sigmanemo");
    load(dens_pt0,"dens_pt0");

    FunctionDefaults<3>::set_k(rho.k());
    FunctionDefaults<3>::set_thresh(rho.thresh());

    rho.reconstruct();
    rhonemo.reconstruct();
    std::vector< std::shared_ptr<Derivative<double,3> > > gradop =
             gradient_operator<double,3>(world);

    // vsig with nemos
    {

    }

    // vsig with mos
    {
        real_function_3d drhoa_x=(*gradop[0])(rho, true).refine();
        real_function_3d drhoa_y=(*gradop[1])(rho, true).refine();
        real_function_3d drhoa_z=(*gradop[2])(rho, true).refine();
        world.gop.fence();

        // assign the reduced densities sigma
        vsigaa=(drhoa_x * drhoa_x + drhoa_y * drhoa_y + drhoa_z * drhoa_z);
    }

    double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
    coord_3d start(0.0); start[0]=-width;
    coord_3d end(0.0); end[0]=width;


    plot_line("line_rho",10000,start,end,rho);
    plot_line("line_vsigaa",10000,start,end,vsigaa);
    plot_line("line_dens_pt0",10000,start,end,dens_pt0);

    SeparatedConvolution<double,3> smooth=SmoothingOperator3D(world,0.005);
    real_function_3d sigmanemo_smooth=smooth(sigmanemo);
    plot_line("line_sigaanemo_smooth",10000,start,end,sigmanemo_smooth);
    real_function_3d dens_pt0_smooth=smooth(dens_pt0);
    plot_line("line_dens_pt0_smooth",10000,start,end,dens_pt0_smooth);

}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);


    // determine the box size L
    double L=-1.0;
    bool no_orient=false;
    std::ifstream f("input");
    position_stream(f, "dft");
    std::string s;
    while (f >> s) {
    	if (s == "end") {
    		break;
    	} else if (s == "L") {
    		f >> L;
        } else if (s == "no_orient") {
                no_orient=true;
        }
    }

    std::string c;
    position_stream(f, "plot");
    while (f >> s) {
    	if (s == "end") {
    		break;
    	} else if (s == "line") {
    		f >> c;
        }
    }

    if (L<0.0) MADNESS_EXCEPTION("box size indetermined",1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<6>::set_cubic_cell(-L,L);
//    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<6>::set_tensor_type(TT_FULL);

    // convert human to mad form (analogous to plot_plane in funcplot.h)
    unsigned int cc;
    if (c == "x1") cc = 0;
    else if (c == "x2") cc = 1;
    else if (c == "x3") cc = 2;
    else {
    	print("plot line axis not defined, plotting x1 axis");
    	c = "x1";
    	cc = 0;
    }

    if (world.rank()==0) {
     	    print("cell size:         ", FunctionDefaults<6>::get_cell_width()[0]);
    }

    // load the function of interest
    std::vector<std::string> filenames;

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="file") {                               // usage: restart=path/to/mo_file
            filenames.push_back(stringify(val));
        }
    }

	// make sure we're doing what we want to do
	if (world.rank()==0) {
		print("polynomial order:  ", FunctionDefaults<6>::get_k());
		print("threshold:         ", FunctionDefaults<6>::get_thresh());
		print("cell size:         ", FunctionDefaults<6>::get_cell()(0,1) - FunctionDefaults<6>::get_cell()(0,0));
		print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
		print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());
		print("");
		print("facReduce          ", GenTensor<double>::fac_reduce());
		print("max displacement   ", Displacements<6>::bmax_default());
		print("apply randomize    ", FunctionDefaults<6>::get_apply_randomize());
		print("world.size()       ", world.size());
		print("no_orient          ", no_orient);
		print("");
	}

//	dostuff(world);

    try {
        static const size_t NDIM=3;
        std::vector<Function<double,NDIM> > vf;
        for (size_t i=0; i<filenames.size(); ++i) {
            real_function_3d tmp;
            try { // load a single function
                load_function(world,tmp,filenames[i]);
                vf.push_back(tmp);
            } catch (...) { // load a vector of functions
                std::vector<Function<double,NDIM> > tmp2;
                load_function(world,tmp2,filenames[i]);
                for (auto& t : tmp2) vf.push_back(t);
            }
        }
		plot_plane(world,vf,filenames[0]);

		double width = FunctionDefaults<3>::get_cell_min_width()/2.0 - 1.e-3;
		coord_3d start(0.0); start[cc]=-width;
		coord_3d end(0.0); end[cc]=width;
		plot_line(("line_"+c+"_"+filenames[0]).c_str(),10000,start,end,vf[0]);

		// plot the Gaussian cube file
		std::vector<std::string> molecular_info=cubefile_header("input",no_orient);
		std::string filename=filenames[0]+".cube";
		plot_cubefile<3>(world,vf[0],filename,molecular_info);

    } catch (...) {
        try {
            static const size_t NDIM=6;
            std::vector<Function<double,NDIM> > vf(filenames.size());
            for (size_t i=0; i<filenames.size(); ++i) load_function(world,vf[i],filenames[i]);
            plot_plane(world,vf,filenames[0]);
        } catch (...) {

        }
    }



    world.gop.fence();
    print("exiting tiny");

    return 0;
}

