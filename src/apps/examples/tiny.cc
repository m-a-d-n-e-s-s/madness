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
  \file helium_mp2.cc
  \brief Solves the Hartree-Fock and MP2 equations for the helium atom
  \defgroup examplehehf Hartree-Fock and MP2 for the helium atom
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/helium_mp2.cc>here</a>.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/operator.h>
#include <mra/funcplot.h>
#include <mra/lbdeux.h>


#include <iostream>

using namespace madness;

template<size_t NDIM>
void load_function(World& world, Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0)  print("loading function ", name);

    archive::ParallelInputArchive ar(world, name.c_str());
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

    archive::ParallelOutputArchive ar(world, name.c_str());
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



template<size_t NDIM>
void draw_plane(World& world, Function<double,NDIM>& function, const std::string restart_name) {

    // determine the plot plane
    std::string c1, c2;	// the coordinates for the two electrons in human form ("x1" or "z2" or so)

    // the coordinates to be plotted
    Vector<double,NDIM> coord(0.0);

    std::ifstream f("input");
    position_stream(f, "plot");
    std::string s;
    while (f >> s) {
    	if (s == "end") {
    		break;
    	} else if (s == "plane") {
    		f >> c1 >> c2;
    	} else if (s == "electron1") {
    		f >> coord[0] >> coord[1] >> coord[2];
    	} else if (s == "electron2") {
    		f >> coord[3] >> coord[4] >> coord[5];
    	}
    }
    // convert human to mad form
    int cc1, cc2;
    if (c1=="x1") cc1=0;
    if (c1=="y1") cc1=1;
    if (c1=="z1") cc1=2;
    if (c1=="x2") cc1=3;
    if (c1=="y2") cc1=4;
    if (c1=="z2") cc1=5;
    if (c2=="x1") cc2=0;
    if (c2=="y1") cc2=1;
    if (c2=="z1") cc2=2;
    if (c2=="x2") cc2=3;
    if (c2=="y2") cc2=4;
    if (c2=="z2") cc2=5;

    std::string filename="plane_"+c1+c2+"_"+restart_name;
    const double scale=0.25;
    // assume a cubic cell
    double lo=-FunctionDefaults<6>::get_cell_width()[0]*0.5;
    lo=lo*scale;
//    const double hi=FunctionDefaults<6>::get_cell_width()[0]*0.5;

    const long nstep=100;
    const double stepsize=FunctionDefaults<6>::get_cell_width()[0]*scale/nstep;

    if(world.rank() == 0) {

    	// plot 3d plot
    	FILE *f =  0;
    	f=fopen(filename.c_str(), "w");
    	if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);

    	for (int i0=0; i0<nstep; i0++) {
    		for (int i1=0; i1<nstep; i1++) {
    			// plot plane
    			coord[cc1]=lo+i0*stepsize;
    			coord[cc2]=lo+i1*stepsize;

    			// other electron
    			fprintf(f,"%12.6f %12.6f %12.6f\n",coord[cc1],coord[cc2],function(coord));

    		}
    		// uncomment for gnuplot-style; leave commented out for mathematica
    		//        fprintf(f,"\n");
    	}
    	fclose(f);

    }

    // plot mra structure
	FILE *file =  0;
	filename="mra_structure_"+c1+c2+"_"+restart_name;
//	file=fopen(filename.c_str(), "w");
	if(!f) MADNESS_EXCEPTION("plot_along: failed to open the plot file", 0);
	function.get_impl()->print_plane(filename.c_str(),cc1,cc2,coord);
//	fclose(file);

}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    static const size_t NDIM=6;

    // determine the box size L
    double L=-1.0;
    std::ifstream f("input");
    position_stream(f, "dft");
    std::string s;
    while (f >> s) {
    	if (s == "end") {
    		break;
    	} else if (s == "L") {
    		f >> L;
    	}
    }
    if (L<0.0) MADNESS_EXCEPTION("box size indetermined",1);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<6>::set_cubic_cell(-L,L);
//    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<6>::set_tensor_type(TT_FULL);


    if (world.rank()==0) {
     	    print("cell size:         ", FunctionDefaults<6>::get_cell_width()[0]);
    }

    // load the function of interest
    std::string restart_name;
    bool restart=false;

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="restart") {                               // usage: restart=path/to/mo_file
            restart_name=stringify(val);
            restart=true;
        }
    }

    Function<double,NDIM> pair;
    if (restart) {
    	Function<double,NDIM> r12phi;
    	load_function(world,pair,restart_name);
//    	load_function(world,r12phi,"r12phi");
//    	pair=pair+r12phi;
    }

    // make sure we're doing what we want to do
    if (world.rank()==0) {
        print("polynomial order:  ", FunctionDefaults<6>::get_k());
        print("threshold:         ", FunctionDefaults<6>::get_thresh());
        print("cell size:         ", FunctionDefaults<6>::get_cell()(0,1) - FunctionDefaults<6>::get_cell()(0,0));
        print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
        print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());
        print("");
        print("orthogonalization  ", OrthoMethod());
        print("facReduce          ", GenTensor<double>::fac_reduce());
        print("max displacement   ", Displacements<6>::bmax_default());
        print("apply randomize    ", FunctionDefaults<6>::get_apply_randomize());
        print("world.size()       ", world.size());
        print("");
    }

//    if (restart) {
//    	draw_line(world,pair,restart_name);
//    	draw_circle(world,pair,restart_name);
//    }
    draw_plane(world,pair,restart_name);
    draw_line(world,pair,restart_name);

    world.gop.fence();
    print("exiting tiny");

    return 0;
}

