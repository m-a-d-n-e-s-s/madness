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
  \file mp2.cc
  \brief Solves molecular MP2 equations
  \defgroup Solves molecular MP2 equations
  \ingroup examples

  The source is
  <a href=http://code.google.com/p/m-a-d-n-e-s-s/source/browse/local/trunk/src/apps/examples/mp2.cc>here</a>.


*/


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <examples/mp2.h>


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

#ifdef MADNESS_HAS_GOOGLE_PERF_MINIMAL
    if (world.rank()==0) print("using gperftools, clearing memory at each fence()");
#endif

    CorrelationFactor f12(world,1.0);

    MP2 mp2(world,f12,"input");

    TensorType tt=TT_2D;
    FunctionDefaults<6>::set_tensor_type(tt);

	// find the pair to compute
    std::ifstream f("input");
    position_stream(f, "mp2");
    std::string s;

    int iteration;

    int i=0,j=0;
    while (f >> s) {
        if (s == "end") break;
        else if (s == "pair") f >> i >> j;
        else if (s == "iteration") f >> iteration;
        else continue;
    }

    if(world.rank() == 0) printf("\ncomputing pair (%d,%d)\n\n", i,j);


    // get command line parameters (overrides input file)
    bool do_test=false;
    for(int ii = 1; ii < argc; ii++) {
        const std::string arg=argv[ii];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="test") do_test=true;
    }


    FunctionDefaults<6>::set_apply_randomize(true);
    if (world.rank()==0) {
        print("polynomial order:  ", FunctionDefaults<6>::get_k());
        print("threshold 3D:      ", FunctionDefaults<3>::get_thresh());
        print("threshold 6D:      ", FunctionDefaults<6>::get_thresh());
        print("cell size:         ", FunctionDefaults<6>::get_cell_width()[0]);
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

    if(world.rank() == 0) printf("\nstarting at time %.1fs\n", wall_time());

    if (do_test) mp2.test(i,j,iteration);
    else mp2.solve_residual_equations(i,j);

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;
}
