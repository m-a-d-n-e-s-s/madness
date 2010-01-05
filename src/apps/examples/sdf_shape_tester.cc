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
  \file examples/sdf_shape_tester.cc
  \brief Demonstrates/tests use of 3D shape functions
  \defgroup shape_tester Demonstrates/tests use of 3D shape functions
  \ingroup examples 

  \par Points of interest
  - Use of shape functions to define interior surfaces
  - Plotting for subsequent visualization with OpenDX
  
  \par Backgroud

  The classes in mra/sdf_shape_3D.h illustrate how to define 
  surfaces as MADNESS functions for subsequent use in calculations
  using them as sources or for enforcing boundary conditions on 
  interior surfaces.  This example instantiates each shape and
  plots it for subsequent visual verification (with OpenDX).

 */


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>

using namespace madness;

int main(int argc, char **argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    static const double L = 2.0;
    
    // Function defaults
    FunctionDefaults<3>::set_k(5);
    FunctionDefaults<3>::set_thresh(1e-3);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    
    Tensor<int> bc(3,2);
    bc(_,0) = 0;          // Dirichlet in all directions
    bc(_,1) = 0;
    FunctionDefaults<3>::set_bc(bc);
    
    // create the shape mask
    coord_3d pt, vec, sides;
    double c;
    pt[0] = 0.0;
    pt[1] = 0.5;
    pt[2] = 0.0;
    vec[0] = 0.0;
    vec[1] = 0.0;
    vec[2] = 1.0;
    sides[0]=0.3; 
    sides[1]=0.6; 
    sides[2]=1.0;
    c = 0.5;

    // for plotting the shapes the old code used VTK, but I dunno how to use that
    // so switched to opendx

    // the following line permutes the "inside" and "outside"
    //mask.unaryop(&mask_complement<double, 3>);

    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Cylinder<double>(0.2, 0.75, 1.0, pt, vec)));
        plotdx(f, "cylinder.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Cube<double>(0.2, sqrt(2.0), pt)));
        plotdx(f, "cube.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Cone<double>(0.2, c, pt, vec)));
        plotdx(f, "cone.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Paraboloid<double>(0.2, c, pt, vec)));
        plotdx(f, "paraboloid.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Plane<double>(0.2, vec, pt)));
        plotdx(f, "plane.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Sphere<double>(0.2, c, pt)));
        plotdx(f, "sphere.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Ellipsoid<double>(0.2, sides, pt)));
        plotdx(f, "ellipsoid.dx");
    }
    {
        real_function_3d f = real_factory_3d(world).functor(real_functor_3d(new SDF_Box<double>(0.2, sides, pt)));
        plotdx(f, "box.dx");
    }

    /*
    char filename[100];
    sprintf(filename, "shape.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -Lplot;
        plothi[i] = Lplot;
        npts[i] = 51;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(mask, "mask", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);
    */
    
    MPI::Finalize();
    
    return 0;
}
