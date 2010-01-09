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
#include <constants.h>

using namespace madness;

int main(int argc, char **argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    static const double L = 2.0;

    double epsilon = 0.2;
    
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
    double c, rad, vol, area;
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
    rad = 0.75;

    // for plotting the shapes the old code used VTK, but I dunno how to use that
    // so switched to opendx

    real_function_3d f;

    // the following line permutes the "inside" and "outside"
    //mask.unaryop(&mask_complement<double, 3>);
    
    f = real_factory_3d(world).functor(shape_mask(epsilon, new SDFCylinder(rad, 1.0, pt, vec)));
    vol = f.trace();
    plotdx(f, "cylinder_mask.dx");
    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFCylinder(rad, 1.0, pt, vec)));
    area = f.trace();
    plotdx(f, "cylinder_surface.dx");
    if (world.rank() == 0)
        print("cylinder: err in volume", vol - constants::pi*rad*rad, "err in area", area-2.0*3.141593*(rad + rad*rad));

    f = real_factory_3d(world).functor(shape_mask(epsilon, new SDFCube(sqrt(2.0), pt)));
    vol = f.trace();
    plotdx(f, "cube_mask.dx");
    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFCube(sqrt(2.0), pt)));
    area = f.trace();
    plotdx(f, "cube_surface.dx");
    if (world.rank() == 0)
        print("cube: err in volume", vol-sqrt(2.0)*sqrt(2.0)*sqrt(2.0), "err in area", area-12.0);

    f = real_factory_3d(world).functor(shape_mask(epsilon, new SDFSphere(c, pt)));
    vol = f.trace();
    plotdx(f, "sphere_mask.dx");
    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFSphere(c, pt)));
    area = f.trace();
    plotdx(f, "sphere_surface.dx");
    if (world.rank() == 0)
        print("sphere: err in volume", vol-4.0*constants::pi/3.0*c*c*c, "err in area", area-4.0*constants::pi*c*c);

    f = real_factory_3d(world).functor(shape_mask(epsilon, new SDFCone(c, pt, vec)));
    plotdx(f, "cone_mask.dx");

    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFCone(c, pt, vec)));
    plotdx(f, "cone_surface.dx");

    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFParaboloid(c, pt, vec)));
    plotdx(f, "paraboloid_surface.dx");

    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFPlane(vec, pt)));
    plotdx(f, "plane_surface.dx");


    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFEllipsoid(sides, pt)));
    plotdx(f, "ellipsoid_surface.dx");

    f = real_factory_3d(world).functor(shape_surface(epsilon, new SDFBox(sides, pt)));
    plotdx(f, "box_surface.dx");
    
    MPI::Finalize();
    
    return 0;
}
