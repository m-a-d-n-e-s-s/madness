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

/** \file tipefield.cc
    \brief 

*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include "tipsurface.h"

using namespace madness;

int main(int argc, char **argv) {
    double eps;
    int k;
    double thresh, phi, d, penalty;
    unsigned int i;
    char funcname[80];

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    
    if (world.rank() == 0) {
        if(argc < 6) {
            print("usage: ./nanophoto k thresh epsilon potential-difference " \
                "tip-surface\n");
            print("potential difference is in mV, and tip-surface distance " \
                "in nm\n");
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        eps = atof(argv[3]) / 0.052918; // convert to a.u.
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        phi = atof(argv[4]) * 3.6749324e-5; // convert to a.u.

        d = atof(argv[5]) / 0.052918; // convert to a.u.
    }
    world.gop.broadcast(phi);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(d);

    // box size
    Tensor<double> cell(3, 2);
    cell(0,0) = cell(1,0) = -300.0 / 0.052918;
    cell(0,1) = cell(1,1) = 300.0 / 0.052918;
    cell(2,0) = 10.0 * (constants::pi - 8.0) / 0.052918;
    cell(2,1) = (10.0 * (constants::pi - 8.0) + 600.0) / 0.052918;
    FunctionDefaults<3>::set_cell(cell);

    penalty = 1.75 / eps;

    // the key data structure: sets up the problem details and projects
    // the initial functions
    std::shared_ptr<TipSurface> tps(new TipSurface(eps, penalty, phi, d));

    if(world.rank() == 0) {
        // print out the arguments
        printf("Tip-Surface Distance: %.6e nm\nPotential Difference: %.6e " \
            "mV\nk: %d\nthresh: %.6e\neps: %.6e nm\n", d * 0.052918,
            phi * 27211.385, k, thresh, eps * 0.052918);

        // project the surface function
        printf("Load Balancing\n   --- Projecting domain mask to low order.\n");
        fflush(stdout);
    }

    real_function_3d dmask;

    // low order defaults
    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_thresh(1.0e-4);

    // domain mask
    tps->fop = DOMAIN_MASK;
    dmask = real_factory_3d(world).functor(tps);

    // do the load balancing
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(dmask, DirichletLBCost<3>(1.0, 1.0));
    // set this map as default and redistribute the existing functions
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0, false));
    dmask.clear();

    // set the defaults to the real deal
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);

#if 0
    // get the domain mask -- this segment can be commented out
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    sprintf(funcname, "domainmask");
    tps->fop = DOMAIN_MASK;
    dmask = real_factory_3d(world).functor(tps);
    vtk_output(world, funcname, dmask);
    dmask.clear();
#endif

    real_function_3d usol;

    // project the surface and rhs functions
    if(world.rank() == 0) {
        printf("Projecting the surface function\n");
        fflush(stdout);
    }
    tps->fop = SURFACE;
    real_function_3d surf = real_factory_3d(world).functor(tps);
    sprintf(funcname, "surface");

    if(world.rank() == 0) {
        printf("Projecting the rhs function\n");
        fflush(stdout);
    }
    tps->fop = DIRICHLET_RHS;
    real_function_3d rhs = real_factory_3d(world).functor(tps);
    sprintf(funcname, "rhs");

    // green's function
    // note that this is really -G...
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // project the r.h.s. function
    // and then convolute with G
    real_function_3d grhs;
    grhs = G(rhs);
    grhs.truncate();
    rhs.clear();

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to rhs in magnitude
    //     for starting in GMRES
    usol = copy(grhs);
    usol.scale(penalty);
    usol.compress();
    DirichletCondIntOp<3> dcio(G, surf);

    // make the operators and prepare GMRES
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-3;
    double update_thresh = 1.0e-3;
    int maxiter = 60;
    GMRES(space, dcio, grhs, usol, maxiter, resid_thresh, update_thresh,
        true);

    // scale back to mV from a.u.
    //usol.scale(27211.385);

    // get the x, y, and z components of the electric field
    real_derivative_3d Dx = free_space_derivative<double, 3>(world, 0);
    real_derivative_3d Dy = free_space_derivative<double, 3>(world, 1);
    real_derivative_3d Dz = free_space_derivative<double, 3>(world, 2);

    real_function_3d ex = Dx(usol);
    real_function_3d ey = Dy(usol);
    real_function_3d ez = Dz(usol);

    // print out the solution
    sprintf(funcname, "solution-d%snm.dat", argv[5]);
    Tensor<double> plotcell(3, 2);
    std::vector<long> npts(3);
    std::vector<double> dpt(3);
    FILE *f;

    plotcell(0, 0) = plotcell(1, 0) = -100.0 / 0.052918;
    plotcell(0, 1) = plotcell(1, 1) = 100.0 / 0.052918;
    plotcell(2, 0) = 0.0 / 0.052918;
    plotcell(2, 1) = 5.0 / 0.052918;
    npts[0] = npts[1] = 241;
    npts[2] = 11;
    for(i = 0; i < 3; ++i)
        dpt[i] = (plotcell(i, 1) - plotcell(i, 0)) / (npts[i] - 1);

    world.gop.barrier();
    usol.verify();
    ex.verify();
    ey.verify();
    ez.verify();

    world.gop.fence();
    Tensor<double> phigrid = usol.eval_cube(plotcell, npts, false);
    Tensor<double> exgrid = ex.eval_cube(plotcell, npts, false);
    Tensor<double> eygrid = ey.eval_cube(plotcell, npts, false);
    Tensor<double> ezgrid = ez.eval_cube(plotcell, npts, false);

    Vector<double, 3> pt;
    pt[0] = pt[1] = -300.0 / 0.052918;
    pt[2] = 3.0 / 0.052918;
    double farfield = -ez(pt);
    pt[0] = pt[1] = 0.0;
    double nearfield = -ez(pt);
    world.gop.fence();

    if(world.rank() == 0) {
        printf("Far from tip, Ez = %.6e a.u.\nNear tip, Ez = %.6e a.u.\n",
            farfield, nearfield);
        f = fopen(funcname, "w");
        if(!f)
            MADNESS_EXCEPTION("failed to open the output file", 0);

        // go through the grid
        for(LowDimIndexIterator it(npts); it; ++it) {
            for(i = 0; i < 3; ++i)
                fprintf(f, "%.6e ", (plotcell(i, 0) + it[i]*dpt[i]));
            fprintf(f, "%.6e ", phigrid(*it));
            fprintf(f, "%.6e ", -exgrid(*it));
            fprintf(f, "%.6e ", -eygrid(*it));
            fprintf(f, "%.6e\n", -ezgrid(*it));
        }

        fclose(f);
    }

    world.gop.fence();

    finalize();
    
    return 0;
}
