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

/** \file nanophoto.cc
    \brief 

*/

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <linalg/gmres.h>
#include "basisfunction.h"
#include "atom.h"
#include "density.h"

int mol_geom(std::vector<Atom*> &atoms);
void read_states(int nstate, int nbasis, Tensor<double> &coeffs);

using namespace madness;

int main(int argc, char **argv) {
    double eps, penalty_prefact;
    int j, k, n, prob, dim;
    double thresh, radius;
    unsigned int i;

    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    FunctionDefaults<3>::set_k(6);
    FunctionDefaults<3>::set_cubic_cell(-15.0, 31.0);
    FunctionDefaults<3>::set_thresh(1.0e-4);
    FunctionDefaults<3>::set_max_refine_level(6);

    std::vector<Atom*> atoms(0);
    std::vector<BasisFunc> basis(0);
    int nstate = mol_geom(atoms);
    int nbasis = 0;

    // make the set of basis functions
    for(i = 0; i < atoms.size(); ++i) {
        j = atoms[i]->dimBasis();
        nbasis += j;

        for(n = 0; n < j; ++n)
            basis.push_back(atoms[i]->getBasisFunc(n));
    }

    // read in the coefficients of the occupied states
    Tensor<double> coeffs(nstate, nbasis);
    if(world.rank() == 0)
        read_states(nstate, nbasis, coeffs);
    world.gop.broadcast(coeffs);

    // make the overall problem-control functor
    SharedPtr<TipMolecule> tpm(new TipMolecule(1.0, 1.0, coeffs, basis));

Vector<double, 3> pt;
print("projecting the density");
    tpm->fop = DENSITY;
pt[0] = pt[1] = pt[2] = 0.0;
print(tpm->operator()(pt));
    real_function_3d density = real_factory_3d(world).functor(tpm);

double norm = density.norm2();
double trace = density.trace();
print(norm, trace);
    char filename[100];
    sprintf(filename, "density.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -10.0;
        plothi[i] = 26.0;
        npts[i] = 101;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(density, "density", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);

#if 0
    // the structures for the problem, unfortunately different DIMs need
    // different structures...
    SharedPtr<EmbeddedDirichlet<2> > functor2;
    SharedPtr<EmbeddedDirichlet<3> > functor3;
    
    if (world.rank() == 0) {
        if(argc < 6) {
            error("bad number of arguments");
        }

        // read in and validate the command-line arguments
        k = atoi(argv[1]);
        if(k < 4) error("cheapskate");

        thresh = atof(argv[2]);
        if(thresh > 1.0e-4) error("use some real thresholds...");

        prob = atoi(argv[3]);
        if(prob < 1 || prob > 6) error("bad problem number");

        eps = atof(argv[4]);
        if(eps <= 0.0) error("eps must be positive, and hopefully small");

        mu::Parser parser;
        try {
            parser.DefineVar("eps", &eps);
            parser.SetExpr(std::string(argv[5]));
            penalty_prefact = parser.Eval();
        }
        catch(mu::Parser::exception_type &e) {
            error(e.GetMsg().c_str());
        }
        if(penalty_prefact <= 0.0) error("penalty_prefact must be positive");

        switch(atoi(argv[6])) {
        case 1:
            mask = LLRV;
            break;
        case 2:
            mask = Gaussian;
            break;
        default:
            error("unknown domain mask type, should be 1 or 2");
            break;
        }

        if(prob == 1 || prob == 2 || prob == 5 || prob == 6) {
            if(argc > 7) {
                radius = atof(argv[7]);
                if(radius <= 0.0) error("radius must be positive");
            }
            else
                radius = 1.0;
        }
    }
    world.gop.broadcast(prob);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(mask);
    world.gop.broadcast(penalty_prefact);
    world.gop.broadcast(radius);

    // do the final problem setup
    switch(prob) {
    case 1:
        functor3 = SharedPtr<EmbeddedDirichlet<3> >(new ConstantSphere(k,
                       thresh, eps, std::string(argv[5]), penalty_prefact,
                       radius, mask));
        dim = 3;
        break;
    case 2:
        functor3 = SharedPtr<EmbeddedDirichlet<3> >(new CosineSphere(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, radius,
                       mask));
        dim = 3;
        break;
    case 3:
        functor3 = SharedPtr<EmbeddedDirichlet<3> >(new Ellipsoid(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, mask));
        dim = 3;
        break;
    case 4:
        functor2 = SharedPtr<EmbeddedDirichlet<2> >(new LLRVCircle(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, mask));
        dim = 2;
        break;
    case 5:
        functor3 = SharedPtr<EmbeddedDirichlet<3> >(new Y20Sphere(k, thresh,
                       eps, std::string(argv[5]), penalty_prefact, radius,
                       mask));
        dim = 3;
        break;
    case 6:
        functor3 = SharedPtr<EmbeddedDirichlet<3> >(new InhomoConstantSphere(k,
                       thresh, eps, std::string(argv[5]), penalty_prefact,
                       radius, mask));
        dim = 3;
        break;
    default:
        dim = 0;
        error("shouldn't be here");
        break;
    }

    if(world.rank() == 0) {
        // print out the arguments
        switch(dim) {
        case 2:
            functor2->printout();
            break;
        case 3:
            functor3->printout();
            break;
        }
    }

    // project the surface function
    if(world.rank() == 0) {
        printf("Projecting the surface function (to low order)\n");
        fflush(stdout);
    }
    real_function_2d surf2;
    real_function_3d surf3;

    switch(dim) {
    case 2:
        functor2->fop = SURFACE;
        surf2 = real_factory_2d(world).k(6).thresh(1.0e-4).functor(functor2);
        break;
    case 3:
        functor3->fop = SURFACE;
        surf3 = real_factory_3d(world).k(6).thresh(1.0e-4).functor(functor3);
        break;
    }

    if(world.rank() == 0) {
        printf("Performing load balancing\n");
        fflush(stdout);
    }
    switch(dim) {
    case 2:
        break;
    case 3:
        functor3->load_balance(world, surf3);
        break;
    }

    // reproject the surface function to the requested threshold / k
    if(k > 6 || thresh < 1.0e-4) {
        if(world.rank() == 0) {
            printf("Reprojecting the surface function to requested order\n");
            fflush(stdout);
        }
        switch(dim) {
        case 2:
            surf2.clear();
            surf2 = real_factory_2d(world).functor(functor2);
            break;
        case 3:
            surf3.clear();
            surf3 = real_factory_3d(world).functor(functor3);
            break;
        }
    }

    // project the domain mask
    real_function_2d phi2;
    real_function_3d phi3;
    if(world.rank() == 0) {
        printf("Projecting the domain mask\n");
        fflush(stdout);
    }
    switch(dim) {
    case 2:
        functor2->fop = DOMAIN_MASK;
        phi2 = real_factory_2d(world).functor(functor2);
        break;
    case 3:
        functor3->fop = DOMAIN_MASK;
        phi3 = real_factory_3d(world).functor(functor3);
        break;
    }

    // print out the errors in perimeter (2-D) and surface area (3-D)
    // these are checks of the diffuse domain approximation
    double surf_integral, anals;
    double vol_integral, analv;
    surf_integral = 0.0;
    anals = 0.0;
    vol_integral = 0.0;
    analv = 0.0;

    switch(dim) {
    case 2:
        surf_integral = surf2.trace(); 
        anals = functor2->SurfaceIntegral();
        vol_integral = phi2.trace();
        analv = functor2->VolumeIntegral();
        break;
    case 3:
        surf_integral = surf3.trace();
        anals = functor3->SurfaceIntegral();
        vol_integral = phi3.trace();
        analv = functor3->VolumeIntegral();
        break;
    }
    if(world.rank() == 0) {
        printf("Error in Surface Integral: %.6e\n",
            fabs(surf_integral/penalty_prefact - anals));
        printf("Error in Volume Integral: %.6e\n",
            fabs(vol_integral - analv));
    }

    // green's function
    // note that this is really -G...
    real_convolution_2d G2 = BSHOperator<2>(world, 0.0, eps*0.1, thresh);
    real_convolution_3d G3 = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // project the r.h.s. function (phi*f - penalty*S*g)
    // and then convolute with G
    real_function_2d usol2, rhs2;
    real_function_3d usol3, rhs3;
    if(world.rank() == 0) {
        printf("Projecting the r.h.s. function\n");
        fflush(stdout);
    }

    switch(dim) {
    case 2:
        functor2->fop = DIRICHLET_RHS;
        usol2 = real_factory_2d(world).functor(functor2);
        rhs2 = G2(usol2);
        rhs2.truncate();
        usol2.clear();
        break;
    case 3:
        functor3->fop = DIRICHLET_RHS;
        usol3 = real_factory_3d(world).functor(functor3);
        rhs3 = G3(usol3);
        rhs3.truncate();
        usol3.clear();
        break;
    }

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to rhs in magnitude for
    //     starting in GMRES
    DirichletCondIntOp<2> *dcio2 = NULL;
    DirichletCondIntOp<3> *dcio3 = NULL;
    switch(dim) {
    case 2:
        usol2 = copy(rhs2);
        usol2.scale(1.0 / penalty_prefact);
        usol2.compress();
        dcio2 = new DirichletCondIntOp<2>(G2, surf2);
        break;
    case 3:
        usol3 = copy(rhs3);
        usol3.scale(1.0 / penalty_prefact);
        usol3.compress();
        dcio3 = new DirichletCondIntOp<3>(G3, surf3);
        break;
    }

    // make the operators and prepare GMRES
    FunctionSpace<double, 2> space2(world);
    FunctionSpace<double, 3> space3(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    switch(dim) {
    case 2:
        GMRES(space2, *dcio2, rhs2, usol2, maxiter, resid_thresh, update_thresh,
              true);
        break;
    case 3:
        GMRES(space3, *dcio3, rhs3, usol3, maxiter, resid_thresh, update_thresh,
              true);
        break;
    }

    // compare to the exact solution
    real_function_2d uexact2, uerror2;
    real_function_3d uexact3, uerror3;
    double error, ratio, exactval;
    error = 0.0;
    ratio = 0.0;
    exactval = 0.0;

    std::vector<Vector<double, 2> > check_pts2;
    std::vector<Vector<double, 3> > check_pts3;

    switch(dim) {
    case 2:
        functor2->fop = EXACT;
        uexact2 = real_factory_2d(world).functor(functor2);
        uerror2 = (usol2 - uexact2)*phi2; // only use interior solution
        error = uerror2.norm2();
        break;
    case 3:
        functor3->fop = EXACT;
        uexact3 = real_factory_3d(world).functor(functor3);
        uerror3 = (usol3 - uexact3)*phi3; // only use interior solution
        error = uerror3.norm2();
        break;
    }

    if(world.rank() == 0) {
        printf("\nu interior error: %.10e\n", error);
        fflush(stdout);
    }

    // check the points prescribed by the problem
    switch(dim) {
    case 2:
        check_pts2 = functor2->check_pts();
        for(std::vector<Vector<double, 2> >::iterator iter =
            check_pts2.begin(); iter != check_pts2.end(); ++iter) {

            ratio = usol2(*iter);
            exactval = functor2->ExactSol(*iter);

            if(world.rank() == 0) {
                printf("u/uexact ratio at (%.2f, %.2f): %.10e\n",
                    (*iter)[0], (*iter)[1], ratio/exactval);
            }
        }
        break;
    case 3:
        check_pts3 = functor3->check_pts();
        for(std::vector<Vector<double, 3> >::iterator iter =
            check_pts3.begin(); iter != check_pts3.end(); ++iter) {

            ratio = usol3(*iter);
            exactval = functor3->ExactSol(*iter);

            if(world.rank() == 0) {
                printf("u/uexact ratio at (%.2f, %.2f, %.2f): %.10e\n",
                    (*iter)[0], (*iter)[1], (*iter)[2], ratio/exactval);
            }
        }
        break;
    }

    // uncomment these lines for various plots
    if(dim == 2) {
        // make a line plot along the positive z axis
        if(world.rank() == 0)
            printf("\n\n");
        double xmin = 0.0;
        double xmax = 2.0;
        int nx = 201;
        coord_2d pt2;
        pt2[1] = 0.0;
        double dx = (xmax - xmin) / (nx - 1);
        for(int i = 0; i < nx; ++i) {
            pt2[0] = xmin + i * dx;
            double uval = usol2(pt2);

            if(world.rank() == 0) {
                printf("%.4e %.4e\n", pt2[0], uval);
            }
        }
    }

    // make a line plot along the positive z axis
    if(dim == 3) {
        if(world.rank() == 0)
            printf("\n\n");
        double zmin = 0.0;
        double zmax = 2.0;
        int nz = 201;
        coord_3d pt;
        pt[0] = pt[1] = 0.0;
        double dz = (zmax - zmin) / (nz - 1);
        for(int i = 0; i < nz; ++i) {
            pt[2] = zmin + i * dz;
            double uval = usol3(pt);

            if(world.rank() == 0) {
                printf("%.4e %.4e\n", pt[2], uval);
            }
        }
    }

    // print out the solution function
    /*char filename[100];
    sprintf(filename, "spheresol.vts");
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = -2.0;
        plothi[i] = 2.0;
        npts[i] = 101;
    }
    plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, "usol", world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);*/
#endif

    finalize();
    
    return 0;
}

/** \brief Make the molecular geometry.

    \return The number of (occupied) states to use in the density. */
int mol_geom(std::vector<Atom*> &atoms) {
    Vector<double, 3> center;

    center[0] = 27.9797743030;
    center[1] = 13.9898852618;
    center[2] = 24.2819925053;
    atoms.push_back(new Hydrogen(center));

    center[0] = 26.4224322192;
    center[1] = 12.9342276332;
    center[2] = 23.5703519339;
    atoms.push_back(new Carbon(center));

    center[0] = 24.6432570915;
    center[1] = 11.7396090042;
    center[2] = 22.7317293351;
    atoms.push_back(new Carbon(center));

    center[0] = 22.6131358013;
    center[1] = 10.2929746277;
    center[2] = 21.7587226521;
    atoms.push_back(new Carbon(center));

    center[0] = 22.5443289883;
    center[1] = 8.3037640154;
    center[2] = 22.2745082426;
    atoms.push_back(new Hydrogen(center));

    center[0] = 20.7886299265;
    center[1] = 11.2477851392;
    center[2] = 20.2100789771;
    atoms.push_back(new Carbon(center));

    center[0] = 20.8714887419;
    center[1] = 13.2423247788;
    center[2] = 19.7126180594;
    atoms.push_back(new Hydrogen(center));

    center[0] = 18.8039395388;
    center[1] = 9.8189122973;
    center[2] = 19.1604854792;
    atoms.push_back(new Carbon(center));

    center[0] = 17.0840677920;
    center[1] = 8.6326481468;
    center[2] = 18.1574737270;
    atoms.push_back(new Carbon(center));

    center[0] = 15.1780391382;
    center[1] = 7.2097203829;
    center[2] = 16.9718785396;
    atoms.push_back(new Carbon(center));

    center[0] = 15.0817292532;
    center[1] = 5.2029315394;
    center[2] = 17.4118672207;
    atoms.push_back(new Hydrogen(center));

    center[0] = 13.4903097407;
    center[1] = 8.1896566911;
    center[2] = 15.2743339053;
    atoms.push_back(new Carbon(center));

    center[0] = 13.5871752051;
    center[1] = 10.1998338132;
    center[2] = 14.8468759972;
    atoms.push_back(new Hydrogen(center));

    center[0] = 11.6540573293;
    center[1] = 6.7673619853;
    center[2] = 13.9877158563;
    atoms.push_back(new Carbon(center));

    center[0] = 10.0735547638;
    center[1] = 5.5816609732;
    center[2] = 12.7727618859;
    atoms.push_back(new Carbon(center));

    center[0] = 8.3174456315;
    center[1] = 4.1700035351;
    center[2] = 11.3647251628;
    atoms.push_back(new Carbon(center));

    center[0] = 8.2357754538;
    center[1] = 2.1445007351;
    center[2] = 11.7122306541;
    atoms.push_back(new Hydrogen(center));

    center[0] = 6.7367106298;
    center[1] = 5.1856670028;
    center[2] = 9.5885263535;
    atoms.push_back(new Carbon(center));

    center[0] = 6.8080553447;
    center[1] = 7.2156427841;
    center[2] = 9.2608535364;
    atoms.push_back(new Hydrogen(center));

    center[0] = 5.0096409167;
    center[1] = 3.7822771158;
    center[2] = 8.1324452368;
    atoms.push_back(new Carbon(center));

    center[0] = 3.4972591937;
    center[1] = 2.5826563821;
    center[2] = 6.8514774685;
    atoms.push_back(new Carbon(center));

    center[0] = 3.0930827097;
    center[1] = -2.2708591530;
    center[2] = 6.4104229818;
    atoms.push_back(new Hydrogen(center));

    center[0] = 1.6664964362;
    center[1] = 1.2864139586;
    center[2] = 5.3697319832;
    atoms.push_back(new Carbon(center));

    center[0] = 1.6410229299;
    center[1] = -1.2984156084;
    center[2] = 5.3279992744;
    atoms.push_back(new Carbon(center));

    center[0] = -0.1947041377;
    center[1] = 2.9238142838;
    center[2] = 3.8989147832;
    atoms.push_back(new Carbon(center));

    center[0] = -0.1803289921;
    center[1] = -2.9329076453;
    center[2] = 3.9001903483;
    atoms.push_back(new Carbon(center));

    center[0] = -2.1520993233;
    center[1] = 2.4815598212;
    center[2] = 4.4429102033;
    atoms.push_back(new Hydrogen(center));

    center[0] = 0.1232724954;
    center[1] = 4.9194253981;
    center[2] = 4.3227897709;
    atoms.push_back(new Hydrogen(center));

    center[0] = -2.1543745534;
    center[1] = -2.5083334590;
    center[2] = 4.4028310048;
    atoms.push_back(new Hydrogen(center));

    center[0] = 0.1458868463;
    center[1] = -4.9299171568;
    center[2] = 4.3174834203;
    atoms.push_back(new Hydrogen(center));

    center[0] = -0.0058222458;
    center[1] = -2.3339911188;
    center[2] = 0.2206576344;
    atoms.push_back(new Silicon(center));

    center[0] = -0.0084073909;
    center[1] = 2.3313058182;
    center[2] = 0.2343033458;
    atoms.push_back(new Silicon(center));

    center[0] = -7.2506688553;
    center[1] = -2.2726789591;
    center[2] = 0.0000000000;
    atoms.push_back(new Silicon(center));

    center[0] = -7.2506688553;
    center[1] = 2.2726789591;
    center[2] = 0.0000000000;
    atoms.push_back(new Silicon(center));

    center[0] = 7.2506669656;
    center[1] = -2.2726789591;
    center[2] = 0.0000000000;
    atoms.push_back(new Silicon(center));

    center[0] = 7.2506669656;
    center[1] = 2.2726789591;
    center[2] = 0.0000000000;
    atoms.push_back(new Silicon(center));

    center[0] = -10.8760023381;
    center[1] = -3.6253334828;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = -3.6253334828;
    center[1] = -3.6253334828;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = 3.6253334828;
    center[1] = -3.6253334828;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = 10.8760023381;
    center[1] = -3.6253334828;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = -10.8760023381;
    center[1] = 3.6253353725;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = -3.6253334828;
    center[1] = 3.6253353725;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = 3.6253334828;
    center[1] = 3.6253353725;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = 10.8760023381;
    center[1] = 3.6253353725;
    center[2] = -2.1790997282;
    atoms.push_back(new Silicon(center));

    center[0] = -3.6253334828;
    center[1] = 0.0000000000;
    center[2] = -4.7425980682;
    atoms.push_back(new Silicon(center));

    center[0] = -5.9250204057;
    center[1] = 0.0000000000;
    center[2] = -6.4544819433;
    atoms.push_back(new Hydrogen(center));

    center[0] = 3.6253334828;
    center[1] = 0.0000000000;
    center[2] = -4.7425980682;
    atoms.push_back(new Silicon(center));

    center[0] = 5.9250204057;
    center[1] = 0.0000000000;
    center[2] = -6.4544819433;
    atoms.push_back(new Hydrogen(center));

    center[0] = 0.0000000000;
    center[1] = 0.0000000000;
    center[2] = -7.3060964083;
    atoms.push_back(new Silicon(center));

    center[0] = 0.0000000000;
    center[1] = -2.2898244430;
    center[2] = -9.0169144778;
    atoms.push_back(new Hydrogen(center));

    center[0] = 0.0000000000;
    center[1] = 2.2898244430;
    center[2] = -9.0169144778;
    atoms.push_back(new Hydrogen(center));

    center[0] = -7.2506688553;
    center[1] = -3.1938825836;
    center[2] = 1.6640927048;
    atoms.push_back(new Hydrogen(center));

    center[0] = -7.2506688553;
    center[1] = 3.1938825836;
    center[2] = 1.6640927048;
    atoms.push_back(new Hydrogen(center));

    center[0] = 7.2506669656;
    center[1] = -3.1938825836;
    center[2] = 1.6640927048;
    atoms.push_back(new Hydrogen(center));

    center[0] = 7.2506669656;
    center[1] = 3.1938825836;
    center[2] = 1.6640927048;
    atoms.push_back(new Hydrogen(center));

    center[0] = -12.4605526966;
    center[1] = -3.1641137301;
    center[2] = -1.2278362324;
    atoms.push_back(new Hydrogen(center));

    center[0] = -10.8760023381;
    center[1] = -5.5020392585;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = -10.8760023381;
    center[1] = -2.6924759181;
    center[2] = -3.8394507756;
    atoms.push_back(new Hydrogen(center));

    center[0] = -3.6253334828;
    center[1] = -5.5020392585;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 3.6253334828;
    center[1] = -5.5020392585;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 10.8760023381;
    center[1] = -5.5020392585;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 10.8760023381;
    center[1] = -2.6924759181;
    center[2] = -3.8394507756;
    atoms.push_back(new Hydrogen(center));

    center[0] = 12.4605526966;
    center[1] = -3.1641137301;
    center[2] = -1.2278362324;
    atoms.push_back(new Hydrogen(center));

    center[0] = -12.4605526966;
    center[1] = 3.1641156199;
    center[2] = -1.2278362324;
    atoms.push_back(new Hydrogen(center));

    center[0] = -10.8760023381;
    center[1] = 5.5020411482;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = -10.8760023381;
    center[1] = 2.6924778079;
    center[2] = -3.8394507756;
    atoms.push_back(new Hydrogen(center));

    center[0] = -3.6253334828;
    center[1] = 5.5020411482;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 3.6253334828;
    center[1] = 5.5020411482;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 10.8760023381;
    center[1] = 5.5020411482;
    center[2] = -2.5045747940;
    atoms.push_back(new Hydrogen(center));

    center[0] = 10.8760023381;
    center[1] = 2.6924778079;
    center[2] = -3.8394507756;
    atoms.push_back(new Hydrogen(center));

    center[0] = 12.4605526966;
    center[1] = 3.1641156199;
    center[2] = -1.2278362324;
    atoms.push_back(new Hydrogen(center));

    return 191;
}

/** \brief Read in the occupied states from file. */
void read_states(int nstate, int nbasis, Tensor<double> &coeffs) {
    // open file dens.dat
    std::ifstream file;
    file.open("dens.dat");
    if(!file.good()) {
        error("Error opening dens.dat");
    }

    int curstate = 0;
    int i;
    char cbuffer[100];
    double dbuffer;

    // coefficients are output in groups of 5 in the file
    for(curstate = 0; curstate < nstate; curstate += 5) {
        // blow past three useless lines
        file.getline(cbuffer, 100);
        file.getline(cbuffer, 100);
        file.getline(cbuffer, 100);

        // start reading real lines
        for(i = 0; i < nbasis; ++i) {
            // some useless stuff at the front of each line
            file.read(cbuffer, 14);

            file >> dbuffer;
            coeffs(curstate, i) = dbuffer;

            file >> dbuffer;
            if(curstate + 1 < nstate)
                coeffs(curstate + 1, i) = dbuffer;

            file >> dbuffer;
            if(curstate + 2 < nstate)
                coeffs(curstate + 2, i) = dbuffer;

            file >> dbuffer;
            if(curstate + 3 < nstate)
                coeffs(curstate + 3, i) = dbuffer;

            file >> dbuffer;
            if(curstate + 4 < nstate)
                coeffs(curstate + 4, i) = dbuffer;

            // flush the rest of the line
            file.getline(cbuffer, 100);
        }

        // there's a blank line in between each block of five
        file.getline(cbuffer, 100);
    }

    // close the file
    file.close();
}
