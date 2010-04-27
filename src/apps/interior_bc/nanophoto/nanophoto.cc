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

void scaled_plotvtk_begin(World &world, const char *filename,
    const Vector<double, 3> &plotlo, const Vector<double, 3> &plothi,
    const Vector<long, 3> &npt, bool binary = false);
int mol_geom(std::vector<Atom*> &atoms);
void read_states(int nstate, int nbasis, Tensor<double> &coeffs);

using namespace madness;

int main(int argc, char **argv) {
    double eps;
    int j, k, n;
    double thresh, phi, d;
    unsigned int i;
    char filename[100];
    char funcname[15];

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

        phi = atof(argv[4]) * 3.6749322e-5; // convert to a.u.

        d = atof(argv[5]) / 0.052918; // convert to a.u.
    }
    world.gop.broadcast(phi);
    world.gop.broadcast(eps);
    world.gop.broadcast(thresh);
    world.gop.broadcast(k);
    world.gop.broadcast(d);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    //FunctionDefaults<3>::set_max_refine_level(6);

    // box size
    Tensor<double> cell(3, 2);
    cell(0,0) = cell(1,0) = -175.0 / 0.052918;
    cell(0,1) = cell(1,1) = 175.0 / 0.052918;
    cell(2,0) = 10.0 * (constants::pi - 5.0) / 0.052918;
    cell(2,1) = 150.0 / 0.052918;
    FunctionDefaults<3>::set_cell(cell);

    // make the basis functions to get the density
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

    // the key data structure: sets up the problem details and projects
    // the initial functions
    SharedPtr<TipMolecule> tpm(new TipMolecule(eps, coeffs, atoms,
        basis, phi, d));

    // ********************** START electrostatics solver
    sprintf(funcname, "potential");
    if(world.rank() == 0) {
        // print out the arguments
        printf("Tip-Surface Distance: %.6e nm\nPotential Difference: %.6e " \
            "mV\nk: %d\nthresh: %.6e\neps: %.6e nm\n", d * 0.052918,
            phi / 3.6749322e-5, k, thresh, eps * 0.052918);

        // project the surface function
        printf("Projecting the surface function (to low order)\n");
        fflush(stdout);
    }

    tpm->fop = DOMAIN_MASK;
    real_function_3d usol = real_factory_3d(world).k(6).thresh(1.0e-4).
        functor(tpm);

#if 0
    tpm->fop = SURFACE;
    real_function_3d surf = real_factory_3d(world).k(6).thresh(1.0e-4).
        functor(tpm);

    if(world.rank() == 0) {
        printf("Performing load balancing\n");
        fflush(stdout);
    }
    tpm->load_balance(world, surf);

    // reproject the surface function to the requested threshold / k
    if(k > 6 || thresh < 1.0e-4) {
        if(world.rank() == 0) {
            printf("Reprojecting the surface function to requested order\n");
            fflush(stdout);
        }
        surf.clear();
        surf = real_factory_3d(world).functor(tpm);
    }

    // green's function
    // note that this is really -G...
    real_convolution_3d G = BSHOperator<3>(world, 0.0, eps*0.1, thresh);

    // project the r.h.s. function (phi*f - penalty*S*g)
    // and then convolute with G
    real_function_3d usol, rhs;
    if(world.rank() == 0) {
        printf("Projecting the r.h.s. function\n");
        fflush(stdout);
    }

    tpm->fop = DIRICHLET_RHS;
    usol = real_factory_3d(world).functor(tpm);
    rhs = G(usol);
    rhs.truncate();
    usol.clear();

    // make an initial guess:
    // uguess = rhs / penalty_prefact
    // the rescaling will make operator(uguess) close to rhs in magnitude for
    //     starting in GMRES
    usol = copy(rhs);
    usol.scale(2.0 / eps);
    usol.compress();
    DirichletCondIntOp<3> dcio(G, surf);

    // make the operators and prepare GMRES
    FunctionSpace<double, 3> space(world);
    double resid_thresh = 1.0e-5;
    double update_thresh = 1.0e-5;
    int maxiter = 30;
    GMRES(space, dcio, rhs, usol, maxiter, resid_thresh, update_thresh, true);
    // ********************** END electrostatics solver
#endif


#if 0
    // ********************** START density plotter
    if(world.rank() == 0) {
        print("projecting the density");
        fflush(stdout);
    }
    tpm->fop = DENSITY;
    real_function_3d usol = real_factory_3d(world).functor(tpm);

    double norm = usol.norm2();
    double trace = usol.trace();
    if(world.rank() == 0)
        print(norm, trace);

    sprintf(funcname, "density");
    // ********************** END density plotter
#endif

    // print out the solution function on the total domain
    sprintf(filename, "%s.vts", funcname);
    Vector<double, 3> plotlo, plothi;
    Vector<long, 3> npts;
    for(int i = 0; i < 3; ++i) {
        plotlo[i] = cell(i, 0);
        plothi[i] = cell(i, 1);
        npts[i] = 71;
    }
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);

    // print out the solution function near the area of interest
    sprintf(filename, "%s-local.vts", funcname);
    plotlo[0] = -20.0 / 0.052918; plothi[0] = 20.0 / 0.052918; npts[0] = 71;
    plotlo[1] = -20.0 / 0.052918; plothi[1] = 20.0 / 0.052918; npts[1] = 71;
    plotlo[2] = -5.0 / 0.052918; plothi[2] = 20.0 / 0.052918; npts[2] = 71;
    scaled_plotvtk_begin(world, filename, plotlo, plothi, npts);
    plotvtk_data(usol, funcname, world, filename, plotlo, plothi, npts);
    plotvtk_end<3>(world, filename);

    finalize();
    
    return 0;
}

/** \brief Same as plotvtk_begin, but scales the coordinates back to nm */
void scaled_plotvtk_begin(World &world, const char *filename,
    const Vector<double, 3> &plotlo, const Vector<double, 3> &plothi,
    const Vector<long, 3> &npt, bool binary) {

    PROFILE_FUNC;

    Tensor<double> cell(3, 2);
    int i;
    for(i = 0; i < 3; ++i) {
        cell(i, 0) = plotlo[i];
        cell(i, 1) = plothi[i];
    }

    FILE *f=0;
    if(world.rank() == 0) {
        f = fopen(filename, "w");
        if(!f)
            MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

        fprintf(f, "<VTKFile type=\"StructuredGrid\" version=\"0.1\"" \
            " byte_order=\"LittleEndian\" compressor=\"" \
            "vtkZLibDataCompressor\">\n");
        fprintf(f, "  <StructuredGrid WholeExtent=\"");
        for(i = 0; i < 3; ++i)
            fprintf(f, "0 %ld ", npt[i]-1);
        for(; i < 3; ++i)
            fprintf(f, "0 0 ");
        fprintf(f, "\">\n");
        fprintf(f, "    <Piece Extent=\"");
        for(i = 0; i < 3; ++i)
            fprintf(f, "0 %ld ", npt[i]-1);
        for(; i < 3; ++i)
            fprintf(f, "0 0 ");
        fprintf(f, "\">\n");
        fprintf(f, "      <Points>\n");
        fprintf(f, "        <DataArray NumberOfComponents=\"3\" " \
            "type=\"Float32\" format=\"ascii\">\n");

        Vector<double, 3> space;
        for(i = 0; i < 3; ++i) {
            if(npt[i] == 1)
                space[i] = 0.0;
            else
                space[i] = (cell(i, 1) - cell(i, 0)) / (npt[i] - 1);
        }

        // go through the grid
        for(LowDimIndexIterator it(npt); it; ++it) {
            for(i = 0; i < 3; ++i)
                fprintf(f, "%f ", (plotlo[i] + it[i]*space[i]) * 0.052918);
            for(; i < 3; ++i)
                fprintf(f, "0.0 ");
            fprintf(f, "\n");
        }

        fprintf(f, "        </DataArray>\n");
        fprintf(f, "      </Points>\n");
        fprintf(f, "      <PointData>\n");
        fclose(f);
    }
    world.gop.fence();
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
