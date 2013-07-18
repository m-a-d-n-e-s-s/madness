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
/*********************
 * Trying to isolate a memory error
 *********************/
#include <fstream>
using std::ofstream;
using std::ios;
#include <iomanip>
#include <string>
using std::string;
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "wavef.h"
#define PRINT(str) if(world.rank()==0) std::cout << str
#define PRINTLINE(str) if(world.rank()==0) std::cout << str << std::endl

using namespace madness;

typedef std::complex<double> complexd;
typedef Vector<double,NDIM> vector3D;
typedef Function<complexd,NDIM> complex_functionT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<complexd,NDIM> complex_factoryT;
typedef FunctionFactory<double,NDIM> factoryT;
typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
const char* wave_function_filename(int step);

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
    return fname;
}
bool wave_function_exists(World& world, int step) {
    return archive::ParallelInputArchive::exists(world, wave_function_filename(step));
}
void wave_function_store(World& world, int step, const complex_functionT& psi) {
    archive::ParallelOutputArchive ar(world, wave_function_filename(step), nIOProcessors);
    ar & psi;
}
complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    archive::ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

// Wrappers & structs
struct PhiKAdaptor : public FunctionFunctorInterface<std::complex<double>,3> {
    PhiK& phik;
    PhiKAdaptor(World& world, PhiK& phik) : phik(phik) {
        phik.Init(world);
    }
    std::complex<double> operator()(const vector3D& x) {
        return phik(x);
    }
};

struct WF {
    int t;
    complex_functionT func;
    WF(const int& timeStep, const complex_functionT& FUNC) : t(timeStep), func(FUNC) {
        func.truncate();
    }
};

void loadParameters2(World& world, int &nGrid, double& th, double& phi, int& wf, double& kMomentum, int& lMAX, int& nPhoton) {
    std::string tag;
    std::ifstream f("input2");
    std::cout << std::scientific;
    if( f.is_open() ) {
        while(f >> tag) {
            if (tag[0] == '#') {
                char ch;
                PRINTLINE("    comment  " << tag.c_str());
                while (f.get(ch)) {
                    PRINTLINE(ch);
                    if (ch == '\n') break;
                }
            }
            else if (tag == "nGrid") {
                f >> nGrid;
                PRINTLINE("nGrid = " << nGrid);
            }
            else if (tag == "th") {
                f >> th;
            }
            else if (tag == "phi") {
                f >> phi;
            }
            else if (tag == "wf") {
                f >> wf;
            }
            else if (tag == "kMomentum") {
                f >> kMomentum;
                PRINTLINE("kMomentum = " << kMomentum);
            }
            else if (tag == "lMAX") {
                f >> lMAX;
                PRINTLINE("lMAX = " << lMAX);
            }
            else if (tag == "nPhoton") {
                f >> nPhoton;
                PRINTLINE("nPhoton = " << nPhoton);
            }
        }
    }
}

void loadParameters(World& world, double& thresh, int& kMAD, double& L, double &Z, int&
        nPhoton, double& cutoff, double& Llarge) {
    std::string tag;
    double Rx, Ry, Rz;
    double omega;
    double tMAX;
    int natom;
    std::ifstream f("input");
    std::cout << std::scientific;
    if( f.is_open() ) {
        while(f >> tag) {
            if (tag[0] == '#') {
                char ch;
                PRINTLINE("    comment  " << tag.c_str());
                while (f.get(ch)) {
                    PRINTLINE(ch);
                    if (ch == '\n') break;
                }
            }
            else if (tag == "L") {
                f >> L;
                PRINTLINE("L = " << L);
            }
            else if (tag == "thresh") {
                f >> thresh;
                PRINTLINE("thresh = " << thresh);
            }
            else if (tag == "k") {
                f >> kMAD;
                PRINTLINE("kMAD = " << kMAD);
            }
            else if (tag == "natom") {
                f >> natom;
                f >> Z >> Rx >> Ry >> Rz;
                PRINTLINE("Z = " << Z);
            }
            else if (tag == "omega") {
                f >> omega;
            }
            else if (tag == "target_time") {
                f >> tMAX;
                PRINTLINE( "tMAX = " << tMAX );
                /// Atomic Units: m = hbar = e = 1
                /// E_n     = n hbar omega
                /// I_p     = Z^2 / 2
                /// KE      = E_n - I_p
                /// 0.5 v^2 = n omega - Z^2/2
                ///      v  = sqrt(2 n omega - Z^2)
                /// dMAX    = v tMAX
                PRINTLINE( "omega = " << omega );
                while( 2*(nPhoton*omega - Z*Z) < 0.0 ) nPhoton++; //just in case nPhoton is too small
                PRINTLINE("nPhoton = " << nPhoton);
                PRINTLINE("The following shoud be positive: 2*nPhoton*omega - Z*Z = " << 2*nPhoton*omega - Z*Z);
                double dMAX = std::sqrt( 2*(nPhoton*omega - Z*Z)) * tMAX;
                PRINTLINE("dMAX = " << dMAX );
                if ( cutoff < dMAX ) cutoff = 0.0;
                while( cutoff < dMAX ) {
                    cutoff += L/128;
                }
                PRINTLINE( "cutoff = " << cutoff );
            }
            else if (tag == "Llarge") {
                f >> Llarge;
                PRINTLINE("Llarge = " << Llarge);
            }
        }
    }
    f.close();
}

string toString(int& i) {
    std::stringstream ss;
    ss << i;
    return ss.str();
}

void doplot(World& world, int step, const complex_functionT& psi, double Lplot, long numpt, const char* fname) {
    double start = wall_time();
    Tensor<double> cell(3,2);
    std::vector<long> npt(3, numpt);
    cell(_,0) = -Lplot;
    cell(_,1) =  Lplot;
    plotdx(psi, fname, cell, npt);
    if (world.rank() == 0) print("plotting used", wall_time()-start);
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);
    // Setup defaults for numerical functions
    double thresh = 1e-3;
    int    k  = 12;
    double L = M_PI;
    double Z = 1.0;
    double cutoff = L;
    double Llarge = L;
    int    nGrid = 5;
    double th = M_PI;
    double phi = 1.0;
    int    wf = 0;
    double kMomentum = 1.0;
    int    lMAX  = 12;
    int    nPhoton = 2;
    bool   asciiFile = false;
    bool   coordEdge = false;
    bool   plotDX    = false;
    loadParameters(world, thresh, k, L, Z, nPhoton, cutoff, Llarge);
    loadParameters2(world, nGrid, th, phi, wf, kMomentum, lMAX, nPhoton);
    FunctionDefaults<NDIM>::set_k(k);               // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(1e-3);       // Accuracy
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_apply_randomize(false);
    FunctionDefaults<NDIM>::set_autorefine(true);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    //FunctionDefaults<NDIM>::set_pmap(pmapT(new LevelPmap(world)));
    FunctionDefaults<NDIM>::set_truncate_on_project(true);
    try{
        //LOAD Psi(t) -> psiIT
        PRINTLINE("cutoff = " << cutoff);
        PRINTLINE("Llarge = " << Llarge);
        std::ifstream f("wf.num");
        if( !f.is_open() ) {
            PRINTLINE("File: wf.num expected to contain a list of integers of loadable wave functions");
        } else {
            std::string step;
            std::vector<WF> psiList;
            while(f >> step) {
                if( !wave_function_exists(world, atoi(step.c_str())) ) {
                    PRINTLINE("Function " << step << " not found");
                } else {
                    PRINTLINE("phi(T=" << step << ",x,z)");
                    int nWF = atoi(step.c_str());
                    const complex_functionT psiT = wave_function_load(world, nWF);
                    double Lplot = 60;
                    long numpt = nGrid;
                    if(plotDX) doplot(world, nWF, psiT, Lplot, numpt,
                            ("wf" + step + ".dx").c_str());
                    FILE* binaryFile = fopen(("wf" + step + ".bin").c_str(), "wb");
                    ofstream asciiFile( ("wf" + step + ".dat").c_str());

                    //Compute grid
                    ///xMax: the largest x & z coordinate
                    ///An odd number of grid points ensures zero is a gird point
                    int n = 2*nGrid + 1;
                    const double xMax = Llarge;
                    const double dr = xMax/nGrid;
                    float buffer[n+1];
                    // coordEdge pads the matrix with the coordinates
                    // otherwise no coordinates will be provided
                    if( coordEdge ) {
                    ///Output for GNUPLOT binary matrix
                    ///<n+1> <y0>   <y1>   <y2>   ...  <yN>
                    ///<x0>  <x0,0> <z0,1> <z0,2> ... <z0,N>
                    ///<x1>  <x1,0> <z1,1> <z1,2> ... <z1,N>
                    /// .      .      .      .    ...   .
                    /// .      .      .      .    ...   .  float buffer[n+1];
                     // initialize first line
                        buffer[0] = (float) n+1;
                        if( asciiFile ) asciiFile << buffer[0] << "\t";
                        for( int i=0; i<n; i++ ) {
                            buffer[i+1] = i*dr - xMax;
                            if( asciiFile ) asciiFile << buffer[i+1] << "\t";
                        }
                        fwrite(buffer, sizeof(float), n+1, binaryFile);
                        if( asciiFile ) asciiFile << std::endl;
                    }
                    // initialize the bulk
                    for( int i=0; i<n; i++ ) {
                        const double x = i*dr - xMax;
                        if( coordEdge ) {
                            buffer[0] = x;
                            if( asciiFile ) asciiFile <<  x << "\t";
                        }
                        for( int j=0; j<n; j++ ) {
                            const double z = j*dr - xMax;
                            const double r[3] = {x, 0, z};
                            const vector3D         rVec(r);
                            if( coordEdge ) {
                                buffer[j+1] =  (float) real(psiT(rVec));
                            } else {
                                buffer[j] =  (float) real(psiT(rVec));
                            }
                            if( asciiFile ) asciiFile <<  real(psiT(rVec)) << "\t";
                        }
                        if( world.rank() == 0 ) {
                            if( coordEdge ) {
                                fwrite(buffer, sizeof(float), n+1, binaryFile);
                            } else {
                                fwrite(buffer, sizeof(float),   n, binaryFile);
                            }
                            if( asciiFile ) asciiFile << std::endl;
                        }
                        world.gop.fence();
                    }
                    fclose(binaryFile);
                    asciiFile.close();
                    world.gop.fence();
                }
            }// done loading wf.num
        }
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
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }
    world.gop.fence();
    finalize();				//FLAG
    return 0;
}
