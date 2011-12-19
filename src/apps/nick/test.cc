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
#include <string>
#include <fstream>
using std::ofstream;
#include <stdlib.h>
#include <iomanip>
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
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);

struct PhiKAdaptor : public FunctionFunctorInterface<std::complex<double>,3> {
    PhiK& phik;
    PhiKAdaptor(World& world, PhiK& phik) : phik(phik) {
        phik.Init(world);
    }
    std::complex<double> operator()(const vector3D& x) {
        return phik(x);
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

void loadParameters(World& world, double& thresh, int& kMAD, double& L, double &Z, int& nPhoton, double& cutoff) {
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
        }
    }
    f.close();
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world, argc, argv);
    // Setup defaults for numerical functions
    int    nPhoton = 5;
    int    k = 12;
    double L = M_PI;
    double Z = 1.0;
    double thresh = 1e-3;
    double cutoff = L;
    loadParameters(world,thresh, k, L, Z, nPhoton, cutoff);
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
        double before = 0.0;
        double after = 0.0;
        if(world.rank()==0) before = wall_time();
        double arr[3] = {0, 0, k};
        const vector3D kVec(arr);
        const double constCutoff = L;
        PRINTLINE("    PhiK phik = PhiK(world, Z, kVec, constCutoff);");
        PhiK phik = PhiK(world, Z, kVec, constCutoff);
        PRINTLINE("    phik.Init(world);");
        phik.Init(world);
        PRINTLINE("    complex_functionT phiK = complex_factoryT(world).functor(functorT( new PhiKAdaptor(phik) ));");
        complex_functionT phiK = complex_factoryT(world).functor(functorT( new PhiKAdaptor(phik) ));
        if(world.rank()==0) after = wall_time();
        PRINT(" took " << after - before << " seconds ");
        // const int n = 9;
        // std::vector<int> a(n);
        // for(int i=0; i<n; i++) {
        //     a[i] = 0;
        // }
        // for(int i=world.rank(); i<n; i+=world.size()) {
        //     a[i] = 2*i + 1;
        // }
        // world.gop.sum(&a[0], n);
        // world.gop.fence();
        // if( world.rank()==0 ){
        //     int total = 0;
        //     for(int i=0; i<n; i++) {
        //         total += a[i];
        //     }
        //     std::cout << total << std::endl;
        // }
    } catch (const MPI::Exception& e) {
        //print(e);
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
