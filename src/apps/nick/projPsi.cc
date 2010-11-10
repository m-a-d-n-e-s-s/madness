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
//\file projPsi.cc
//\brief Projects a time evolved wave function onto an arbitrary number of bound states
/***************************************************************************************
 * By: Nick Vence
 * This code must handled with care for the following reasons:
 * 1) It uses the following Libraries:
 *    GNU Scientific Library      http://www.gnu.org/software/gsl/
 *    GNU Multiprecision Library  http://www.gmplib.org/
 *    MultiPrecision Floating point (with correct Rounding)  http://www.mpfr.org/
 *    MADNESS needs to be configured with
 *    ./configure LIBS="-lgsl -lgslblas -lmpfr -lgmp"
 * 2) Is designed modularly. Uncomment the desired functions as desired
 *    projectPsi: Loads wave functions from disk, projects them onto an arbitrary basis
 *    projectZdip:For perturbation calculations
 *    printBasis: For debugging purposes
 *    belkic:     Reproduces an analytic integral
 * 3) projectPsi needs the following files
 *    wf.num                      a list of wave function ID numbers
 *    bound.num OR unbound.num    a list of states to project on
 ***************************************************************************************/

#include "wavef.h"
#include "extra.h"
#include <mra/lbdeux.h>
#include <string>
#include <fstream>
using std::ofstream;
using std::ofstream;
#include <stdlib.h>

using namespace madness;

const double PI = M_PI;
const complexd I(0,1);
const complexd one(1,0);

struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0) 
        : leaf_value(leaf_value)
        , parent_value(parent_value) 
    {}

    double operator()(const Key<3>& key, const FunctionNode<double_complex,3>& node) const {
        if (key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};


/******************************************************
 * <Yl0|Psi(t)>
 * Needs: input2 
 * Parsing is set up to do only one time step at a time
 ******************************************************/
void projectL(World& world, const double L, const int wf, const int n, const int lMAX) {
    complexd output;
    PRINTLINE("\t\t\t\t\t\t|<Yl0|Psi(t)>|^2 ");
    PRINT("\t\t\t\t\t\t");
    //LOAD Psi(t)
    complex_functionT psi;
    if( !wave_function_exists(world, wf) ) {
        PRINTLINE("Function " << wf << " not found");
        exit(1);
    } else {
        psi = wave_function_load(world, wf);
        psi.reconstruct();
        LoadBalanceDeux<3> lb(world);
        lb.add_tree(psi, LBCost(1.0,1.0));
        FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0,false));
        PRINTLINE("|" << wf << ">\t\t");
    }
    PRINTLINE("");
    const double PI = M_PI;
    const double dr = L/n;
    const int smallN = n;
    const double dTH = PI/smallN;
    const double dPHI = 2*PI/smallN;
    const bool printR = true;
    clock_t before=0, after=0, middle=0;
    std::vector<WF>::iterator psiT;
    std::vector<complexd> YlPsi(n);
    for( int l=0; l<=lMAX; l++) {
        PRINT("Y"<< l << "0: \t\t\t\t\t\t");
        psi.reconstruct();
        if(world.rank()==0) before = clock();
        Yl0 yl0(L, l);
        for( int i=0; i<n; i++ ) {
            YlPsi[i] = 0.0;
        }
        for( int i=world.rank(); i<n; i+=world.size() ) {
            const double r = (0.5 + i)*dr;
            complexd Rl = 0.0;
            for( int j=0; j<smallN ; j++ ) {
                const double th = (0.5 + j)*dTH;
                const double sinTH = std::sin(th);
                for( int k=0; k<smallN; k++ ) {
                    const double phi = k*dPHI;
                    const double a[3] = {r*sinTH*std::cos(phi), r*sinTH*std::sin(phi), r*std::cos(th)};
                    const vector3D rVec(a);
                    Rl += psi.eval(rVec).get() * yl0(rVec) * sinTH*dTH*dPHI;
                }
            }
            YlPsi[i] = conj(Rl)*Rl * r*r*dr;
            if(printR) PRINT(std::real(YlPsi[i]) << "\t");
        }
        if(printR) PRINTLINE("");
        if(world.rank()==0) middle = clock();
        world.gop.sum(&YlPsi[0], n);
        world.gop.fence();
        double Pl = 0.0;
        for( int i=0; i<n; i++ ) {
            Pl += real( YlPsi[i] );
        }
        if(world.rank()==0) after = clock();
        //PRINT(" Integration took " << (middle - before)/CLOCKS_PER_SEC << " seconds ");
        PRINTLINE(std::setprecision(6) << std::scientific << Pl);
    }
    PRINTLINE("");
}


/****************************************
 * Needs: input input2
 ****************************************/
void zSlice(World& world, const int n, double L, double th, double phi, const int wf) {
    complex_functionT psiT;
    PRINTLINE(std::setprecision(2) << std::fixed);
    if( !wave_function_exists(world, wf) ) {
        PRINTLINE("Function " << wf << " not found");
        exit(1);
    } else {
        psiT = wave_function_load(world, wf);
        PRINTLINE("phi(T=" << wf << ",r) =\t th=0 \t\t\t th=" << th << "  phi = " << phi);
    }// done loading wf
    complexd output;
    const double dr = L/n;
    for( int i=0; i<n; i++ ) {
        const double r = i*dr;
        const double a[3] = {0, 0, r};
        const double b[3] = {r*std::sin(th)*std::sin(phi), r*std::sin(th)*std::cos(phi), r};
        const vector3D aVec(a);
        const vector3D bVec(b);
        double psiA = real(psiT(aVec));
        double psiB = real(psiT(bVec));
        PRINT(std::fixed << std::setprecision(2));
        PRINTLINE(r << " \t\t " << std::scientific << std::setprecision(6) << psiA << " \t " << psiB);
    }
}


//Clunky code! Design out of this!
struct PhiKAdaptor : public FunctionFunctorInterface<std::complex<double>,3> {
    PhiK& phik;
    PhiKAdaptor(PhiK& phik) : phik(phik) {}
    std::complex<double> operator()(const vector3D& x) const {
        return phik(x);
    }
};

//Testing projPsi against mathematica via an integral that is very local and analytic
void testIntegral(World& world, double L, const double Z, double k) {
    double arr[3] = {0, 0, k};
    const vector3D kVec(arr);
    const double constCutoff = L;
    PRINTLINE("    PhiK phik = PhiK(world, Z, kVec, constCutoff);");
    PhiK phik = PhiK(world, Z, kVec, constCutoff);
    PRINTLINE("    phik.Init(world);");
    phik.Init(world);
    PRINTLINE("    complex_functionT phiK = complex_factoryT(world).functor(functorT( new PhiKAdaptor(phik) ));");
    complex_functionT phiK = complex_factoryT(world).functor(functorT( new PhiKAdaptor(phik) ));
    complexd output;
    for( int i=1; i<=5; i++ ) {
        double a = 0.1*i*i*i*i;
        complex_functionT guass = complex_factoryT(world).functor(functorT( new Gaussian(a) ));
        output = inner(phiK, guass);
        PRINT(std::setprecision(1) << std::fixed << "<k=" << k << "|exp(-" << a << "r^2)> = ");
        PRINTLINE( std::setprecision(8) << output);
    }
}

void debugSlice(World& world, const int n, double L, double Z, double k) {
    const double dr = L/n;
    const double a[3] = {0, 0, k};
    const vector3D kVec(a);
    ofstream fout;
    fout.open("phi0.5.dat");
    PhiK phi = PhiK(world, Z, kVec, L);
    phi.Init(world);
    for( int i=1; i<n; i++ ) {
        double r = i*dr;
        const double b[3] = {0, 0, r};
        const vector3D rVec(b);
        fout << std::fixed << std::setprecision(2)
             << r << " \t " << std::scientific << std::setprecision(16)
             << real(phi(rVec)) << " \t " << imag(phi(rVec)) << std::endl;
    }
}


 
/************************************************************************************
 * The correlation amplitude |<Psi(+)|basis>|^2 are dependent on the following files:
 * wf.num                  Integer time step of the Psi(+) to be loaded
 * bound.num               Integer triplets of quantum numbers   2  1  0 
 * unbound.num             Double triplets of momentum kx ky kz  0  0  0.5
 ************************************************************************************/
void projectPsi(World& world, std::vector<std::string> boundList, std::vector<std::string> unboundList, const double Z, double cutoff) {
    PRINTLINE("\t\t\t\t\t\t|<basis|Psi(t)>|^2 ");
    std::ifstream f("wf.num");
    if( !f.is_open() ) {
        PRINTLINE("File: wf.num expected to contain a list of integers of loadable wave functions");
    } else {
        if(boundList.empty() && unboundList.empty()) {
            boundList.push_back("1 0 0");
            boundList.push_back("2 1 0");
        }
        //LOAD Psi(t)
        std::string tag;
        std::vector<WF> psiList;
        complexd output;
        PRINT("\t\t\t\t\t\t");
        while(f >> tag) {
            if( !wave_function_exists(world, atoi(tag.c_str())) ) {
                PRINTLINE("Function " << tag << " not found");
            } else {
                WF psi_t = WF(tag, wave_function_load(world, atoi(tag.c_str())));
                psiList.push_back(WF(tag, psi_t.func));
                PRINT("|" << tag << ">\t\t");
            }
        }// done loading wf.num
        PRINT("\n");
        //psiIT holds the time evolved wave functions
        //LOAD bound states
        complex_functionT psi0;
        if( wave_function_exists(world, 0) ) {
            psi0 = wave_function_load(world,0);
        } else {
            PRINTLINE("psi0 must be present in this directory i.e.  data-00000.0000*");
            PRINTLINE("Consider changing the restart file to 0 and rerunning tdse");
            exit(1);
        }
        std::vector<WF>::const_iterator psiIT; 
        if( !boundList.empty() ) {
            // <phi_bound|Psi(t)>
            std::vector<std::string>::const_iterator boundIT;
            int N, L, M;
            std::cout.precision(6);
            for(boundIT = boundList.begin(); boundIT !=boundList.end(); boundIT++ ) {
                std::stringstream ss(*boundIT);
                ss >> N >> L >> M;
                //PROJECT Phi_nlm into MADNESS
                complex_functionT phi_nlm = complex_factoryT(world).
                    functor(functorT( new BoundWF(Z, N, L, M)));
                complexd n_overlap_0 = inner(phi_nlm,psi0);
                PRINT(*boundIT << "\t\t\t\t\t");
                //loop through time steps
                for( psiIT=psiList.begin(); psiIT !=  psiList.end(); psiIT++ ) {
                    //|PSI(t)> = |Psi(t)> - <phi_k|Psi(0)>|Psi(0)>
                    //<phi_nl|PSI(t)> = <phi_nl|Psi(t)> - <phi_nl||Psi(0)> <Psi(0)|Psi(t)>
                    output =  inner(phi_nlm, psiIT->func) - n_overlap_0  * inner(psi0,psiIT->func);
                    PRINT(std::scientific <<"\t" << real(output*conj(output)));
                }
                PRINT("\n");
            }            
        }
        clock_t before=0, after=0;
        //LOAD unbound states
        if( !unboundList.empty() ) {
            std::vector<std::string>::const_iterator unboundIT;
            for( unboundIT=unboundList.begin(); unboundIT !=  unboundList.end(); unboundIT++ ) {
                //parsing unboundList
                double KX, KY, KZ;
                std::stringstream ss(*unboundIT);
                ss >> KX >> KY >> KZ;
                double dArr[3] = {KX, KY, KZ};
                const vector3D kVec(dArr);
                std::cout.precision(12);
                PRINT( std::fixed << KX << " " << KY << " " << KZ << "  ");
                if((dArr[1]>0.0 || dArr[1]<0.0) || (dArr[2]>0.0 || dArr[2]<0.0)) { //removing k={0,0,0}
                    //PROJECT Psi_k into MADNESS
                    if(world.rank()==0) before = clock();
                    const double constcutoff = cutoff;
                    PhiK phik = PhiK(world, Z, kVec, constcutoff);
                    phik.Init(world);
                    complex_functionT phi_k = complex_factoryT(world).functor(functorT( new PhiKAdaptor(phik) ));
                    // W/O timing
                    //complex_functionT phi_k = 
                    //complex_factoryT(world).functor(functorT( new PhiK(Z, kVec, cutoff) ));
                    if(world.rank()==0) after = clock();
                    std::cout.precision( 8 );
                     //<phi_k|Psi(0)>
                    complexd k_overlap_0 = inner(phi_k,psi0);
                    //loop through time steps
                    for( psiIT=psiList.begin(); psiIT !=  psiList.end(); psiIT++ ) {
                        //|PSI(t)> = |Psi(t)> - <phi_k|Psi(0)>|Psi(0)>
                        //<phi_k|PSI(t)> = <phi_k|Psi(t)>   - <phi_k||Psi(0)> <Psi(0)|Psi(t)>
                        output =  inner(phi_k, psiIT->func) - k_overlap_0  * inner(psi0,psiIT->func);
                        PRINT( std::scientific << "\t" << real(output*conj(output)) );
                    }
                    PRINT(" took " << (after - before)/CLOCKS_PER_SEC << " seconds ");
                    int WFsize = phi_k.size();
                    PRINT("and has " << WFsize << " coefficients.\n");
                }
            }
        }
    }
}

void loadParameters2(World& world, int &n, double& th, double& phi, int& wf, double& k, int& lMAX) {
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
            else if (tag == "n") {
                f >> n;
                PRINTLINE("n = " << n);
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
            else if (tag == "k") {
                f >> k;
                PRINTLINE("k = " << k);
            }
            else if (tag == "lMAX") {
                f >> lMAX;
            }
        }
    }
    f.close();
}

void loadParameters(World& world, double& thresh, int& k, double& L, double &Z, double &cutoff) {
    std::string tag;
    int natom;
    double Rx, Ry, Rz;
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
                f >> k;
                PRINTLINE("k = " << k);
            }
            else if (tag == "natom") {
                f >> natom;
                f >> Z >> Rx >> Ry >> Rz;
                PRINTLINE("Z = " << Z);
            }
            else if (tag == "omega") {
                double omega;
                f >> omega;
                ///E_n = n hbar omega
                ///I_p = 0.5 Z^2
                ///KE = E_n - I_p
                ///Atomic Units: hbar = m = 1
                ///v = sqrt( 2n omega - Z^2)
                ///cutoff > dMAX = v t
                double dMAX = std::sqrt(2*3*omega - Z*Z) * 10;
                cutoff = 0.0;
                while( cutoff < dMAX ) { cutoff += L/16; }
                PRINTLINE("dMAX = " << dMAX);
                PRINTLINE("cutoff = " << cutoff);
            }
        }
    }
    f.close();
}

int main(int argc, char**argv) {
    // INITIALIZE the parallel programming environment
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    PRINTLINE("After initialize");
    // Load info for MADNESS numerical routines
    startup(world,argc,argv);
    PRINTLINE("world.size() = " << world.size());
    // Setup defaults for numerical functions
    int    k   = 12;
    double L   = 1.0;
    double Z   = 1.0;
    double thresh = 1e-6;
    double cutoff = L;
    double th  = 0.0;
    double phi = 0.0;
    double kMomentum   = 1.0;
    int    n   = 10;
    int    wf  = 0;
    int   lMAX = 0;
    loadParameters(world, thresh, k, L, Z, cutoff);
    loadParameters2(world, n, th, phi, wf, kMomentum, lMAX);
    FunctionDefaults<NDIM>::set_k(k);                 // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(thresh);       // Accuracy
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_apply_randomize(false);
    FunctionDefaults<NDIM>::set_autorefine(false);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    FunctionDefaults<NDIM>::set_pmap(pmapT(new LevelPmap(world)));
    FunctionDefaults<NDIM>::set_truncate_on_project(true);
    try {
        std::vector<std::string> boundList;
        std::vector<std::string> unboundList;
        projectL(world, L, wf, n, lMAX);
        //zSlice(world, n1, L, th, phi);
        //testIntegral(world, L, Z, kMomentum);
        //debugSlice(world, n1, L, Z, kMomentum);
        loadList(world, boundList, unboundList);
        projectPsi(world, boundList, unboundList, Z, cutoff);
        //PRINTLINE("Z = " << Z);
        //std::vector<WF> boundList;
        //std::vector<WF> unboundList;
        //compareGroundState(world, Z);
        //compare1F1(world, cutoff);
        //printBasis(world, Z, cutoff);
        //loadBasis(world,  boundList, unboundList, Z, cutoff);
        //belkic(world, cutoff);
        //projectZdip(world, unboundList);
        //groundOverlap(world, boundList, unboundList, Z, cutoff);
        world.gop.fence();
//         if (world.rank() == 0) {
//             world.am.print_stats();
//             world.taskq.print_stats();
//             world_mem_info()->print();
//             WorldProfile::print(world);
//         }
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
