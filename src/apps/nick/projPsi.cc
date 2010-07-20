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
 * 3) k (The wavelet order) must be the same as the projected functions: see main()
 *    12 has been the default
 ***************************************************************************************/

#include <string>
#include <fstream>
using std::ofstream;
using std::ofstream;
#include <stdlib.h>
#include <iomanip>
#include <time.h>
#include "wavef.h"
#define PRINT(str) if(world.rank()==0) std::cout << str 
#define PRINTLINE(str) if(world.rank()==0) std::cout << str << std::endl

using namespace madness;

const int nIOProcessors =1;
const std::string prefix = "data";
typedef std::complex<double> complexd;
typedef Vector<double,NDIM> vector3D;
typedef Function<complexd,NDIM> complex_functionT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<complexd,NDIM> complex_factoryT;
typedef FunctionFactory<double,NDIM> factoryT;
typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;
const char* wave_function_filename(int step);
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);

struct InputParameters {
  static const int MAXNATOM=99;

    // IF YOU ADD A NEW PARAMETER DON'T FORGET TO INCLUDE IT IN
    // a) read()
    // b) serialize()
    // c) operator<<()
  
  double L;           // Box size for the simulation
  double Lsmall;      // Box size for small (near nucleus) plots
  double Llarge;      // Box size for large (far from nucleus) plots
  double F;           // Laser field strength
  double omega;       // Laser frequency
  double ncycle;      // Number of laser cycles in envelope
  int natom;          // Number of atoms
  double Z[MAXNATOM]; // Nuclear charge of atoms
  double R[MAXNATOM][3]; // Coordinates of atoms
  int k;              // wavelet order
  double thresh;      // precision for truncating wave function
  double safety;      // additional precision (thresh*safety) for operators and potential
  double cut;         // smoothing parameter for 1/r (same for all atoms for now)
  std::string iState ; // initial state = "1s" or "2s"
  std::string prefix; // Prefix for filenames
  int ndump;          // dump wave function to disk every ndump steps
  int nplot;          // dump opendx plot to disk every nplot steps
  int nprint;         // print stats every nprint steps
  int nloadbal;       // load balance every nloadbal steps
  int nio;            // Number of IO nodes 
  double tScale;      // Scaling parameter for optimization
  double target_time; // Target end-time for the simulation
  
  void read(const char* filename) {
    std::ifstream f(filename);
    std::string tag;
    iState = "1s";
    printf("\n");
    printf("       Simulation parameters\n");
    printf("       ---------------------\n");
    while(f >> tag) {
        if (tag[0] == '#') {
            char ch;
            printf("    comment  %s ",tag.c_str());
            while (f.get(ch)) {
                printf("%c",ch);
                if (ch == '\n') break;
            }
        }
        else if (tag == "L") {
            f >> L;
            printf("             L = %.1f\n", L);
        }
        else if (tag == "Lsmall") {
            f >> Lsmall;
            printf("        Lsmall = %.1f\n", Lsmall);
        }
        else if (tag == "Llarge") {
            f >> Llarge;
            printf("        Llarge = %.1f\n", Llarge);
        }
        else if (tag == "F") {
            f >> F;
            printf("             F = %.6f\n", F);
        }
        else if (tag == "omega") {
            f >> omega;
            printf("         omega = %.6f\n", omega);
        }
        else if (tag == "ncycle") {
            f >> ncycle;
            printf("         ncycle = %.6f\n", ncycle);
        }
        else if (tag == "natom") {
            f >> natom;
            printf("         natom = %d\n", natom);
            for (int i=0; i<natom; i++) {
                f >> Z[i] >> R[i][0] >> R[i][1] >> R[i][2];
                printf("           atom %2d   %.1f  %10.6f  %10.6f  %10.6f\n", i, Z[i], R[i][0], R[i][1], R[i][2]);
            }
        }
        else if (tag == "k") {
            f >> k;
            printf("             k = %d\n", k);
        }
        else if (tag == "thresh") {
            f >> thresh;
            printf("        thresh = %.1e\n", thresh);
        }
        else if (tag == "safety") {
            f >> safety;
            printf("        safety = %.1e\n", safety);
        }
        else if (tag == "cut") {
            f >> cut;
            printf("           cut = %.2f\n", cut);
        }
        else if (tag == "iState") {
            f >> iState;
            printf("        iState = %s\n", iState.c_str());
        }
        else if (tag == "prefix") {
            f >> prefix;
            printf("        prefix = %s\n", prefix.c_str());
        }
        else if (tag == "ndump") {
            f >> ndump;
            printf("         ndump = %d\n", ndump);
        }
        else if (tag == "nplot") {
            f >> nplot;
            printf("         nplot = %d\n", nplot);
        }
        else if (tag == "nprint") {
            f >> nprint;
            printf("         nprint = %d\n", nprint);
        }
        else if (tag == "nloadbal") {
            f >> nloadbal;
            printf("       nloadbal = %d\n", nloadbal);
        }
        else if (tag == "nio") {
            f >> nio;
            printf("            nio = %d\n", nio);
        }
        else if (tag == "target_time") {
            f >> target_time;
            printf("    target_time = %.3f\n", target_time);
        }
        else if (tag == "tScale") {
            f >> tScale;
            printf("         tScale = %.5f\n", tScale);
        }
        else {
            MADNESS_EXCEPTION("unknown input option", 0);
        }
    }
  }
    
  template <typename Archive>
  void serialize(Archive & ar) {
    ar & L & Lsmall & Llarge & F & omega & ncycle & natom & Z;
    ar & archive::wrap(&(R[0][0]), 3*MAXNATOM);
    ar & k & thresh & safety & cut & iState & prefix & ndump & nplot & nprint & nloadbal & nio;
    ar & target_time & tScale;
  }
};

std::ostream& operator<<(std::ostream& s, const InputParameters& p) {
    s << p.L<< " " << p.Lsmall<< " " << p.Llarge<< " " << p.F << " " << p.omega <<
        " " << p.ncycle << " " << p.Z << " " << p.R[0]<< " " << p.k<< " " <<
        p.thresh<< " " << p.cut<< " " << p.iState << " " << p.prefix<< " " << p.ndump<< " " <<
        p.nplot << " " << p.nprint << " "  << p.nloadbal << " " << p.nio << p.tScale << std::endl;
return s;
}

InputParameters param;


struct WF {
    std::string str;
    complex_functionT func;
    WF(const std::string& STR, const complex_functionT& FUNC) 
        : str(STR)
        , func(FUNC) 
    {
        func.truncate();
    }
};

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
    return fname;
}
bool wave_function_exists(World& world, int step) {
    return ParallelInputArchive::exists(world, wave_function_filename(step));
}
void wave_function_store(World& world, int step, const complex_functionT& psi) {
    ParallelOutputArchive ar(world, wave_function_filename(step), nIOProcessors);
    ar & psi;
}
complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}
// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    // Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;

        // This randomly hashes levels 0-2 and then
        // hashes nodes by their grand-parent key so as
        // to increase locality separately on each level.
        //if (n <= 2) hash = key.hash();
        //else hash = key.parent(2).hash();

        // This randomly hashes levels 0-3 and then 
        // maps nodes on even levels to the same
        // random node as their parent.
        // if (n <= 3 || (n&0x1)) hash = key.hash();
        // else hash = key.parent().hash();

        // This randomly hashes each key
        hash = key.hash();

        return hash%nproc;
    }
};

template<class T>
std::string toString( const T& a ) {
    std::ostringstream o;
    o << a[0] << ", " << a[1] << ", " << a[2];
    return o.str();
}
void loadDefaultBasis(World& world, std::vector<WF>& boundList, double Z) {
    PRINT("Loading the default basis");
    const int NBSt = 3;    //Number of Bound States
    const int bSt[][3] = { {1,0,0},
                           {2,0,0},
                           {2,1,0} };
    for( int i=0; i<NBSt; i++ ) {
       boundList.push_back( WF(toString(bSt[i]), 
                 FunctionFactory<complexd,NDIM>(world).functor(functorT(
                 new BoundWF(Z , bSt[i][0], bSt[i][1], bSt[i][2]) ))));
    }
    PRINTLINE("Done loading the standard basis");
}

void loadList(World& world, std::vector<std::string>& boundList, std::vector<std::string>& unboundList) {
    std::ifstream bound("bound.num");
    std::ifstream unbound("unbound.num");
    if( ! bound.is_open() && ! unbound.is_open() ) {
        PRINTLINE("bound.num and unbound.num not found");
        boundList.push_back("1 0 0");
        boundList.push_back("2 1 0");
    } else {
        if(bound.is_open()) {
            int n,l,m;
            std::ostringstream nlm;
            while(bound >> n) {
                bound >> l;
                bound >> m;
                nlm << n << " " << l << " " << m;
                boundList.push_back(nlm.str());
                nlm.str("");
            }
            bound.close();
        } else PRINTLINE("bound.num not found");
        if(unbound.is_open()) {
            double kx, ky, kz;
            std::ostringstream kxyz;
            while(unbound >> kx) {
                unbound >> ky;
                unbound >> kz;
                kxyz.precision( 8 );
                kxyz << std::fixed;
                kxyz << kx << " " << ky << " " << kz;
                unboundList.push_back(kxyz.str());
                kxyz.str("");
            }
            unbound.close();
        } else PRINTLINE("unbound.num not found");
    }
}

// void loadBasis(World& world, std::vector<WF>& boundList, std::vector<WF>& unboundList, const double Z, const double cutoff) {
//     std::ifstream bound("bound.num");
//     std::ifstream unbound("unbound.num");
//     if( ! bound.is_open() && ! unbound.is_open() ) {
//         PRINTLINE("bound.num and unbound.num not found");
//         loadDefaultBasis(world,boundList,Z);
//     } else {
//         if(bound.is_open()) {
//             PRINTLINE("Calculating bound quantum states");
//             int n,l,m;
//             std::ostringstream nlm;
//             while(bound >> n) {
//                 bound >> l;
//                 bound >> m;
//                 nlm << n << l << m;
//                 double start = wall_time();
//                 PRINT(nlm.str());                
//                 boundList.push_back(WF(nlm.str()+"           ",
//                                        FunctionFactory<complexd,NDIM>(world).
//                                        functor(functorT(new BoundWF(Z,n,l,m)))));
//                 double used = wall_time() - start;
//                 PRINTLINE("\t" << used << " sec" );
//                 nlm.str("");
//             }
//             bound.close();
//         } else PRINTLINE("bound.num not found");
//         if(unbound.is_open()) {
//             PRINTLINE("Calculating unbound quantum states");
//             double kx, ky, kz;
//             std::ostringstream kxyz;
//             while(unbound >> kx) {
//                 unbound >> ky;
//                 unbound >> kz;
//                 kxyz.precision( 2 );
//                 kxyz << std::fixed;
//                 kxyz << kx << " " << ky << " " << kz;
//                 PRINT(kxyz.str());
//                 double start = wall_time();
//                 const double kvec[] = {kx, ky, kz};
//                 unboundList.push_back(WF(kxyz.str(),
//                                        FunctionFactory<complexd,NDIM>(world).
//                                        functor(functorT(new PhiK(Z,kvec,cutoff)))));
//                 double used = wall_time() - start;
//                 PRINTLINE("\t" << used << " sec");
//                 kxyz.str("");
//             }
//             unbound.close();
//         } else PRINTLINE("unbound.num not found");
//     }
// }
complexd zdipole( const vector3D& r) {
    return complexd(r[2],0.0);
}

/*****************************************************************
 * Dipole matrix elements available for perturbation calculations
 * |<phi_A|z|phi_B>|^2
 *****************************************************************/
void projectZdip(World& world, std::vector<WF> stateList) {
    std::cout.precision(8);
    complex_functionT z = complex_factoryT(world).f(zdipole);
    complexd output;
    std::vector<WF>::iterator basisI;
    std::vector<WF>::iterator basisII;
    PRINT(std::endl << "\t\t|<basis_m|z|basis_n>|^2 " << std::endl << "\t\t");
    for(basisII=stateList.begin(); basisII != stateList.end(); basisII++) {
        PRINT("|" << basisII->str << ">");
    }
    PRINT("\n");
    for(basisI = stateList.begin(); basisI !=stateList.end(); basisI++ ) {
        PRINT("<" << basisI->str << "|" );
        for(basisII = stateList.begin(); basisII != stateList.end(); basisII++) {
            output = inner(basisII->func,z*basisI->func); 
            PRINT(" " << real(output*conj(output)) <<"\t");
        }
        PRINT("\n");
    }
    PRINT("\n");
}

void csToFile(World& world, std::vector<WF> basisList, WF psi_t, std::string suffix) {
    std::vector<WF>::iterator basisI;
    complexd output;
    ofstream file;
    file.open((psi_t.str+"."+suffix).c_str());
    if( file.is_open() ) {
        for(basisI = basisList.begin(); basisI !=basisList.end(); basisI++ ) {
            output = psi_t.func.inner( basisI->func ); 
            file << basisI->str << "\t" << real(output*conj(output)) << "\n";
        }
        file.close();
    }
}

void displayToScreen(World& world, std::vector<WF> basisList, std::vector<WF> psiList, std::string header) {
    PRINTLINE(header);
    complexd output;
    std::vector<WF>::iterator basisI;
    std::vector<WF>::iterator psiPlusI;
    //Display table header
    PRINT("\t\t");
    for(psiPlusI=psiList.begin(); psiPlusI != psiList.end(); psiPlusI++) {
        PRINT("|" << psiPlusI->str << ">\t\t");
    }
    //Display table row
    PRINT("\n");
    for(basisI = basisList.begin(); basisI !=basisList.end(); basisI++ ) {
        PRINT("<" << basisI->str << "|" );
        for(psiPlusI = psiList.begin(); psiPlusI != psiList.end(); psiPlusI++) {
            output = psiPlusI->func.inner(basisI->func); 
            PRINT(" " << real(output*conj(output)) << "\t");
        }
        PRINT("\n");
    }
}

// ///What is the overlap of the Coulombic eigenfunctions on the ground state?
// void groundOverlap(World& world, std::vector<std::string> boundList, std::vector<std::string> unboundList, double Z, double cutoff) {
//     //LOAD Psi(0)
//     if( wave_function_exists(world,'0') ) {
//         complex_functionT psi0 = wave_function_load(world,0);
//         PRINTLINE("\t\t|<basis|Psi(0)>|^2 ");
//         if(boundList.empty() && unboundList.empty()) {
//             boundList.push_back("1 0 0");
//             boundList.push_back("2 1 0");
//         }
//         std::string tag;
//         std::vector<WF> psiList;
//         complexd output;
//         PRINTLINE("|t=0>");
//         //phi_nlm
//         if( !boundList.empty() ) {
//             std::vector<std::string>::const_iterator boundIT;
//             int N, L, M;
//             for(boundIT = boundList.begin(); boundIT !=boundList.end(); boundIT++ ) {
//                 std::stringstream ss(*boundIT);
//                 ss >> N >> L >> M;
//                 //PROJECT Phi_nlm into MADNESS
//                 complex_functionT phi_nlm = complex_factoryT(world).
//                     functor(functorT( new BoundWF(Z, N, L, M)));
//                 PRINT(*boundIT << "   ");
//                 // <phi_bound|Psi(0)>
//                 output = phi_nlm.inner( psi0 ); 
//                 PRINTLINE(std::scientific <<"\t" << real(output*conj(output)));
//             }
//         }
//         if( !unboundList.empty() ) {
//             std::vector<std::string>::const_iterator unboundIT;
//             for( unboundIT=unboundList.begin(); unboundIT !=  unboundList.end(); unboundIT++ ) {
//                 //parsing unboundList
//                 double KX, KY, KZ;
//                 std::stringstream ss(*unboundIT);
//                 ss >> KX >> KY >> KZ;
//                 double dArr[3] = {KX, KY, KZ};
//                 const vector3D kVec(dArr);
//                 //screening out the zero vector
//                 if((dArr[1]>0.0 || dArr[1]<0.0) || (dArr[2]>0.0 || dArr[2]<0.0)) {
//                     //PROJECT Psi_k into MADNESS
//                     complex_functionT phi_k = 
//                         complex_factoryT(world).functor(functorT( new PhiK(world, Z, kVec, cutoff) ));
//                     phi_k.init();
//                     output =  inner(psi0,phi_k);
//                     std::cout.precision( 8 );
//                     PRINT( std::fixed << KX << " " << KY << " " << KZ << "\t" <<
//                            std::scientific << "\t" << real(output*conj(output)) );
//                 }
//             }
//         }else {
//             PRINTLINE("psi0 must be present in this directory i.e.  data-00000.0000*");
//         }
//     }
// } 

/************************************************************************************
 * The correlation amplitude |<Psi(+)|basis>|^2 are dependent on the following files:
 * wf.num                  Integer time step of the Psi(+) to be loaded
 * bound.num               Integer triplets of quantum numbers   2  1  0 
 * unbound.num             Double triplets of momentum kx ky kz  0  0  0.5
 ************************************************************************************/

void projectPsi(World& world, std::vector<std::string> boundList, std::vector<std::string> unboundList, const double Z, double cutoff) {
    PRINTLINE("\t\t|<basis|Psi(t)>|^2 ");
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
        PRINT("\t\t");
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
            for(boundIT = boundList.begin(); boundIT !=boundList.end(); boundIT++ ) {
                std::stringstream ss(*boundIT);
                ss >> N >> L >> M;
                //PROJECT Phi_nlm into MADNESS
                complex_functionT phi_nlm = complex_factoryT(world).
                    functor(functorT( new BoundWF(Z, N, L, M)));
                complexd n_overlap_0 = inner(phi_nlm,psi0);
                PRINT(*boundIT << "   ");
                //loop through time steps
                for( psiIT=psiList.begin(); psiIT !=  psiList.end(); psiIT++ ) {
                    //|PSI(t)> = |Psi(t)> - <phi_k|Psi(0)>|Psi(0)>
                    //<phi_k|PSI(t)> = <phi_k|Psi(t)> - <phi_k||Psi(0)> <Psi(0)|Psi(t)>
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
                if((dArr[1]>0.0 || dArr[1]<0.0) || (dArr[2]>0.0 || dArr[2]<0.0)) { //removing k={0,0,0}
                    //PROJECT Psi_k into MADNESS
                    if(world.rank()==0) before = clock();
                    const double constcutoff = cutoff;
                    PhiK phik = PhiK(world, Z, kVec, constcutoff);
                    phik.init(world);
                    complex_functionT phi_k = complex_factoryT(world).functor(functorT( &phik ));
                    // W/O timing
                    //complex_functionT phi_k = 
                    //complex_factoryT(world).functor(functorT( new PhiK(Z, kVec, cutoff) ));
                    if(world.rank()==0) after = clock();
                    std::cout.precision( 8 );
                    PRINT( std::fixed << KX << " " << KY << " " << KZ << "  ");
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

void compare1F1(World& world, double cutoff) {
    //load param
    std::string tag;
    double rMIN = 0.0;
    double rMAX = 10.0;
    double dr   = 1.0;
    double k    = 1.0;
    double Z    = 1.0;
    /***************************************
     *Load graphing parameters from the file: param
     * rMIN 0.0
     * rMAX 10.0
     * dr   1.0
     * TH   0.0
     ****************************************/
    std::ifstream f("param");
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
            else if (tag == "rMIN") {
                f >> rMIN;
                PRINTLINE("rMIN = " << rMIN);
            }
            else if (tag == "rMAX") {
                f >> rMAX;
                PRINTLINE("rMAX = " << rMAX);
            }
            else if (tag == "dr") {
                f >> dr;
                PRINTLINE("dr = " << dr);
            }
            else if (tag == "k") {
                f >> k;
                PRINTLINE("k = " << k);
            }
        }
    }
    //make functor
    const double kvec[3] = {0, 0, k};
    PhiK phi_k =  PhiK(Z, kvec, cutoff);
    complexd ONE(1.0,0.0);
    complexd I(0.0,1.0);
    std::cout << std::fixed;
    for(double r=rMIN; r<rMAX; r+=dr) {
        complexd ZZ(0.0,-r);
        std::cout.precision(2);
        PRINT(r                         << "\t");
        std::cout.precision(8);
        PRINT(real(conhyp(-I/k,ONE,ZZ)) << "\t");
        PRINT(imag(conhyp(-I/k,ONE,ZZ)) << "\t");
        PRINT(real(phi_k.aForm(ZZ))     << "\t");
        PRINT(imag(phi_k.aForm(ZZ))     << "\n");
    }
}



// Invoke as \c u(r/c)/c where \c c is the radius of the smoothed volume.  
static double smoothed_potential(double r) {
    double r2 = r*r;
    double pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2) {
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }
    return pot;
}

// Nuclear attraction potential
complexd V(const vector3D& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      double rr = sqrt(xx*xx+yy*yy+zz*zz);
      sum +=  -param.Z[i]*smoothed_potential(rr/param.cut)/param.cut;
    }
    return complexd(sum,0.0);
}
double myreal(double t) {return t;}
double myreal(const double_complex& t) {return real(t);}
// Given psi and V evaluate the energy ... leaves psi compressed, potn reconstructed
template <typename T>
double energy(World& world, const Function<T,3>& psi, const Function<T,3>& potn) {
    // First do all work in the scaling function basis
    psi.reconstruct();
    bool DOFENCE = false;
    Derivative<T,3> Dx(world,0), Dy(world,1), Dz(world,2);
    Function<T,3> dx = Dx(psi,DOFENCE);
    Function<T,3> dy = Dy(psi,DOFENCE);
    Function<T,3> dz = Dz(psi,DOFENCE);
    Function<T,3> Vpsi = psi*potn;
    // Now do all work in the wavelet basis
    psi.compress(DOFENCE); Vpsi.compress(DOFENCE); dx.compress(DOFENCE); dy.compress(DOFENCE); dz.compress(true);
    T S = psi.inner(psi);
    T PE = psi.inner(Vpsi);
    T KE = 0.5*(inner(dx,dx) + inner(dy,dy) + inner(dz,dz));
    T E = (KE+PE)/S;
    dx.clear(); dy.clear(); dz.clear(); Vpsi.clear(); // To free memory on return
    world.gop.fence();
    return myreal(E);
}
void converge(World& world, complex_functionT& potn, complex_functionT& psi, double& eps) {
     for (int iter=0; iter<30; iter++) {
         SeparatedConvolution<double,NDIM> op =
             BSHOperator3D(world, sqrt(-2*eps), param.cut, param.thresh);
         complex_functionT Vpsi = (potn*psi);
         Vpsi.scale(-2.0).truncate();
         complex_functionT tmp = apply(op,Vpsi).truncate(param.thresh);
         double norm = tmp.norm2();
         complex_functionT r = tmp-psi;
         double rnorm = r.norm2();
         double eps_new = eps - 0.5*real(inner(Vpsi,r))/(norm*norm);
         // if (world.rank() == 0) {
         //     print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
         // }
         psi = tmp.scale(1.0/norm);
         eps = eps_new;
         if (rnorm < std::max(1e-5,param.thresh)) break;
     }
     psi.truncate(param.thresh);     
 }


void compareGroundState(World& world, double Z) {
    //Create the softened ground state for a list of cuts
    const int nCUT = 5;
    double cut[nCUT] = { 0.3, 0.2, 0.1, 0.05, 0.03 };
    complexd output;
    if (world.rank() == 0) param.read("input");
    PRINTLINE("            |<1s|psi0>|^2 \t |<2s|psi0>|^2 \t |<3s|psi0>|^2 \t |<4s|psi0>|^2");
    //make |ns> states
    complex_functionT oneS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 1, 0, 0)));
    complex_functionT twoS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 2, 0, 0)));
    complex_functionT threeS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 3, 0, 0)));
    complex_functionT fourS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 4, 0, 0)));
    for(int i=0; i<nCUT; i++) {
        //Generate psi0 cut
        complex_functionT psi0 = FunctionFactory<complexd,NDIM>(world).
            functor(functorT(new BoundWF(Z, 1, 0, 0)));
        psi0.scale(1/psi0.norm2());
        psi0.truncate();//Does this throw away small coefficents?
        psi0.scale(1/psi0.norm2());
        //Make potn
        param.cut = cut[i]; //The function V calls param.cut
        complex_functionT potn = complex_factoryT(world).f(V);
        potn.truncate(param.thresh);
        double eps = energy(world, psi0, potn);
        converge(world, potn, psi0, eps);
        output = inner(oneS,psi0);
        PRINT(std::setprecision(2) << "cut = " << cut[i] << "  " <<std::setprecision(12) 
              << real(conj(output)*output) << "\t");
        output = inner(twoS,psi0);
        PRINT(real(conj(output)*output) << "\t");
        output = inner(threeS,psi0);
        PRINT(real(conj(output)*output) << "\t");
        output = inner(fourS,psi0);
        PRINTLINE(real(conj(output)*output) << "\t");
    }
}








/*************************************************
 * If you're curious about a wave function's value
 *************************************************/
void printBasis(World& world, double Z, double cutoff = 10.0) {
    complexd output, output2;
    double sinTH, cosTH, sinPHI, cosPHI;
    std::string tag;
    double rMIN = 0.0;
    double rMAX = 10.0;
    double dr = 1.0;
    double TH = 0.0;
    double PHI = 0.0;
    double k = 1.0;
    /***************************************
     *Load graphing parameters from the file: param
     * rMIN 0.0
     * rMAX 10.0
     * dr   1.0
     * TH   0.0
     ****************************************/
    std::ifstream f("param");
    std::cout << std::fixed;
    std::cout.precision(1);
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
            else if (tag == "rMIN") {
                f >> rMIN;
                PRINTLINE("rMIN   = " << rMIN);
            }
            else if (tag == "rMAX") {
                f >> rMAX;
                PRINTLINE("rMAX   = " << rMAX);
            }
            else if (tag == "dr") {
                f >> dr;
                PRINTLINE("dr     = " << dr);
            }
            else if (tag == "TH") {
                f >> TH;
                PRINTLINE("TH     = " << TH);
            }
            else if (tag == "k") {
                f >> k;
                PRINTLINE("k      = " << k);
            }
            else if (tag == "cutoff") {
                f >> cutoff;
                PRINTLINE("cutoff = " << cutoff);
            }
        }
    }
    //make functions
    double kvec[3] = {0, 0, k};
    //    vector3D kVec(dARR);
    PhiK phi_k(Z, kvec, cutoff);
    //for(double TH=0; TH<3.14; TH+=0.3 ) {
    //    for(double r=0; r<sqrt(3)*phi_k.domain*phi_k.k; r+=1.0 ) {
    PRINT("r \tRe:phi(rVec) \t Im:phi(rVec) \t");
    PRINTLINE("Re:diff(rVec) \t Im:diff(rVec) \t");
    for(double r=rMIN; r<rMAX; r+=dr ) {
        std::cout.precision(3);
        std::cout << std::fixed;
        cosTH =  std::cos(TH);
        sinTH =  std::sin(TH);
        cosPHI = std::cos(PHI);
        sinPHI = std::sin(PHI);
        double rvec[3] = {r*sinTH*cosPHI, r*sinTH*sinPHI, r*cosTH};
        output = phi_k(rvec);
        PRINT(r);
        std::cout.precision(7);
        PRINT("\t" << real(output) << "\t" << imag(output));
//         std::cout << std::scientific;
//         PRINT("\t" << phi_k.diffR(2*k*r) << "\t" << phi_k.diffI(2*k*r));
//      PRINT("\t" << real(phi_k.f11(2*k*r)) << "\t" << imag(phi_k.f11(2*k*r)));
//      PRINT("\t" << real(conhyp(-I/k,ONE,-I*k*r)) << "\t" << imag(conhyp(-I/k,ONE,-I*k*r)));
        PRINTLINE(" ");
    }
    //    use sed to make the complexd output standard
    //    system("sed -i '' -e's/\\+/, /' -e's/j//' f11.out");
}

/****************************************************************************
 * Reproduces atomic form integrals given by Dz Belkic's analytic expression
 ****************************************************************************/
void belkic(World& world, double cutoff) {
    /************************************************************************
     * qVec is the momentum transfered from the laser field
     * kVec is the momentum of the ejected electron
     * I will take these to be colinear
     ************************************************************************/
    PRINTLINE("1 0 0");
    double Z = 1.0;
    complex_functionT b1s = complex_factoryT(world).functor(functorT( 
                                                      new BoundWF(Z, 1, 0, 0) ));
    double dARR[3] = {0, 0, 0.5};
    const vector3D kVec(dARR);
    PRINTLINE("|" << kVec << ">");
    PhiK phik = PhiK(Z, kVec, cutoff);
    phik.init(world);
    complex_functionT phi_k = complex_factoryT(world).functor(functorT( &phik ));
    dARR[2] =  1.5; // {0, 0, 1.5}
    const vector3D qVec(dARR);
    PRINTLINE("Exp[I" << qVec << ".r>");
    complex_functionT expikDOTr = complex_factoryT(world).functor(functorT(
                                                           new Expikr(qVec) ));
    PRINTLINE("<k=0.5| Exp[iqVec.r] |100>");
    complexd output = inner(phi_k, expikDOTr*b1s);
    PRINTLINE(output);
}

void loadParameters(World& world, int& k, double& L, double &Z, double &cutoff) {
    std::string tag;
    int natom;
    double Rx, Ry, Rz;
    std::ifstream f("input");
    std::cout << std::fixed;
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
                //E_n = n hbar omega
                //I_p = 0.5 Z^2
                //KE = E_n - I_p
                //Atomic Units: hbar = m = 1
                //v = sqrt( 2n omega - Z^2)
                //cutoff > dMAX = v t
                double dMAX = std::sqrt(2*3*omega - Z*Z) * 10;
                cutoff = 0.0;
                while( cutoff < dMAX ) { cutoff += L/16; }
                PRINTLINE("dMAX = " << dMAX);
                PRINTLINE("cutoff = " << cutoff);
            }
        }
    }
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
    int    k = 12;
    double L = 10.0;
    double Z = 1.0;
    double cutoff = L;
    loadParameters(world, k, L, Z, cutoff);
    FunctionDefaults<NDIM>::set_k(k);               // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(1e-6);       // Accuracy
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
