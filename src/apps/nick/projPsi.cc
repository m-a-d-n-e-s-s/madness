//\file projPsi.cc
//\brief Projects a time evolved wave function onto an arbitrary number of bound states
/***************************************************************************************
 * By: Nick Vence
 * This code must handled with care for the following reasons:
 * 1) It uses the Gnu Scientific Library      http://www.gnu.org/software/gsl/
 * MADNESS should be reconfigured to include these libraries
 *  ./configure LIBS="-lgsl -lgslblas"
 * 2) A list of wave functions from tdse is expected in this directory
 *    List the number from the above wave functions in wf.num
 * 3) It is reading files from the disk
 *    wf.num                  Integer of the WF to be loaded
 *    bound.num               Integer triplets of quantum numbers   2    1  0 
 *    unbound.num             Double triplets momentum kx ky kz     0.5  0  0
 * 4) k (The wavelet order) must be the same as the projected functions: see main()
 *    12 has been the default
 *****************************************/
#include <mra/mra.h>
#include <complex>
#include <string>
#include <fstream>
using std::ofstream;
#include "wavef.h"
#define PRINT(str) if(world.rank()==0) cout << str 
#define PRINTLINE(str) if(world.rank()==0) cout << str << endl

using namespace madness;

const int nIOProcessors =1;
const string prefix = "data";
typedef std::complex<double> complexd;
typedef Vector<double,NDIM> vector3D;
typedef Function<complexd,NDIM> complex_functionT;
typedef Function<double,NDIM> functionT;
typedef FunctionFactory<complexd,NDIM> complex_factoryT;
typedef FunctionFactory<double,NDIM> factoryT;
const char* wave_function_filename(int step);
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);

struct WF {
    WF(const string& str, const complex_functionT& func) : str(str), func(func) {}
    string str;
    complex_functionT func;
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

template<class T>
string toString( const T& a ) {
    ostringstream o;
    o << a[0] << ", " << a[1] << ", " << a[2];
    return o.str();
}

void loadDefaultBasis(World& world, std::vector<WF>& stateList) {
    /******
     * For now I'm going to use a 3xN array to house my coordinates so they 
     * are easy to see, and a vector because it easy to manipulate generally.
     * NEXT: load from a file
     ******/
    const int NBSt = 3;    //Number of Bound States
    const int bSt[][3] = { {1,0,0},
                           {2,0,0},
                           {2,1,0}};
    const int NkSt = 3;    //Number of k States
    const double kSt[][3]  = { {0.52, 0.0, 0.0}
                             };
    double Z = 1.0;
    for( int i=0; i<NBSt; i++ ) {
       stateList.push_back( WF(toString(bSt[i]), 
                 FunctionFactory<complexd,NDIM>(world).functor(functorT(
                 new BoundWF(Z , bSt[i][0], bSt[i][1], bSt[i][2]) ))));
    }
    for( int i=0; i<NkSt; i++ ) {
        stateList.push_back(WF(toString(kSt[i]), 
                 FunctionFactory<complexd,NDIM>(world).functor(functorT(
                 new ScatteringWF(Z, kSt[i]) ))));
    }
    PRINTLINE("Done loading the standard basis");
}
void loadBasis(World& world, std::vector<WF>& stateList) {
    ifstream bound("bound.num");
    ifstream unbound("unbound.num");
    if( ! bound.is_open() && ! unbound.is_open() ) {
        PRINTLINE("bound.num and unbound.num not found");
        loadDefaultBasis(world,stateList);
    } else {
        double Z = 1.0;
        if(bound.is_open()) {
            PRINTLINE("Loading bound quantum states");
            int n,l,m;
            std::ostringstream nlm;
            while(bound >> n) {
                bound >> l;
                bound >> m;
                nlm << n << l << m;
                PRINTLINE(nlm.str());
                stateList.push_back(WF(nlm.str(),FunctionFactory<complexd,NDIM>(world).
                                       functor(functorT(new BoundWF(Z,n,l,m)))));
                nlm.str("");
            }
            bound.close();
        } else PRINT("bound.num not found");
        if(unbound.is_open()) {
            PRINTLINE("Loading unbound quantum states");
            double kx, ky, kz;
            std::ostringstream kxyz;
            while(unbound >> kx) {
                unbound >> ky;
                unbound >> kz;
                kxyz << kx << " " << ky << " " << kz;
                const double kvec[] = {kx, ky, kz};
                stateList.push_back(WF(kxyz.str(),FunctionFactory<complexd,NDIM>(world).
                                       functor(functorT(new ScatteringWF(Z,kvec)))));
                kxyz.str("");
            }
            unbound.close();
        } else PRINTLINE("unbound.num not found");
    }
}
/*
void projectOp(World& world, vector<WF>& aVec, , vector<WF>& bVec) {
void testBasis(World& world, vector<WF>& basisF) {
    PRINTLINE("Testing basis function orthonormality");
    vector<WF>::iterator psiA, psiB;
    complexd output;
    for(psiA=basisF.begin(); psiA!=basisF.end(); psiA++) {
            PRINT("\t|" << psiA->str << ">\t\t");
    }
    PRINT("\n");
    for(psiA=basisF.begin(); psiA!=basisF.end(); psiA++) {
        PRINT("<" << basisI->str << "|\t");
        for(psiB=basisF.begin(); psiB <= psiA; psiB++) {
            output = psiA->func.inner(psiB.func);
            PRINT("\n");
            void display(
                         */
complexd zdipole( const vector3D& r) {
    return complexd(r[2],0.0);
}

void doWork(World& world) {
    std::vector<WF> stateList;    
    std::vector<WF> psiList;
    loadBasis(world, stateList);
    ifstream f("wf.num");
    if(f.is_open()) {
        string tag;
        while(f >> tag) {
            //LOAD Psi(+)
            if(wave_function_exists(world, atoi(tag.c_str())) ) {
                psiList.push_back(WF(tag, wave_function_load(world, atoi(tag.c_str()))));
            } else {
                PRINT("Function: " << tag << " not found"<< endl);
            }
            f.close();
        }
//         PRINTLINE("testing matrix_inner_op()");
//         complex_functionT z = complex_factoryT(world).f(zdipole);
//         PRINT(matrix_inner_op(world, stateList, z, psiList) << endl);
        PRINTLINE("The Psi(+) are loaded");
        complexd output;
        std::vector<WF>::iterator basisI;
        std::vector<WF>::iterator psiPlusI;
        //DISPLAY <basis| (x) |psi+> 
        for(psiPlusI=psiList.begin(); psiPlusI != psiList.end(); psiPlusI++) {
//        for(psiPlusI=stateList.begin(); psiPlusI != stateList.end(); psiPlusI++) {
            PRINT("\t|" << psiPlusI->str << ">\t\t");
        }
        PRINT("\n");
        for(basisI = stateList.begin(); basisI !=stateList.end(); basisI++ ) {
            PRINT("<" << basisI->str << "|" );
            for(psiPlusI = psiList.begin(); psiPlusI != psiList.end(); psiPlusI++) {
//            for(psiPlusI = stateList.begin(); psiPlusI != stateList.end(); psiPlusI++) {
                output = psiPlusI->func.inner(basisI->func);
                PRINT("\t" << output*conj(output));
            }
            PRINT("\n");
        }
     } else {
        PRINTLINE("File: wf.num expected to contain a " 
               << "list of integers of loadable wave functions");
    }
    PRINTLINE("End of Job");
}


int main(int argc, char**argv) {
    // Initialize the parallel programming environment
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    // Load info for MADNESS numerical routines
    startup(world,argc,argv);
    // Setup defaults for numerical functions
    FunctionDefaults<NDIM>::set_k(12);              // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(1e-5);       // Accuracy
    FunctionDefaults<NDIM>::set_cubic_cell(-1000.0, 1000.0);
    FunctionDefaults<NDIM>::set_initial_level(4);
    FunctionDefaults<NDIM>::set_cubic_cell(-1000,1000);
    FunctionDefaults<NDIM>::set_apply_randomize(false);
    FunctionDefaults<NDIM>::set_autorefine(false);
    FunctionDefaults<NDIM>::set_truncate_mode(1);
    try {
        doWork(world);
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

    MPI::Finalize();				//FLAG
    return 0;
}
