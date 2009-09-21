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
#include "wavef.h"
#include <string>
#include <fstream>
using std::ofstream;
using std::ofstream;
#include <stdlib.h>
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
    string str;
    complex_functionT func;

    WF(const string& STR, const complex_functionT& FUNC) 
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
template<class T>
string toString( const T& a ) {
    ostringstream o;
    o << a[0] << ", " << a[1] << ", " << a[2];
    return o.str();
}
void loadDefaultBasis(World& world, std::vector<WF>& boundList) {
    PRINT("Loading the default basis");
    const int NBSt = 3;    //Number of Bound States
    const int bSt[][3] = { {1,0,0},
                           {2,0,0},
                           {2,1,0} };
    double Z = 1.0;
    for( int i=0; i<NBSt; i++ ) {
       boundList.push_back( WF(toString(bSt[i]), 
                 FunctionFactory<complexd,NDIM>(world).functor(functorT(
                 new BoundWF(Z , bSt[i][0], bSt[i][1], bSt[i][2]) ))));
    }
    PRINTLINE("Done loading the standard basis");
}

void loadBasis(World& world, double Z, std::vector<WF>& boundList, std::vector<WF>& unboundList) {
    ifstream bound("bound.num");
    ifstream unbound("unbound.num");
    if( ! bound.is_open() && ! unbound.is_open() ) {
        PRINTLINE("bound.num and unbound.num not found");
        loadDefaultBasis(world,boundList);//HERE;
    } else {
        double Z = 1.0;
        if(bound.is_open()) {
            PRINTLINE("Calculating bound quantum states");
            int n,l,m;
            std::ostringstream nlm;
            while(bound >> n) {
                bound >> l;
                bound >> m;
                nlm << n << l << m;
                double start = wall_time();
                PRINT(nlm.str());                
                boundList.push_back(WF(nlm.str()+"           ",
                                       FunctionFactory<complexd,NDIM>(world).
                                       functor(functorT(new BoundWF(Z,n,l,m)))));
                double used = wall_time() - start;
                PRINTLINE("\t" << used << " sec" );
                nlm.str("");
            }
            bound.close();
        } else PRINTLINE("bound.num not found");
        if(unbound.is_open()) {
            PRINTLINE("Calculating unbound quantum states");
            double kx, ky, kz;
            std::ostringstream kxyz;
            while(unbound >> kx) {
                unbound >> ky;
                unbound >> kz;
                kxyz.precision( 2 );
                kxyz << fixed;
                kxyz << kx << " " << ky << " " << kz;
                PRINT(kxyz.str());
                double start = wall_time();
                const double kvec[] = {kx, ky, kz};
                unboundList.push_back(WF(kxyz.str(),
                                       FunctionFactory<complexd,NDIM>(world).
                                       functor(functorT(new ScatteringWF(Z,kvec)))));
                double used = wall_time() - start;
                PRINTLINE("\t" << used << " sec");
                kxyz.str("");
            }
            unbound.close();
        } else PRINTLINE("unbound.num not found");
    }
}
complexd zdipole( const vector3D& r) {
    return complexd(r[2],0.0);
}

/*****************************************************************
 * Dipole matrix elements available for perturbation calculations
 * |<phi_A|z|phi_B>|^2
 *****************************************************************/
void projectZdip(World& world, std::vector<WF> stateList) {
    cout.precision(8);
    complex_functionT z = complex_factoryT(world).f(zdipole);
    complexd output;
    std::vector<WF>::iterator basisI;
    std::vector<WF>::iterator basisII;
    PRINT(endl << "\t\t|<basis_m|z|basis_n>|^2 " << endl << "\t\t");
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

void csToFile(World& world, std::vector<WF> basisList, WF psi_t, string suffix) {
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

void displayToScreen(World& world, std::vector<WF> basisList, std::vector<WF> psiList, string header) {
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

/****************************************************************************
 * The correlation amplitude |<Psi(+)|basis>|^2
 * wf.num                  Integer time step of the Psi(+) to be loaded
 * bound.num               Integer triplets of quantum numbers   2  1  0 
 * unbound.num             Double triplets of momentum kx ky kz  0  0  0.5
 ****************************************************************************/
void projectPsi(World& world, std::vector<WF> boundList, std::vector<WF> unboundList) {
    std::vector<WF> psiList;
    if(boundList.empty() && unboundList.empty())
        loadDefaultBasis(world, boundList);
    ifstream f("wf.num");
    if(f.is_open()) {
        string tag;
        //LOAD Psi(+)
        while(f >> tag) {
            if(wave_function_exists(world, atoi(tag.c_str())) ) {
                WF psi_t = WF(tag, wave_function_load(world, atoi(tag.c_str())));
                psiList.push_back(WF(tag, psi_t.func));
                if( !boundList.empty() )
                    csToFile(world, boundList, psi_t, "bnd");
                if( !unboundList.empty() )
                    csToFile(world, unboundList, psi_t, "unb");
            } else {
                PRINT("Function: " << tag << " not found"<< endl);
            }
        }
        //join boundList -> unboundList
        std::vector<WF> basisList;
        std::vector<WF>::iterator basisI;
        for(basisI = boundList.begin(); basisI != boundList.end(); basisI++) {
            basisList.push_back(*basisI);
        } for(basisI = unboundList.begin(); basisI != unboundList.end(); basisI++) {
            basisList.push_back(*basisI);
        }
        displayToScreen(world, basisList, psiList,"\t\t|<Psi(+)|basis>|^2 ");
        f.close();
    } else {
        PRINTLINE("File: wf.num expected to contain a " 
                  << "list of integers of loadable wave functions");
    }
}

void compareGroundState(World& world, double Z) {
    //import psi0
    complex_functionT psi0;
    if(wave_function_exists(world, 0) ) {
        psi0 = wave_function_load(world, 0);
    } else {
        PRINTLINE("Psi( t=0 ) must be present");
    }
    //make 1s
    complex_functionT oneS = FunctionFactory<complexd,NDIM>(world).
                 functor(functorT(new BoundWF(Z, 1, 0, 0)));
    //Read in Psi(+)
    ifstream f("wf.num");
    if(f.is_open()) {
        string tag;
        complexd output;
        //LOAD Psi(+)
        output = oneS.inner(oneS);
        PRINTLINE(      "|<1s|1s>|^2 = " << real(output*conj(output)));
        while(f >> tag) {
            if(wave_function_exists(world, atoi(tag.c_str())) ) {
                complex_functionT psi_t = wave_function_load(world, atoi(tag.c_str()));
                output = psi0.inner(psi_t);
                PRINT(      "|<psi0|" << tag << ">|^2 = " << real(output*conj(output)));
                output = oneS.inner(psi_t);
                PRINTLINE("\t|< 1s |" << tag << ">|^2 = " << real(output*conj(output)));
            } else {
                PRINT("Function: " << tag << " not found"<< endl);
            }
        }
    }
}


/*************************************************
 * If you're curious about a wave function's value
 *************************************************/
void printBasis(World& world, double Z) {
    complexd output, output2;
    double sinTH, cosTH, sinPHI, cosPHI;
    //make functions
    complex_functionT psi0;
    if(wave_function_exists(world, 0) ) {
        psi0 = wave_function_load(world, 0);
    } else {
        PRINTLINE("Psi( t=0 ) must be present");
    }

    double dARR[3] = {0, 0, 0.5};
    vector3D kVec(dARR);
    BoundWF psi_100(Z,1,0,0);
    //    ScatteringWF psi_k(Z, kVec);
    double PHI = 0.0;
    double TH = 0.0;
    //for(double TH=0; TH<3.14; TH+=0.3 ) {
    PRINTLINE("k = {" << kVec );
    cout.precision(2);
    //    for(double r=0; r<sqrt(3)*psi_k.domain*psi_k.k; r+=1.0 ) {
    for(double r=0; r<sqrt(3)*10; r+=1.0 ) {
        cout << scientific;
        cosTH =  std::cos(TH);
        sinTH =  std::sin(TH);
        cosPHI = std::cos(PHI);
        sinPHI = std::sin(PHI);
        double dARR[3] = {r*sinTH*cosPHI, r*sinTH*sinPHI, r*cosTH};        
        //        PRINTLINE(r << "\t" << psi_k.diffR(r) << " + I" << psi_k.diffI(r));
        //output = psi_k(dARR);
        output = psi_100(dARR);
        output2 = psi0(dARR);
        PRINTLINE(r << "\t" << real(output) << "\t" << real(output2) << "\t" << dARR);
    }
    //    use sed to make the complexd output standard
    //    system("sed -i '' -e's/\\+/, /' -e's/j//' f11.out");
}

/****************************************************************************
 * Reproduces atomic form integrals given by Dz Belkic's analytic expression
 ****************************************************************************/
void belkic(World& world) {
    /************************************************************************
     * qVec is the momentum transfered from the laser field
     * kVec is the momentum of the ejected electron
     * I will take these to be colinear
     ************************************************************************/
    PRINTLINE("1 0 0");
    complex_functionT b1s = complex_factoryT(world).functor(functorT( 
                                                      new BoundWF(1.0, 1, 0, 0) ));
    double dARR[3] = {0, 0, 0.5};
    const vector3D kVec(dARR);
    PRINTLINE("|" << kVec << ">");
    complex_functionT psi_k = complex_factoryT(world).functor(functorT(
                                                      new ScatteringWF(1.0, kVec) ));
    dARR[2] =  1.5;
    const vector3D qVec(dARR);
    PRINTLINE("Exp[I" << qVec << ".r>");
    complex_functionT expikDOTr = complex_factoryT(world).functor(functorT(
                                                      new Expikr(qVec) ));
    PRINTLINE("<k=0.5| Exp[iqVec.r] |100>");
    complexd output = inner(psi_k, expikDOTr*b1s);
    PRINTLINE(output);
}

void loadParameters(World& world, int& k, double& L, double &Z) {
    string tag;
    int natom;
    double Rx, Ry, Rz;
    ifstream f("input");
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
        }
    }
    //fflush(stdout);
}

int main(int argc, char**argv) {
    // INITIALIZE the parallel programming environment
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    PRINTLINE("After initialize");
    // Load info for MADNESS numerical routines
    startup(world,argc,argv);
    // Setup defaults for numerical functions
    int    k = 12;
    double L = 10.0;
    double Z = 1.0;
    loadParameters(world, k, L, Z);
    FunctionDefaults<NDIM>::set_k(k);               // Wavelet order
    FunctionDefaults<NDIM>::set_thresh(1e-4);       // Accuracy
    FunctionDefaults<NDIM>::set_cubic_cell(-L, L);
    FunctionDefaults<NDIM>::set_initial_level(3);
    FunctionDefaults<NDIM>::set_apply_randomize(false);
    FunctionDefaults<NDIM>::set_autorefine(false);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_truncate_mode(0);
    FunctionDefaults<NDIM>::set_truncate_on_project(true);
    try {
        std::vector<WF> boundList;
        std::vector<WF> unboundList;
        double dARR[3] = {0, 0, 0.5};
        vector3D rVec(dARR);
        //compareGroundState(world, Z);
        //printBasis(world,Z);
        //loadBasis(world, Z, boundList,unboundList);
        //belkic(world);
        //projectZdip(world, unboundList);
        projectPsi(world, boundList, unboundList);         
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
