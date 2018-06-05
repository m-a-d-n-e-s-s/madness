#ifndef INPUT_H
#define INPUT_H

#define NO_GENTENSOR
#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/constants.h>
#include <madness/mra/lbdeux.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace madness;

// Definitions of functions, function-vectors, operators
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > real_functorT;
typedef std::shared_ptr< FunctionFunctorInterface<double_complex,3> > comp_functorT;
typedef Function<double,3> real_functionT;
typedef Function<double_complex,3> comp_functionT;
typedef std::vector<real_functionT> real_vecfuncT;
typedef std::vector<comp_functionT> comp_vecfuncT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;
typedef Tensor<double> real_tensorT;
typedef Tensor<double_complex> comp_tensorT;
typedef FunctionFactory<double,3> real_factoryT;
typedef FunctionFactory<double_complex,3> comp_factoryT;
typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;


// Basic simulation setup parameters 
struct Setup
{
    Setup() : A(208), 
              Z(82), 
              length(200.0), 
              initial(3), 
              boundary(1), 
              knumber(5.0),
              project(1)
              {}
    Setup(const std::string file): A(208), 
                                   Z(82), 
                                   length(200), 
                                   initial(3), 
                                   boundary(1),
                                   knumber(5.0),
                                   thresh(1.e-3),
                                   IO_nodes(20),
                                   project(1) {
        std::ifstream f(file.c_str());
        position_stream(f, "setup:");
        std::string s;
        while (f >> s) {
            if (s=="end_setup") break;
            else if (s=="A") {f>>A;}
            else if (s=="Z") {f>>Z;}
            else if (s=="length") {f>>length;}
            else if (s=="initial") {f>>initial;}
            else if (s=="boundary") {f>>boundary;}
            else if (s=="knumber") {f>>knumber;}
            else if (s=="thresh") {f>>thresh;}
            else if (s=="IO_nodes") {f>>IO_nodes;}
            else if (s=="project") {f>>project;}
            else {std::cout << "mshf.input: unrecognized input for setup " << s << std::endl;}
        }
    }       

    int A;
    int Z;
    double length;
    int    initial;
    int    boundary;
    double knumber;
    double thresh;
    int    IO_nodes;
    int    project;
};


// Nuclear matter parameters
struct Nuclear 
{
    Nuclear() : spinorbit(1),
                meff(1),
                jellium(0),
                screening(0),
                screenl(3.0),
                lap_comp(2)
                {}
    Nuclear(const std::string file): spinorbit(1),
                                     meff(1),
                                     jellium(0),
                                     screening(0),
                                     screenl(3.0),
                                     lap_comp(2) {
        std::ifstream f(file.c_str());

        position_stream(f, "nuclear:");
        std::string n;
        while (f >> n) {
            if (n=="end_nuclear") break;
            else if (n=="spinorbit") {f>>spinorbit;}
            else if (n=="meff") {f>>meff;}
            else if (n=="jellium") {f>>jellium;}
            else if (n=="screening") {f>>screening;}
            else if (n=="screenl") {f>>screenl;}
            else if (n=="lap_comp") {f>>lap_comp;}
            else {std::cout << "mshf.input: Unrecognized input for nuclear" << n << std::endl;}
        }
    }       

    int    spinorbit;
    int    meff;
    int    jellium;
    int    screening;
    double screenl;
    int    lap_comp;
};


// Mixing parameters
struct Mixing
{
    Mixing() : avg_pot(0),
               avg_lap(1),
               avg_wav(1)
               {}
    Mixing(const std::string file): avg_pot(0),
                                    avg_lap(1),
                                    avg_wav(1) {
        std::ifstream f(file.c_str());

        position_stream(f, "mixing:");
        std::string m;
        while (f >> m) {
            if (m == "end_mixing") break;
            else if (m == "avg_pot") {f>>avg_pot;}
            else if (m == "avg_lap") {f>>avg_lap;}
            else if (m == "avg_wav") {f>>avg_wav;}
            else {std::cout << "mshf.input: Unrecognized input for mixing" << m << std::endl;}
        }
    }       

    int avg_pot;
    int avg_lap;
    int avg_wav;
};


// Output parameters
struct Output
{
    Output() : vtk_output(1),
               txt_output(1)
               {}
    Output(const std::string file): vtk_output(1),
                                    txt_output(1) {
        std::ifstream f(file.c_str());

        position_stream(f, "output:");
        std::string o;
        while (f >> o) {
            if (o == "end_output") break;
            else if (o == "vtk_output") {f>>vtk_output;}
            else if (o == "txt_output") {f>>txt_output;}
            else {std::cout << "mshf.input: Unrecognized input for output" << o << std::endl;}
        }
    }       

    int vtk_output;
    int txt_output;
};


// Additional parameters
struct Additional
{
    Additional() : use_infile(0),
                   prec(0.1),
                   chi(0.6),
                   tol(1.e-5),
                   brad(0.25),
                   timing(1)
                   {}
    Additional(const std::string file): use_infile(0),
                                        prec(0.1),
                                        chi(0.6),
                                        tol(1.e-5),
                                        brad(0.25),
                                        timing(1) {
        std::ifstream f(file.c_str());

        position_stream(f, "additional:");
        std::string a;
        while (f >> a) {
            if (a == "end_additional") break;
            else if (a == "use_infile") {f>>use_infile;}
            else if (a == "prec") {f>>prec;}
            else if (a == "chi") {f>>chi;}
            else if (a == "tol") {f>>tol;}
            else if (a == "brad") {f>>brad;}
            else if (a == "timing") {f>>timing;}
            else {std::cout << "mshf.input: Unrecognized input for additional" << a << std::endl;}
        }      
    }       

    int    use_infile;
    double prec;
    double chi;
    double tol;
    double brad;
    int    timing;
};


// Parameters for the used Skyrme force 
struct Skyrme
{
    Skyrme() : t0(-1879.64001), 
               t1(313.74933), 
               t2(112.67627), 
               t3(12527.38965), 
               t4(124.63330), 
               x0(0.25855), 
               x1(-0.38169), 
               x2(-2.82364), 
               x3(0.12323), 
               bp4(34.11170),
               alpha(0.30), 
               k_fn((20.72126e0 + 20.74982e0)/2.0), 
               k_fp((20.72126e0 + 20.74982e0)/2.0)         
    {}
    Skyrme(const std::string file): t0(-1879.64001), 
                                    t1(313.74933), 
                                    t2(112.67627), 
                                    t3(12527.38965), 
                                    t4(124.63330),
                                    x0(0.25855), 
                                    x1(-0.38169), 
                                    x2(-2.82364), 
                                    x3(0.12323), 
                                    bp4(34.11170),
                                    alpha(0.30), 
                                    k_fn((20.72126e0 + 20.74982e0)/2.0), 
                                    k_fp((20.72126e0 + 20.74982e0)/2.0) {
        std::ifstream f(file.c_str());
        position_stream(f, "skyrme");
        std::string s;
        while (f >> s) {
            if(s == "end") break;
            else if(s == "t0")f>>t0;
            else if(s == "t1")f>>t1;
            else if(s == "t2")f>>t2;
            else if(s == "t3")f>>t3;
            else if(s == "t4")f>>t4;
            else if(s == "x0")f>>x0;
            else if(s == "x1")f>>x1;
            else if(s == "x2")f>>x2;
            else if(s == "x3")f>>x3;
            else if(s == "bp4")f>>bp4;
            else if(s == "alpha")f>>alpha;
            else if(s == "k_fn")f>>k_fn;
            else if(s == "k_fp")f>>k_fp;
            else {std::cout << "mshf.input: Unrecognized input for skyrme" << s << std::endl;}
        }
    }
    double t0;
    double t1;
    double t2;
    double t3;
    double t4;
    double x0;
    double x1;
    double x2;
    double x3;
    double bp4;
    double alpha;
    double k_fn;
    double k_fp;
};




// Physical parameters
//------------------
const double e2        = 1.43989e0;                         // square of elementary charge [MeV*fm]
const double_complex I = double_complex(0.0, 1.0);
static const double PI = 3.141592653589793238;

#endif
