#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include "madness/mra/commandlineparser.h"
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/molecule.h>

using namespace madness;

int main(int argc, char** argv) {

    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);


    commandlineparser parser(argc,argv);
    bool silent=false;
    if (parser.key_exists("silent")){
    	silent = (bool) std::stoi(parser.value("silent"));
    }

    if (world.rank()==0){
    	std::cout << "--------------------------\n";
    	std::cout << "Madness Orbital to Plane\n\n";
    	std::cout << "usage:\n";
    	std::cout << "plot2cube file=name\n\n";
    	std::cout << "with e.g. name=my_mra_function.00000\n";
    	std::cout << "\n\narguments:\n";
    	std::cout << std::setw(25) << std::left << "\tno_orient=<bool> " <<  "orient the molecule according to madness defaults or not - default=true\n";
    	std::cout << std::setw(25) << std::left << "\tinput<string> " <<  "madness input file - default='input'\n";
    	std::cout << std::setw(25) << std::left << "\toutfile=<string> " << "name of outputfile\n";
    	std::cout << std::setw(25) << std::left << "\tsilent<bool> " << "silence output - default=false\n";
    	std::cout << "--------------------------\n";
    	std::cout << "necessary files:\n";
    	std::cout << std::setw(25) << std::left << "\tinputfile" << "the inputfile for madnesss that should hold plot parameters.\nE\t\texample:\n";
    	std::cout << "\n";
    	std::cout << "\t\tplot\n";
    	std::cout << "\t\tplane=x1,x2\n";
    	std::cout << "\t\tpoints=100\n";
    	std::cout << "\t\tzoom=1.5\n";
    	std::cout << "\t\toutput=gnuplot\n";
    	std::cout << "\t\torigin=0.0 0.0 0.0\n";
    	std::cout << "\t\tend\n\n";
    	std::cout << std::setw(25) << std::left << "\tname.00000" << "the madness function file - filename given by plot2cube file=name\n";
    	std::cout << "--------------------------\n";
    }

    // get names of input files
    // input: the name of the file with the madness input (mol geometry is needed)
    // file:  the name of the file with the MRA data of the madness orbital
    std::string inputfile = "input";
    std::string filename = "";
    std::string outfile="";
    double L=-1.0;
    int npoints=200;
    double zoom=2.0;
    auto no_orient=false;
    if (parser.key_exists("input")){
    	inputfile = parser.value("input");
    }
    if (parser.key_exists("file")){
    	filename = parser.value("file");
    }else{
    	throw std::invalid_argument("file=name not provided");
    }
    if (parser.key_exists("outfile")){
    	outfile = parser.value("outfile");
    }else outfile = filename;

    outfile += ".plane";


    // check if files exist
    CalculationParameters dummy;
    auto file_ok=dummy.file_exists(world,parser.value("input"));
    if (not file_ok){
        std::string msg="could not find user-defined input file: "+parser.value("input")+"\n";
        throw std::invalid_argument(msg);
    }
    file_ok=dummy.file_exists(world,parser.value("file")+".00000");
    if (not file_ok){
        std::string msg="could not find data file: "+parser.value("file")+"\n";
        throw std::invalid_argument(msg);
    }

    // load functions
    Function<double,3> f = FunctionFactory<double,3>(world);
    if(world.rank()==0) std::cout << "load function " << filename << "\n";
    load(f,filename);
    if(world.rank()==0) std::cout << "... success\n";

    // fetch defaults from f;

    FunctionDefaults<3>::set_k(f.k());
    FunctionDefaults<3>::set_thresh(f.thresh());

    // read in molecule
    if(world.rank()==0) std::cout << "initializing molecule from input=" << inputfile << "\n";
    Molecule molecule;
    molecule.read_file(inputfile);

    if(no_orient==false){
    	if(world.rank()) std::cout << "re-orienting\n";
    	molecule.orient();
    }
    if(world.rank()==0) std::cout << "... success\n";

    // if L could not determined until this point
    // guess it with the madness defaults
    if(L<0.0){
    	dummy = CalculationParameters(world, parser);
    	AtomicBasisSet aobasis("sto-3g");
    	dummy.set_derived_values(molecule, aobasis, parser);
    	L = molecule.bounding_cube();
    }
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    if(world.rank()==0){
    	std::cout << "using madness parameters:\n";
    	std::cout << "\tthresh=" << f.thresh() << "\n";
    	std::cout << "\tk=" << f.k() << "\n";
    	std::cout << "\tL=" << L << "\n";
    }

	// plot the Gaussian cube file

    std::vector<Function<double,3>> vf;
    vf.push_back(f);
    if(world.rank()==0) std::cout << "creating plot data " << outfile << "\n";
    plot_plane(world,vf,outfile,inputfile);
    if(world.rank()==0) std::cout << "... success\n";

    return 0;

}
