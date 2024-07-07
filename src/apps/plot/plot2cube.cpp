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
    	std::cout << "Madness Orbital to Cube\n\n";
    	std::cout << "usage:\n";
    	std::cout << "plot2cube file=name\n\n";
    	std::cout << "with e.g. name=my_mra_function.00000\n";
    	std::cout << "\n\narguments:\n";
    	std::cout << std::setw(25) << std::left << "\tno_orient=<bool> " << "orient the molecule according to madness defaults or not - default=true\n";
    	std::cout << std::setw(25) << std::left << "\tinput<string> " << "madness input file - default='input'\n";
    	std::cout << std::setw(25) << std::left << "\toutfile=<string> " << "name of outputfile\n";
    	std::cout << std::setw(25) << std::left << "\tzoom<double> " << "zoom factor for the cube file - default=1.0\n";
    	std::cout << std::setw(25) << std::left << "\tnpoints<int> " << "number of points for each dimension in the cube file -- default=150\n";
    	std::cout << std::setw(25) << std::left << "\tsilent<bool> " << "silence output - default=false\n";
    	std::cout << "--------------------------\n";
    	std::cout << "necessary files:\n";
    	std::cout << "\tinputfile" << std::setw(25) << "the inputfile for madnesss that should hold the molecular geometry\n";
    	std::cout << "\tname.00000" << std::setw(25) << "the madness function file - filename given by plot2cube file=name\n";
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
    double zoom=1.0;
    auto no_orient=false;
    if (parser.key_exists("input")){
    	inputfile = parser.value("input");
    }
    if (parser.key_exists("file")){
    	filename = parser.value("file");
    }else{
    	throw std::invalid_argument("file=name not provided");
    }
    if (parser.key_exists("zoom")){
    	zoom = std::stod(parser.value("zoom"));
    }
    if (parser.key_exists("no_orient")){
    	no_orient = (bool) std::stoi(parser.value("no_orient"));
    }
    if (parser.key_exists("l")){
    	L = std::stod(parser.value("l"));
    }
    if (parser.key_exists("npoints")){
    	npoints = std::stoi(parser.value("npoints"));
    }
    if (parser.key_exists("outfile")){
    	outfile = parser.value("outfile");
    }else outfile = filename;

    outfile += ".cube";


    // if L, no_orient are not given:
    // check if they are in the input file
    // if not guess it with current madness defaults
    std::ifstream filestream(inputfile);
    position_stream(filestream, "dft");
    std::string s;
    while (filestream >> s) {
    	if (s == "end") {
    		break;
    	} else if (L < 0.0 and (s == "L" or s == "l")) {
    		filestream >> L;
        } else if (s == "no_orient") {
            filestream >> no_orient;
        }
    }

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
    Molecule molecule(world, parser);

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
    	L = dummy.L();
    }
    FunctionDefaults<3>::set_cubic_cell(-L,L);

    if(world.rank()==0){
    	std::cout << "using madness parameters:\n";
    	std::cout << "\tthresh=" << f.thresh() << "\n";
    	std::cout << "\tk=" << f.k() << "\n";
    	std::cout << "\tL=" << L << "\n";
    	std::cout << "\tnpoints=" << npoints << "\n";
    	std::cout << "\tzoom=" << zoom << "\n";
    }

	// plot the Gaussian cube file

    if(world.rank()==0) std::cout << "creating cubefile " << outfile << "\n";
	plot_cubefile<3>(world,f,outfile,molecule.cubefile_header(), npoints,zoom);
    if(world.rank()==0) std::cout << "... success\n";

    return 0;

}
