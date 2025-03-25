#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include "madness/mra/commandlineparser.h"
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/molecule.h>

using namespace madness;

template<std::size_t NDIM>
void do_plot (World& world, const std::vector<std::string>& filenames, const PlotParameters& pparam) {

	std::vector<Function<double,NDIM>> vf;
	for (const auto& filename : filenames) {
		Function<double,NDIM> f = FunctionFactory<double,NDIM>(world);
		if(world.rank()==0) std::cout << "load function " << filename << "\n";
		load(f,filename);
		f.get_impl()->verify_tree_state_local();

		// annoying temporary fix
		TensorArgs targs(FunctionDefaults<NDIM>::get_tensor_type(),f.get_impl()->get_thresh());
		f.change_tensor_type(targs);

		if(world.rank()==0) std::cout << "... success\n";
		vf.push_back(f);
	}

	std::string outfile=filenames[0]+".plane";
	if(world.rank()==0) std::cout << "creating plot data in file " << outfile << "\n";
	plot_plane(world,vf,outfile,pparam);
	if(world.rank()==0) std::cout << "... success\n";

}




int main(int argc, char** argv) {
	World& world=initialize(argc, argv);
	startup(world,argc,argv);
	std::cout.precision(6);

	commandlineparser parser(argc,argv);
	bool printhelp=parser.key_exists("help");

	// if no filename is given print help
	std::vector<std::string> filenames;
	if (parser.key_exists("file")) {
		filenames = parser.split(parser.value("file"),",");
	} else {
		std::cout << "\ninvalid usage: no filename(s) given\n";
		printhelp=true;
	}

	// check if file exists using c++ filesystem
	try {
		for (const auto& filename : filenames) {
			bool file_ok=std::filesystem::exists(filename+".00000");
			// if (not file_ok) throw std::runtime_error("file not found");
			if (not file_ok) {
				std::cout << "\ncould not find data file: "+filename+"\n";
				printhelp=true;
			}
		}
	} catch (...) {
		std::cout << "\ncould not find data file(s): ";
		printhelp=true;
	}


	if (printhelp) {
		std::cout << "\nUSAGE: plot2plane file=name1[,name2][,name3] [options]\n\n";
		std::cout << "name1,name2,..  names of the files to plot, each must contain one 3d function only \n";
		std::cout << "                filename without extension, e.g. nemo4 instead of nemo4.00000\n\n";
		std::cout << "OPTIONS: command line options similar to the electronic structure calculations \n";
		std::cout << "         will read the input file if it exists, command line options will override the input file\n";
		std::cout << std::endl;
		std::cout << "           e.g. plot2plane file=nemo4,nemo0 --geometry=h2o\n";
		std::cout << "           e.g. plot2plane file=nemo4 --geometry=\"no_orient=true; source_name=h2o\"\n";
		std::cout << "           e.g. plot2plane file=nemo4 --geometry=h2o --plot=\"zoom 4;  origin=[3, 3, 3]\"\n";
		std::cout << "           e.g. plot2plane file=nemo4 --geometry=h2o --plot=\"plane=x1x2;  origin=[3, 3, 3]\"\n\n";
		std::cout << std::endl;
		std::cout << "default plot options are: \n";
		PlotParameters pparam;
		pparam.print("plot","end");
		finalize();
		return 0;
	}


	// get molecule and calculation parameters
	Molecule molecule(world,parser);
	CalculationParameters param(world,parser);
	AtomicBasisSet aobasis(param.aobasis()); // not actually needed here..
	param.set_derived_values(molecule,aobasis,parser);
	PlotParameters pparam(world,parser);

	if (world.rank()==0) {
		print_header2("plotting parameters");
		param.print("dft","end");
		molecule.print();
		pparam.print("plot","end");
	}

	// load functions
	FunctionDefaults<3>::set_cubic_cell(-param.L(),param.L());
	FunctionDefaults<6>::set_cubic_cell(-param.L(),param.L());


	try {
		do_plot<3>(world,filenames,pparam);
	} catch (...) {
		std::cout << "plotting in 3d failed\n";
	}

	try {
		do_plot<6>(world,filenames,pparam);
	} catch (...) {
		std::cout << "plotting in 6d failed\n";
	}

	finalize();
	return 0;
}
