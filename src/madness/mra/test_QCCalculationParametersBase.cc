/*
 * test_QCCalculationParametersBase.cc
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */



#include"QCCalculationParametersBase.h"
#include"parallel_archive.h"

using namespace madness;

class Parameters : public QCCalculationParametersBase {
public:
	Parameters() : QCCalculationParametersBase() {

		// initialize with: key, value, comment (optional), allowed values (optional)
		initialize<double>("econv",1.e-5,"recommended values: 1.e-4 < econv < 1.e-8");
		initialize<double>("dconv",3.e-4,"recommended values: 1.e-4 < econv < 1.e-8");
		initialize<int>("maxiter",15);
		initialize<bool>("localize",true);
		initialize<std::string>("local","BOYS","localization method",{"boys","pm","canon"});
		initialize<std::string>("xc","hf","free-form string value, blanks and all");
		initialize<std::vector<double> >("proto",std::vector<double>{1,2});
		initialize<std::pair<std::string,double> >("ncf",{"slater",2.0});
	}

	std::string get_tag() const override {
		return std::string("test");
	}


	void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
		read_input_and_commandline_options(world,parser,tag);
		set_derived_value("dconv",sqrt(get<double>("econv"))*0.1);
	}

	// convenience getters
	double econv() const {return get<double>("econv");}
	double dconv() const {return get<double>("dconv");}
	bool localize() const {return get<bool>("localize");}
	std::string local() const {return get<std::string>("local");}
	std::string xc() const {return get<std::string>("xc");}
	std::pair<std::string,double> ncf() const {return get<std::pair<std::string,double> >("ncf");}
	int maxiter() const {return get<int>("maxiter");}


};

template<typename T>
void test_same(const T& t1, const T& t2) {
	if (t1!=t2) {
		print("t1, t2", t1,t2);
		using madness::operators::operator<<;
		std::cout << "++"<< t1 << "++"<< std::endl;
		std::cout << "++" << t2 << "++" <<std::endl;

		throw std::runtime_error("failure in test");;
	}
}

/// check that an error message carries \p snippet, wherever in the message it sits
void test_contains(const std::string& msg, const std::string& snippet) {
	if (msg.find(snippet)==std::string::npos) {
		std::cout << "++" << msg << "++" << std::endl;
		std::cout << "does not contain ++" << snippet << "++" << std::endl;

		throw std::runtime_error("failure in test");
	}
}

struct inputfile {
	std::string fname;
	inputfile(const std::string filename, const std::string lines) {
		fname=filename;
		std::ofstream myfile;
		myfile.open (fname);
		myfile << lines << std::endl;
		myfile.close();
	}

	~inputfile() {
		remove(fname.c_str());
	}
};

/// a deck that cannot be parsed must be rejected, naming the offending key

/// these used to pass vacuously: read_input_and_commandline_options printed the
/// error and swallowed it, so the catch block below was never entered
void test_deck_is_rejected(World& world, const std::string& inputlines, const std::string& key) {
	inputfile ifile("input1",inputlines);

	bool found_exception=false;
	try {
		Parameters param;
		commandlineparser parser;
		parser.set_keyval("input","input1");
		param.read_and_set_derived_values(world,parser,"mp3");
	} catch (std::exception& err) {
		found_exception=true;
		test_contains(std::string(err.what()),"found an error for key >> "+key);
	}
	if (not found_exception) throw std::runtime_error("expected exception for key "+key);
}

bool test_serialize(World& world) {

	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
	param.read_and_set_derived_values(world,parser,"mp3");


	param.print("defined parameters","foot\n\n");

	const std::string name = "test.dat";
	madness::archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> oar(world, name.c_str(), 1);
	oar & param;
	oar.close();

	Parameters param1;
	param1.print("default parameters","foot\n\n");

	madness::archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> iar(world, name.c_str(), 1);
	iar & param1;
	iar.close();

	param1.print("serialized parameters","foot\n\n");

	test_same(param.econv(),param1.econv());
	test_same(param.dconv(),param1.dconv());
	test_same(param.maxiter(),param1.maxiter());
	test_same(param.ncf(),param1.ncf());

	return true;
}

bool test_type_conversion1(World& world) {
	print("entering test_type_conversion1");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.ncf(),std::pair<std::string,double>("slater",1.2));
	return true;
}


bool test_type_conversion2(World& world) {
	print("entering test_type_conversion2");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.ncf(),std::pair<std::string,double>("slater",1.2));
	return true;
}


bool test_type_conversion3(World& world) {
	print("entering test_type_conversion3");
	std::string inputlines=R"input(mp3 
			econv 1.d-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2
			end)input";
	test_deck_is_rejected(world,inputlines,"econv");
	return true;
}

bool test_type_conversion4(World& world) {
	print("entering test_type_conversion4");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2 
			localize tru
			end)input";
	test_deck_is_rejected(world,inputlines,"localize");
	return true;
}


bool test_type_conversion5(World& world) {
	print("entering test_type_conversion5");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			localize true
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.localize(),true);
	return true;
}

bool test_type_conversion6(World& world) {
	print("entering test_type_conversion6");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			localize 1
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.localize(),true);
	return true;
}


bool test_type_conversion7(World& world) {
	print("entering test_type_conversion7");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			localize 0
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.localize(),false);
	return true;
}

bool test_type_conversion8(World& world) {
	print("entering test_type_conversion8");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.localize(),false);
	return true;
}

/// a value that does not convert must abort the run, not be printed and ignored

/// swallowing it left the calculation running with the default value, so a deck
/// that asked for something the parser choked on quietly computed something else
bool test_trailing_characters(World& world) {
	print("entering test_trailing_characters");
	std::string inputlines=R"input(mp3 
			econv 1.e-4a   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2
			end)input";
	test_deck_is_rejected(world,inputlines,"econv");
	return true;
}

/// blank-separated values for free-form string keys, e.g. `xc GGA_X_PBE 0.75 ...`

/// double quotes used to be mandatory: an unquoted line kept only the first word
/// and then tripped the trailing-character check, so the value was discarded
bool test_blank_separated_string(World& world) {
	print("entering test_blank_separated_string");
	const std::string expected("gga_x_pbe 0.75 gga_c_pbe 1.0 hf_x 0.25");

	std::vector<std::string> variants={
		"xc GGA_X_PBE 0.75 GGA_C_PBE 1.0 HF_X 0.25   # spelled-out PBE0",
		"xc \"GGA_X_PBE 0.75 GGA_C_PBE 1.0 HF_X 0.25\"    # the historical spelling",
		"xc=GGA_X_PBE 0.75 GGA_C_PBE 1.0 HF_X 0.25",
		"xc GGA_X_PBE 0.75 GGA_C_PBE 1.0 HF_X 0.25\t  ",
		// the line as print_to_string() writes it: the printed deck must be re-readable
		"                  xc  gga_x_pbe 0.75 gga_c_pbe 1.0 hf_x 0.25 # defined   free-form"
	};
	for (const auto& variant : variants) {
		inputfile ifile("input1","mp3\n"+variant+"\nend");
		Parameters param;
		commandlineparser parser;
		parser.set_keyval("input","input1");
		param.read_and_set_derived_values(world,parser,"mp3");
		test_same(param.xc(),expected);
	}

	// .. but a key with no value at all is still an error
	{
		inputfile ifile("input1","mp3\nxc\nend");
		bool found_exception=false;
		try {
			Parameters param;
			commandlineparser parser;
			parser.set_keyval("input","input1");
			param.read_and_set_derived_values(world,parser,"mp3");
		} catch (std::exception&) {
			found_exception=true;
		}
		if (not found_exception) throw std::runtime_error("expected exception for a key without a value");
	}
	return true;
}

/// the same values, coming in through --<tag>="key=val; key=val"
bool test_blank_separated_string_commandline(World& world) {
	print("entering test_blank_separated_string_commandline");

	Parameters param;
	commandlineparser parser;
	parser.set_keyval("input","input1");   // does not exist -- defaults only
	parser.set_keyval("mp3","xc=GGA_X_PBE 0.75 GGA_C_PBE 1.0 HF_X 0.25; maxiter=7");
	param.read_and_set_derived_values(world,parser,"mp3");

	test_same(param.xc(),std::string("gga_x_pbe 0.75 gga_c_pbe 1.0 hf_x 0.25"));
	test_same(param.maxiter(),7);
	return true;
}

bool test_capitalization(World& world) {
	print("entering test_capitalization");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
 			# econv 1.e-3
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			LocAl CanON
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.local(),std::string("canon"));
	return true;
}

bool test_not_allowed(World& world) {
	print("entering test_not_allowed");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
 			# econv 1.e-3
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			LocAl canoni
			end)input";
	inputfile ifile("input1",inputlines);

	bool found_exception=true;
	try {
		Parameters param;
        commandlineparser parser;
        parser.set_keyval("input","input1");
        param.read_and_set_derived_values(world,parser,"mp3");


		found_exception=false;
    } catch (std::invalid_argument& err) {
        // the message names the offending keyword before explaining the problem
        std::string errmsg(err.what());
        test_contains(errmsg,std::string("in keyword `local`"));
        test_contains(errmsg,std::string("trying to assign a value that's not allowed"));
	} catch (std::runtime_error& err) {
		std::string errmsg=std::string(err.what()).substr(0,30);
		test_same(errmsg,std::string("found an error for key >> loca"));
	}
	if (not found_exception) throw std::runtime_error("expected exception");
	return true;
}


bool test_comment_lines(World& world) {
	print("entering test_commentlines");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
 			# econv 1.e-3
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.econv(),1.e-4);
	return true;
}
bool test_empty_lines(World& world) {
	print("entering test_empty_lines");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			dconv 1.e-4
 				
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.econv(),1.e-4);
	test_same(param.dconv(),1.e-4);
	return true;
}



bool test_derived(World& world) {
	print("entering test_derived");
	std::string inputlines=R"input(mp3 
			econv 1.e-4   
			#dconv 1.e-4
			maxiter 12# asd 
			ncf (slater,1.2)
			localize no
			end)input";
	inputfile ifile("input1",inputlines);

	Parameters param;
    commandlineparser parser;
    parser.set_keyval("input","input1");
    param.read_and_set_derived_values(world,parser,"mp3");


	test_same(param.econv(),1.e-4);
	test_same(param.dconv(),sqrt(param.econv())*0.1);
	return true;
}




int main(int argc, char** argv) {
	initialize(argc, argv);
	World world(SafeMPI::COMM_WORLD);
//	startup(world,argc,argv);
	std::cout.precision(6);

	int success=0;
	try {
		test_serialize(world);
		test_type_conversion1(world);
		test_type_conversion2(world);
		test_type_conversion3(world);
		test_type_conversion4(world);
		test_type_conversion5(world);
		test_type_conversion6(world);
		test_type_conversion7(world);
		test_type_conversion8(world);
		test_capitalization(world);
		test_not_allowed(world);
		test_trailing_characters(world);
		test_blank_separated_string(world);
		test_blank_separated_string_commandline(world);
		test_comment_lines(world);
		test_empty_lines(world);
		test_derived(world);


	} catch (std::exception& e) {
		print("\n\tan error occurred .. ");
		print(e.what());
		success=1;
	} catch (...) {
		print("\n\tan unknown error occurred .. ");
		success=1;
	}


	if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
	world.gop.fence();
	finalize();

	return success;
}
