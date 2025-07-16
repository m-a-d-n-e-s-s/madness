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
	inputfile ifile("input1",inputlines);

	try {
		Parameters param;
        commandlineparser parser;
        parser.set_keyval("input","input1");
        param.read_and_set_derived_values(world,parser,"mp3");


	} catch (std::runtime_error& err) {
		std::string errmsg=std::string(err.what()).substr(0,30);
		test_same(errmsg,std::string("found an error for key >> econ"));
	}
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
	inputfile ifile("input1",inputlines);

	try {
		Parameters param;
        commandlineparser parser;
        parser.set_keyval("input","input1");
        param.read_and_set_derived_values(world,parser,"mp3");


	} catch (std::runtime_error& err) {
		std::string errmsg=std::string(err.what()).substr(0,30);
		test_same(errmsg,std::string("found an error for key >> loca"));
	}
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

template<typename T>
bool test_trailing_characters(World& world) {
	print("entering test_trailing_characters");
	std::string inputlines=R"input(mp3 
			econv 1.e-4a   
			dconv 1.e-4
			maxiter 12# asd 
			ncf slater 1.2
			end)input";
	inputfile ifile("input1",inputlines);

	try {
		Parameters param;
        commandlineparser parser;
        parser.set_keyval("input","input1");
        param.read_and_set_derived_values(world,parser,"mp3");


	} catch (std::runtime_error& err) {
		std::string errmsg=std::string(err.what()).substr(0,30);
		test_same(errmsg,std::string("found an error for key >> econ"));
	}
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
        std::string errmsg=std::string(err.what()).substr(0,30);
        test_same(errmsg,std::string("\ntrying to assign a value that"));
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
		test_type_conversion5(world);
		test_type_conversion6(world);
		test_type_conversion7(world);
		test_type_conversion8(world);
		test_capitalization(world);
		test_not_allowed(world);
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
