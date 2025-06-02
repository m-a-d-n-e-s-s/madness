/*
 * test_qc.cc
 *
 *  Created on: Aug 13, 2019
 *      Author: fbischoff
 */


#include <madness/mra/mra.h>
#include<madness/chem/CalculationParameters.h>
#include<madness/chem/SCF.h>

using namespace madness;


/// helper class to set up a vector of possible values for a given key

/// derived from the CalculationParameters class; parameters member variable is kept constant,
/// the variations are kept in the parameter_variations variable.
class TestCalculationParameters : public CalculationParameters {
	typedef std::map<std::string,std::vector<QCParameter> > ParameterVectorContainerT;
	ParameterVectorContainerT parameter_variations;

public:

	TestCalculationParameters(const CalculationParameters& cp) : CalculationParameters(cp) {}

	template<typename T>
	void extend_parameters(const std::string& key, std::vector<T> values) {

		// check if key exists at all
		QCParameter parameter=get_parameter(key);
		if (not check_type_silent<T>(parameter)) {
			throw std::runtime_error("type error in set_user_defined_value");
		}

		// set up the different values for this key
		std::vector<QCParameter> parametervalues(values.size(),parameter);
		for (std::size_t i=0; i<values.size(); i++) {
			parametervalues[i].set_user_defined_value(tostring(values[i]));
		}
		parameter_variations.insert(std::make_pair(key,parametervalues));
	}

	CalculationParameters copy_and_replace_key_in_parameters(const CalculationParameters& ref_parameters, const std::string& key, int index) const {

		std::size_t size=parameter_variations.find(key)->second.size();
		if (index>size) throw std::runtime_error("asdf");

		const QCParameter p=parameter_variations.find(key)->second[index];
		CalculationParameters cp(ref_parameters);
		cp.get_parameter(key)=p;
		return cp;
	}

	std::size_t get_parameter_range(const std::string& key) const {
		return parameter_variations.find(key)->second.size();
	}

	std::vector<CalculationParameters> make_all_parameters_for_one_key(const std::string& key) const {
		std::vector<CalculationParameters> vcp;
		for (std::size_t i=0; i<get_parameter_range(key); ++i) {
			vcp.push_back(copy_and_replace_key_in_parameters(*this,key,i));
		}
		return vcp;
	}

	std::vector<CalculationParameters> make_all_parameters_for_two_keys(const std::string& key1, const std::string& key2) const {
		std::vector<CalculationParameters> vcp;
		for (std::size_t i=0; i<get_parameter_range(key1); ++i) {
			for (std::size_t j=0; j<get_parameter_range(key2); ++j) {
				CalculationParameters cp1=copy_and_replace_key_in_parameters(*this,key1,i);
				vcp.push_back(copy_and_replace_key_in_parameters(cp1,key2,j));
			}
		}
		return vcp;
	}

	std::vector<CalculationParameters> make_all_parameter_singles() const {
		std::vector<CalculationParameters> vcp;
		for (auto pv : parameter_variations) {
			std::string key=pv.first;
			madness::print("make all single variations for",key);
			auto vcp1=make_all_parameters_for_one_key(key);
			vcp.insert(vcp.end(),vcp1.begin(),vcp1.end());
		}
		madness::print("made ",vcp.size()," singles permutations");
		return vcp;
	}

	std::vector<CalculationParameters> make_all_parameter_doubles() const {
		std::vector<CalculationParameters> vcp;
		for (ParameterVectorContainerT::const_iterator pv1 = parameter_variations.begin(); pv1 != parameter_variations.end(); ++pv1) {
			for (ParameterVectorContainerT::const_iterator pv2 = pv1; pv2 != parameter_variations.end(); ++pv2) {
				if (pv1==pv2) continue;
				std::string key1=pv1->first;
				std::string key2=pv2->first;
				madness::print("make all double variations for keys ",key1, key2);
				auto vcp1=make_all_parameters_for_two_keys(key1,key2);
				vcp.insert(vcp.end(),vcp1.begin(),vcp1.end());
			}
		}
		madness::print("made ",vcp.size()," doubles permutations");
		return vcp;
	}

};



/// will write a test input and remove it from disk upon destruction
struct write_test_input {

	double eprec=1.e-6;

	std::string filename_;
	write_test_input(const TestCalculationParameters& param, const std::string& mol="lih") : filename_("test_input") {
		std::ofstream of(filename_);
		of << "dft\n";
		of << param.print_to_string(true);
		of << "end\n";

		if (mol=="lih") {
			// of << "geometry\n";
			of << "molecule\n";
			of << "eprec " << eprec << std::endl;
			of << "Li 0.0    0.0 0.0\n";
			of << "H  1.4375 0.0 0.0\n";
			of << "end\n";
		} else if (mol=="hf") {
			double eprec=1.e-5;
			// of << "geometry\n";
			of << "molecule\n";
			of << "eprec " << eprec << std::endl;
			of << "F  0.1    0.0 0.2\n";
			of << "H  1.4375 0.0 0.0\n";
			of << "end\n";
		}
		of.close();
	}

	~write_test_input() {
		std::remove(filename_.c_str());
	}

	std::string filename() const {return filename_;}
};

int run_all_calculations(World& world, const std::vector<CalculationParameters>& all_parameters,
		bool test_result=false) {
	int success=0;
	for (auto cp : all_parameters) {

		print(cp.print_to_string(true));

		write_test_input test_input(cp,"lih");
        commandlineparser parser;
        parser.set_keyval("input",test_input.filename());
        SCF calc(world,parser);
		calc.set_protocol<3>(world, 1e-4);
		MolecularEnergy ME(world, calc);
		double energy=ME.value(calc.molecule.get_all_coords().flat()); // ugh!
		print("energy(LiH)",energy);
		if (test_result) {
			double thresh=cp.econv();
			print("energy, hard-wire, diff",energy,7.703833,energy+7.703833e+00);
			if (std::abs(energy+7.703833e+00)>thresh) success+=1;
		}
	}
	return success;
}

int main(int argc, char** argv) {


    madness::World& world=madness::initialize(argc, argv);
    startup(world,argc,argv);

    int result=0;
    bool small = true;
    

    // default set of parameters for closed shell
    CalculationParameters cparam;
    cparam.set_user_defined_value("print_level",10);
    cparam.set_user_defined_value("save",false);

    try {
        // check for correct result
    	{
    		result+=run_all_calculations(world, {cparam},true);
    		if (result>0) throw std::runtime_error("moldft returns incorrect result");
    	}

    	// check for no error (run 1 iteration only)
    	{
    		// store the variations of the default set
    		TestCalculationParameters tparam(cparam);
			tparam.set_user_defined_value("maxiter",1);

			if (small) {
				tparam.extend_parameters<double>("econv",{1.e-4});
				//tparam.extend_parameters<bool>("derivatives",{true}); // ,false
				//tparam.extend_parameters<int>("k",{6});
				//tparam.extend_parameters<double>("l",{25.0});
			}
			else {
				tparam.extend_parameters<double>("econv",{1.e-4,1.e-5}); // default and higher accuracy
				tparam.extend_parameters<std::string>("localize",{"canon","boys","new"}); //
				tparam.extend_parameters<bool>("spin_restricted",{false});
				tparam.extend_parameters<bool>("no_orient",{true});
				tparam.extend_parameters<bool>("derivatives",{true,false}); //
				tparam.extend_parameters<bool>("dipole",{true});
				tparam.extend_parameters<std::string>("xc",{"hf","lda"}); //
				tparam.extend_parameters<int>("k",{6});
				tparam.extend_parameters<double>("l",{25.0});
			}

			std::vector<CalculationParameters> all_singles=tparam.make_all_parameter_singles();
			run_all_calculations(world, all_singles);
		}

		// need to modify two parameters for this test
		if (! small)
		{
			// store the variations of the default set
			TestCalculationParameters tparam(cparam);
			tparam.set_user_defined_value("charge",1.0);
			tparam.set_user_defined_value("nopen",1);
			run_all_calculations(world, {tparam});
    	}
    } catch (std::exception& e) {
    	print("an error occurred, moldft tests failed");
    	print(e.what());
    	result=1;
    }


//    std::vector<CalculationParameters> all_doubles=tparam.make_all_parameter_doubles();
//    run_all_calculations(world, all_doubles);
    madness::finalize();
    return result;
}
