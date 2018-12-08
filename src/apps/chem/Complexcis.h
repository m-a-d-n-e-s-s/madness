/*
 * Complexcis.h
 *
 *  Created on: 21 Nov 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_COMPLEXCIS_H_
#define SRC_APPS_CHEM_COMPLEXCIS_H_

#include <chem/Nemocomplex.h>


namespace madness {



class Complex_CIS_Parameters : public CalculationParametersBase {
public:
	enum parameterenum {guess_excitation_operators_,exops_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<std::string>(guess_excitation_operators_,{"guess_excitation_operators",{"dipole+"}}),
        		init<std::vector<std::string> >(exops_,{"exops",{"x 1.0","y 1.0","z 1.0","x 2.0 , y 2.0 , z 2.0"}})
    };

	/// ctor reading out the input file
	Complex_CIS_Parameters(World& world) {

		// read input file
		read(world,"input","response",params);

		// set derived values
//		params[param2_].set_derived_value(this->get<int>(param1_)*10.0);

		// print final parameters
//		if (world.rank()==0) print(params,"Our parameters");
	}

	std::string guess_excitation_operators() const {return get<std::string>(guess_excitation_operators_);};
	std::vector<std::string> exops() const {return get<std::vector<std::string> >(exops_);};

	/// return the value of the parameter
	template<typename T>
	T get(parameterenum k) const {
		if (params.find(int(k))!=params.end()) {
			return params.find(int(k))->second.get_parameter<T>().get();
		} else {
			MADNESS_EXCEPTION("could not fine parameter ",1);
		}
	}
};


class Complex_cis {
public:
	Complex_cis(World& w, Nemo_complex& n) : world(w), cis_param(world), nemo(n) {};
	virtual ~Complex_cis() {};

	double value();

	std::vector<std::vector<complex_function_3d> > read_guess(const std::string spin) const;

	std::vector<std::vector<complex_function_3d> > make_guess(const std::string spin) const;

	Tensor<double_complex> make_CIS_matrix(std::vector<complex_function_3d> virtuals) const;

	/// the world
	World& world;

	/// the parameters
	Complex_CIS_Parameters cis_param;

	/// the reference
	Nemo_complex nemo;

	/// the x vectors
	std::vector<std::vector<complex_function_3d> > roots;
};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_COMPLEXCIS_H_ */
