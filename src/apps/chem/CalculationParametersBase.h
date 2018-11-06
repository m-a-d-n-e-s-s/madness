/*
 * CalculationParametersBase.h
 *
 *  Created on: 29 Oct 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_
#define SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_

#include<string>
#include<iomanip>
#include <algorithm>
#include <typeindex>
#include <typeinfo>
#include <madness/misc/misc.h>


namespace madness {

/**

	There are three main classes defined:
	 - Parameter, templated with the data type (e.g. double, int,..) holding the actual value of a parameter
	 - CalculationParametersBase, a class holding basic input/output functionality

 **/

/// inverting the print method from print.h for std::vector
template <typename Q>
std::istream& operator>>(std::istream& is, std::vector<Q>& v) {

	// get the full line from opening to closing brackets [ .. ]
	std::string word, line="";
	while (is >> word) {
		line+=word;
		if (word.find(']')!=std::string::npos) break;
	}

	// remove enclosing brackets and commas
	auto find_c = [](char& c){ return ((c==',') or (c=='[') or (c==']')); };
	std::replace_if (line.begin(), line.end(), find_c, ' '); // 0 2 0 4 0 6 0 8 0

	// stream the values into the container
	std::stringstream sline(line);
	Q tmp;
	while (sline >> tmp) v.push_back(tmp);
	if (sline.bad()) {
		madness::print("error while reading vector from istream: ");
		madness::print(line,"\n");
		MADNESS_EXCEPTION("IO error",1);
	}

    return is;
}


/// parameter class
template<typename T>
struct Parameter {

	/// the constructor takes the default value
	Parameter(std::string k, const T d) : key(k), value(d), default_value(d) {
		set_all();
	}

	/// set the user defined value
	void set_user_defined_value(const T u) {
		user_defined_value=u;
		set_all();
	}

	/// set a derived value
	void set_derived_value(const T d) {
		derived_value=d;
		set_all();
	}

	const std::string& get_key() const {return key;}
	/// return the value of the parameter
	T get() const {return value;}

	T operator()() const {return value;}

    const std::type_info& value_type_info() const {return typeid(value);}

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & value & default_value & user_defined_value & derived_value & null;
    }

private:
    const std::string key;		///< the key for the key/value pair
	T value;					///< the actual value
	const T default_value;		///< the default value
	T user_defined_value = T();	///< the value if defined by the user
	T derived_value = T();		///< the value if derived from another parameter
	const T null = T();			///< for comparison

	// set the final value in order: user > derived > default;
	void set_all() {
		value=default_value;
		if (derived_value!=null) value=derived_value;
		if (user_defined_value!=null) value=user_defined_value;
	}

    static bool stringtobool(std::string str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str=="true" or str=="1" or str=="yes") return true;
        if (str=="false" or str=="0" or str=="no") return false;
        throw("unknown boolean ");
        return 0;
    }

public:
	template<class Q = T>
	friend typename std::enable_if<std::is_same<Q,bool>::value, std::istream&>::type
	operator>>(std::istream& is, Parameter<T>& param) {
		std::string word;
		is >> word;
		param.set_user_defined_value(stringtobool(word));
		return is;
	}


	template<class Q = T>
	friend typename std::enable_if<!std::is_same<Q,bool>::value, std::istream&>::type
	operator>>(std::istream& is, Parameter<T>& param) {
		T word;
		is >> word;
		param.set_user_defined_value(word);
		return is;
	}

	friend std::ostream& operator<<(std::ostream& os, const Parameter<T>& param) {
		os << param();
		return os;
	}
};

/// base class for calculation parameters
class CalculationParametersBase {
public:

	/// ctor for testing
	CalculationParametersBase() {}

	/// ctor for testing
	virtual ~CalculationParametersBase() {}

	/// print all parameters in the original ordering
    virtual void print() const = 0;

    /// read the stream, starting from tag

    /// only parameters that are defined in the constructor will be processed,
    /// all others will be discarded.
    virtual void read(const std::string filename, const std::string tag) = 0;

    template<typename T>
    void print_twocolumn_centered(const std::string s, const T t, const std::string s2="") const {
    	int len=20-s.size();
    	for (int i=0; i<len; ++i) std::cout << " ";
    	std::ostringstream os;
    	os << std::scientific << std::boolalpha << t;

    	std::cout << s <<  "  "<<std::setw(15) << os.str() << "  " << s2<<std::endl;
    }
};


} // namespace madness


#endif /* SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_ */
