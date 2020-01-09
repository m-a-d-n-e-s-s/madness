/*
 * QCCalculationParametersBase.h
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_
#define SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_

#include<string>
#include <algorithm>
#include<iomanip>
#include <typeindex>
#include <map>
#include <typeinfo>
#include <madness/misc/misc.h>
#include <madness/world/archive.h>
#include <madness/world/world.h>


namespace madness {


template<typename T>
static typename std::enable_if<std::is_floating_point<T>::value, void>::type
check_for_inf(const std::string& str, T& arg) {
	std::string sinf;
	std::stringstream ss;
	ss << std::numeric_limits<T>::infinity();
	sinf=ss.str();

	if (sinf==str) arg=std::numeric_limits<T>::infinity();
}

template<typename T>
static typename std::enable_if<!std::is_floating_point<T>::value, void>::type
check_for_inf(const std::string& str, T& arg) {
	return;
}

/// inverting the print method from print.h for std::vector
/// TODO: move this where it belongs (into print.h ??)
template <typename T, typename A=std::allocator<T> >
std::istream& operator>>(std::istream& is, std::vector<T,A>& v) {

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
	T tmp;
	while (sline >> word) {
		std::stringstream sword(word);
		sword >> tmp;
		check_for_inf(word,tmp);
		v.push_back(tmp);
	}
	if (sline.bad()) {
		madness::print("error while reading vector from istream: ");
		madness::print(line,"\n");
		MADNESS_EXCEPTION("IO error",1);
	}

	return is;
}




/// inverting the print method from print.h for std::vector
/// TODO: move this where it belongs (into print.h ??)
template <typename Q, typename T>
std::istream& operator>>(std::istream& is, std::pair<T,Q>& p) {

	// get all words from the line
	std::string word, line="";
	while (is >> word) {
		line+=word + " ";
	};
	// have to set this here to account for the error handling later..
	is.clear();

	// remove enclosing brackets and commas
	auto find_c = [](char& c){ return ((c==',') or (c=='(') or (c==')')); };
	std::replace_if (line.begin(), line.end(), find_c, ' '); // 0 2 0 4 0 6 0 8 0

	// stream the values into the container
	std::stringstream sline(line);
	T tmp1;
	Q tmp2;
	sline >> tmp1 >> tmp2;
	if (sline.bad() or sline.fail()) {
		madness::print("error while reading vector from istream: ");
		madness::print(line,"\n");
		MADNESS_EXCEPTION("IO error",1);
	}
	p=std::pair<T,Q>(tmp1,tmp2);

	return is;
}


/// structure holding the value for a given parameter

/// keeps logic about default, derived and user-defined settings (with increasing priority),
/// as well as comments for the user (e.g. a recommended range for a given parameter).
///
/// might be extended to hold allowed values for certain parameters, e.g. localization procedures.
///
/// all values are stored as strings and must be converted to their respective types by the
/// QCCalculationParametersBase class (see below)
struct QCParameter {
public:
	QCParameter() {};

	QCParameter(const std::string v, const std::string t, const std::string comment="",
			const std::vector<std::string> allowed_values1={})
	: default_value(v), type(t), comment(comment), allowed_values(allowed_values1) {
		static int i=0;
		print_order=i++;
		set_all();
	}

	void set_derived_value(const std::string val) {
		precedence=std::max(derived,precedence);
		derived_value=val;
		set_all();
	}
	void set_user_defined_value(const std::string val) {
		precedence=std::max(defined,precedence);
		user_defined_value=val;
		set_all();
	}

	std::string get_value() const {return value;}
	std::string get_type() const {return type;}
	std::string get_comment() const {return comment;}
	std::string print_precedence() const {
		if (precedence==def) return "default";
		if (precedence==derived) return "derived";
		if (precedence==defined) return "defined";
		std::stringstream sprecedence;
		sprecedence << precedence;
		throw std::runtime_error("unknown precedence in QCParameter"+sprecedence.str());
		return "darn";
	}

	int get_print_order() const {return print_order;}

	std::string print_line(const std::string& key) const {

		auto fill_left = [](const int size, const std::string word) {
			int nspaces=std::max(int(0),size-int(word.length()));
			return std::string(nspaces, ' ')+word;
		};
		auto fill_right = [](const int size, const std::string word) {
			int nspaces=std::max(int(0),size-int(word.length()));
			return word+std::string(nspaces, ' ');
		};

		std::string result=fill_left(20,key)+"  "+fill_right(10,get_value()) + " # "
				+fill_right(10,print_precedence())
				//				+fill_right(5,get_type())
				+ fill_right(45,get_comment());
		if (allowed_values.size()>0) {
			std::stringstream ss;
			ss << allowed_values;
			result+=ss.str();
		}

		// trim result
		std::size_t last = result.find_last_not_of(' ');
		return result.substr(0, last+1);
	}

	template <typename Archive> void serialize (Archive& ar) {
		ar & value & default_value & derived_value & user_defined_value & type & null &
		comment & allowed_values & print_order & precedence;
	}

	enum {def, derived, defined} precedence=def;

private:

	void set_all() {
		value=default_value;
		if (derived_value!=null) value=derived_value;
		if (user_defined_value!=null) value=user_defined_value;
		if (not check_allowed()) throw std::runtime_error(not_allowed_errmsg());
	}

	bool check_allowed() {
		if (allowed_values.size()==0)  return true;
		auto it = std::find(allowed_values.begin(), allowed_values.end(), value);
		return (it!=allowed_values.end());
	}

	std::string not_allowed_errmsg() const {
		std::stringstream ss;
		ss<< allowed_values;
		std::string errmsg="\ntrying to assign a value that's not allowed\n\n";
		errmsg+="\tuser-defined value: " + value + "\n";
		errmsg+="\tallowed values:     " + ss.str() + "\n\n";
		return errmsg;
	}


	std::string value="";
	std::string default_value="";
	std::string derived_value="";
	std::string user_defined_value="";
	std::string type="";
	std::string null="";
	std::string comment="";
	std::vector<std::string> allowed_values=std::vector<std::string>();
	int print_order=0;		// use this for printing the parameters in the same order as they are defined

};


/// class for holding the parameters for calculation

/// Actual parameter classes will be derived from this class with a simple constructor
/// (see test_QCCalculationParametersBase.cc for an example) and convenience
/// getters for the parameters of each parameter class.
/// Having the base class will allow consistent parameter input/output handling
/// for all madness programs and reuse of the parsing methods.
/// Even if the same parameter is used in different programs, default might differ
/// (i.e. econv for 3D/6D calculations). The parameter class will effectively serve as
/// a factory for the calculation classes (SCF, nemo, mp2, etc)
///
/// parameters are kept in a map with key (string) and value (QCParameter),
/// types are converted whenever a parameter is accessed (i.e. should not be
/// done too frequently in an inner loop)
class QCCalculationParametersBase {

public:
	/// print all parameters
	void print(const std::string header="", const std::string footer="") const;

	std::string print_to_string(bool non_defaults_only=false) const;

	template<typename T>
	T get(const std::string key) const {
		const QCParameter& parameter=get_parameter(key);
		MADNESS_ASSERT(check_type<T>(key,parameter));
		return fromstring<T>(parameter.get_value());
	}

	template <typename Archive> void serialize (Archive& ar) {
		ar & parameters & print_debug;
	}

protected:

	typedef std::map<std::string,QCParameter> ParameterContainerT;
	ParameterContainerT parameters;

	/// read the parameters from file

	/// only world.rank()==0 reads the input file and broadcasts to all other nodes,
	/// so we don't need to serialize the ParameterMap
	virtual void read(World& world, const std::string filename, const std::string tag);


	bool print_debug=false;

	/// ctor for testing
	QCCalculationParametersBase() {}

	/// copy ctor
	QCCalculationParametersBase(const QCCalculationParametersBase& other)
		: parameters(other.parameters)
		, print_debug(other.print_debug) {
	}

	/// destructor
	virtual ~QCCalculationParametersBase() {}

	template<typename T>
	void initialize(const std::string& key, const T& value, const std::string comment="",
			const std::vector<T> allowed_values={}) {

		if (parameters.find(key)!=parameters.end()) {
			madness::print("you cannot initialize a parameter twice: ",key);
			throw std::runtime_error("initialization error");
		}

		std::string svalue=tostring(value);
		std::string type = std::type_index(typeid(T)).name();

		// transform everything to lower case
		std::string key_lower=key;
		std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(), ::tolower);
		std::transform(svalue.begin(), svalue.end(), svalue.begin(), ::tolower);
		std::vector<std::string> av_lower_vec;
		for (auto av : allowed_values) {
			std::string av_lower=tostring(av);
			std::transform(av_lower.begin(), av_lower.end(), av_lower.begin(), ::tolower);
			av_lower_vec.push_back(av_lower);
		}

		parameters.insert(std::make_pair<std::string, QCParameter>
		(std::string(key_lower),QCParameter(svalue,type,comment,av_lower_vec)));
	}

public:
	template<typename T>
	void set_derived_value(const std::string& key, const T& value) {

		QCParameter& parameter=get_parameter(key);
		if (not check_type_silent<T>(parameter)) {
			throw std::runtime_error("type error in set_derived_value for key "+key);
		}
		parameter.set_derived_value(tostring(value));
	}

protected:

	template<typename T>
	bool try_setting_user_defined_value(const std::string& key, const std::string& val) {

		if (not check_type_silent<T>(get_parameter(key))) return false;

		if (print_debug) ::madness::print("key:",key,"will set type" ,std::type_index(typeid(T)).name());
		T value=fromstring<T>(val);
		set_user_defined_value<T>(key,value);
		return true;
	}

public:
	template<typename T>
	void set_user_defined_value(const std::string& key, const T& value) {

		QCParameter& parameter=get_parameter(key);
		if (not check_type_silent<T>(parameter)) {
			throw std::runtime_error("type error in set_user_defined_value");
		}

		parameter.set_user_defined_value(tostring(value));
	}
protected:

	const QCParameter& get_parameter(const std::string key) const {
		if (not parameter_exists(key)) {
			throw std::runtime_error("could not find parameter for key "+key);
		}
		const QCParameter& parameter=parameters.find(key)->second;
		return parameter;
	}

public:
	QCParameter& get_parameter(const std::string key) {
		if (not parameter_exists(key)) {
			throw std::runtime_error("could not find parameter for key "+key);
		}
		QCParameter& parameter=parameters.find(key)->second;
		return parameter;
	}

	bool parameter_exists(const std::string& key) const {
		return (parameters.find(key)!=parameters.end());
	}


	template<typename T>
	static bool check_type(const std::string key, const QCParameter& parameter) {
		if (check_type_silent<T>(parameter)) return true;

		madness::print("trying to get the wrong type in QCCalculationParametersBase");
		madness::print("key             ",key);
		madness::print("parameter type  ",parameter.get_type());
		madness::print("setting type    ",std::type_index(typeid(T)).name());
		madness::print("value           ",parameter.get_value());
		return false;
	}

	template<typename T>
	static bool check_type_silent(const QCParameter& parameter) {
		return (parameter.get_type()==std::type_index(typeid(T)).name());
	}

	/// read the stream, starting from tag

	/// only parameters that are defined in the constructor will be processed,
	/// all others will be discarded.
	virtual void read_internal(World& world, std::string& filecontents, std::string tag);


	static std::string tostring(const bool& arg) {
		std::ostringstream ss;
		ss << std::boolalpha << arg;
		return ss.str();
	}

	template<typename T>
	static std::string tostring(const T& arg) {
		std::ostringstream ss;

		ss<<std::scientific  << std::setprecision(4) << arg;
		std::string str=ss.str();
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);

		overwrite_if_inf(str,arg);
		return str;
	}

	template<typename T>
	static typename std::enable_if<std::is_floating_point<T>::value, void>::type
	overwrite_if_inf(std::string& str, const T& arg) {
		if (std::isinf(arg)) {
			std::stringstream ss;
			ss << std::numeric_limits<T>::infinity();
			str=ss.str();
		}
	}

	template<typename T>
	static typename std::enable_if<!std::is_floating_point<T>::value, void>::type
	overwrite_if_inf(std::string& str, const T& arg) {
		return;
	}

	template<typename T>
	static typename std::enable_if<!std::is_same<T,bool>::value, T>::type
	fromstring(const std::string& arg) {

		std::stringstream ssvalue(arg);
		T result=T();
		ssvalue >> result;

		bool type_conversion_failed=ssvalue.fail();

		// check for infinity in floating point conversions
		if (type_conversion_failed and (std::is_floating_point<T>::value)) {

			const static T inf=std::numeric_limits<T>::infinity();
			std::string sinf=tostring(inf);         // repeat type conversion from above
			if (sinf==arg) result=inf;
			type_conversion_failed=false;
		}

		if (type_conversion_failed) {

			std::string errmsg="error in type conversion for argument >> " + arg
					+ " << to type " + std::type_index(typeid(T)).name();
			throw std::runtime_error(errmsg);
		}

		// check for trailing characters
		std::string word;
		while (ssvalue >> word) {
			std::string errmsg="trailing characters in arguement >> " + arg + " <<";
			throw std::runtime_error(errmsg);
		}
		return result;
	}


	template<typename T>
	static typename std::enable_if<std::is_same<T,bool>::value, T>::type
	fromstring(const std::string& arg) {
		std::string str=arg;
		std::transform(str.begin(), str.end(), str.begin(), ::tolower);
		if (str=="true" or str=="1" or str=="yes") return true;
		if (str=="false" or str=="0" or str=="no") return false;
		std::string errmsg="error in type conversion for argument >> " + arg
				+ " << to type " + std::type_index(typeid(T)).name();
		throw std::runtime_error(errmsg);
		return 0;
	}

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_QCCALCULATIONPARAMETERSBASE_H_ */
