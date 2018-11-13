/*
 * CalculationParametersBase.h
 *
 *  Created on: 29 Oct 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_
#define SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_

#include<string>
#include <algorithm>
#include<iomanip>
#include <typeindex>
#include <typeinfo>
#include <madness/misc/misc.h>
#include <madness/world/archive.h>


namespace madness {

/**

	There are three main classes defined:
	 - Parameter, templated with the data type (e.g. double, int,..) holding the actual value of a parameter
	 - ParameterObject, not templated, implementing the type erasure idiom for parameter
	 - CalculationParametersBase, a class implementing basic input/output functionality

	The design idea is to have a map of parameters. Since parameters have different types,
	and there is no heterogeneous map the type erasure idiom has been implemented, hiding
	the actual parameter type from the map. This way all parameters are consistently set and
	are consistently printed.

	Parameter values can be set by different mechanisms, namely default, derived, and user-defined,
	with increasing priority.

	All parameters are defined and set to their default values in the constructor of the derived
	parameter class. User-defined values (and derived values, if applicable) are set through the
	input file, or from the command line (not yet implemented). An example for a derived parameter
	class is given below. A life example can be found in mp2.h/mp2.cc

	class DerivedParamterClass : public CalculationParametersBase {
		enum parameterenum {param1, param2};

		/// the parameters with the enum key, the constructor taking the input file key and a default value
		ParameterMap params={
            		init<int>(param1,{"param1_key",2}),
            		init<double>(param2,{"param2_key",1.e-3})
        };

		/// ctor reading out the input file
		Parameters(World& world) {

			// read input file
			read(world,"input","parameter_tag",params);

			// set derived values
			params[param2].set_derived_value(this->get<int>(param1)*10.0);

			// print final parameters
			if (world.rank()==0) print(params,"Our parameters");
		}

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


	To access a parameter value the type of the value has to be known at compile time. Since the
	syntax is a bit cumbersome shortcuts may be defined.

	..
	int param1() const {return param.get<int>(Parameters::param1);}
	..
	int use_this_parameter = param.param1();
	..

 **/


/// inverting the print method from print.h for std::vector
/// TODO: move this where it belongs (into print.h ??)
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

	/// @param[in]	k	the key for this parameter
	/// @param[in]	d	the default and current value
	Parameter(const std::string k, const T d) : key(k), value(d), default_value(d) {
		set_all();
	}

	/// set the user defined value

	/// @param[in]	u	a new, user-defined value
	void set_user_defined_value(const T u) {
		user_defined_value=u;
		set_all();
	}

	/// set a derived value

	/// certain parameters will influence others, e.g. econv will determine dconv
	/// @param[in]	u	a new, from other parameters derived value
	void set_derived_value(const T d) {
		derived_value=d;
		set_all();
	}

	/// return the key of this parameter
	const std::string& get_key() const {return key;}

	/// return the value of the parameter
	T get() const {return value;}

	/// return the value of the parameter
	T operator()() const {return value;}

    /// return the type of the parameter value
    const std::type_info& value_type_info() const {return typeid(value);}

private:
    std::string key;			///< the key in the input file
	T value;					///< the actual value
	T default_value;			///< the default value
	T user_defined_value = T();	///< the value if defined by the user
	T derived_value = T();		///< the value if derived from another parameter
	T null = T();				///< for comparison

	/// set the final value in order: user > derived > default;

	/// Use the default value if no other information is available.
	/// Always respect the user request!
	void set_all() {
		value=default_value;
		if (derived_value!=null) value=derived_value;
		if (user_defined_value!=null) value=user_defined_value;
	}

	/// helper function to convert alphanumerical bool into actual bool
    static bool stringtobool(std::string str) {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str=="true" or str=="1" or str=="yes") return true;
        if (str=="false" or str=="0" or str=="no") return false;
        throw("unknown boolean ");
        return 0;
    }

public:
    /// read the value for this parameter from a string -- specialization for bool

    /// @param[inout]	is	the input stream to read the value from (key is already read away)
    /// @param[inout]	param	the parameter into which the value is read
	/// @return			the input stream for chaining
	template<class Q = T>
	friend typename std::enable_if<std::is_same<Q,bool>::value, std::istream&>::type
	operator>>(std::istream& is, Parameter<T>& param) {
		std::string word;
		if (is >> word) param.set_user_defined_value(stringtobool(word));
		if (is.fail()) throw std::runtime_error("I/O error");
		return is;
	}

    /// read the value for this parameter from a string -- general implementation

    /// @param[inout]	is	the input stream to read the value from (key is already read away)
    /// @param[inout]	param	the parameter into which the value is read
	/// @return			the input stream for chaining
	template<class Q = T>
	friend typename std::enable_if<!std::is_same<Q,bool>::value, std::istream&>::type
	operator>>(std::istream& is, Parameter<T>& param) {
		T word;
		if (is >> word) param.set_user_defined_value(word);
		if (is.fail()) throw std::runtime_error("I/O error");
		return is;
	}

	/// user-friendly output of the parameter value

	/// @param[inout]	os	the output stream
	/// @param[in]		param the parameter to print
	/// @return			the output stream for chaining
	friend std::ostream& operator<<(std::ostream& os, const Parameter<T>& param) {
		os << param();
		return os;
	}
};


/// implementation of the type erasure idiom
class ParameterObject {

public:

	template <typename T>
	ParameterObject(const T& obj): object(std::make_shared<Model<T>>(std::move(obj))){}

	/// default constructor
	ParameterObject(): object(0){}

	/// shallow copy constructor
	ParameterObject(const ParameterObject& other) = default;

	/// assignment constructor
	ParameterObject& operator=(const ParameterObject& other) = default;

	/// given an input stream assign the parameter value to this

	/// @param[in]	key	assign only if key equals this' key
	/// @param[in]	is	input stream with the parameter data
	/// @return		true if key matches and assignment was successful, false otherwise
	bool parse(const std::string& key, std::istream& is) {
		bool success=false;
		if (key==object->get_key()) {
		    try {
		    	is >> *this;
		    	success=true;
		    } catch (...) {
		    	std::string errmsg="error in reading parameter line for key: "+key;
		    	throw std::runtime_error(errmsg);
		    }
		}
		return success;
	}

	/// return a reference to the parameter
	template<typename T>
	Parameter<T>& get_parameter() {
		if (std::type_index(typeid(T)) == std::type_index(object->value_type_info())) {
			return static_cast<Model<Parameter<T> >&>(*object).Mobject;
		} else {
			MADNESS_EXCEPTION("faulty type in set_user_defined_value ",1);
		}
	}

	/// return a reference to the parameter
	template<typename T>
	const Parameter<T>& get_parameter() const {
		if (std::type_index(typeid(T)) == std::type_index(object->value_type_info())) {
			return static_cast<Model<Parameter<T> >&>(*object).Mobject;
		} else {
			MADNESS_EXCEPTION("faulty type in set_user_defined_value ",1);
		}
	}

	/// set a user-defined value by an istream
	friend std::istream& operator>>(std::istream& is, ParameterObject& po) {
		po.object->set_value_by_istream(is);
		return is;
	}

	/// set a user-defined value by assignment
	template<typename T>
	ParameterObject& operator=(const T u) {
		set_user_defined_value(u);
		return *this;
	}

	/// set a derived value
	template<typename T>
	void set_derived_value(const T d) {
		get_parameter<T>().set_derived_value(d);
	}

	/// return the key of this parameter
	std::string get_key() const {
		return object->get_key();
	}

	/// return the value of the parameter -- you need to know the type at compile time
	template<typename T>
	T get() const {
		return get_parameter<T>().get();
	}

	/// print the value of this parameter -- you don't need to know the type at compile time
	friend std::ostream& operator<<(std::ostream& os, const ParameterObject& po) {
		return po.object->get_ostream(os);
	}

private:

	template<typename T>
	void set_user_defined_value(const T u) {
		get_parameter<T>().set_user_defined_value(u);
	}

	struct Concept {
		virtual ~Concept() {}
		/// return the key
		virtual std::string get_key() const = 0;
		/// return the type of the parameter, e.g. Parameter<double>
		virtual const std::type_info& type_info() const = 0;
		/// return the type of the parameter value, e.g. double
		virtual const std::type_info& value_type_info() const = 0;
		/// set a parameter value using an input stream
		virtual void set_value_by_istream(std::istream& is) = 0;
		/// print the parameter value into an output stream
		virtual std::ostream& get_ostream(std::ostream& os) const = 0 ;

	};

	template< typename T >
	struct Model : Concept {
		/// constructor
		Model(const T& t) : Mobject(t) {}
		/// return the key
		std::string get_key() const override {return Mobject.get_key();}
		/// return the type of the parameter, e.g. Parameter<double>
		const std::type_info& type_info() const override { return typeid(T); }
		/// return the type of the parameter value, e.g. double
		const std::type_info& value_type_info() const override {return Mobject.value_type_info();}
		/// set a parameter value using an input stream
		void set_value_by_istream(std::istream& is) override {is >> Mobject;}
		/// print the parameter value into an output stream
		virtual std::ostream& get_ostream(std::ostream& os) const override {
			os << Mobject;
			return os;
		}

		/// the parameter object
		T Mobject;
	};

	std::shared_ptr<Concept> object;
};




/// base class for calculation parameters
class CalculationParametersBase {
public:

	typedef  std::map<int,ParameterObject> ParameterMap;

	/// ctor for testing
	CalculationParametersBase() {
	}

	/// destructor
	virtual ~CalculationParametersBase() {}

	/// print all parameters
	virtual void print(const ParameterMap params, const std::string header="",
			const std::string footer="") const {
		if (header.size()>0) madness::print(header);
		for (auto& p : params) print_twocolumn_centered(p.second);
		if (footer.size()>0) madness::print(footer);
	}


    /// read the parameters from file

	/// only world.rank()==0 reads the input file and broadcasts to all other nodes,
	/// so we don't need to serialize the ParameterMap
    virtual void read(World& world, const std::string filename, const std::string tag,
    		ParameterMap& param) {

    	std::string filecontents, line;
    	if (world.rank()==0) {
    		std::ifstream f(filename.c_str());
    		while (std::getline(f,line)) filecontents+=line+"\n";
    	}

    	// broadcast the input file to all nodes
    	world.gop.broadcast_serializable(filecontents, 0);
    	read_internal(world, filecontents,tag,param);
    }

    /// read the stream, starting from tag

    /// only parameters that are defined in the constructor will be processed,
    /// all others will be discarded.
    virtual void read_internal(World& world, std::string& filecontents, std::string tag, ParameterMap& params) {
    	std::stringstream f(filecontents);
		position_stream(f, tag);
		std::string line, word;

		// read input lines
		while (std::getline(f,line)) {

			std::stringstream sline(line);
			sline >> word;

			// skip comment line
			if (word[0]=='#') continue;
			if (word=="end") break;

			bool found=false;
			for (auto& p : params) {
				try {
					found=p.second.parse(word,sline);
				} catch (std::exception& e) {
					// the key was found but the parsing failed
					if (world.rank()==0) madness::print(e.what());
					found=true;
				}
				if (found) break;
			}

			if (not found) {
				if (world.rank()==0) madness::print("\nignoring unknown keyword in group ",tag,":",word);
			}
		}
    };

    /// pretty print the parameters with optional comment

    /// @param[in]	t	the parameter
    /// @param[in]	s2	optional comment line
    template<typename T>
    void print_twocolumn_centered(const T t, const std::string s2="") const {
    	std::string key=t.get_key();
    	std::ostringstream os1, os2;
    	os1 << std::scientific << std::boolalpha << t;
    	int len=20-key.size();
    	for (int i=0; i<len; ++i) os2 << " ";
//    	os2 << key <<  "  "<<std::setw(15) << os1.str();
    	os2 << key <<  "  " << os1.str();
    	if (s2.size()>0) os2 << "  #  " << s2;
    	madness::print(os2.str());
    }

protected:

    /// initialize a key/value pair with an integerized enum and a parameter

    /// see mp2.h for an example
    /// @param[in]	i	enum as integer
    /// @param[in]	p	templated parameter
    template<typename T>
    static std::pair<int,ParameterObject> init(int i, Parameter<T> p) {
    	return std::make_pair(i,ParameterObject(p));
    }

};


} // namespace madness


#endif /* SRC_APPS_CHEM_CALCULATIONPARAMETERSBASE_H_ */
