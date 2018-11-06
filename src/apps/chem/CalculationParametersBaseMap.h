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
	Parameter(const std::string k, const T d) : key(k), value(d), default_value(d) {
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

    std::string getName() const {                           // (8)
    	std::string name="parameter with type "+std::string(typeid(T).name());
        return name;;
    }

    /// return the type of the parameter value
    const std::type_info& value_type_info() const {return typeid(value);}

private:
    std::string key;			///< the key in the input file
	T value;					///< the actual value
	T default_value;		///< the default value
	T user_defined_value = T();	///< the value if defined by the user
	T derived_value = T();		///< the value if derived from another parameter
	T null = T();			///< for comparison

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
		if (is >> word) param.set_user_defined_value(stringtobool(word));
		if (is.fail()) throw std::runtime_error("I/O error");
		return is;
	}


	template<class Q = T>
	friend typename std::enable_if<!std::is_same<Q,bool>::value, std::istream&>::type
	operator>>(std::istream& is, Parameter<T>& param) {
		T word;
		if (is >> word) param.set_user_defined_value(word);
		if (is.fail()) throw std::runtime_error("I/O error");
		return is;
	}

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
		    	std::cout << "error in reading parameter line for key: " << key<<std::endl;
		    	throw std::runtime_error("I/O error");
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
			std::cout << "faulty type in set_user_defined_value " << std::endl;
			throw;
		}
	}

	/// return a reference to the parameter
	template<typename T>
	const Parameter<T>& get_parameter() const {
		if (std::type_index(typeid(T)) == std::type_index(object->value_type_info())) {
			return static_cast<Model<Parameter<T> >&>(*object).Mobject;
		} else {
			std::cout << "faulty type in set_user_defined_value " << std::endl;
			throw;
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

	std::string getName() const {
		return object->getName();
	}

	std::string get_key() const {
		return object->get_key();
	}

	/// return the value of the parameter
	template<typename T>
	T get() const {
		return get_parameter<T>().get();
	}

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
		virtual std::string getName() const = 0;
		virtual std::string get_key() const = 0;
		virtual const std::type_info& type_info() const = 0;
		virtual const std::type_info& value_type_info() const = 0;
		virtual void set_value_by_istream(std::istream& is) = 0;
		virtual std::ostream& get_ostream(std::ostream& os) const = 0 ;

	};

	template< typename T >
	struct Model : Concept {
		Model(const T& t) : Mobject(t) {}
		std::string getName() const override {
			return Mobject.getName();
		}
		std::string get_key() const override {
			return Mobject.get_key();
		}
		const std::type_info& type_info() const override { return typeid(T); }
		const std::type_info& value_type_info() const override {return Mobject.value_type_info();}

		void set_value_by_istream(std::istream& is) override {
			is >> Mobject;
		}

		virtual std::ostream& get_ostream(std::ostream& os) const override {
			os << Mobject;
			return os;
		}

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
		/*
		add the default parameters in the derived constructor, e.g.

		params[ParameterKey(0,"econv")]=ParameterObject(Parameter<double>(1.e-5));
		params[ParameterKey(1,"dconv")]=ParameterObject(Parameter<double>(1.e-4));

		read("input","mp3");
		print();

		*/
	}

	/// destructor
	virtual ~CalculationParametersBase() {}

	/// return the parameter map of the derived class
	virtual ParameterMap& get_ParameterMap() = 0;

	/// print all parameters
	virtual void print(const ParameterMap params, const std::string header="",
			const std::string footer="") const {
		if (header.size()>0) madness::print(header);
		for (auto& p : params) print_twocolumn_centered(p.second.get_key(), p.second);
		if (footer.size()>0) madness::print(footer);
	}


    /// read the parameters from file -- no need to reimplement this
    virtual void read(World& world, const std::string filename, const std::string tag,
    		ParameterMap& param) {
    	std::ifstream f(filename.c_str());
    	read(f,tag,param);
    }

    /// read the stream, starting from tag

    /// only parameters that are defined in the constructor will be processed,
    /// all others will be discarded.
    virtual void read(std::istream& f, std::string tag, ParameterMap& params) {
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
//				std::cout << " parsing " << sline.str() <<"|---" <<  std::endl;
				try {
					found=p.second.parse(word,sline);
				} catch (...) {
					// the key was found but the parsing failed
					std::cout << "could not read properly from line\n" << line << std::endl;
					found=true;
				}
				if (found) break;
			}

			if (not found) {
				madness::print("\nignoring unknown keyword in group ",tag,":",word);
			}
		}
    };

    /// pretty print the parameters with optional comment
    template<typename T>
    void print_twocolumn_centered(const std::string s, const T t, const std::string s2="") const {
    	int len=20-s.size();
    	for (int i=0; i<len; ++i) std::cout << " ";
    	std::ostringstream os;
    	os << std::scientific << std::boolalpha << t;
    	std::cout << s <<  "  "<<std::setw(15) << os.str() << "  " << s2<<std::endl;
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
