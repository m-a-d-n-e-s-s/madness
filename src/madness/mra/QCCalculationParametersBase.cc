/*
 * QCCalculationParametersBase.cpp
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */

#include"QCCalculationParametersBase.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include"worldgop.h"



namespace madness {




/// print all parameters
void QCCalculationParametersBase::print(const std::string header,
		const std::string footer) const {

	std::string body=print_to_string({"all"});
	if (header.size()>0) madness::print(header);
	::madness::print(body);
	if (footer.size()>0) madness::print(footer);
}

std::string QCCalculationParametersBase::print_to_string(const std::list<std::string> precedence) const {

	// sort parameters according to increasing print_order
	typedef std::tuple<int,std::string,QCParameter> keyvalT;
	std::list<keyvalT> list;
	for (auto& p : parameters) list.push_back(keyvalT(p.second.get_print_order(),p.first,p.second));
	list.sort([](const keyvalT& first, const keyvalT& second) {return std::get<0>(first) < std::get<0>(second);});

	// check if we have to print only non-default parameters
	bool print_all = std::find(precedence.begin(), precedence.end(), "all") != precedence.end();
	bool print_default = std::find(precedence.begin(), precedence.end(), "default") != precedence.end();
	bool print_derived = std::find(precedence.begin(), precedence.end(), "derived") != precedence.end();
	bool print_defined = std::find(precedence.begin(), precedence.end(), "defined") != precedence.end();

	std::stringstream ss;
    int counter=0;
	for (auto& p : list) {
		const QCParameter& param=std::get<2>(p);
		bool skip=
			(not print_default and param.precedence == QCParameter::def) or
			(not print_derived and param.precedence == QCParameter::derived) or
			(not print_defined and param.precedence == QCParameter::defined);
		if (print_all) skip=false;
		if (skip) continue;
		if ((counter++)>0) ss << std::endl; // no newline at the very end
		ss << param.print_line(std::get<1>(p));
	}
	return ss.str();
}


bool QCCalculationParametersBase::file_exists(World& world, std::string filename) const {
    bool file_exists = true;
    if (world.rank() == 0) {
        std::ifstream ifs(filename);
        if (not ifs.is_open()) file_exists=false;
        ifs.close();
    }
    world.gop.broadcast_serializable(file_exists, 0);
    return file_exists;
}

/// read the parameters from file and broadcast
void QCCalculationParametersBase::read_input(World& world, const std::string filename, const std::string tag) {

	std::string filecontents, line;
    std::string errmsg;
	if (world.rank()==0) {
        try {
            std::ifstream f(filename.c_str());
            while (std::getline(f, line)) filecontents += line + "\n";
            read_internal(world, filecontents, tag);
        } catch (std::invalid_argument& e) {
            errmsg=e.what();
            throw;
        } catch (std::exception& e) {
            std::stringstream ss;
            ss << "could not read data group >>" << tag << "<< in file " << filename;
            errmsg=ss.str();
        }
    }
	world.gop.broadcast_serializable(*this, 0);
    if (errmsg.size()>0) throw std::runtime_error(errmsg);
}

/// read the parameters from the command line and broadcast

/// syntax is: qcprogram --mp2='maxiter 10; freeze 1' --dft:maxiter=20 --Xmpi:debug=true
/// the argument in quotes is the value of the parser keys
void QCCalculationParametersBase::read_commandline_options(World& world, const commandlineparser& parser,
														   const std::string tag) {
	if (not parser.key_exists(tag)) return;
	std::string value=parser.value(tag);
	// turn this into a fake input file
	if (world.rank()==0) {
        std::string q="'";
		std::replace_copy(value.begin(), value.end(), value.begin(), q.c_str()[0] , ' ');
		std::replace_copy(value.begin(), value.end(), value.begin(), ';', '\n');
		value=tag+"\n"+value+"\nend";
//        print("value",value);
		read_internal(world, value,tag);
	}
	world.gop.broadcast_serializable(*this, 0);
}



/// read the stream, starting from tag

/// only parameters that are defined in the constructor will be processed,
/// all others will be discarded.
void QCCalculationParametersBase::read_internal(World& world, std::string& filecontents, std::string tag) {
	std::stringstream f(filecontents);
	position_stream_to_word(f, tag, '#', true, true);
	std::string line, key,value;

	// read input lines
	while (std::getline(f,line)) {

		// all in lower case
		std::transform(line.begin(), line.end(), line.begin(), ::tolower);

		// remove comments from line
		std::size_t last = line.find_first_of('#');
		line=line.substr(0,last);
        std::replace_copy(line.begin(), line.end(), line.begin(),'=', ' ');

		std::stringstream sline(line);

		// sline might be empty by now
		if (not (sline >> key)) continue;
		if (print_debug) ::madness::print("reading key ",key);

		// skip comment line
		if (key[0]=='#') continue;
		if (key=="end") break;

		// check if key exists in the initialized parameter list
		if (not (parameter_exists(key))) {
            if (not ignore_unknown_keys) {
                if (world.rank()==0) {
                    ::madness::print("found unknown key: ",key);
                    ::madness::print("in datagroup:      ",tag);
                }
                throw std::runtime_error("input error");
            }
            if ((not ignore_unknown_keys_silently)
                and (world.rank()==0)) madness::print("ignoring unknown parameter in input file: ",key);
			continue;
		}

		std::string word,line1;
		while (sline >> word) {line1+=word+" ";}
		// trim result
		last = line1.find_last_not_of(' ');
		line1=line1.substr(0, last+1);


		// check which type to expect from the given key
		bool success=false;
		try {
			success = try_setting_any<all_parameter_types>(key, line1);
        } catch (std::invalid_argument& e) {
            throw;
		} catch (std::exception& e) {
			std::string errmsg="found an error for key >> "+key+" << \n" +e.what();
            throw std::runtime_error(errmsg);
		}
		if (not success) {
			madness::print("\n\ncould not assign the input parameter for key ",key);
			std::string requested_type=get_parameter(key).get_type();
			madness::print("\trequested type: ",requested_type,"\n");
			madness::print("add the corresponding type to QCCalculationsParametersBase.h\n\n");
			throw std::invalid_argument("add the corresponding type to QCCalculationsParametersBase.h");
		}
	}
};


bool operator==(const QCCalculationParametersBase& p1,
			   const QCCalculationParametersBase& p2) {
	if (p1.parameters.size() != p2.parameters.size()) return false;
	for (const auto& [key, value] : p1.parameters) {
		if (not p2.parameter_exists(key)) return false;
		if (not (value.get_value() == p2.get_parameter(key).get_value())) return false;
	}
	return true;
}


bool operator!=(const QCCalculationParametersBase& p1,
                const QCCalculationParametersBase& p2) {
    return !(p1 == p2);
}
} /* namespace madness */
