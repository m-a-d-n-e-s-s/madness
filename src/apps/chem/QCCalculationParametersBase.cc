/*
 * QCCalculationParametersBase.cpp
 *
 *  Created on: 27 Jun 2019
 *      Author: fbischoff
 */

#include <chem/QCCalculationParametersBase.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include<madness/world/worldgop.h>



namespace madness {




/// print all parameters
void QCCalculationParametersBase::print(const std::string header,
		const std::string footer) const {

	std::string body=print_to_string();
	if (header.size()>0) madness::print(header);
	::madness::print(body);
	if (footer.size()>0) madness::print(footer);
}

std::string QCCalculationParametersBase::print_to_string(bool non_defaults_only) const {

	// sort parameters according to increasing print_order
	typedef std::tuple<int,std::string,QCParameter> keyvalT;
	std::list<keyvalT> list;
	for (auto& p : parameters) list.push_back(keyvalT(p.second.get_print_order(),p.first,p.second));
	list.sort([](const keyvalT& first, const keyvalT& second) {return std::get<0>(first) < std::get<0>(second);});

	std::stringstream ss;
	for (auto& p : list) {
		const QCParameter& param=std::get<2>(p);
		if (non_defaults_only and (param.precedence==QCParameter::def)) continue;
		ss << param.print_line(std::get<1>(p));
		ss << std::endl;
	}
	return ss.str();
}

/// read the parameters from file and broadcast
void QCCalculationParametersBase::read(World& world, const std::string filename, const std::string tag) {

	std::string filecontents, line;
//	if (world.rank()==0) {
		std::ifstream f(filename.c_str());
		while (std::getline(f,line)) filecontents+=line+"\n";
//	}
	read_internal(world, filecontents,tag);
//	world.gop.broadcast_serializable(*this, 0);

}

/// read the stream, starting from tag

/// only parameters that are defined in the constructor will be processed,
/// all others will be discarded.
void QCCalculationParametersBase::read_internal(World& world, std::string& filecontents, std::string tag) {
	std::stringstream f(filecontents);
	position_stream_to_word(f, tag);
	std::string line, key,value;

	// read input lines
	while (std::getline(f,line)) {

		// all in lower case
		std::transform(line.begin(), line.end(), line.begin(), ::tolower);

		// remove comments from line
		std::size_t last = line.find_first_of('#');
		line=line.substr(0,last);

		std::stringstream sline(line);

		// sline might be empty by now
		if (not (sline >> key)) continue;
		if (print_debug) ::madness::print("reading key ",key);

		// skip comment line
		if (key[0]=='#') continue;
		if (key=="end") break;

		// check if key exists in the initialized parameter list
		if (not (parameter_exists(key))) {
			if (world.rank()==0) madness::print("ignoring unknown parameter in input file: ",key);
			continue;
		}

		std::string word,line1;
		while (sline >> word) {line1+=word+" ";}
		// trim result
		last = line1.find_last_not_of(' ');
		line1=line1.substr(0, last+1);


		// check which type to expect from the given key
		// this look clumsy, maybe there are more elegant solutions?
		bool success=false;
		try {
			success=try_setting_user_defined_value<double>(key,line1) or success;
			success=try_setting_user_defined_value<int>(key,line1) or success;
			success=try_setting_user_defined_value<unsigned int>(key,line1) or success;
			success=try_setting_user_defined_value<long>(key,line1) or success;
			success=try_setting_user_defined_value<std::size_t>(key,line1) or success;
			success=try_setting_user_defined_value<bool>(key,line1) or success;
			success=try_setting_user_defined_value<std::string>(key,line1) or success;
			success=try_setting_user_defined_value<std::vector<double> >(key,line1) or success;
			success=try_setting_user_defined_value<std::vector<int> >(key,line1) or success;
            success=try_setting_user_defined_value<std::vector<std::size_t> >(key,line1) or success;
			success=try_setting_user_defined_value<std::vector<std::string> >(key,line1) or success;
			success=try_setting_user_defined_value<std::pair<std::string,double> >(key,line1) or success;

		} catch (std::runtime_error& e) {
			std::string errmsg="found an error for key >> "+key+" << \n" +e.what();
			throw std::runtime_error(errmsg);
		}
		if (not success) {
			madness::print("\n\ncould not assign the input parameter for key ",key);
			std::string requested_type=get_parameter(key).get_type();
			madness::print("\trequested type: ",requested_type,"\n");
			madness::print("add the corresponding type to QCCalculationsParametersBase.h\n\n");
			MADNESS_EXCEPTION("add the corresponding type to QCCalculationsParametersBase.h",1);
		}
	}
};


} /* namespace madness */
