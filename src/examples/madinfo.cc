/*
 * cheminfo_prog.cc
 *
 *  Created on: Jan 20, 2020
 *      Author: fbischoff
 */

#include <string>
#include <iostream>

#include <madness/madness_config.h>
#include "madness/misc/gitinfo.h"
int main() {

	std::string line1(MADNESS_REVISION);
	std::cout << "The git revision at configure time is: " << line1 << std::endl;
	std::string line(madness::info::mad_git_commit());
	std::cout << "The git revision at compile time is:   " << line << std::endl;
	std::string line2(madness::info::build_time());
	std::string line3(madness::info::build_date());
	std::cout << "MADNESS built on:                      " << line3  << " at " << line2 << std::endl;

	return 0;

}
