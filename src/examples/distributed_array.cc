/*
 * cheminfo_prog.cc
 *
 *  Created on: Jan 20, 2020
 *      Author: fbischoff
 */

#include <string>
#include <iostream>

#include "madness/misc/info.h"
#include <madness/madness_config.h>
int main() {

        std::cout << "The git source tree description at compile time is:   " << madness::info::git_source_description() << std::endl;
	std::string line2(madness::info::build_time());
	std::string line3(madness::info::build_date());
	std::cout << "MADNESS built on:                      " << line3  << " at " << line2 << std::endl;

	return 0;

}
