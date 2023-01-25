/*
 * gitversion.cc
 *
 *  Created on: Jan 21, 2020
 *      Author: fbischoff
 */


#include <madness/madness_config.h>
#include <sstream>
#include <string>
#include <iostream>

namespace madness {
    namespace info {

        const char* version() {
#ifdef MADNESS_VERSION // from madness_config.h
          return MADNESS_VERSION;
#else
          return "unavailable";
#endif
        }


        const char* git_commit() {
            return MADNESS_GIT_REVISION;
        }

        const char* git_source_description() {
          return MADNESS_GIT_DESCRIPTION;
        }

        const char* build_time() {
            return __TIME__;
        }

        const char* build_date() {
            return __DATE__;
        }

        std::string print_revision_information() {
        	std::stringstream ss;
                ss << "    git source description ...     " << info::git_source_description() << std::endl;
        	const std::string time(build_time());
        	const std::string date(build_date());
        	ss << "    built date/time ...     " << date << "/" << time << std::endl;
        	return ss.str();
        }

    } // namespace info
} // namespace madness

