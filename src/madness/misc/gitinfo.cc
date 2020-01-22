/*
 * gitversion.cc
 *
 *  Created on: Jan 21, 2020
 *      Author: fbischoff
 */


#include <madness/madness_config.h>
#include <madness/misc/gitversion.h>
#include <sstream>
#include <string>
#include <iostream>

namespace madness {
    namespace info {

        const char* mad_git_commit() {
            return MADNESS_GITREVISION;
        }

        const char* build_time() {
            return __TIME__;
        }

        const char* build_date() {
            return __DATE__;
        }

        std::string print_revision_information() {
        	std::stringstream ss;
#ifdef MADNESS_REVISION
        	const  char* gitrev =  MADNESS_REVISION;
        	const std::string gitrevision(gitrev);
        	ss << "    git revision at configure time ... " << gitrevision  << std::endl;
#endif
        	const std::string gitrevision1(info::mad_git_commit());
        	ss << "    git revision at build time ...     " << gitrevision1 << std::endl;
        	const std::string time(build_time());
        	const std::string date(build_date());
        	ss << "    build time at ...                  " << time << " on " << date << std::endl;
        	return ss.str();
        }

    } // namespace info
} // namespace madness

