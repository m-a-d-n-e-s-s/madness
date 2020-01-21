/*
 * gitversion.cc
 *
 *  Created on: Jan 21, 2020
 *      Author: fbischoff
 */


#include <madness/madness_config.h>
#include <madness/misc/gitversion.h>

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

    } // namespace info
} // namespace madness

