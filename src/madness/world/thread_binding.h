/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

/**
 \file thread_binding.h
 \brief Startup check that the MPI launcher did not pin all threads of a rank
        onto too few hardware threads.
 \ingroup world

 Most MPI launchers bind a rank to a single core by default (e.g. Open MPI's
 `--bind-to core` when there are more cores than ranks on a node).  A
 multi-threaded MADNESS rank inherits that affinity mask for *every* thread it
 spawns, so the whole rank time-shares one hardware thread.  The job still runs
 and produces correct results, but one to two orders of magnitude slower, which
 is easy to miss.  \c check_thread_binding() detects this at startup and throws
 a \c ThreadBindingException carrying an actionable error message.
*/

#ifndef MADNESS_WORLD_THREAD_BINDING_H__INCLUDED
#define MADNESS_WORLD_THREAD_BINDING_H__INCLUDED

#include <madness/world/madness_exception.h>
#include <madness/world/safempi.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace madness {

    /// Description of the CPU affinity of the calling process.

    /// Filled in by get_thread_binding_info(); all information is obtained from
    /// the OS (sched_getaffinity(2) and sysfs), no external library is needed.
    struct ThreadBindingInfo {
        /// True if the platform exposes CPU affinity (Linux).  On macOS, or if
        /// the affinity query failed, this is false and all other members are
        /// meaningless.
        bool available = false;
        /// Logical (hardware thread) ids this process is allowed to run on.
        std::vector<int> cpus;
        /// Distinct NUMA nodes spanned by \c cpus, sorted; empty if the NUMA
        /// topology could not be determined.
        std::vector<int> numa_nodes;
        /// Number of logical CPUs configured on this host (not just allowed).
        int ncpu_configured = 0;

        /// Number of hardware threads this process may run on.
        std::size_t ncpu_allowed() const { return cpus.size(); }
    };

    namespace detail {
        /// Owns the (long) diagnostic text of a ThreadBindingException.

        /// \c MadnessException only stores a <tt>const char*</tt>, so the text has to
        /// live somewhere else.  It is held through a \c shared_ptr so that the copies
        /// C++ is allowed to make while the exception is in flight keep pointing at
        /// the same, still-alive buffer.
        struct ThreadBindingMessage {
            std::shared_ptr<const std::string> text;
            explicit ThreadBindingMessage(std::string t)
                : text(std::make_shared<const std::string>(std::move(t))) {}
        };
    } // namespace detail

    /// Thrown by check_thread_binding() when the CPU binding cannot support the thread pool.

    /// Derives from MadnessException, so existing `catch(const MadnessException&)`
    /// handlers pick it up; \c what() returns the full multi-line diagnostic.
    /// \note The private base must stay first in the base list -- bases are
    ///       initialized in declaration order, and MadnessException's constructor
    ///       needs the text to exist already.
    class ThreadBindingException : private detail::ThreadBindingMessage,
                                   public MadnessException {
    public:
        ThreadBindingException(std::string text, const int line, const char* function,
                               const char* file)
            : detail::ThreadBindingMessage(std::move(text))
            , MadnessException(detail::ThreadBindingMessage::text->c_str(),
                               "CPU binding over-subscribes the hardware threads", 1,
                               line, function, file) {}
    };

    /// Query the CPU affinity mask and NUMA placement of the calling process.
    ThreadBindingInfo get_thread_binding_info();

    /// Render a list of cpu ids in compact form, e.g. {0,1,2,3,8} -> "0-3,8".
    std::string cpu_list_to_string(const std::vector<int>& cpus);

    /// Throw if the CPU binding cannot accommodate the thread pool.

    /// Two conditions are checked, both of which make the calculation
    /// pathologically slow:
    ///  - this rank runs more threads than there are hardware threads in its
    ///    affinity mask (the classic "`--bind-to core` clogs the hwthreads"
    ///    failure), and
    ///  - all ranks sharing this host together request more threads than there
    ///    are hardware threads available to them (over-subscription that
    ///    `--bind-to none` does not cure).
    ///
    /// In addition a non-fatal warning is issued if a rank's affinity mask
    /// straddles several NUMA domains while several ranks share the host, since
    /// that costs memory locality.
    ///
    /// The check is collective over \c comm and must be called after the thread
    /// pool has been created.  It is a no-op on platforms without an affinity
    /// API, and can be disabled with `MAD_CHECK_BINDING=OFF`.
    ///
    /// The verdict is agreed on collectively, so either every rank throws or none
    /// does -- no rank is left waiting in a collective for a rank that has bailed
    /// out.  The offending rank has already written the full diagnostic to both
    /// \c stdout and \c stderr by the time it throws, and carries that same text
    /// in \c what(); the remaining ranks throw a one-line message pointing at it,
    /// so the diagnostic is not repeated once per rank.
    ///
    /// \param[in] comm the communicator MADNESS was initialized with
    /// \param[in] nthread_app number of application threads per rank, i.e. the
    ///            thread pool plus the main thread (the MPI communication thread
    ///            is deliberately not counted, see the implementation)
    /// \throws ThreadBindingException if the binding over-subscribes the hardware
    void check_thread_binding(const SafeMPI::Intracomm& comm, int nthread_app);

} // namespace madness

#endif // MADNESS_WORLD_THREAD_BINDING_H__INCLUDED
