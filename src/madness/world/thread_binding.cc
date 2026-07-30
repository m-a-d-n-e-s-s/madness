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
 \file thread_binding.cc
 \brief Implementation of the startup CPU-binding sanity check.
 \ingroup world
*/

#include <madness/madness_config.h>
#include <madness/world/thread_binding.h>
#include <madness/world/worldinit.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#include <unistd.h>

#ifndef ON_A_MAC
#  include <sched.h>
#endif

namespace madness {

    namespace {

        /// Parse a Linux cpu list such as "0-3,8,10-11" into individual ids.
        std::vector<int> parse_cpu_list(const std::string& s) {
            std::vector<int> cpus;
            std::stringstream ss(s);
            std::string item;
            while (std::getline(ss, item, ',')) {
                if (item.empty()) continue;
                const std::size_t dash = item.find('-');
                try {
                    if (dash == std::string::npos) {
                        cpus.push_back(std::stoi(item));
                    }
                    else {
                        const int lo = std::stoi(item.substr(0, dash));
                        const int hi = std::stoi(item.substr(dash + 1));
                        for (int i = lo; i <= hi; ++i) cpus.push_back(i);
                    }
                }
                catch (...) {  // malformed sysfs content -- just ignore it
                    continue;
                }
            }
            return cpus;
        }

        /// Read a one-line file; returns false if it does not exist.
        bool read_line_from_file(const std::string& path, std::string& line) {
            std::ifstream f(path.c_str());
            if (!f) return false;
            if (!std::getline(f, line)) return false;
            // strip trailing whitespace
            while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
                line.pop_back();
            return true;
        }

        /// Determine which NUMA nodes the given cpus belong to.

        /// Uses sysfs only, so no dependency on libnuma; returns an empty vector
        /// if the topology is not exposed (non-Linux, containers without sysfs).
        std::vector<int> numa_nodes_of(const std::vector<int>& cpus) {
            std::vector<int> nodes;
            std::string online;
            if (!read_line_from_file("/sys/devices/system/node/online", online)) return nodes;
            for (int node : parse_cpu_list(online)) {
                std::string cpulist;
                std::ostringstream path;
                path << "/sys/devices/system/node/node" << node << "/cpulist";
                if (!read_line_from_file(path.str(), cpulist)) continue;
                const std::vector<int> node_cpus = parse_cpu_list(cpulist);
                for (int cpu : cpus) {
                    if (std::find(node_cpus.begin(), node_cpus.end(), cpu) != node_cpus.end()) {
                        nodes.push_back(node);
                        break;
                    }
                }
            }
            std::sort(nodes.begin(), nodes.end());
            return nodes;
        }

        /// Number of NUMA domains physically present (online) on this host.

        /// Note this is the *host* topology, not what the affinity mask allows;
        /// returns 0 if sysfs does not expose it.
        int numa_nodes_on_host() {
            std::string online;
            if (!read_line_from_file("/sys/devices/system/node/online", online)) return 0;
            return static_cast<int>(parse_cpu_list(online).size());
        }

        std::string hostname() {
            char buf[256];
            if (gethostname(buf, sizeof(buf)) != 0) return "<unknown>";
            buf[sizeof(buf) - 1] = '\0';
            return std::string(buf);
        }

        /// Recipe block appended to every binding complaint.

        /// \param[in] nthread_app threads per rank that the user asked for
        /// \param[in] nnuma number of NUMA domains on this host (0 if unknown)
        std::string binding_advice(const int nthread_app, const int nnuma) {
            std::ostringstream o;
            o << "  HOW TO FIX THIS\n"
              << "  Give every rank as many hardware threads as it runs threads.  Prefer a\n"
              << "  binding that keeps a rank inside a single NUMA domain: switching binding\n"
              << "  off entirely (--bind-to none) removes the over-subscription, but then the\n"
              << "  OS is free to migrate threads across sockets, which typically costs another\n"
              << "  20-50% on the memory-bound parts of MADNESS.\n"
              << "\n";
            if (nnuma > 0) {
                o << "  This host has " << nnuma << " NUMA domain(s), so a good starting layout is one\n"
                  << "  rank per NUMA domain -- " << nnuma << " rank(s) per host with "
                  << nthread_app << " thread(s) each.\n\n";
            }
            o << "  Open MPI -- NUMA-aware, recommended (one rank per NUMA domain):\n"
              << "      mpirun --map-by ppr:1:numa:PE=" << nthread_app << " --bind-to core \\\n"
              << "             -x MAD_NUM_THREADS=" << nthread_app << " ...\n"
              << "  Open MPI -- several ranks per NUMA domain (N ranks of " << nthread_app << " threads):\n"
              << "      mpirun --map-by ppr:N:numa:PE=" << nthread_app << " --bind-to core \\\n"
              << "             -x MAD_NUM_THREADS=" << nthread_app << " ...\n"
              << "  Open MPI -- quick workaround, no NUMA locality:\n"
              << "      mpirun --bind-to none -x MAD_NUM_THREADS=" << nthread_app << " ...\n"
              << "  Slurm:\n"
              << "      srun --cpus-per-task=" << nthread_app << " --cpu-bind=cores ...\n"
              << "      export MAD_NUM_THREADS=$SLURM_CPUS_PER_TASK\n"
              << "      (--hint=nomultithread keeps a rank off the sibling hyperthreads)\n"
              << "  MPICH:\n"
              << "      mpiexec -bind-to numa ...        or   mpiexec -bind-to none ...\n"
              << "  Intel MPI / MPICH derivatives:\n"
              << "      export I_MPI_PIN_DOMAIN=numa     or   export I_MPI_PIN=0\n"
              << "\n"
              << "  Alternatively reduce the thread count to match the binding you have:\n"
              << "      export MAD_NUM_THREADS=<hardware threads available per rank>\n"
              << "\n"
              << "  To run anyway (deliberate over-subscription, e.g. in a small test):\n"
              << "      export MAD_CHECK_BINDING=OFF\n";
            return o.str();
        }

        /// Emit to both stdout and stderr.

        /// Batch schedulers routinely capture only one of the two streams, and a
        /// fatal start-up diagnostic must not be the thing that gets lost.
        void print_to_both(const std::string& msg) {
            std::cout << msg << std::flush;
            std::cerr << msg << std::flush;
        }

        bool env_says_off(const char* value) {
            if (!value) return false;
            std::string v(value);
            for (char& c : v) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            return v == "off" || v == "0" || v == "false" || v == "no" || v == "none";
        }

    } // anonymous namespace


    std::string cpu_list_to_string(const std::vector<int>& cpus) {
        if (cpus.empty()) return "<none>";
        std::ostringstream o;
        std::size_t i = 0;
        bool first = true;
        while (i < cpus.size()) {
            std::size_t j = i;
            while (j + 1 < cpus.size() && cpus[j + 1] == cpus[j] + 1) ++j;
            if (!first) o << ",";
            first = false;
            if (j == i)          o << cpus[i];
            else if (j == i + 1) o << cpus[i] << "," << cpus[j];
            else                 o << cpus[i] << "-" << cpus[j];
            i = j + 1;
        }
        return o.str();
    }


    ThreadBindingInfo get_thread_binding_info() {
        ThreadBindingInfo info;
#ifdef ON_A_MAC
        // macOS exposes no CPU affinity API; nothing to check.
        return info;
#else
        int nconf = static_cast<int>(sysconf(_SC_NPROCESSORS_CONF));
        if (nconf <= 0) return info;
        info.ncpu_configured = nconf;

        // N.B. CPU_ALLOC rather than a plain cpu_set_t: the latter is fixed at
        // CPU_SETSIZE (1024) and silently truncates on very large nodes.
        cpu_set_t* set = CPU_ALLOC(nconf);
        if (!set) return info;
        const std::size_t setsize = CPU_ALLOC_SIZE(nconf);
        CPU_ZERO_S(setsize, set);
        if (sched_getaffinity(0, setsize, set) != 0) {
            CPU_FREE(set);
            return info;
        }
        for (int i = 0; i < nconf; ++i)
            if (CPU_ISSET_S(i, setsize, set)) info.cpus.push_back(i);
        CPU_FREE(set);

        info.available = !info.cpus.empty();
        if (info.available) info.numa_nodes = numa_nodes_of(info.cpus);
        return info;
#endif
    }


    void check_thread_binding(const SafeMPI::Intracomm& comm, const int nthread_app) {
        if (env_says_off(getenv("MAD_CHECK_BINDING"))) return;

        const ThreadBindingInfo info = get_thread_binding_info();
        if (!info.available) return;  // no affinity API, or query failed -- nothing to say

        const int rank = comm.Get_rank();
        const int nproc = comm.Get_size();
        const int ncpu_rank = static_cast<int>(info.ncpu_allowed());

        // ---- aggregate over the ranks that share this host --------------------
        // Blocking SafeMPI collectives rather than World::gop: this runs inside
        // madness::initialize() before RMI::begin(), so there is no task scheduler
        // yet that a blocking call could park a thread outside of.
        SafeMPI::Intracomm node = comm.Split_type(SafeMPI::Intracomm::SHARED_SPLIT_TYPE);
        const int nrank_on_host = node.Get_size();

        int demand_host = 0;
        node.Allreduce(&nthread_app, &demand_host, 1, MPI_INT, MPI_SUM);

        // union of the affinity masks of all ranks on this host; _SC_NPROCESSORS_CONF
        // is identical for all of them, but reduce it anyway so the buffer sizes agree
        int nconf_host = 0;
        node.Allreduce(&info.ncpu_configured, &nconf_host, 1, MPI_INT, MPI_MAX);
        std::vector<unsigned char> mine(nconf_host, 0), anyone(nconf_host, 0);
        for (int cpu : info.cpus)
            if (cpu < nconf_host) mine[cpu] = 1;
        node.Allreduce(mine.data(), anyone.data(), nconf_host, MPI_UNSIGNED_CHAR, MPI_MAX);

        std::vector<int> host_cpus;
        for (int i = 0; i < nconf_host; ++i)
            if (anyone[i]) host_cpus.push_back(i);
        const int ncpu_host = static_cast<int>(host_cpus.size());
        const int nnuma_host = numa_nodes_on_host();

        // ---- the two fatal conditions ----------------------------------------
        // The MPI communication thread is deliberately not counted: it spends
        // essentially all of its time blocked in MPI, and counting it would make
        // the natural choice MAD_NUM_THREADS == cores-per-rank trip the check.
        const bool rank_oversubscribed = nthread_app > ncpu_rank;
        const bool host_oversubscribed = demand_host > ncpu_host;
        const bool failed = rank_oversubscribed || host_oversubscribed;

        // let the lowest-numbered offending rank do the talking
        int mine_rank = failed ? rank : std::numeric_limits<int>::max();
        int culprit = std::numeric_limits<int>::max();
        comm.Allreduce(&mine_rank, &culprit, 1, MPI_INT, MPI_MIN);

        if (culprit != std::numeric_limits<int>::max()) {
            // Every rank agrees on `culprit`, so every rank throws -- nobody is left
            // waiting in a later collective for a rank that has already bailed out.
            std::string message;
            if (rank == culprit) {
                std::ostringstream o;
                o << "\n"
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                  << "!! MADNESS FATAL ERROR: the CPU binding of this job over-subscribes the    !!\n"
                  << "!! hardware threads -- the calculation would run, but crippled.           !!\n"
                  << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
                  << "\n"
                  << "  host                                : " << hostname() << "\n"
                  << "  MPI rank                            : " << rank << " of " << nproc << "\n"
                  << "  MPI ranks on this host              : " << nrank_on_host << "\n"
                  << "  threads per rank (pool + main)      : " << nthread_app
                  << "   [MAD_NUM_THREADS]\n"
                  << "  hw threads this rank may use        : " << ncpu_rank
                  << "   [cpus " << cpu_list_to_string(info.cpus) << "]\n"
                  << "  hw threads usable on this host      : " << ncpu_host
                  << "   [cpus " << cpu_list_to_string(host_cpus) << "]\n"
                  << "  hw threads present on this host     : " << nconf_host << "\n"
                  << "  threads requested on this host      : " << demand_host << "\n";
                if (nnuma_host > 0)
                    o << "  NUMA domains on this host           : " << nnuma_host << "\n";
                if (!info.numa_nodes.empty())
                    o << "  NUMA domains in this rank's mask    : " << info.numa_nodes.size()
                      << "   [node(s) " << cpu_list_to_string(info.numa_nodes) << "]\n";
                o << "\n";

                if (rank_oversubscribed) {
                    o << "  WHAT WENT WRONG\n"
                      << "  This rank is pinned to " << ncpu_rank << " hardware thread(s) but runs "
                      << nthread_app << " thread(s), so\n"
                      << "  all of its threads time-share the same core(s).  This is almost always\n"
                      << "  the MPI launcher binding each rank to a single core by default (Open MPI\n"
                      << "  does this whenever a node holds more cores than ranks).  Expect a slowdown\n"
                      << "  of one to two orders of magnitude, not a crash -- which is exactly why\n"
                      << "  MADNESS refuses to continue rather than let you burn the allocation.\n\n";
                }
                else {
                    o << "  WHAT WENT WRONG\n"
                      << "  Each rank on its own has enough hardware threads, but the " << nrank_on_host
                      << " rank(s) on\n"
                      << "  this host together ask for " << demand_host << " thread(s) while only " << ncpu_host
                      << " hardware thread(s) are\n"
                      << "  available to them -- their CPU sets overlap, so the threads time-share cores.\n";
                    if (ncpu_host < nconf_host) {
                        o << "  Note the host has " << nconf_host << " hardware thread(s) in total; the launcher or the\n"
                          << "  batch system has confined this job to " << ncpu_host
                          << " of them.  Check your cpus-per-task\n"
                          << "  / cgroup allocation as well as the binding options below.\n\n";
                    }
                    else {
                        o << "  The host is simply too small for ranks x threads as requested: reduce\n"
                          << "  MAD_NUM_THREADS or the number of ranks per host.  Unbinding\n"
                          << "  (--bind-to none) will not help.\n\n";
                    }
                }

                o << binding_advice(nthread_app, nnuma_host);
                o << "\n";
                message = o.str();
                print_to_both(message);
            }
            else {
                // Short text for everyone else: the exception surfaces on every rank,
                // and repeating the full diagnostic once per rank would bury it.
                std::ostringstream o;
                o << "MADNESS: aborting because the CPU binding of rank " << culprit
                  << " cannot support " << nthread_app
                  << " thread(s) per rank -- see the diagnostic printed by that rank."
                  << "  Set MAD_CHECK_BINDING=OFF to run anyway.";
                message = o.str();
            }
            comm.Barrier();  // let the diagnostic reach the output before anyone unwinds
            exception_break(MADNESS_DISPLAY_EXCEPTION_BREAK_MESSAGE);
            throw ThreadBindingException(message, __LINE__, __FUNCTION__, __FILE__);
        }

        // ---- non-fatal: mask straddles several NUMA domains -------------------
        // Typical after "--bind-to none" on a multi-socket node: nothing is
        // over-subscribed, but threads of a rank wander across NUMA domains.
        if (!quiet() && nrank_on_host > 1 && info.numa_nodes.size() > 1) {
            int mine_warn = rank;
            int first_warn = 0;
            comm.Allreduce(&mine_warn, &first_warn, 1, MPI_INT, MPI_MIN);
            if (rank == first_warn) {
                std::ostringstream o;
                o << "!!MADNESS WARNING: rank " << rank << " on " << hostname() << " may run on "
                  << info.numa_nodes.size() << " NUMA domains\n"
                  << "!!MADNESS WARNING: (nodes " << cpu_list_to_string(info.numa_nodes) << ", cpus "
                  << cpu_list_to_string(info.cpus) << ") while " << nrank_on_host
                  << " ranks share this host.\n"
                  << "!!MADNESS WARNING: Threads will migrate between sockets and lose memory\n"
                  << "!!MADNESS WARNING: locality.  Consider confining each rank to one NUMA domain:\n"
                  << "!!MADNESS WARNING:     mpirun --map-by ppr:1:numa:PE=" << nthread_app
                  << " --bind-to core ...\n"
                  << "!!MADNESS WARNING:     srun --cpus-per-task=" << nthread_app
                  << " --cpu-bind=cores --hint=nomultithread ...\n"
                  << "!!MADNESS WARNING: Set MAD_CHECK_BINDING=OFF to silence this.\n";
                print_to_both(o.str());
            }
        }
    }

} // namespace madness
