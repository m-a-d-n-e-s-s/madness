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
*/

/// \file test_rmi_huge.cc
/// \brief Messages too large for an MPI count of bytes, sent as the argument of a remote task.
///
/// A message longer than \c MAD_BUFFER_SIZE takes RMI's huge-message protocol, and one of
/// more than \c INT_MAX bytes must further be described to MPI by a derived datatype
/// ( \c SafeMPI::ByteCount ). Compiled as is, the test lowers the size at which the
/// derived datatype takes over, so that messages of a few MiB cover that path: one that
/// leaves a remainder after whole chunks, one that does not, one that stays a plain count,
/// and several at once. test_rmi_huge_large.cc compiles it with
/// \c MADNESS_TEST_RMI_HUGE_LARGE , which sends one message of more than \c INT_MAX bytes
/// at the default setting; that takes about 9 GiB, so it carries the "large" label.
///
/// Needs at least 2 ranks to send anything; with 1 the tasks are local.

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <madness/world/MADworld.h>
#include <madness/world/safempi.h>

using namespace madness;

namespace {

/// the byte at position \c i of the message with the given \c seed
inline unsigned char pattern(std::size_t i, std::uint64_t seed) {
    std::uint64_t x = (i + seed) * 0x9E3779B97F4A7C15ull;
    return static_cast<unsigned char>(x >> 56);
}

std::vector<unsigned char> make_message(std::size_t n, std::uint64_t seed) {
    std::vector<unsigned char> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = pattern(i, seed);
    return v;
}

/// runs on the receiver: the number of bytes that differ from the pattern, or
/// \c size_t(-1) if the size is wrong
std::size_t count_errors(const std::vector<unsigned char>& v, std::size_t n, std::uint64_t seed) {
    if (v.size() != n) return std::size_t(-1);
    std::size_t nerr = 0;
    for (std::size_t i = 0; i < n; ++i) nerr += (v[i] != pattern(i, seed));
    return nerr;
}

/// sends messages of the given sizes to the next rank, all pending at once
bool send_and_check(World& world, const std::vector<std::size_t>& sizes, const char* what) {
    const ProcessID dest = (world.rank() + 1) % world.size();
    std::vector<Future<std::size_t>> results;
    for (std::size_t k = 0; k < sizes.size(); ++k) {
        const std::uint64_t seed = 1000 * world.rank() + k;
        results.push_back(world.taskq.add(dest, &count_errors, make_message(sizes[k], seed), sizes[k], seed));
    }
    world.gop.fence();
    bool ok = true;
    for (std::size_t k = 0; k < sizes.size(); ++k) {
        const std::size_t nerr = results[k].get();
        if (nerr != 0) {
            ok = false;
            std::cout << "rank " << world.rank() << ": " << what << ": message of " << sizes[k]
                      << " bytes arrived with " << (nerr == std::size_t(-1) ? "the wrong size" : "wrong bytes")
                      << std::endl;
        }
    }
    int all_ok = ok;
    world.gop.logic_and(&all_ok, 1);
    if (world.rank() == 0) std::cout << what << ": " << (all_ok ? "ok" : "FAILED") << std::endl;
    return all_ok;
}

/// SafeMPI::ByteCount itself, with no message: which sizes take a derived datatype
bool test_byte_count() {
    bool ok = true;
    {
        SafeMPI::ByteCount bc(1000, 1000);
        ok = ok && !bc.derived() && bc.count() == 1000 && bc.datatype() == MPI_BYTE;
    }
    if (SafeMPI::Is_initialized()) {
        for (std::size_t n : {std::size_t(1001), std::size_t(2000), std::size_t(123457)}) {
            SafeMPI::ByteCount bc(n, 1000);
            ok = ok && bc.derived() && bc.count() == 1;
#ifndef STUBOUTMPI
            MPI_Count size = 0;
            MPI_Type_size_x(bc.datatype(), &size);
            ok = ok && static_cast<std::size_t>(size) == n;
#endif
        }
    }
    return ok;
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    bool ok = true;

#ifdef MADNESS_TEST_RMI_HUGE_LARGE
    // one message of more than INT_MAX bytes, from rank 0 only: each one is held four
    // times over (the vector, its serialized form, the receive buffer, the vector again)
    const std::size_t n = (std::size_t(1) << 31) + (std::size_t(1) << 20) + 7;
    ok = send_and_check(world, world.rank() == 0 ? std::vector<std::size_t>{n} : std::vector<std::size_t>{},
                        "a message of more than INT_MAX bytes") && ok;
#else
    const std::size_t MiB = std::size_t(1) << 20;
    // a huge message in RMI's terms, whatever MAD_BUFFER_SIZE is; with one rank there is
    // no RMI to ask, and nothing is sent
    const std::size_t huge = (world.size() > 1 ? RMI::max_msg_len() : MiB) + 3 * MiB;

#ifndef STUBOUTMPI
    if (world.rank() == 0) {
        const bool bc_ok = test_byte_count();
        std::cout << "ByteCount: " << (bc_ok ? "ok" : "FAILED") << std::endl;
        ok = bc_ok && ok;
    }
#endif

    ok = send_and_check(world, {huge + 12345}, "a plain count") && ok;

    SafeMPI::ByteCount::set_max_plain_count(MiB);
    ok = send_and_check(world, {huge + 12345}, "whole chunks and a remainder") && ok;
    // the serialized task is larger than its payload, so pick the datatype's chunk to
    // divide the message rather than the payload: a chunk of 1 byte always does
    SafeMPI::ByteCount::set_max_plain_count(1);
    ok = send_and_check(world, {huge / 8}, "whole chunks only") && ok;
    SafeMPI::ByteCount::set_max_plain_count(MiB);
    ok = send_and_check(world, {huge + 1, 2 * huge + 77, huge + MiB, 100, huge + 5}, "several at once") && ok;
    SafeMPI::ByteCount::set_max_plain_count(0);
#endif

    world.gop.fence();
    finalize();
    if (!ok) return 1;
    return 0;
}
