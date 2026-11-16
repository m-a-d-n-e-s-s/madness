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

/// \file test_rmi_huge.h
/// \brief Messages too large for an MPI count of bytes, sent as the argument of a remote task.
///
/// A message longer than \c MAD_BUFFER_SIZE takes RMI's huge-message protocol, and one of
/// more than \c INT_MAX bytes must further be described to MPI by a derived datatype
/// ( \c SafeMPI::ByteCount ). This is the body of two tests, each a source file that
/// includes it:
///   - test_rmi_huge.cc lowers the size at which the derived datatype takes over, so that
///     messages of a few MiB cover that path: one that leaves a remainder after whole
///     chunks, one that does not, one that stays a plain count, and several at once. It
///     also checks the datatypes ByteCount builds for sizes beyond \c INT_MAX , which
///     takes no memory.
///   - test_rmi_huge_large.cc defines \c MADNESS_TEST_RMI_HUGE_LARGE and sends one message
///     of more than \c INT_MAX bytes at the default setting; that takes about 9 GiB, so it
///     carries the "large" label.
///
/// Needs at least 2 ranks to send anything; with 1 the tasks are local.

#ifndef MADNESS_WORLD_TEST_RMI_HUGE_H__INCLUDED
#define MADNESS_WORLD_TEST_RMI_HUGE_H__INCLUDED

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <madness/world/MADworld.h>
#include <madness/world/safempi.h>

using namespace madness;

namespace {

/// The payload is words rather than bytes: serialization counts the elements of an array
/// in an unsigned int, so a vector of bytes could not test a message of 4 GiB or more.
using wordT = std::uint64_t;

/// the word at position \c i of the message with the given \c seed
inline wordT pattern(std::size_t i, std::uint64_t seed) {
    return (i + seed) * 0x9E3779B97F4A7C15ull;
}

/// \param nbyte the size of the message, rounded up to whole words
std::vector<wordT> make_message(std::size_t nbyte, std::uint64_t seed) {
    std::vector<wordT> v((nbyte + sizeof(wordT) - 1) / sizeof(wordT));
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = pattern(i, seed);
    return v;
}

/// runs on the receiver: the number of words that differ from the pattern, or
/// \c size_t(-1) if the size is wrong
std::size_t count_errors(const std::vector<wordT>& v, std::size_t nbyte, std::uint64_t seed) {
    if (v.size() != (nbyte + sizeof(wordT) - 1) / sizeof(wordT)) return std::size_t(-1);
    std::size_t nerr = 0;
    for (std::size_t i = 0; i < v.size(); ++i) nerr += (v[i] != pattern(i, seed));
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
                      << " bytes arrived with " << (nerr == std::size_t(-1) ? "the wrong size" : "wrong contents")
                      << std::endl;
        }
    }
    int all_ok = ok;
    world.gop.logic_and(&all_ok, 1);
    if (world.rank() == 0) std::cout << what << ": " << (all_ok ? "ok" : "FAILED") << std::endl;
    return all_ok;
}

#ifndef STUBOUTMPI
/// whether \c bc describes \c n contiguous bytes by a derived datatype
bool describes(const SafeMPI::ByteCount& bc, std::size_t n) {
    if (!bc.derived() || bc.count() != 1) return false;
    // contiguous bytes from 0: the extent, with and without padding, is the size
    MPI_Aint lb = -1, extent = -1, true_lb = -1, true_extent = -1;
    MPI_Type_get_extent(bc.datatype(), &lb, &extent);
    MPI_Type_get_true_extent(bc.datatype(), &true_lb, &true_extent);
    return lb == 0 && true_lb == 0 && static_cast<std::size_t>(extent) == n &&
           static_cast<std::size_t>(true_extent) == n;
}

/// SafeMPI::ByteCount itself, with no message and no buffer: which sizes take a derived
/// datatype, and what it describes, sizes beyond \c INT_MAX at the default chunk included
bool test_byte_count() {
    const std::size_t int_max = static_cast<std::size_t>(std::numeric_limits<int>::max());
    bool ok = true;
    {
        SafeMPI::ByteCount bc(1000, 1000);
        ok = ok && !bc.derived() && bc.count() == 1000 && bc.datatype() == MPI_BYTE;
    }
    {
        SafeMPI::ByteCount bc(int_max);
        ok = ok && !bc.derived() && bc.count() == std::numeric_limits<int>::max();
    }
    for (std::size_t n : {std::size_t(1001), std::size_t(2000), std::size_t(123457)}) {
        SafeMPI::ByteCount bc(n, 1000);
        ok = ok && describes(bc, n);
    }
    for (std::size_t n : {int_max + 1, 2 * int_max, 2 * int_max + 1, (std::size_t(1) << 31) + (std::size_t(1) << 20) + 7,
                          std::size_t(5) << 32, (std::size_t(1) << 40) + 3}) {
        SafeMPI::ByteCount bc(n);
        ok = ok && describes(bc, n);
    }
    return ok;
}
#endif

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    bool ok = true;

#ifdef MADNESS_TEST_RMI_HUGE_LARGE
    if (world.size() == 1) {
        std::cout << "a message of more than INT_MAX bytes: skipped, needs 2 ranks" << std::endl;
        finalize();
        return 0;
    }
    // one message of more than INT_MAX bytes, from rank 0 only: each one is held four
    // times over (the vector, its serialized form, the receive buffer, the vector again)
    // MAD_TEST_RMI_HUGE_BYTES overrides the size, to try others by hand
    std::size_t n = (std::size_t(1) << 31) + (std::size_t(1) << 20) + 7;
    if (const char* env = std::getenv("MAD_TEST_RMI_HUGE_BYTES")) n = std::strtoull(env, nullptr, 10);
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

#endif // MADNESS_WORLD_TEST_RMI_HUGE_H__INCLUDED
