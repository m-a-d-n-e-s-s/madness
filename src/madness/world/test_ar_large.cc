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

/// \file test_ar_large.cc
/// \brief An array of more than 2^32 elements through the buffer and vector archives.
///
/// The count of an array used to be an unsigned int, so such an array was cut to its
/// size modulo 2^32 when stored, with no error. Takes about 13 GiB, so the test carries
/// the "large" label; test_ar.cc checks the count alone, which takes no memory.

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <madness/world/MADworld.h>
#include <madness/world/buffer_archive.h>
#include <madness/world/vector_archive.h>

using namespace madness;

namespace {

inline unsigned char pattern(std::size_t i) {
    return static_cast<unsigned char>(((i + 17) * 0x9E3779B97F4A7C15ull) >> 56);
}

std::size_t count_errors(const std::vector<unsigned char>& v, std::size_t n) {
    if (v.size() != n) return std::size_t(-1);
    std::size_t nerr = 0;
    for (std::size_t i = 0; i < n; ++i) nerr += (v[i] != pattern(i));
    return nerr;
}

bool report(const char* what, std::size_t nerr) {
    std::cout << what << ": " << (nerr == 0 ? "ok" : nerr == std::size_t(-1) ? "FAILED, wrong size" : "FAILED, wrong contents")
              << std::endl;
    return nerr == 0;
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    bool ok = true;

    if (world.rank() == 0) {
        const std::size_t n = (std::size_t(1) << 32) + 5;
        std::vector<unsigned char> v(n);
        for (std::size_t i = 0; i < n; ++i) v[i] = pattern(i);

        // a std::vector through a buffer archive, as the argument of a remote task goes
        {
            archive::BufferOutputArchive count;
            count & v;
            const std::size_t nbyte = count.size();
            ok = report("size of the stored vector", nbyte >= n ? 0 : std::size_t(-1)) && ok;
            std::vector<unsigned char> buf(nbyte);
            archive::BufferOutputArchive out(buf.data(), nbyte);
            out & v;
            std::vector<unsigned char> w;
            archive::BufferInputArchive in(buf.data(), nbyte);
            in & w;
            ok = report("vector through a buffer archive", count_errors(w, n)) && ok;
        }

        // a bare array through a vector archive
        {
            std::vector<unsigned char> buf;
            archive::VectorOutputArchive out(buf);
            out << archive::wrap(v.data(), n);
            ok = report("size of the stored array", buf.size() >= n ? 0 : std::size_t(-1)) && ok;
            std::vector<unsigned char> w(n);
            archive::VectorInputArchive in(buf);
            in >> archive::wrap(w.data(), n);
            ok = report("array through a vector archive", count_errors(w, n)) && ok;
        }
    }

    world.gop.fence();
    finalize();
    return ok ? 0 : 1;
}
