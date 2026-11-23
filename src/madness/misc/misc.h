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

#ifndef MADNESS_MISC_MISC_H__INCLUDED
#define MADNESS_MISC_MISC_H__INCLUDED

/// \file misc.h
/// \brief Header to declare stuff which has not yet found a home

#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>

namespace madness {
    /// True if `x` is finite (neither NaN nor +/-Inf).
    /// Inspects IEEE-754 bits directly so `-ffast-math` / `-ffinite-math-only`
    /// cannot fold the check to `true`.
    inline constexpr bool is_finite(double x) noexcept {
        static_assert(std::numeric_limits<double>::is_iec559 &&
                      sizeof(double) == sizeof(std::uint64_t));
        constexpr std::uint64_t exp_mask = 0x7ff0000000000000ULL;
        return (std::bit_cast<std::uint64_t>(x) & exp_mask) != exp_mask;
    }

    inline constexpr bool is_finite(float x) noexcept {
        static_assert(std::numeric_limits<float>::is_iec559 &&
                      sizeof(float) == sizeof(std::uint32_t));
        constexpr std::uint32_t exp_mask = 0x7f800000U;
        return (std::bit_cast<std::uint32_t>(x) & exp_mask) != exp_mask;
    }

    /// True if `x` is NaN, even under `-ffast-math` / `-ffinite-math-only`.
    inline constexpr bool is_nan(double x) noexcept {
        static_assert(std::numeric_limits<double>::is_iec559 &&
                      sizeof(double) == sizeof(std::uint64_t));
        constexpr std::uint64_t exp_mask = 0x7ff0000000000000ULL;
        constexpr std::uint64_t mag_mask = 0x7fffffffffffffffULL;
        return (std::bit_cast<std::uint64_t>(x) & mag_mask) > exp_mask;
    }

    inline constexpr bool is_nan(float x) noexcept {
        static_assert(std::numeric_limits<float>::is_iec559 &&
                      sizeof(float) == sizeof(std::uint32_t));
        constexpr std::uint32_t exp_mask = 0x7f800000U;
        constexpr std::uint32_t mag_mask = 0x7fffffffU;
        return (std::bit_cast<std::uint32_t>(x) & mag_mask) > exp_mask;
    }
    unsigned long checksum_file(const char* filename);
    std::istream& position_stream(std::istream& f, const std::string& tag, bool rewind=true);
    std::istream& position_stream_to_word(std::istream& f, const std::string& tag,
                                          const char comment='#', bool rewind=true, bool silent=false,
                                          bool first_word_only=false);
    std::string lowercase(const std::string& s);
    void gprofexit(int id, int nproc);
    /// creates a unique filename, using PBS ID if available
    std::string unique_fileid();

    // get the memory usage of the current process in byte
    long get_memory_usage();

    bool nearlyEqual(double a, double b, double epsilon = 1e-9);
}

#endif // MADNESS_MISC_MISC_H__INCLUDED



