/*
  This file is part of MADNESS.

  Copyright (C) 2019 Virginia Tech

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
 \file archive.cc
 \brief Definitions of serialization functions
 \ingroup serialization
*/

#include <cstring>
#include <cstddef>

#include <madness/world/archive.h>

namespace madness {
namespace archive {

namespace detail {
struct Ref {
  void fn() {}
};
}  // namespace detail

std::ptrdiff_t fn_ptr_origin() {
  static const std::ptrdiff_t result = []() {
    std::ptrdiff_t ptr;
    const auto ref_fn_ptr = &detail::Ref::fn;
    std::memcpy(&ptr, &ref_fn_ptr, sizeof(std::ptrdiff_t));
    return ptr;
  }();
  return result;
}

}  // namespace archive
}  // namespace madness
