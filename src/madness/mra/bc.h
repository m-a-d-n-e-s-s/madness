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

#ifndef MADNESS_MRA_BC_H__INCLUDED
#define MADNESS_MRA_BC_H__INCLUDED


/// \file bc.h
/// \brief Provides BoundaryConditions
/// \ingroup mrabcext

#include <madness/world/madness_exception.h>

#include <cstddef>
#include <iostream>
#include <vector>

namespace madness {

enum BCType {
  BC_ZERO,
  BC_PERIODIC,
  BC_FREE,
  BC_DIRICHLET,
  BC_ZERONEUMANN,
  BC_NEUMANN
};

/*!
  \brief This class is used to specify boundary conditions for all operators
  \ingroup mrabcext

  Exterior boundary conditions (i.e., on the simulation domain)
  are associated with operators (not functions).  The types of
  boundary conditions available are in the enum BCType.

  The default boundary conditions are obtained from the FunctionDefaults.
  For non-zero Dirichlet and Neumann conditions additional information
  must be provided when derivative operators are constructed. For integral
  operators, only periodic and free space are supported.
*/
template <std::size_t NDIM> class BoundaryConditions {
private:
  // Used to use STL vector but static data on  a MAC was
  // causing problems.
  BCType bc[NDIM * 2];

public:
  /// Constructor. Default boundary condition set to free space
  BoundaryConditions(BCType code = BC_FREE) {
    for (std::size_t i = 0; i < NDIM * 2; ++i)
      bc[i] = code;
  }

  /// Copy constructor is deep
  BoundaryConditions(const BoundaryConditions<NDIM> &other) { *this = other; }

  /// Assignment makes deep copy
  BoundaryConditions<NDIM> &operator=(const BoundaryConditions<NDIM> &other) {
    if (&other != this) {
      for (std::size_t i = 0; i < NDIM * 2; ++i)
        bc[i] = other.bc[i];
    }
    return *this;
  }

  /// Returns value of boundary condition

  /// @param d Dimension (0,...,NDIM-1) for boundary condition
  /// @param i Side (0=left, 1=right) for boundary condition
  /// @return Value of boundary condition
  BCType operator()(std::size_t d, int i) const {
    MADNESS_ASSERT(d < NDIM && i >= 0 && i < 2);
    return bc[2 * d + i];
  }

  /// Returns reference to boundary condition

  /// @param d Dimension (0,...,NDIM-1) for boundary condition
  /// @param i Side (0=left, 1=right) for boundary condition
  /// @return Value of boundary condition
  BCType &operator()(std::size_t d, int i) {
    MADNESS_ASSERT(d < NDIM && i >= 0 && i < 2);
    return bc[2 * d + i];
  }

  template <typename Archive> void serialize(const Archive &ar) { ar & bc; }

  /// Translates code into human readable string

  /// @param code Code for boundary condition
  /// @return String describing boundary condition code
  static const char *code_as_string(BCType code) {
    static const char *codes[] = {"zero",      "periodic",     "free",
                                  "Dirichlet", "zero Neumann", "Neumann"};
    return codes[code];
  }

  /// Convenience for application of integral operators

  /// @return Returns a vector indicating if each dimension is periodic
  std::vector<bool> is_periodic() const {
    std::vector<bool> v(NDIM);
    for (std::size_t d = 0; d < NDIM; ++d)
      v[d] = (bc[2 * d] == BC_PERIODIC);
    return v;
  }

  /// Checks whether the boundary condition along any axis is periodic

  /// @return Returns true if any dimension is periodic
  bool is_periodic_any() const {
    for (std::size_t d = 0; d < NDIM; ++d)
      if (bc[2 * d] == BC_PERIODIC)
        return true;
    return false;
  }
};

template <std::size_t NDIM>
static inline std::ostream &operator<<(std::ostream &s,
                                       const BoundaryConditions<NDIM> &bc) {
  s << "BoundaryConditions(";
  for (unsigned int d = 0; d < NDIM; ++d) {
    s << bc.code_as_string(bc(d, 0)) << ":" << bc.code_as_string(bc(d, 1));
    if (d == NDIM - 1)
      s << ")";
    else
      s << ", ";
  }
  return s;
}

}  // namespace madness

#endif // MADNESS_MRA_BC_H__INCLUDED
