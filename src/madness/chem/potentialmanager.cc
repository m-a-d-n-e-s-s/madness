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

  $Id$
*/

/// \file moldft/potentialmanager.cc
/// \brief Definition of commonly used density/potential related classes and functions

#include <madness/chem/potentialmanager.h>

namespace madness {

NuclearDensityFunctor::NuclearDensityFunctor(
    const Molecule &atoms,
    const BoundaryConditions<3> &bc,
    const Tensor<double> &cell,
    int special_level, double rscale)
    : atoms(atoms), bc_(bc), cell(cell), special_level_(special_level),
      special_points_(this->atoms.get_all_coords_vec()), rscale(rscale) {
  if (bc_.is_periodic_any()) {
    MADNESS_ASSERT(cell.ndim() == 2 && cell.dim(0) == 3 && cell.dim(1) == 2);
    this->maxR = 1;
  } else {
    this->maxR = 0;
  }
}

double NuclearDensityFunctor::operator()(const coord_3d &x) const {
  const auto tol = + 6.0*atoms.smallest_length_scale();
  double sum = 0.0;
  for (int i=-maxR; i<=+maxR; i++) {
    const double ext_x = cell(0,1) - cell(0,0);
    double xx = x[0]+i*ext_x;
    if (xx < cell(0,1) + tol && xx > cell(0,0) - tol) {
      for (int j=-maxR; j<=+maxR; j++) {
        const double ext_y = cell(1,1) - cell(1,0);
        double yy = x[1]+j*ext_y;
        if (yy < cell(1,1) + tol && yy > cell(1,0) - tol) {
          for (int k=-maxR; k<=+maxR; k++) {
            const double ext_z = cell(2,1) - cell(2,0);
            double zz = x[2]+k*ext_z;
            if (zz < cell(2,1) + tol && zz > cell(2,0) - tol) {
              sum += atoms.nuclear_charge_density(
                  xx, yy, zz, rscale);
            }
          }
        }
      }
    }
  }
  return sum;
}

std::vector<coord_3d> NuclearDensityFunctor::special_points() const {return special_points_;}

Level NuclearDensityFunctor::special_level() {
  return special_level_;
}

NuclearDensityFunctor &NuclearDensityFunctor::set_rscale(double rscale) {
  this->rscale = rscale;
  return *this;
}

}  // namespace madness
