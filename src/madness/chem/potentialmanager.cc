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
#include <madness/chem/SAP_interpolators.h>

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

Level NuclearDensityFunctor::special_level() const {
  return special_level_;
}

NuclearDensityFunctor &NuclearDensityFunctor::set_rscale(double rscale) {
  this->rscale = rscale;
  return *this;
}

double WignerSeitzPotentialFunctor::operator()(const coord_3d &r) const {
  enum { x = 0, y = 1, z = 2 };

  const auto natoms = atoms.natom();
  double E = 0.;

  const auto eval = [&](const auto &r1, const auto &r2, const auto rcut) {
    double restrictor_factor = 1.;
    double r12_sq = 0.;
    for (int d = 0; d != 3; ++d) {
      const auto x1 = r1[d];
      const auto x2 = r2[d];
      const auto x12 = x1 - x2;
      r12_sq += x12 * x12;
      restrictor_factor *= range[d].value(std::abs(x1 - x2) * rcell_width[d]);
    }
    const auto r12 = std::sqrt(r12_sq);
    return restrictor_factor * smoothed_potential(r12 * rcut) * rcut;
  };

  double result = 0.0;
  for (int cx = -lattice_sum_range[x]; cx <= lattice_sum_range[x]; ++cx) {
    for (int cy = -lattice_sum_range[y]; cy <= lattice_sum_range[y]; ++cy) {
      for (int cz = -lattice_sum_range[z]; cz <= lattice_sum_range[z]; ++cz) {
        // const auto intracell = cx==0 && cy==0 && cz==0;

        // r - (A + C) = (r - C) - A
        const std::array<double, 3> rC{r[x] - cx * cell_width[0],
                                       r[y] - cy * cell_width[1],
                                       r[z] - cz * cell_width[2]};

        for (std::size_t a = 0; a != natoms; ++a) {
          const auto &atom = atoms.get_atom(a);

          // make sure this isn't a pseudo-atom
          if (atom.pseudo_atom)
            continue;

          result -=
              atom.q * eval(rC, std::array<double, 3>{atom.x, atom.y, atom.z},
                            atoms.get_rcut()[a]);
        }
      }
    }
  }
  return result;
}

SAPFunctor::SAPFunctor(
    const Atom &atom,
    const BoundaryConditions<3> &bc,
    const Tensor<double> &cell,
    int special_level)
    : atom(atom), bc_(bc), cell(cell), special_level_(special_level),
      special_points_({atom.get_coords()}) {
}

double SAPFunctor::operator()(const coord_3d &x) const {
  const double zero_cutoff = 1e-10;
  double sum = 0.0;
  const auto ext_x = cell(0,1) - cell(0,0);
  const auto ext_y = cell(1,1) - cell(1,0);
  const auto ext_z = cell(2,1) - cell(2,0);
    const auto& interpolator = SAPCharges[atom.q - 1];
    auto maxX = bc_.is_periodic()[0] ? static_cast<int>(ceil(interpolator.get_hi() / ext_x)) : 0;
    auto maxY = bc_.is_periodic()[1] ? static_cast<int>(ceil(interpolator.get_hi() / ext_y)) : 0;
    auto maxZ = bc_.is_periodic()[2] ? static_cast<int>(ceil(interpolator.get_hi() / ext_z)) : 0;
    for (int dx = -maxX; dx <= maxX; dx++) {
      double xx = atom.x + dx*ext_x - x[0];
      for (int dy = -maxY; dy <= maxY; dy++) {
        double yy = atom.y + dy*ext_y - x[1];
        for (int dz = -maxZ; dz <= maxZ; dz++) {
          double zz = atom.z + dz*ext_z - x[2];
          auto r = std::sqrt(xx * xx + yy * yy + zz * zz);
          if (r > zero_cutoff and r < interpolator.get_hi()) sum += interpolator(r) / r;
        }
      }
    }

  return sum;
}

}  // namespace madness
