/*
  This file is part of MADNESS.

  Copyright (C) 2025 Virginia Tech

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


#ifndef MADNESS_MRA_MW_H__INCLUDED
#define MADNESS_MRA_MW_H__INCLUDED

#include <madness/madness_config.h>
#include <madness/mra/functypedefs.h>
#include <madness/mra/function_interface.h>
#include <madness/mra/legendre.h>

#include <array>
#include <functional>
#include <vector>

namespace madness {

/// this FunctionFunctorInterface evaluates Legendre scaling functions
template <std::size_t NDIM, typename Enabler = void>
struct ScalingFunctionFunctor;

// 1-d scaling function in {n,l} box on [-L,L]
template <>
struct ScalingFunctionFunctor<1> : public FunctionFunctorInterface<double, 1> {
  double L = -1;
  Level n = -1;
  Translation l;
  Translation twon;
  int k = -1;
  double one_over_twon;
  double sqrt_twon;
  double one_over_twoL;
  double one_over_sqrttwoL;

  ScalingFunctionFunctor() = default;

  ScalingFunctionFunctor(double L, Level n, Translation l, int k)
      : L(L), n(n), l(l), k(k), twon(1 << n), one_over_twon(1. / twon),
        sqrt_twon(sqrt(twon)), one_over_twoL(0.5/L), one_over_sqrttwoL(sqrt(one_over_twoL)) {
    MADNESS_ASSERT(L>0);
    MADNESS_ASSERT(n>=0);
    MADNESS_ASSERT(l >= 0 && l < twon);
    MADNESS_ASSERT(k>=0);
  }
  double operator()(const coord_1d &r) const final {
    const auto x = (r[0] + L)/(2 * L);
    double values[50];
    if (x < l * one_over_twon || x > (l + 1) * one_over_twon)
      return 0.;
    else {
      MADNESS_ASSERT(k < sizeof(values) / sizeof(double));
      legendre_scaling_functions(twon * x - l, k+1, values);
    }
    return values[k] * sqrt_twon * one_over_sqrttwoL;
  }

  virtual Level special_level() const final {
    return n;
  }

  virtual std::vector< Vector<double,1> > special_points() const final {
    return std::vector< Vector<double,1> >(1, Vector<double,1>{-L + 2*L*std::max(l+0.5,0.) * one_over_twon});
  }
};

template <std::size_t NDIM>
struct ScalingFunctionFunctor<NDIM, std::enable_if_t<std::greater{}(NDIM,1)>> : public FunctionFunctorInterface<double, NDIM> {
  std::array<double, NDIM> L;
  Key<NDIM> key;
  std::array<int, NDIM> k;
  std::array<ScalingFunctionFunctor<1>, NDIM> sf1;

  ScalingFunctionFunctor(const std::array<double, NDIM>& L, const Key<NDIM>& key, const std::array<int, NDIM>& k)
      : L(L), key(key), k(k) {
    for (int d = 0; d != NDIM; ++d)
      sf1[d] = ScalingFunctionFunctor<1>(L[d], key.level(), key[d], k.at(d));
  }
  double operator()(const Vector<double, NDIM> &r) const final {
    double result = 1.0;
    int d = 0;
    while (result != 0. && d < NDIM) {
      result *= sf1[d]({r[d]});
      ++d;
    }
    return result;
  }

  virtual Level special_level() const final {
    return key.level();
  }

  virtual std::vector< Vector<double,NDIM> > special_points() const final {
    Vector<double, NDIM> r;
    for (int d = 0; d != NDIM; ++d) {
      r[d] = sf1[d].special_points().at(0)[0];
    }
    return std::vector<Vector<double, NDIM>>(1, r);
  }
};

}

#endif // MADNESS_MRA_MW_H__INCLUDED
