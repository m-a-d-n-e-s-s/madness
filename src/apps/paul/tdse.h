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
#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <tensor/vmath.h>
#include <mra/lbdeux.h>

using namespace madness;

template <int NDIM>
class MyCostFn {
private:
  static SharedPtr<ApplyTime<NDIM> > apply_time;
  static double wt_factor;
public:
//    MyCostFn(SharedPtr<ApplyTime<NDIM> > app_ptr = SharedPtr<ApplyTime<NDIM> > ()) :
//      apply_time(app_ptr)
//      , wt_factor(100.0) {
//      if (apply_time) {
//        // find smallest non-zero entry in apply_time
//        // set wt_factor to order of magnitude of that entry
//      }
//    }

  static void set_defaults() {
    apply_time = SharedPtr<ApplyTime<NDIM> > ();
    wt_factor = 100.0;
  }

  template <typename T>
  static Cost cost_fun(const Key<NDIM> & key, const FunctionNode<T,NDIM>& node) {
    return actual_cost_fun(key);
  }
  static Cost actual_cost_fun(const Key<NDIM> & key) {
    if (!(apply_time)) return 0;
    double f = apply_time->get(key);
    return (Cost) (f*wt_factor);
  }
//   static long cost_fun(const Key<NDIM> & key, const FunctionNode<double_complex,NDIM>& node) {
//     return cost_fun(key);
//   }
//   static long cost_fun(const Key<NDIM> & key, const FunctionNode<double,NDIM>& node) {
//     return cost_fun(key);
//   }

  static void set_apply_time(SharedPtr<ApplyTime<NDIM> > apt) {
    apply_time = apt;
  }

  static void set_wt_factor() {
    // find smallest non-zero entry in apply_time
    // set wt_factor to order of magnitude of that entry
  }
  static void set_wt_factor(double w) {
    wt_factor = w;
  }

};
