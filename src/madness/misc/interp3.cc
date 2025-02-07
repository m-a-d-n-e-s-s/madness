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
#include <iostream>
#include <cmath>
#include <vector>
#include <madness/world/worldinit.h>
#include <madness/misc/interpolation_1d.h>

using namespace std;

/// A simple program for testing the CubicInterpolationTable class.

double func(double x) {
    return sin(x);
}

int main(int argc, char* argv[]) {
  auto& world = madness::initialize(argc, argv);
  cout.precision(12);

  // Uniform mesh for sin(x)
  madness::CubicInterpolationTable<double> fit(world, -10.0, 30.0, 1000001, func);
  cout << "maxerr " << fit.err(func) << endl;

  std::vector<double> xs{0.1, 0.25, 0.3, 0.42, 0.5, 0.7, 0.8};
  std::vector<double> ys(xs.size());

  for (size_t i = 0; i < ys.size(); i++) {
    auto& x = xs[i];
    ys[i] = 1 + 2 * x + 3 * x*x + 4.5 * x*x*x;
  }
  madness::CubicInterpolationTable<double> points(xs, ys);

  for (size_t i = 0; i < ys.size(); i++) {
    if(abs(points(xs[i]) - ys[i]) > 1E-12) {
      throw std::runtime_error("Point based cubic table fails");
    }
  }

  madness::finalize();
  return 0;
}
