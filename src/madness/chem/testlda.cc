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

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <fstream>
#include "xcfunctional.h"

using namespace madness;

struct xcfunc_data_point
{
  double rhoa, rhob;
  double sigmaaa, sigmaab, sigmabb;
  double zk;
  double vrhoa, vrhob;
  double vsigmaaa, vsigmaab, vsigmabb;
  double v2rhoa2, v2rhoab, v2rhob2;
};

void test_lda(World& world)
{
    XCfunctional xcfunc;
    xcfunc.initialize("LDA RHOTOL 1e-12 RHOMIN 1e-12", false, world);

    /*

      generated WITHOUT libxc

100 -887.752 -5.88809
10 -42.0542 -2.78016
1 -2.01597 -1.32674
0.1 -0.0981113 -0.64225
0.01 -0.00485142 -0.316006
0.001 -0.000242858 -0.157725
0.0001 -1.22083e-05 -0.0793057
1e-05 -6.10474e-07 -0.0397809
1e-06 -3.01609e-08 -0.0197451
1e-07 -1.46896e-09 -0.00966166
1e-08 -7.06144e-11 -0.00466261
1e-09 -3.35942e-12 -0.00222473
1e-10 -1.58607e-13 -0.00105256
1e-11 -7.44811e-15 -0.000494995
1e-12 -3.48465e-16 -0.000231817
1e-13 -1.38544e-16 -0.000184377
1e-14 -1.38544e-16 -0.000184377
1e-15 -1.38544e-16 -0.000184377

     */

    const long N = 18;
    Tensor<double> rho(N);
    double x = 100.0;
    for (int i=0; i<18; i++) {
        rho[i] = x;
        x *= 0.1;
    }

    std::vector<Tensor<double>> t {rho};

    Tensor<double> e = xcfunc.exc(t);
    std::vector<Tensor<double>> v = xcfunc.vxc(t, 0);

    for (int i=0; i<N; i++) {
        print(rho[i], e[i], v[0][i]);
    }
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();

    test_lda(world);

    madness::finalize();
    return 0;
}
