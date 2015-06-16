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
#include <cmath>
#include <cstdio>
using namespace madness;

class A {
public:
    void member() {
        PROFILE_MEMBER_FUNC(A);
        double sum = 0.0;
        for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
        std::printf("from A::member sum=%.6f\n", sum);
    }
};

void b() {
    PROFILE_FUNC;

    double sum = 0.0;
    for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
    std::printf("from B before sum=%.6f\n", sum);

    A a;
    a.member();

    for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
    std::printf("from B after  sum=%.6f\n", sum);

}

void a() {
    PROFILE_FUNC;

    double sum = 0.0;
    for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
    std::printf("from A before sum=%.6f\n", sum);

    b();

    for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
    std::printf("from A after  sum=%.6f\n", sum);
}

int main(int argc, char** argv) {
    initialize(argc, argv);

    World world(SafeMPI::COMM_WORLD);
    world.args(argc,argv);

    {
        PROFILE_BLOCK(main);

        double sum = 0.0;
        for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
        std::printf("from main before sum=%.6f\n", sum);

        a();

        for (int i=0; i<100000.0; i++) sum += sin(i*0.001);
        std::printf("from main after  sum=%.6f\n", sum);
    }

    print_stats(world);
    finalize();
    return 0;
}
