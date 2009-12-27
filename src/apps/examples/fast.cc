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
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>
#include <mra/lbdeux.h>

using namespace madness;
using namespace std;

Void doit(int niter) {
    Tensor<double> c(20,20,20), d(20,20,20), e(20,20,20);
    Tensor<double> a(20,20);
    for (int i=0; i<niter; i++) {
        fast_transform(c,a,d,e);
    }
    return None;
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    try {
        World world(MPI::COMM_WORLD);
        cout.precision(8);

        int niter = 100000;

        world.gop.fence();
        double start = wall_time();
        if (world.rank() == 0) print("starting at", start);
        world.gop.fence();
        for (unsigned int i=0; i<=ThreadPool::size(); i++) world.taskq.add(doit, niter); // Note +1
        world.taskq.fence();
        world.gop.fence();
        double end = wall_time();
        if (world.rank() == 0) print("starting at", end);
        world.gop.fence();

        ThreadPool::end();
        print_stats(world);

        if (world.rank() == 0) {
            double nflop = (ThreadPool::size()+1.0)*20.0*20.0*20.0*20.0*2.0*3.0*niter;
            print("");
            print("NFLOPS ", nflop);
            print("  USED ", end-start);
            print("  RATE ", nflop/(end-start));
            print("");
        }
    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }
    
    finalize();
    return 0;
}
