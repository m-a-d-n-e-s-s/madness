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

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>
#include <madness/mra/lbdeux.h>
#include <madness/constants.h>
#include <stdlib.h>
using namespace madness;

// A class that behaves like a function to compute a Gaussian of given origin and exponent
class Gaussian : public FunctionFunctorInterface<double,3> {
public:
    const coord_3d center;
    const double exponent;
    const double coefficient;
    std::vector<coord_3d> specialpt;
    
    Gaussian(const coord_3d& center, double exponent, double coefficient)
        : center(center), exponent(exponent), coefficient(coefficient), specialpt(1)
    {
        specialpt[0][0] = center[0];
        specialpt[0][1] = center[1];
        specialpt[0][2] = center[2];
    }
    
    // MADNESS will call this interface
    double operator()(const coord_3d& x) const {
        double sum = 0.0;
        for (int i=0; i<3; i++) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        return coefficient*exp(-exponent*sum);
    }
    
    // // By default, adaptive projection into the spectral element basis
    // // starts uniformly distributed at the initial level.  However, if
    // // a function is "spiky" it may be necessary to project at a finer
    // // level but doing this uniformly is expensive.  This method
    // // enables us to tell MADNESS about points/areas needing deep
    // // refinement (the default is no special points).
    // std::vector<coord_3d> special_points() const {
    //     return specialpt;
    // }
};

/// A pmap that spatially decomposes the domain and by default slightly overdcomposes to attempt to load balance
class PartitionPmap : public WorldDCPmapInterface<Key<3>> {
private:
    const int nproc;
    Level target_level;
public:
    PartitionPmap()
        : nproc(1)
        , target_level(3)
    {};

    // Default is to try to optimize the target_level, but you can specify any value > 0
    PartitionPmap(size_t nproc, const Level target_level=0)
        : nproc(nproc)
    {
        if (target_level > 0) {
            this->target_level = target_level;
        }
        else {
            this->target_level = 1;
            int p = nproc-1;
            while (p) {
                p >>= 3;
                this->target_level++;
            }
        }            
    }

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        hashT hash;
        if (key.level() <= target_level) {
            hash = key.hash();
        }
        else {
            hash = key.parent(target_level - key.level()).hash();
        }
        return hash%nproc;
    }
};


void test2(World& world, size_t nfunc, size_t k, double thresh) {
    // Extra scope level to be sure all global data is freed before leaving
    {
        double expnt = 30000.0;
        double fac = std::pow(2.0*expnt/M_PI,0.25*3);
        
        FunctionDefaults<3>::set_cubic_cell(-6.0,6.0);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(false);
        FunctionDefaults<3>::set_truncate_on_project(true);
        FunctionDefaults<3>::set_truncate_mode(0); // Or should be 1 ????????????????????????????????????????????
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_pmap(std::shared_ptr<WorldDCPmapInterface<Key<3>>>(static_cast<WorldDCPmapInterface<Key<3>>*>(new PartitionPmap(world.size()))));
        {
            const int N = 6; // looking for where exp(-a*x^2) < 10**-N
            const int K = 6; // typically the lowest order of the polyn
            const double log10 = std::log(10.0);
            const double log2 = std::log(2.0);
            const double L = 12.0;
            const double a = expnt*L*L;
            double n = std::log(a/(4*K*K*(N*log10+std::log(fac))))/(2*log2);
            //std::cout << expnt << " " << a << " " << n << std::endl;
            Level initlev = Level(n<2 ? 2.0 : std::ceil(n)); // ????????????? should subtract one ???????????
            FunctionDefaults<3>::set_initial_level(initlev);
        }
        
        srand48(5551212); // for reproducible results
        for (size_t i=0; i<10000; i++) drand48(); // warmup generator
        
        double start_total = wall_time();
        double start;
        
        start = wall_time();
        std::vector<real_function_3d> funcs;
        for (size_t i=0; i<nfunc; i++) {
            coord_3d r;
            for (size_t d=0; d<3; d++) {
                r[d] = -6.0 + 12.0*drand48();
            }
            print(r, expnt, fac);
            funcs.push_back(real_factory_3d(world).functor(Gaussian(r, expnt, fac)).nofence());
        }
        world.gop.fence();
        double used_project = wall_time() - start;
        
        start = wall_time();
        std::vector<double> norms = norm2s(world, funcs);
        double used_norms_project = wall_time() - start;
        if (world.rank() == 0) print("norms after projection", norms);
        
        start = wall_time();
        compress(world, funcs);
        double used_compress = wall_time() - start;
        
        start = wall_time();
        norms = norm2s(world, funcs);
        double used_norms_compress = wall_time() - start;
        if (world.rank() == 0) print("norms after compress", norms);
        
        start = wall_time();
        reconstruct(world, funcs);
        double used_reconstruct = wall_time() - start;
        
        start = wall_time();
        norms = norm2s(world, funcs);
        double used_norms_reconstruct = wall_time() - start;
        if (world.rank() == 0) print("norms after reconstruct", norms);
        double used_total = wall_time() - start_total;
        
        if (world.rank() == 0) {
            print("  #processes : ", world.size());
            print("    #threads : ", ThreadPool::size()+1); // does NOT include the communication thread and assumes main thread is computing
            print("     project : ", used_project);
            print("    compress : ", used_compress);
            print(" reconstruct : ", used_reconstruct);
            print("      norm-1 : ", used_norms_project);
            print("      norm-2 : ", used_norms_compress);
            print("      norm-3 : ", used_norms_reconstruct);
            print("       total : ", used_total);
        }
    }
    // Extra scope level to be sure all global data is freed before leaving
    world.gop.fence();
    world.gop.fence();
}

int main(int argc, char** argv) {
    // Initialize the parallel programming environment
    initialize(argc,argv);
    // Extra scope to make sure world is destroyed before finalizing
    {
        World world(SafeMPI::COMM_WORLD);
        startup(world,argc,argv);
        test2(world, 1, 10, 1e-8);
    }
    finalize();
    return 0;
}
