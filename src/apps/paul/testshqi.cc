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
/// \file tdse.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <tensor/vmath.h>

#include <mra/loadbal.h>
//#include <mra/lbdeux.h>
#include <mra/funcimpl.h>

using namespace madness;

const int LBMODE=1; // 1 for constant-cost, 2 for profiled-cost
const double PI = 3.1415926535897932384;

// typedefs to make life less verbose
typedef Singleton< LBTimer<3> > CostFun;
typedef Vector<double,3> coordT;
typedef SharedPtr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef SharedPtr< FunctionFunctorInterface<double_complex,3> > complex_functorT;
typedef Function<double_complex,3> complex_functionT;
typedef FunctionFactory<double_complex,3> complex_factoryT;
typedef SeparatedConvolution<double_complex,3> complex_operatorT;
typedef SharedPtr< WorldDCPmapInterface< Key<3> > > pmapT;

// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    // Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;

        // This randomly hashes levels 0-2 and then
        // hashes nodes by their grand-parent key so as
        // to increase locality separately on each level.
        //if (n <= 2) hash = key.hash();
        //else hash = key.parent(2).hash();

        // This randomly hashes levels 0-3 and then 
        // maps nodes on even levels to the same
        // random node as their parent.
        // if (n <= 3 || (n&0x1)) hash = key.hash();
        // else hash = key.parent().hash();

        // This randomly hashes each key
        hash = key.hash();

        return hash%nproc;
    }
};

double myreal(double t) {return t;}

double myreal(const double_complex& t) {return real(t);}

template <typename T, int NDIM>
class GaussianFunctor : public FunctionFunctorInterface<T,NDIM> {
private:
    typedef Vector<double,NDIM> coordT;
    const std::vector<coordT> center;
    const std::vector<double> exponent;
    const std::vector<T> coefficient;

public:
    GaussianFunctor(const std::vector<coordT>& center, std::vector<double> exponent, std::vector<T> coefficient)
            : center(center), exponent(exponent), coefficient(coefficient) {};

    GaussianFunctor(const coordT& center, double exponent, T coefficient)
            : center(std::vector<coordT>(1,center)), exponent(std::vector<double>(1,exponent)), coefficient(std::vector<T>(1,coefficient)) {};

    virtual T operator()(const coordT& x) const {
        T retval = 0;
        for (unsigned int j=0; j<center.size(); j++) {
            double sum = 0.0;
            for (int i=0; i<NDIM; i++) {
                double xx = center[j][i]-x[i];
                sum += xx*xx;
            }
            retval += coefficient[j]*exp(-exponent[j]*sum);
        }
        return retval;
    };
  GaussianFunctor operator+ (const GaussianFunctor& other) const {
        std::vector<coordT> newcenter = this->center;
        std::vector<double> newexponent = this->exponent;
        std::vector<T> newcoefficient = this->coefficient;
        newcenter.insert(newcenter.end(), other.center.begin(), other.center.end());
        newexponent.insert(newexponent.end(), other.exponent.begin(), other.exponent.end());
        newcoefficient.insert(newcoefficient.end(), other.coefficient.begin(), other.coefficient.end());
        return (GaussianFunctor(newcenter, newexponent, newcoefficient));
    };

    void print() const {
        madness::print("Sum of", center.size(), "gaussians:");
        for (unsigned int i = 0; i < center.size(); i++) {
            madness::print("   g[", i, "] : =", coefficient[i], "* exp(", -exponent[i], "(", center[i], "- x )^2 )");
        }
    };
};


void doit(World& world) {
   if (world.rank() == 0) print("at beginning of test_loadbal");
   if (world.rank() == 0) print("lbmode = ", LBMODE);

    //    for (int i=0; i<3; i++) {
    //        FunctionDefaults<3>::cell(i,0) = -10.0;
    //        FunctionDefaults<3>::cell(i,1) =  10.0;
    //    }

    FunctionDefaults<3>::set_cubic_cell(-10.0, 10.0);
 int nspikes = 2;
    //int nspikes = 10;
    std::vector<coordT> vorigin(nspikes);
    const double expnt = 64.0;
    const double expnt1 = 4096;
    std::vector<double> vexpnt(nspikes);
    Vector<double, 3> dcell, avgcell;
    for (int i = 0; i < 3; i++) {
        double cell0 = FunctionDefaults<3>::get_cell()(i,0);
        double cell1 = FunctionDefaults<3>::get_cell()(i,1);
        dcell[i] = cell0 - cell1;
        avgcell[i] = (cell0 + cell1)/2;
    }
    for (int i = 0; i < nspikes; i++) {
        Vector<double, 3> v(0);
        for (int j = 0; j < 3; j++) {
            v[j] = 0.2;
            //v[j] = ((double) rand()/(double) RAND_MAX)*dcell[j] - FunctionDefaults<3>::cell(j,1);
        }
        if (i%2) {
            vexpnt[i] = expnt1;
        }
        else {
            vexpnt[i] = expnt;
        }
        vorigin[i] = v;
    }
  const double coeff = pow(2.0/PI,0.25*3);
    std::vector<double> vcoeff(nspikes,coeff);
    if (world.rank() == 0) print("about to make gf");
    GaussianFunctor<double,3> *gf = new GaussianFunctor<double,3>(vorigin, vexpnt, vcoeff);
    if (world.rank() == 0) print("about to make functor");
    functorT functor(gf);
    if (world.rank() == 0) {
        gf->print();
    }

    for (int k=4; k >=4; k-=2) {
        //MyPmap<3> temppmap(world);
        FunctionDefaults<3>::set_pmap(SharedPtr<MyPmap<3> >(new MyPmap<3>(world)));
        //if (LBMODE == 2) CostFun::Instance(world).set_master_status(true);
        CostFun::Instance(world).init(FunctionDefaults<3>::get_pmap());
        CostFun::Instance(world).set_master_status(true);
        CostFun::Instance(world).reset();
        /*
        if (world.rank() == 0) {
            FunctionDefaults<3>::pmap->print();
        }
        */
        //int n = 2;
        int n = 4;
     //   int n = 7;
   //double thresh = 1e-12;
        //double thresh = 5e-4;
        //double thresh = 1e-3;
        //double thresh = 1e-6;
        double thresh = 5e-4;
        if (world.rank() == 0) {
            printf("k=%d:\n", k);
        }
        double t0 = MPI::Wtime();
        //Function<double,3> f = FunctionFactory<double,3>(world).functor(functor).norefine().initial_level(n).k(k);
        Function<double,3> f = FunctionFactory<double,3>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        if (world.rank() == 0) {
            print("just made function f");
        }
        double t1 = MPI::Wtime();
        double err2 = f.err(*functor);
        std::size_t size = f.size();
        std::size_t tree_size = f.tree_size();
        std::size_t maxsize = f.max_nodes();
        std::size_t minsize = f.min_nodes();
        std::size_t maxdepth = f.max_depth();
        double t2 = MPI::Wtime();
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err2, double(size), double(tree_size), double(maxdepth), double(maxsize), double(minsize), abs(log(err2)/n/k));
        }
   world.gop.fence();
        double t3 = MPI::Wtime();
        Function<double,3> g = copy(f);
        world.gop.fence();
        double t4 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO COMPRESS");
        f.compress(true);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t5 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO MULTIPLY");
        Function<double,3> ff = f*f;
        world.gop.fence();
        double t51 = MPI::Wtime();
        ff = f*f;
        world.gop.fence();
        double t52 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO RECON");
        f.reconstruct(true);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t6 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO CREATE LOADBAL");
        //  typedef Cost (*my_default_fun<3,double>) (const Key<3>&, const FunctionNode<double,3>&);
        //LoadBalImpl<3> lb(g, lbcost<double,3>(1.0, 1.0));
        if (LBMODE == 1) CostFun::Instance(world).set_master_status(false);
        CostFun::Instance(world).compute_min_cost();
        LoadBalImpl<3> lb(g, CostFun::Instance(world));
        //        LoadBalImpl<3> lb(g, &my_default_fun<double,3>(Key<3>, FunctionNode<double,3>));
        //  LoadBalImpl<3> lb(g, &f_ptr);
        if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        double t7 = MPI::Wtime();
        if (world.rank() == 0) print("ABOUT TO DO LOADBAL");
        FunctionDefaults<3>::set_pmap(lb.load_balance());
        CostFun::Instance(world).remap(FunctionDefaults<3>::get_pmap(),true,true);

        /*
        if (world.rank() == 0) {
            FunctionDefaults<3>::pmap->print();
        }
        */
        if (LBMODE == 1) CostFun::Instance(world).set_master_status(true);
     if (world.rank() == 0) print("ABOUT TO FENCE");
        world.gop.fence();
        Function<double,3> h = FunctionFactory<double,3>(world).functor(functor).thresh(thresh).initial_level(n).k(k);
        double t8 = MPI::Wtime();
        f = copy(g,FunctionDefaults<3>::get_pmap());
        double t9 = MPI::Wtime();
        double err21 = h.err(*functor);
        std::size_t size1 = h.size();
        std::size_t tree_size1 = h.tree_size();
        std::size_t maxsize1 = h.max_nodes();
        std::size_t minsize1 = h.min_nodes();
        std::size_t maxdepth1 = h.max_depth();
        double t10 = MPI::Wtime();
        if (world.rank() == 0) {
            printf("   n=%d err=%.2e #coeff=%.2e tree_size=%.2e max_depth=%.2e max_nodes=%.2e min_nodes=%.2e log(err)/(n*k)=%.2e\n",
                   n, err21, double(size1), double(tree_size1), double(maxdepth1), double(maxsize1), double(minsize1), abs(log(err21)/n/k));
        }
        world.gop.fence();
        double t11 = MPI::Wtime();
        h.compress(true);
        world.gop.fence();
        double t12 = MPI::Wtime();
        Function<double,3> hh = h*h;
        world.gop.fence();
        double t121 = MPI::Wtime();
        h.reconstruct(true);
        world.gop.fence();
        double t13 = MPI::Wtime();


    if (world.rank() == 0) {
            madness::print("for f with k =", k, "n =", n, "and thresh =",
                           thresh, "running on", world.nproc(), "processors:");
            madness::print("Routine            |  Time");
            madness::print("-------------------+--------------");
            madness::print("Init f             | ", t1-t0);
            madness::print("Diagnostics for f  | ", t2-t1);
            madness::print("Copy f to g        | ", t4-t3);
            madness::print("Compress f         | ", t5-t4);
            madness::print("Multiply f         | ", t51-t5);
            madness::print("Multiply f 2       | ", t52-t51);
            madness::print("Reconstruct f      | ", t6-t52);
            madness::print("Create LoadBalImpl | ", t7-t6);
            madness::print("Load Balance       | ", t8-t7);
            madness::print("Copy g to f        | ", t9-t8);
            madness::print("Diagnostics for h  | ", t10-t9);
            madness::print("Compress h         | ", t12-t11);
            madness::print("Multiply h         | ", t121-t12);
            madness::print("Reconstruct h      | ", t13-t121);
            madness::print("-------------------+--------------");
            madness::print("Total Time         | ", t13-t0);
            madness::print("");
        }
    }

    if (world.rank() == 0) print("test loadbal OK\n\n");
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    
    startup(world,argc,argv);

    try {
        doit(world);
    } catch (const MPI::Exception& e) {
        //print(e); std::cout.flush();
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e); std::cout.flush();
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e); std::cout.flush();
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s); std::cout.flush();
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what()); std::cout.flush();
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }


    world.gop.fence();
    
    ThreadPool::end();
    print_stats(world);
    finalize();
    return 0;
}


