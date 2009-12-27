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
#include <world/world.h>

#include <utility>
#include <tensor/tensor.h>

#include <ii/systolic.h>

using namespace madness;


template <typename T>
class TestSystolicMatrixAlgorithm : public SystolicMatrixAlgorithm<T> {
    volatile int niter;
    DistributedMatrix<T>& A;
    World& world;
public:
    TestSystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag) 
        : SystolicMatrixAlgorithm<T>(A, tag) 
        , niter(0)
        , A(A)
        , world(A.get_world())
    {
        madness::print("Testing SystolicMatrixAlgorithm ", 
                       SystolicMatrixAlgorithm<T>::get_coldim(), 
                       SystolicMatrixAlgorithm<T>::get_rowdim());
    }
    
    /* impliment kernel */
    /* this time, check row i,j element for all rows */

    void kernel(int i, int j, T* rowi, T* rowj) {
        for (int k=0; k < SystolicMatrixAlgorithm<T>::get_rowdim(); k++) {
            MADNESS_ASSERT(rowi[k] == i);
            MADNESS_ASSERT(rowj[k] == j);
        }
        print("In kernel column", i, ",", j);
    }

    void start_iteration_hook(const TaskThreadEnv& env) {
        int id = env.id();
        if (id == 0) {
            niter++;
            world.gop.max(A.data()(0,0));
        }
        env.barrier();
    }

    bool converged(const TaskThreadEnv& env) const {
        if (niter >= 3) { /* except first 3 times */
            if (env.id() == 0) {
                madness::print("    done!");
            }
            return true;
        }
        else {
            return false;
        }
    }
};

class local_iterator{
private:
    int64_t first;
    int64_t last;
    int64_t size;
    int64_t rest;
public:
    local_iterator(int64_t first, int64_t last_):
        first(first),
        last(last_ +1),
        size(last-first),
        rest(size)
    {}

    template <class T>
    local_iterator(const DistributedMatrix<T> &A):
        size(A.local_coldim()),
        rest(size)
    {
        A.local_colrange(first, last);
        ++last;
    }

    ~local_iterator(){}

    int64_t begin() { return first; } // index of first element in global tensor
    int64_t end() { return last; } // index of end element in global tensor
    bool has_next() { return rest; }
    void next() { --rest; }

    // return current index in local processor
    int64_t local() { return size-rest; } 

    // return current index in global tensor 
    int64_t global() { return first + size-rest; } 

    // given a global index, return local index
    /* 7/Aug/09 Takahiro
        what do processors return if they don't have the element of given index
    */
    //int64_t local_at(int64_t i) { return i-first; }

    //
    void reset() { rest = size; }

};

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    madness::World world(MPI::COMM_WORLD);

    redirectio(world); /* redirect a results to file "log.<number of processor>" */
    
    try {
        for (int64_t n=1; n<10; n++) {
            int64_t m = 2*n;
            DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, m);
            /*
            for(int64_t i=ilo; i<=ihi; i++){
                print("global, local", i, i-ilo);
                //A.data()(i-ilo, _) = i;
            }
            */  
            for (local_iterator index(A); index.has_next(); index.next()){
                //print("iterator:global, local", index.global(), index.local());
                A.data()(index.local(), _) = index.global();
            }

            //TestSystolicMatrixAlgorithm<double> a(A, 3333);
            //a.solve();

            world.taskq.add(new TestSystolicMatrixAlgorithm<double>(A, 3333));
            world.taskq.fence();

            int64_t ilo, ihi;
            A.local_colrange(ilo, ihi); // get column range of A 
            if(A.local_size() > 0){
                for (int i=ilo; i<=ihi; i++) {
                    for (int k=0; k<m; k++) {
                        MADNESS_ASSERT(A.data()(i-ilo,k) == i);
                    }
                }
                for(local_iterator i(A); i.has_next(); i.next()) {
                    for (int k=0; k<m; k++){
                        MADNESS_ASSERT(A.data()(i.local(),k) == i.global());
                    }
                }
            }
        }
    }
    catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
    catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    }
    catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }
    
    MPI::Finalize();
}

