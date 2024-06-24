#include <madness/world/MADworld.h>
#include <utility>
#include <madness/tensor/tensor.h>
#include <madness/tensor/systolic.h>

using namespace madness;

template <typename T>
class TestSystolicMatrixAlgorithm : public SystolicMatrixAlgorithm<T> {
    volatile int niter;
public:
    TestSystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag)
        : SystolicMatrixAlgorithm<T>(A, tag)
        , niter(0)
    {
        madness::print("Testing SystolicMatrixAlgorithm ",
                       SystolicMatrixAlgorithm<T>::get_coldim(),
                       SystolicMatrixAlgorithm<T>::get_rowdim());
    }
    
    void kernel(int i, int j, T* rowi, T* rowj) {
        for (int k=0; k < SystolicMatrixAlgorithm<T>::get_rowdim(); ++k) {
            MADNESS_CHECK(rowi[k] == i);
            MADNESS_CHECK(rowj[k] == j);
        }
    }
    
    void start_iteration_hook(const TaskThreadEnv& env) {
        int id = env.id();
        if (id == 0) {
            ++niter;
        }
    }
    
    bool converged(const TaskThreadEnv& env) const {
        if (niter >= 3) {
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


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    redirectio(world);

    try {
        for (int64_t n=1; n<100; ++n) {
            int64_t m = 2*n;
            DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, m);
            int64_t ilo, ihi;
            A.local_colrange(ilo, ihi);
            for (int i=ilo; i<=ihi; ++i) A.data()(i-ilo,_) = i;

            world.taskq.add(new TestSystolicMatrixAlgorithm<double>(A, 3333));
            world.taskq.fence();

            for (int i=ilo; i<=ihi; ++i) {
                for (int k=0; k<m; ++k) {
                    MADNESS_CHECK(A.data()(i-ilo,k) == i);
                }
            }
        }
    }
    catch (const SafeMPI::Exception& e) {
        print(e);
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
    catch (char* s) {
        print(s);
        error("caught a c-string exception");
    }
//    catch (const char* s) {
//        print(s);
//        error("caught a c-string exception");
//    }
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

    world.gop.fence();
    finalize();
}
