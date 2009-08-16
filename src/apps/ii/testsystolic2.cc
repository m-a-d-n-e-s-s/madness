#include <world/world.h>
#include <utility>
#include <tensor/tensor.h>
#include <mra/systolic.h>

using namespace madness;


template <typename T>
class TestSystolicMatrixAlgorithm : public SystolicMatrixAlgorithm<T> { 
    volatile int niter; 
    
    public:
    TestSystolicMatrixAlgorithm(DistributedMatrix<T>& A, int tag) 
        : SystolicMatrixAlgorithm<T>(A, tag) 
        , niter(0)
    {
        madness::print("Testing SystolicMatrixAlgorithm: col ", 
                       SystolicMatrixAlgorithm<T>::get_coldim(), 
                       ", row ",
                       SystolicMatrixAlgorithm<T>::get_rowdim());
    }
    
    /* impliment kernel */
    /* this time, check row i,j element for all rows */

    void kernel(int i, int j, T* rowi, T* rowj) {
        for (int k=0; k < SystolicMatrixAlgorithm<T>::get_rowdim(); k++) {
            MADNESS_ASSERT(rowi[k] == i);
            MADNESS_ASSERT(rowj[k] == j);
        }
        print("columns:", i, j);
    }

    void start_iteration_hook(const TaskThreadEnv& env) {
        int id = env.id();
        if (id == 0) {
            niter++;
        }
    }

    bool converged(const TaskThreadEnv& env) const {
        if (niter >= 1) { 
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

/// make identity matrix with same propaties of A
template <typename T>
DistributedMatrix<T> idMatrix(const DistributedMatrix<T>& A){
    int64_t n, m, coltile, rowtile;
    n = A.coldim();
    m = A.rowdim();
    coltile = A.coltile();
    rowtile = A.rowtile();
    MADNESS_ASSERT(n==m);
    DistributedMatrix<T> result(A.get_world(), n, m, coltile, rowtile );
    
    int64_t ilo, ihi;
    result.local_colrange(ilo,ihi);
    for(int64_t i=0; i<=(ihi-ilo); ++i) {
        result.data()(i, i+ilo) = 1;
    }
    
    return result;
}

template <typename T>
class SystolicEigensolver : public SystolicMatrixAlgorithm<T> {
private:
    DistributedMatrix<T>& AV;
    World& world;
    volatile int niter;
    int nrot, nrotsum; 
    static const T tolmin = (T)5.0e-16; //threshld
    T tol, maxd, maxdaij;
    T *ai, *aj, *vi, *vj;
    int size;
    
public:
    SystolicEigensolver<T>(DistributedMatrix<T>& AV, int tag):
        SystolicMatrixAlgorithm<T>( AV, tag ),
        AV(AV),
        world(AV.get_world()),
        niter(0),
        nrot(0),
        nrotsum(0),
        tol((T)1e-2),
        maxd(0),
        maxdaij(1e3), /// just I want a very big value
        size(AV.rowdim()/2)
    {
        MADNESS_ASSERT(AV.is_column_distributed());
        MADNESS_ASSERT(AV.coldim()*2 == AV.rowdim());
        print("One-sided Jacobi start");
        
    }

    void kernel(int i, int j, T* rowi, T* rowj) {
        /// get elements of A and V from concatenated row
        //print("at iteration kernel", i,j);
        ai = rowi;
        vi = rowi+size;
        aj = rowj;
        vj = rowj+size;

        T aii = inner(vi, ai);
        T ajj = inner(vj, aj);
        T aij = inner(vi, aj);
        T daij = fabs(aij);

        maxd = std::max<T>( std::max<T>(fabs(aii),fabs(ajj)), maxd );
        maxdaij = std::max<T>( maxdaij, daij/maxd );

        if( daij < maxd*tol ) return;

        T s = ajj-aii;
        T ds = fabs(s);
          
        if( daij > ds*tolmin ){
            nrot++;
            T c,t,u;
            /// make prameters of rotation matrix
            if( tolmin*daij > ds ) c = s = 1/sqrt(2.0);
            else{
                t = aij/s;
                u = 0.25/sqrt(0.25+t*t);
                c = sqrt(0.5+u);
                s = 2.0*t*u/c;
            }

            /// update all elements
            for (int k=0; k < size; k++) {
                t = ai[k];
                u = aj[k];
                ai[k] = c*t - s*u;
                aj[k] = c*u + s*t;

                t = vi[k];
                u = vj[k];
                vi[k] = c*t - s*u;
                vj[k] = c*u + s*t;
            }
        }
    }

    void start_iteration_hook(const TaskThreadEnv& env) {
        int id = env.id();
        //print("check1");

        if (id == 0) world.gop.max(maxdaij);
        env.barrier();
        //print("check2");

        if (id == 0){
            tol = std::min( tol, std::min(maxdaij*0.1, maxdaij*maxdaij) ); 
            tol = std::max( tol, tolmin );
            niter++;
            //world.gop.max(tol);
            //env.barrier();

            maxdaij = 0;
            nrot = 0; // don't move this line to above
        }
    }

    void end_iteration_hook(const TaskThreadEnv& env) {

        int id = env.id();

        if (id == 0) world.gop.sum(nrot);
        env.barrier();
        //env.barrier();

        nrotsum += nrot;
    }

    bool converged(const TaskThreadEnv& env) const {
        int id = env.id();
        //if (id == 0) world.gop.sum(nrot);
        //env.barrier();
        if ( id == 0 ){
            print("nrot: ", nrot);
            print("tol: ", tol);
        }
        
        if (niter >= 50) {
            if (id == 0) {
                madness::print("    Did not converged in 50 iteration!");
            }
            print("");
            return true;
        }
        else if(nrot == 0 && tol <= tolmin){
            if (id == 0) {
                madness::print("    Converged! ", AV.rowdim()/2);
            }
            print("");
            return true;
        }
        else{
            print("");
            return false;
        }
    }

    //DistributedMatrix<T>
    void get_eval() const{
        for(int64_t i=0; i<size; i++){
            for(int64_t j=0; j<size; j++){
                Tensor<T> ai= AV.data()(i, Slice(0, size-1)); 
                Tensor<T> vj= AV.data()(j, Slice(size, -1));
                T aij = inner(vj, ai);
    
                if(i!=j){
                    print(inner(vj,ai));
                    MADNESS_ASSERT(inner(vj,ai)<=tolmin*10e2);
                }
                else print("eigen value: ", aij);
            }
        }
    }

    DistributedMatrix<T> get_evec() const{
        int64_t m = AV.local_coldim();
        int64_t n = AV.local_rowdim()/2;
        DistributedMatrix<T> result = column_distributed_matrix<T>(world, m, n);
        result.data() = AV.data()(_,Slice(size,-1));

        return result;
    }

private:
    T inner(const T* a, const T* b ) const{
        T s=0;
        for(int64_t i=0; i<size; i++){
            s += a[i] * b[i];
        }
        return s;
    }

    T inner(const Tensor<T>& a, const Tensor<T>& b) const{
        T s=0;
        for(int64_t i=0; i<size; i++) {
            s += a[i] * b[i];
        }

        return s;
    } 
};
/* trial function.
template <typename T>
//SystolicEigensolver<T> 
void systolic_eigensolver (DistributedMatrix<T>& A, int tag )
{
    MADNESS_ASSERT(A.is_column_distributed() == true);
    /// initialize V as identity matrix of size(A)
    DistributedMatrix<T> V = column_distributed_matrix<T>( A.get_world(), A.coldim(), A.rowdim() );

    int64_t ilo, ihi;
    V.local_colrange(ilo, ihi);
    for(int i=ilo; i<=ihi; i++){
        V.data()(i-ilo,i) = 1.0;
    }
  
    DistributedMatrix<T> A_V = concatenate_rows(A,V);

    //print("matrix A_V");
    //print(A_V.data());

    A.get_world().taskq.add(new SystolicEigensolver<T>(A_V, tag));
    A.get_world().taskq.fence();

}
*/
int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    madness::World world(MPI::COMM_WORLD);

    redirectio(world); /* redirect a results to file "log.<number of processor>" */
    
    try {
        print("Test of testsystolic2.cc\n");
        for (int64_t n=8; n>1; n-=1) {
            DistributedMatrix<double> A = column_distributed_matrix<double>(world, n, n);

            int64_t ilo, ihi, jlo, jhi;
            A.local_rowrange(ilo, ihi); // local row range is equal to global row range
            A.local_colrange(jlo, jhi); /* get column range of A */
            for (int j=jlo; j<=jhi; j++) {
                for (int i=ilo; i<=ihi; i++) {
                    A.data()(j-jlo, i-ilo) = (i + j) * 0.1 ; //in this way, matrix is symmetric
                }
                A.data()(j-jlo,j) = ( A.data()(j-jlo,j)+0.5 ) * n;  
            }

            DistributedMatrix<double>  V = idMatrix(A);
            DistributedMatrix<double> AV = concatenate_rows(A, V);
            SystolicEigensolver<double> sA(AV, 3334);

            sA.solve();

            /*
            print("matrix A");
            print(A.data());
            
            DistributedMatrix<double> B = column_distributed_matrix<double>(world, n, l);
            B.local_rowrange(ilo, ihi);
            B.local_colrange(jlo, jhi); // get column range of B
            for (int i=ilo; i<=ihi; i++) {
                for (int j=jlo; j<=jhi; j++) {
                    B.data()(j-jlo,i-ilo) = j + i*100 ;
                }
            }
            print("matrix B");
            print(B.data());

            DistributedMatrix<double> D = concatenate_rows(A, B);

            print("matrix D");
            print(D.data()); //... O.K

            world.taskq.add(new TestSystolicMatrixAlgorithm<double>(C, 3333));
            world.taskq.fence();

            world.taskq.add(new SystolicEigensolver<double>(A, 3334));
            world.taskq.fence();
            */
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

