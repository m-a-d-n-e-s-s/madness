#include <world/world.h>
#include <mra/mra.h>
#include <utility>
#include <ctime>
#include <cmath> 
#include <tensor/tensor.h>
#include <ii/systolic.h>

using namespace madness;

template <typename T>
class LocalizeBoys : public SystolicMatrixAlgorithm<T>
{
public:
    LocalizeBoys<T>( DistributedMatrix<T> &M, const std::vector<int>& set, long nmo, int tag,
        const double threash = 1e-9, const double thetamax = 0.5, const bool randomize = true
    ):
        SystolicMatrixAlgorithm<T>(M, tag),
        M(M),
        world(M.get_world()),
        set(set),
        nmo(nmo),
        ndone_iter(0),
        threash(threash),
        thetamax(thetamax),
        randomize(randomize),
        tol(thetamax),
        niter(0),
        nrot(0)
    {
        madness::print("Start boys localization\n");
    }

    DistributedMatrix<T> get_U(){ return M.data()(_, Slice(0, nmo)); }

    void start_iteration_hook(const TaskThreadEnv& env);
    void kernel(int i, int j, T* rowi, T* rowj);
    void end_iteration_hook(const TaskThreadEnv& env);
    bool converged(const TaskThreadEnv& env) const;

private:
    World& world;
    DistributedMatrix<T> M;
    std::vector<int>& set;
    long nmo, ndone_iter;
    volatile int64_t niter;
    int64_t nrot;
    double threash, thetamax, tol, maxtheta;
    bool randomize;
    void drot(T* restrict a, T* restrict b, double sin, double cos); 
    inline T DIP(T* e_ij, T* e_kl);
    inline T inner(const T* a, const T* b ) const;
};
template <typename T>
void LocalizeBoys<T>::start_iteration_hook(const TaskThreadEnv& env)
{
    T sum = 0.0;
    T* tmp_ptr = M.data().ptr();
    for(long i = 0;i < M.local_coldim();i++){
        sum += DIP(tmp_ptr, tmp_ptr);
    }
    if (env.id() == 0) world.gop.sum(sum);


    long ndone_iter = 0; /// number of iteration
    maxtheta = 0.0; /// maximum rotation angle

    madness::print("\t", "iteration %ld, sum=%.4f ndone=%.2e\n", niter, sum, ndone_iter, tol);
}
template <typename T>
void LocalizeBoys<T>::kernel(int i, int j, T* rowi, T* rowj)
{
    if(set[i] != set[j]) return;

    // make rowi and rowj since we're using one-sided jacobi
    T *ui = rowi + 3*nmo;
    T *uj = rowj + 3*nmo;
    T *xi = rowi;
    T *xj = rowj; 
    T *yi = rowi + nmo;
    T *yj = rowj + nmo;
    T *zi = rowi + 2*nmo;
    T *zj = rowj + 2*nmo;
    T ii[] = { inner(ui, xi), inner(ui, yi), inner(ui, zi) };
    T ij[] = { inner(ui, xj), inner(ui, yj), inner(ui, zj) };
    T jj[] = { inner(uj, xj), inner(uj, yj), inner(uj, zj) };

    double g = DIP(ij, jj) - DIP(ij, ii);
    double h = 4.0 * DIP(ij, ij) + 2.0 * DIP(ii, jj) - DIP(ii, ii) - DIP(jj, jj);
    bool doit = false;

    if (h >= 0.0) {
        doit = true;
        h = -1.0;
    }
    double theta = -g / h;

    maxtheta = std::max<double>(std::abs(theta), maxtheta);

    /// restriction
    if (fabs(theta) > thetamax){
        doit = true;
        if (g < 0) theta = -thetamax;
        else theta = thetamax * 0.8;
    }

    // randomize will be implemented here
    // double sij = DIP(ij, ij); // this line will be used by randomize

    if (fabs(theta) >= tol || doit || randomize){
        ndone_iter++;
        
        double c = cos(theta);
        double s = sin(theta);
        drot (&xi, &xj, s, c);
        drot (&yi, &yj, s, c);
        drot (&zi, &zj, s, c);
        drot (&ui, &uj, s, c);
    }
}
template <typename T>
void LocalizeBoys<T>::end_iteration_hook(const TaskThreadEnv& env)
{
    int id = env.id();

    if (id == 0) world.gop.max(maxtheta);
}
template <typename T>
bool LocalizeBoys<T>::converged(const TaskThreadEnv& env) const
{
    int id = env.id();
    if( ndone_iter == 0 && tol == threash){
        if( id == 0) madness::print("\tBoys localization converged in", ndone_iter, "\n");
        return true;
    }
    else if(niter >= 300){
        if( id == 0) madness::print("\tDid not converged in 300 iteration!\n");
        return true;
    }
    else 
        return true;
}
/// rotate matrix using sin and cos
template <typename T>
void LocalizeBoys<T>::drot(T* restrict a, T* restrict b, double sin, double cos)
{
    for ( long i=0; i<nmo; i++ ) {
        T aa = a[i]*cos - b[i]*sin;
        T bb = b[i]*cos + a[i]*sin;
        a[i] = aa;
        b[i] = bb;
    }
}
template <typename T>
inline T LocalizeBoys<T>::DIP(T* ij, T* kl)
{
    return ij[0] * kl[0] + ij[1] * kl[1] + ij[2] * kl[2];
}
template <typename T>
inline T LocalizeBoys<T>::inner(const T* a, const T* b ) const {
    T s=0;
    for(int64_t i=0; i<nmo; i++){
        s += a[i] * b[i];
    }
    return s;
}
        
/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const Vector<double,3>& x) const {
        return x[axis];
    }
};

/// This is an adapter function. returns a distributed and localized set
/// \param mo molecules
/// \param set basis set
template <typename T>
DistributedMatrix<T> plocalize_boys( World & world, const std::vector< Function<T,3> > & mo, const std::vector<int> & set, const double vtol,
        const double thresh = 1e-9,
        const double thetamax = 0.5,
        const bool randomize = true 
        )
{
    //START_TIMER(world);
    const bool doprint = false; /// do print
    long nmo = mo.size(); /// number of molecule

    /** make big matrix such like ...
        +---+-------+-------+-------+
        | I | dip_x | dip_y | dip_z |
        +---+-------+-------+-------+
    */
    DistributedMatrix<T> dip[3] = {
        column_distributed_matrix<T>(world, nmo, nmo), 
        column_distributed_matrix<T>(world, nmo, nmo), 
        column_distributed_matrix<T>(world, nmo, nmo) 
    };

    /// for all axis, make dipole function
    for(int axis=0;axis < 3;axis++){  
        Function<T,3> fdip = FunctionFactory<T,3>(world).functor(SharedPtr< FunctionFunctorInterface<T,3> >(new DipoleFunctor(axis))).initial_level(4);
        dip[axis].copyin( matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true) );
    }

    DistributedMatrix<T> id = idMatrix(dip[0]);
    DistributedMatrix<double> M = concatenate_rows(id, dip[0], dip[1], dip[2]);

    LocalizeBoys<T> lb(M, set, nmo, 3335, thresh, thetamax, randomize );
    lb.solve();

    /// image of usage
    return lb.get_U();
}

int main(int argc, char** argv) {
    initialize(argc, argv);

    madness::World world(MPI::COMM_WORLD);

    redirectio(world); /* redirect a results to file "log.<number of processor>" */
    
    try {
        print("Test of localize_boys.cc\n");
        print("label size time eig_val");
        for (int64_t n=10; n>1; n-=2) {
            
            /// image of usage
            /*
            Tensor<double>U = plocalize_boys(world, mo, set);
            double t = cpu_time();
            */

            //print("result:", n, cpu_time()-t, M.localized_set());
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
    
    finalize();
}

