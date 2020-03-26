
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/nonlinsol.h>

using namespace madness;

static const double R =1.4;  // bond length
static const double L = 64.0*R; //box size
static const long k = 8;    //wavelet order
static const double thresh = 1e-5; //precision

static double guess(const coord_3d& r){
    const double x=r[0], y=r[1], z=r[2];// from const coord_3d& r
    return exp(-sqrt(x*x +y*y + (z-R/2)*(z-R/2)+1e-8))+
            exp(-sqrt(x*x +y*y + (z+R/2)*(z+R/2)+1e-8));
}

static double V(const coord_3d& r){
    const double x=r[0], y=r[1], z=r[2];// from const coord_3d& r
    return -1.0/sqrt(x*x +y*y + (z-R/2)*(z-R/2)+1e-8)+
            -1.0/sqrt(x*x +y*y + (z+R/2)*(z+R/2)+1e-8);
}

double rifunction(const coord_3d& r){
    return r[2]; //z
}

double iterate_ground(World& world, NonlinearSolver& solver, real_function_3d V, real_function_3d psi, double& eps){




}

