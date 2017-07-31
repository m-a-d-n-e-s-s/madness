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


This proram simulates the effect of surface solute interaction between a colloid and a hydrogen atom
*/


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
//#include <madness/mra/operator.h>
#include "molecularmask.h"
#include <madness/mra/nonlinsol.h>
#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include <madness/misc/ran.h>
#include <madness/tensor/solvers.h>
#include <ctime>
#include <list>
#include <jacob/molecule.h>
#include <madness/mra/sdf_shape_3D.h>
#include <madness/mra/funcplot.h>
#include <madness/constants.h>
#include <cmath>
#include <vector>
using namespace madness;
using namespace std;
typedef real_function_3d realfunc;
//typedef SeparatedConvolution<double,3> operatorT;
//typedef std::shared_ptr<operatorT> poperatorT;

//compute the distance between two points
inline double distance1(const coord_3d& r, const coord_3d& center){
    double x1 = r[0], y1 = r[1], z1 = r[2];
    double x2 = center[0], y2 = center[1], z2 = center[2];
    double xx = x1-x2;
    double yy = y1-y2;
    double zz = z1-z2;
    return sqrt(xx*xx + yy*yy + zz*zz);
}

double nuclear_charge_function(const coord_3d& r) {
    const double expnt = 100.0;
    const double coeff = pow(1.0/constants::pi*expnt,0.5*3);
    return coeff*exp(-expnt*(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double electronic_charge_function(const coord_3d& r) {
    const double coeff = 1.0/constants::pi;
    return coeff*exp(-2.0*sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double charge_function(const coord_3d& r) {
    return nuclear_charge_function(r) - electronic_charge_function(r);
}


class SurfaceMoleculeInteraction {
private:
    const double& d; //distance between reaction surface and molecular center
    const double& R; //radius of the spheres forming the colloid surface
    const std::vector<double> charge_center; //defines the center of  nuclear charge
    realfunc& rhot; //total charge density
    const double& sigma;  //surface width
    World& world;
    const int& maxiter;
    //compute the reciprocal of a madness function
    template <typename T, int NDIM>
    struct Reciprocal {
        void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
            UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = 1.0/(*_p0));
        }
        template <typename Archive> void serialize(Archive& ar) {}
    };

    //this is used to perform a binary operation
    struct Bop {
        void operator()(const Key<3>& key,
                        real_tensor rfunc,
                        const real_tensor& func1,
                        const real_tensor& func2) const {
            ITERATOR(rfunc,
                     double d = func1(IND);
                     double p = func2(IND);
                     rfunc(IND) = d*p;
                     );
        }

        template <typename Archive>
        void serialize(Archive& ar) {}
    };
public:
    //set coordinates of colloid atoms
    std::vector< madness::Vector<double,3> > colloid_coords()const{
        // const std::vector<double> cc(0.0);// = charge_center;
        std::vector< madness::Vector<double,3> > c(6); //6 spheres on the colloid surface
        double sqrttwo = std::sqrt(2.0);
        double dist= (sqrttwo/2.0)*R;
        //double x = cc[0], y =  cc[1], z = cc[2];
        double x = 0.0, y =  0.0, z = 0.0;
        c[0][0]= x - dist, c[0][1]= y - d - dist, c[0][2] = z;
        c[1][0]= x - dist, c[1][1]= y - d + dist, c[1][2] = z;
        c[2][0]= x + dist, c[2][1]= y - d + dist, c[2][2] = z;
        c[3][0]= x + dist, c[3][1]= y - d - dist, c[3][2] = z;
        c[4][0]= x , c[4][1]= y - d , c[4][2] = z + R;
        c[5][0]= x , c[5][1]= y - d, c[5][2] = z - R;
        return c;
    }
    std::vector<double> colloid_radii()const {
        int nsphere = colloid_coords().size();
        std::vector<double> c(nsphere);
        for(int i=0; i<nsphere; i++)
            c[i] = R;
        return c;
    }

    //make surface charge
    realfunc make_surfcharge(const realfunc& u,const realfunc& surface,const realfunc& volume) const {
        real_derivative_3d Dx = free_space_derivative<double,3>(rhot.world(), 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(rhot.world(), 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(rhot.world(), 2);
        realfunc Revolume = copy(volume);
        Revolume.unaryop(Reciprocal<double,3>());
        coord_3d lo(0.0),hi(0.0);
        lo[0] = -20.0;
        hi[0] = 20.0;
        plot_line("Revolume.dat",401,lo,hi,Revolume);
        realfunc gradu = Dx(u) + Dy(u) + Dz(u);
        const double rfourpi = 1.0/(4.0*constants::pi);
        realfunc ratio_surface_vol = binary_op(Revolume,surface, Bop());
        realfunc scharge = binary_op(gradu,ratio_surface_vol, Bop());
        return scharge.scale(rfourpi);
    }
    //The perturbed potential near the colloid as a result of the presence of the molecular charge distribution
    realfunc perturbed_molecular_pot(const realfunc& surface,const realfunc& volume) const {
        const bool USE_SOLVER = true;
        double tol = std::max(1e-4,FunctionDefaults<3>::get_thresh());
        real_convolution_3d op = madness::CoulombOperator(world, tol*10.0, tol*0.1);
        realfunc U0 = op(rhot);
        coord_3d lo(0.0),hi(0.0);
        lo[0] = -20.0;
        hi[0] = 20.0;
        plot_line("initial_pot.dat",401,lo,hi,U0);
        realfunc charge = make_surfcharge(U0,surface,volume);
        realfunc U = op(charge);
        double unorm = U.norm2();
        if (USE_SOLVER) {
            madness::NonlinearSolver solver;//(5);
            realfunc uvec, rvec;
            if (world.rank() == 0){
                print("\n\n");//for formating output
                madness::print("      Computing the perturbed solute potential         ");
                madness::print("              ______________________           \n ");

                madness::print("iteration ","    "," residue norm2\n");
            }
            realfunc charg = make_surfcharge(U,surface,volume);
            for (int iter=0; iter<maxiter; iter++) {
                uvec = U;
                //coord_3d lo(0.0),hi(0.0);
                //lo[0] = -20.0;
                // hi[0] = 20.0;
                // plot_line("imp_Surfacepot.dat",1001,lo,hi,U);

                rvec = (U -op(charg)).truncate() ;
                realfunc U_new = solver.update(uvec,rvec);
                double err = rvec.norm2();
                if (world.rank()==0)
                    madness::print("  ", iter,"             " , err);
                if (err >0.3*unorm){ U = 0.5*U + 0.5*U_new;
                }
                else
                    U = U_new;
                if (err < 10.0*tol) break;
            }
        }
        return U;
    }

    SurfaceMoleculeInteraction(const double& d, const double& R,const
                               std::vector<double>& charge_center,
                               realfunc& rhot, const double& sigma, World& world,const int& maxiter)
        :d(d),R(R),                      //d and R are in a.u
         charge_center(0.0),
         rhot(rhot),sigma(sigma),world(world),maxiter(maxiter){
        MADNESS_ASSERT(colloid_coords().size()==colloid_radii().size());
    }
};

int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    const int k = 6; // wavelet order

    const double thresh = 1e-6; // truncation threshold

    const double L = 50; // box is [-L,L]

    //const int natom = 6; // number of atoms

    double sigma = 0.5; // Surface width

    const double R = 3.2503287; // radius of colloid sphere

    const double d = 6.61404096; // distance between center of charge and coilloid center

    const int maxiter = 25; // maximum iteration


    // Function defaults

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(2);
    //Creat an object
    //Molecule molecule;
    const std::vector<double> charge_center(0.0); //defines the center of  nuclear charge
    real_function_3d rhot   = real_factory_3d(world).f(charge_function);
    SurfaceMoleculeInteraction SMI(d,R,charge_center,rhot,sigma,world,maxiter);

    real_functor_3d volume_functor(new MolecularVolumeMask(sigma, SMI.colloid_radii(), SMI.colloid_coords()));
    real_functor_3d surface_functor(new MolecularSurface(sigma, SMI.colloid_radii(), SMI.colloid_coords()));
    real_function_3d vol = real_factory_3d(world).functor(volume_functor);
    real_function_3d surface = real_factory_3d(world).functor(surface_functor).truncate_on_project();
    realfunc pot = SMI.perturbed_molecular_pot(surface,vol);
     coord_3d lo(0.0),hi(0.0);
     lo[0] = -20.0;
     hi[0] = 20.0;
     plot_line("colloid_vol.dat",1001,lo,hi,vol);
     plot_line("colloid_surface.dat",1001,lo,hi,surface);
     plot_line("colloid_pot.dat",1001,lo,hi,pot);

    print("the volume is", vol.trace());
    print("the area   is", surface.trace());
    finalize();
    return 0;
}
