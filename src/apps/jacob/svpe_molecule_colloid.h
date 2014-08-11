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
#ifndef SVPE_MOLECULE_COLLOID_H
#define SVPE_MOLECULE_COLLOID_H

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/tensor/solvers.h>
#include <examples/molecularmask.h>
#include <examples/nonlinsol.h>
#include <madness/constants.h>
#include <vector>

using namespace madness;
using namespace std;

//set coordinates of colloid atoms                      
std::vector< madness::Vector<double,3> > colloid_coords(const double&d,const double R,const std::vector<double> cc) {
    int nsphere = 1;
    std::vector< madness::Vector<double,3> > c(nsphere, madness::Vector<double,3>(0.0)); //6 spheres on the colloid surface                                 
    /* double sqrttwo = std::sqrt(2.0);
       double dist= (sqrttwo/2.0)*R */
    double x = cc[0], y =  cc[1], z = cc[2];
/*  c[0][0]= x - R,  c[0][1]= y - d - 2.0*R,   c[0][2] = z; //A  base sphere
    c[1][0]= x ,     c[1][1]= y - d - R,       c[1][2] = z;     //B base sphere
    c[2][0]= x + R,  c[2][1]= y - d + 2.0*R,   c[2][2] = z;  //C  base sphere
    c[3][0]= x ,     c[3][1]= y - d + 3.0*R,   c[3][2] = z;     //D base spheres
    c[4][0]= x ,     c[4][1]= y - d - 2.0*R ,  c[4][2] = z - R; //E bottom sphere
    c[5][0]= x ,     c[5][1]= y - d - 2.0*R,   c[5][2] = z + R;*/ // top sphere
    c[0][0]= x, c[0][1]= y - R - d ,    c[0][2] = z;// for a single sphere
    return c;
}
/*std::vector< madness::Vector<double,3> > colloid_coordss(const double&d,const double R,std::vector<double> cc) {
    std::vector< madness::Vector<double,3> > c(16,madness::Vector<double,3>(0,0)); //16 spheres on the colloid surface 
    double sqrttwo = std::sqrt(2.0);
    double dist= (R/sqrttwo);
    double x = cc[0], y =  cc[1], z = cc[2];
    c[0][0]=  x - dist,           c[0][1]= y - dist,           c[0][2] = z-d-R; //A
    c[1][0]=  x - dist,           c[1][1]= y + dist,           c[1][2] = z-d-R; //B
    c[2][0]=  x + dist,           c[2][1]= y + dist,           c[2][2] = z-d-R;  //C
    c[3][0]=  x + dist,           c[3][1]= y - dist,           c[3][2] = z-d-R;   //D
    c[4][0]=  x + dist,           c[4][1]= y + 2.0*R - dist,   c[4][2] = z-d-R; //C'
    c[5][0]=  x + dist + 2.0*R,   c[5][1]= y + dist,           c[5][2] = z-d-R; //C"
    c[6][0]=  x + 2.0*R + dist ,  c[6][1]= y + 2.0*R + dist,   c[6][2] = z-d-R;//C'" 
    c[7][0]=  x - dist ,          c[7][1]= y + 2.0*R + dist,   c[7][2] = z-d-R; //B'
    c[8][0]=  x - 2.0*R - dist ,  c[8][1]= y + dist,           c[8][2] = z-d-R; //B"  
    c[9][0]=  x - 2.0*R - dist ,  c[9][1]= y + 2.0*R + dist,   c[9][2] = z-d-R; //B"'
    c[10][0]=  x - dist ,         c[10][1]= y - 2.0*R - dist,  c[10][2] = z-d-R; //A'  
    c[11][0]=  x - 2.0*R - dist , c[11][1]= y  - 2.0*R - dist, c[11][2] = z-d-R; //A"
    c[12][0]=  x - 2.0*R - dist , c[12][1]= y - dist,          c[12][2] = z-d-R; //A"' 
    c[13][0]=  x + 2.0*R + dist , c[13][1]= y - dist,          c[13][2] = z-d-R; //D'
    c[14][0]=  x + 2.0*R + dist , c[14][1]= y  - 2.0*R - dist, c[14][2] = z-d-R; //D"
    c[15][0]=  x + dist ,         c[15][1]= y - 2.0*R - dist,  c[15][2] = z-d-R; //D'" 
    return c;
}*/

//colloid radii
//if number of colloid spheres changes don't forget to change it here
std::vector<double> colloid_radii(const double& R) {
    int nsphere = 1; //number of colloid spheres
    std::vector<double> c(nsphere,0.0);
    for(int i=0; i<nsphere; i++)
        c[i] = R;
    return c;
}

class SVPEColloidSolver {
    World& world;
    const double thresh;
    const double minlen;
    const double sigma;                 //< Width of surface layer
    const double epsilon_0;             //< Interior dielectric
    const double epsilon_1;             //< Exterior dielectric
    real_convolution_3d op;             //< Coulomb operator (1/r ... no 4pi)
    //std::vector<real_convolution_3d_ptr> gop; //< Gradient of the Coulomb operator
    vector_real_function_3d dlog; //< Log-derivative of the dielectric
    real_function_3d rdielectric; //< Reciprocal of the dielectric
    real_function_3d volume; //< volume function of the sphere/colloid
    double L;                            //<half lenght of the simulation volume
    //    vector_real_function_3d Eind; //<surface induced electric 
    const int maxiter;
    static const double cutrho = 1e-12;  //cutoff value of the surface charge
    //this is used to perform a pointwise multiplication
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
    real_function_3d func_pdt(const real_function_3d& func1,real_function_3d& func2) const {
        return binary_op(func1,func2, Bop());
    }
public:
    SVPEColloidSolver(World& world,
                      double sigma, double epsilon_0, double epsilon_1, 
                      const vector_real& atomic_radii, const vector_coord_3d& atomic_coords,
                      const double minlen,double L, const int maxiter)
        :world(world) 
        ,thresh(FunctionDefaults<3>::get_thresh())
        , minlen(minlen)
        , sigma(sigma)
        , epsilon_0(epsilon_0)
        , epsilon_1(epsilon_1)
        , op(CoulombOperator(world, minlen, thresh))
        , dlog(3)
        ,L(L)
        , maxiter(maxiter)
    {
        MADNESS_ASSERT(atomic_radii.size() == atomic_coords.size());//check on the consistency
        // Functors for mask related quantities
        real_functor_3d rdielectric_functor(new MolecularVolumeExponentialSwitchReciprocal(sigma,epsilon_0,epsilon_1,atomic_radii,atomic_coords));
        real_functor_3d gradx_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,0));
        real_functor_3d grady_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,1));
        real_functor_3d gradz_functor(new MolecularVolumeExponentialSwitchLogGrad(sigma, epsilon_0, epsilon_1, atomic_radii, atomic_coords,2));
        real_functor_3d volume_functor(new MolecularVolumeMask(sigma,atomic_radii, atomic_coords));
        // Make the actual functions
        // const double rfourpi = 1.0/(4.0*constants::pi);
        rdielectric = real_factory_3d(world).functor(rdielectric_functor).nofence();
        volume = real_factory_3d(world).functor(volume_functor).nofence();
        dlog[0] = real_factory_3d(world).functor(gradx_functor).nofence().truncate_on_project();
        dlog[1] = real_factory_3d(world).functor(grady_functor).nofence().truncate_on_project();
        dlog[2] = real_factory_3d(world).functor(gradz_functor).truncate_on_project(); // FENCE
        // scale(world, dlog, rfourpi);
        rdielectric.truncate(false);
        truncate(world, dlog);
    }

    // Given the full Coulomb potential computes the surface charge
    real_function_3d make_surface_charge(const real_function_3d& u) const {
        real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);
        //double fac = -1.0/(4.0*constants::pi);
        real_function_3d dx =  Dx(u);
        real_function_3d dy =  Dy(u);
        real_function_3d dz =  Dz(u);
        real_function_3d sc = func_pdt(dlog[0],dx).truncate() + func_pdt(dlog[1],dy).truncate() + func_pdt(dlog[2],dz).truncate();
        //coord_3d lo,hi;
        //lo[1]=-L, hi[1]=L;
        //plot_line("colloid_surf_charge.dat", 1001, lo, hi, sc);
        return sc; //.scale(fac);
    }
    // Given the full Laplace potential compute the surface charge
    real_function_3d make_Laplace_surface_charge(const real_function_3d& u,std::vector<double> E=std::vector<double>(3,0.0)) const {
        real_derivative_3d Dx = free_space_derivative<double,3>(world, 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(world, 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(world, 2);
        //double fac = -1.0/(4.0*constants::pi);
        real_function_3d dx = -E[0] + Dx(u);
        real_function_3d dy = -E[1] + Dy(u);
        real_function_3d dz = -E[2] + Dz(u);
        real_function_3d sc = func_pdt(dlog[0],dx).truncate() + func_pdt(dlog[1],dy).truncate() + func_pdt(dlog[2],dz).truncate();
        return sc;
    }
    //computes components of the the electric field due to the surface charge 
    vector_real_function_3d make_electric_field(const real_function_3d& u) const {
        // throw "not correct";
        vector_real_function_3d E(3);// = madness::vector_factory(0.0, 0.0 ,0.0);
        real_derivative_3d Dx = free_space_derivative<double,3>(u.world(), 0);
        real_derivative_3d Dy = free_space_derivative<double,3>(u.world(), 1);
        real_derivative_3d Dz = free_space_derivative<double,3>(u.world(), 2);
        double fac = -1.0/(4.0*constants::pi);
        real_function_3d Sigmax =(Dx(u)).scale(fac), Sigmay = (Dy(u)).scale(fac),Sigmaz = (Dz(u)).scale(fac); //excess charge on colloid surface
        coord_3d lo(0.0),hi(0.0);
        lo[1]=-L, hi[1]=L;
        plot_line("ab_sigma.dat", 1001, hi, lo, Sigmax, Sigmay, Sigmaz);
        E[0] = Dx(op(Sigmax)), E[1] = Dy(op(Sigmay)), E[2] = Dz(op(Sigmaz)) ;
        return E;
    }

    /// Solve for the full Coulomb potential using the free-particle GF
    real_function_3d solve(const real_function_3d& rho) const {
        const double fac = -1.0/(4.0*constants::pi);
        real_function_3d charge = (rdielectric*rho).scale(-1.0);

        // Initial guess is constant dielectric        
        real_function_3d u0 = op(charge);//.truncate();
        //real_function_3d u = uguess.is_initialized() ? uguess : u0;
        real_function_3d u = u0; // op(rho);  ?? RJH
        double unorm = u.norm2();///volume.trace();
        if (world.rank()==0)
            print("UNORM IS", unorm);
        NonlinearSolver solver;
        if (world.rank()==0){
            print("\n\n");//for formating output 
            madness::print("              Computing the Perturbed Potential (due to molecule)        ");
            madness::print("                         Near the Colloid Surface                         ");
            madness::print("                          ______________________                       \n ");
        }
        for (int iter=0; iter<maxiter; iter++) {
            double start = wall_time();
            real_function_3d surface_charge = make_surface_charge(u);
            real_function_3d r = (u - u0 - op(surface_charge).scale(fac)).truncate(.032*FunctionDefaults<3>::get_thresh());
            double sigtot = surface_charge.trace()*fac;
            //surface_charge.clear();
            // u0.clear();
            real_function_3d unew = solver.update(u,r);
            double change = (unew-u).norm2();//volume.trace(); //times 1000 when not in solvent
            if(world.rank()==0){
                print("iter", iter, "change", change,
                      "soln(10.0)", u(coord_3d(10.0)),
                      "surface charge", sigtot,"used",wall_time()-start);
               
            }
            // Step restriction 
            if (change > 0.3*unorm) 
                u = 0.5*unew + 0.5*u;
            else 
                u = unew;
            //if (change < std::max(1e-4,10.0*thresh/volume.trace())) break;
            if (change < 10.0*thresh) break;
        }
        if(world.rank()==0)
            print("\n\n");//format printing
        real_function_3d urxn = op(make_surface_charge(u)).scale(fac);//u - op(rho);
        coord_3d lo,hi;
        lo[1]=-L, hi[1]=L;
        plot_line("colloid_rxn_pot.dat", 1001, lo, hi, urxn);
        plot_line("colloid_surfcharge.dat", 1001, lo, hi, make_surface_charge(urxn).scale(fac));
        //real_tensor cell(3,2);
        // cell(_,0) = -L;
        //cell(_,1) =  L;
        // plotdx(volume, "testpot.dx", cell);
        return urxn;
    }
    /// Solve for the full Coulomb potential using the free-particle GF
    real_function_3d solve_Laplace(std::vector<double>E) const {
        // Initial guess is constant dielectric        
        const double fac = -1.0/(4.0*constants::pi);
        real_function_3d u(world);//  guess pot is zero;
        double unorm = 0.0;//u.norm2();
        NonlinearSolver solver;
        //print for formating
        if (world.rank()==0){
            print("\n\n");//for formating output 
            madness::print("            Computing the Perturbed Potential (due to external field)           ");
            madness::print("                           ______________________                            \n ");
        }
        for (int iter=0; iter<maxiter; iter++) {
            double start = wall_time();
            real_function_3d surface_charge = make_Laplace_surface_charge(u,E);
            real_function_3d r = (u - op(surface_charge).scale(fac)).truncate(.032*FunctionDefaults<3>::get_thresh());
            double sigtot = surface_charge.trace()*fac;
            //surface_charge.clear();
            real_function_3d unew = solver.update(u, r);
            double change = (unew-u).norm2();//volume.trace();///(8.0*std::pow(L,3.0));
            if (world.rank()==0){
                print("iter", iter, "change", change,
                      "soln(10.0)", u(coord_3d(10.0)),
                      "surface charge", sigtot,"used",wall_time()-start);
            }
            // Step restriction 
            if (change > 0.3*unorm) 
                u = 0.5*unew + 0.5*u;
            else 
                u = unew;
            
            if (change < std::min(1e-4,10.0*thresh)) break;
        }
        if (world.rank()==0)
            print("\n\n");
        coord_3d lo,hi;
        lo[1]=-L, hi[1]= L;
        plot_line("laplace_surfcharge.dat", 1001, lo, hi, make_Laplace_surface_charge(u,E).scale(fac));
        plot_line("laplace_pot.dat", 1001, lo, hi, u);
        return u;
    }
    //calculate the average reaction field(\int C(r)F_r(r)d \tau/\int C(r)d\tau
    //the mask is that of the molecule because the average field is that felt by the molecule
    double ave_rxn_field(const real_function_3d& u,const real_function_3d& mask)const {
        real_function_3d  pdtx = mask*make_electric_field(u)[0]; 
        real_function_3d  pdty = mask*make_electric_field(u)[1]; 
        real_function_3d  pdtz = mask*make_electric_field(u)[2]; 
        double numx = pdtx.trace(), numy = pdty.trace(), numz = pdtz.trace();
        double denominator = mask.trace();
        double Favx = numx/denominator, Favy = numy/denominator,Favz = numz/denominator;
        return std::sqrt(std::pow(Favx,2.0) + std::pow(Favy,2.0) + std::pow(Favz,2.0));
    }
};

#endif

