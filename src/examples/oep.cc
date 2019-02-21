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
*/

//#define WORLD_INSTANTIATE_STATIC_TEMPLATES

/*!
  \file examples/oep.cc
  \brief optimized effective potentials for DFT
*/

#include <chem/nemo.h>
#include <chem/cheminfo.h>
#include <chem/SCFOperators.h>
#include <chem/projector.h>

using namespace madness;

struct dens_inv{

    /// @param[out] U   result
    /// @param[in]  t   numerator
    /// @param[in]  inv density to be inverted >0.0
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& t,
            const Tensor<double>& inv) const {
        ITERATOR(
            U, double d = t(IND);
            double p = std::max(inv(IND),1.e-8);
            U(IND) = d/p;
        );
   }
    template <typename Archive>
    void serialize(Archive& ar) {}

};

struct binary_munge{

    /// @param[out] U   result
    /// @param[in]  r   refdensity
    /// @param[in]  f   function to be mungend
    void operator()(const Key<3>& key, Tensor<double>& U, const Tensor<double>& refdens,
            const Tensor<double>& f) const {
        ITERATOR(
            U, double r = refdens(IND);
            double ff = f(IND);
            U(IND) = (r>1.e-8) ? ff : 0.0;
        );
    }

    template <typename Archive>
    void serialize(Archive& ar) {}

};

/// simple structure to take the pointwise logarithm of a function, shifted by +14
struct logme{
    typedef double resultT;
    struct logme1 {
        double operator()(const double& val) {return log(std::max(1.e-14,val))+14.0;}
    };
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& val) const {
        Tensor<double> result=copy(val);
        logme1 op;
        return result.unaryop(op);
    }

    template <typename Archive>
    void serialize(Archive& ar) {}
};

class OEP : public Nemo {

public:
    OEP(World& world, const std::shared_ptr<SCF> calc) : Nemo(world,calc) {}

    real_function_3d slater_potential(const vecfuncT& nemo) const {
        Exchange K(world,this,0);
        vecfuncT Knemo=K(nemo);
        real_function_3d numerator=R_square*dot(world,nemo,Knemo);
        real_function_3d rho=R_square*dot(world,nemo,nemo);
        real_function_3d slater=-1.0*binary_op(numerator,rho,dens_inv());
        save(slater,"Slater");
        return slater;
    }


    /// compute the kinetic energy potential using the Kohut trick Eq. (30)
    real_function_3d kinetic_energy_potential2(const vecfuncT& nemo) const {

        // the density
        real_function_3d rhonemo=dot(world,nemo,nemo);
        real_function_3d rho=R_square*rhonemo;

        // compute tauL: Eq. (20) of Kohut
        vecfuncT Rnemo=R*nemo;
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_Rnemo=Laplace(Rnemo);
        real_function_3d tauL=-0.5*dot(world,laplace_Rnemo,Rnemo);
        save(tauL,"tauL");
        real_function_3d tauL_div_rho=binary_op(tauL,rho,dens_inv());
        save(tauL_div_rho,"tauL_div_rho");

        // compute tau = | grad(mo) |^2
        //  = dR.dR * |nemo|^2 + 2*dR.grad(nemo) + R^2*|grad(nemo)|^2
        real_function_3d tau=real_factory_3d(world).compressed();
        vecfuncT dR=this->nuclear_correlation->U1vec();

        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();

        for (const real_function_3d& n : nemo) {
            vecfuncT gradnemo=grad(n);
            tau+=U1dot*n*n;
            tau-=2.0*n*dot(world,dR,gradnemo);  // grad(R) = -U1
            tau+=dot(world,gradnemo,gradnemo);
        }
        tau=0.5*tau*R_square;
        real_function_3d tau_div_rho=binary_op(tau,rho,dens_inv());
        save(tau_div_rho,"tau_div_rho");

        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rho,logme());
        save(logdensa,"logdensa");
        vecfuncT gradzeta=grad(logdensa);
        madness::save_function(gradzeta,"gradzeta");
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        save(d2rho_div_rho,"d2rho_div_rho");
        real_function_3d result=tau_div_rho-0.25*d2rho_div_rho;
        save(result,"kinetic2");
        return result;
    }


    real_function_3d kinetic_energy_potential(const vecfuncT& nemo) const {

        const Nuclear U_op(world,this->nuclear_correlation);
        const Nuclear V_op(world,this->get_calc().get());

        const vecfuncT Vnemo=V_op(nemo);  // eprec is important here!
        const vecfuncT Unemo=U_op(nemo);

        // nabla^2 nemo
        Laplacian<double,3> Laplace(world,0.0);
        vecfuncT laplace_nemo=Laplace(nemo);

        vecfuncT tmp=Unemo-this->nuclear_correlation->U2()*nemo;
        laplace_nemo-=2.0*tmp;
        vecfuncT D2Rnemo=R*laplace_nemo;

//        // double check result: recompute the density from its laplacian
//        vecfuncT nemo_rec=apply(world,*poisson,D2Rnemo);
//        scale(world,nemo_rec,-1./(4.*constants::pi));
        vecfuncT Rnemo=mul(world,R,nemo);
//        vecfuncT diff=sub(world,Rnemo,nemo_rec);
//        double dnorm=norm2(world,diff);
//        print("dnorm of laplacian phi ",dnorm);

        // compute \sum_i \phi_i \Delta \phi_i
        real_function_3d phiD2phi=dot(world,Rnemo,D2Rnemo);
        save(phiD2phi,"phiD2phi");

        // compute \sum_i \phi_i \epsilon_i \phi_i
        vecfuncT R2nemo=R_square*nemo;
        const real_function_3d rho=2.0*dot(world,nemo,R2nemo);
        const real_function_3d rhonemo=2.0*dot(world,nemo,nemo);

        // compute T1
        double T1 = 0.0;
        for (int axis = 0; axis < 3; axis++) {
            real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
            const vecfuncT dnemo = apply(world, D, Rnemo);
            T1 += 2.0* 0.5 * (inner(world, dnemo, dnemo)).sum();
        }
        printf("T1 %16.8f \n",T1);

        std::vector<double> eps(nemo.size());
        for (std::size_t i=0; i<eps.size(); ++i) eps[i]=calc->aeps(i);
        scale(world,R2nemo,eps);
        real_function_3d phiepsilonphi=dot(world,R2nemo,nemo);

        // divide by the density
        real_function_3d numerator=2.0*(-0.5*phiD2phi);   // fac 2 for sum over spin orbitals
//        real_function_3d numerator=2.0*(-0.5*phiD2phi-phiepsilonphi);   // fac 2 for sum over spin orbitals
        real_function_3d kinetic1=binary_op(numerator,rho,dens_inv());
        save(kinetic1,"kinetic1");
        real_function_3d nu_bar=kinetic1 + 2.0*(-0.5)*binary_op(phiepsilonphi,rho,dens_inv());

        // reintroduce the nuclear potential *after* smoothing
        real_function_3d uvnuc=calc->potentialmanager->vnuclear()-nuclear_correlation->U2();
        nu_bar=nu_bar-uvnuc;

        // compute T2
        vecfuncT dipole(3), drho(3);
        dipole[0]=real_factory_3d(world).functor(MomentFunctor(1,0,0));
        dipole[1]=real_factory_3d(world).functor(MomentFunctor(0,1,0));
        dipole[2]=real_factory_3d(world).functor(MomentFunctor(0,0,1));
        drho[0]=make_ddensity(rhonemo,0);
        drho[1]=make_ddensity(rhonemo,1);
        drho[2]=make_ddensity(rhonemo,2);
        real_function_3d one=real_factory_3d(world).functor(MomentFunctor(0,0,0));
        real_function_3d arg=3.0*rho+dot(world,dipole,drho);
        double T2=0.5*inner(nu_bar,arg);
        printf("T2 %16.8f \n",T2);

        Coulomb J(world,this);
        XCOperator xc(world,this);

        real_function_3d vne=calc->potentialmanager->vnuclear();
        real_function_3d vj=J.potential();
        real_function_3d vxc=xc.make_xc_potential();

        real_function_3d euler=nu_bar + vne + vj + vxc;
        double eulernorm=euler.norm2();
        print("eulernorm ",eulernorm);
        save(euler,"euler");
        save(vne,"vne");
        save(vj,"vj");
        save(vxc,"vxc");

        return nu_bar;
    }

    /// compute the laplacian of the density using the log trick
    real_function_3d make_laplacian_density_oep(const real_function_3d& rhonemo) const {


        // U1^2 operator
        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nuclear_correlation.get());
        const real_function_3d U1dot=real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        real_function_3d result=(2.0*U1dot).truncate();

        // U2 operator
        const real_function_3d V=calc->potentialmanager->vnuclear();
        const real_function_3d U2=nuclear_correlation->U2();

        real_function_3d term2=2.0*(U2-V).truncate();
        result-=term2;


        // compute the laplacian of the density with the log trick
        real_function_3d logdensa=unary_op(rhonemo,logme());
        vecfuncT gradzeta=grad(logdensa);

        // zeta contribution
        real_function_3d d2rho_div_rho=div(gradzeta) + dot(world,gradzeta,gradzeta);
        result+=d2rho_div_rho;

        real_function_3d term3=4.0*dot(world,nuclear_correlation->U1vec(),gradzeta);
        result-=term3;

        result=(R_square*rhonemo*result).truncate();
        save(result,"d2rho");

        // double check result: recompute the density from its laplacian
        real_function_3d rho_rec=-1./(4.*constants::pi)*(*poisson)(result);
        save(rho_rec,"rho_reconstructed");

        real_function_3d rho=rhonemo*R_square;
        save(rho,"rho");
        real_function_3d laplace_rho=div(grad(rho));
        save(laplace_rho,"d2rho_direct");


        return result;
    }


};


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    if (world.rank() == 0) {
    	print("\n  OEP -- optimized effective potentials for DFT  \n");
    	printf("starting at time %.1f\n", wall_time());
    }
    startup(world,argc,argv);
    std::cout.precision(6);

    const std::string input="input";
    std::shared_ptr<SCF> calc(new SCF(world,input.c_str()));
    if (world.rank()==0) {
        calc->molecule.print();
        print("\n");
        calc->param.print(world);
    }

    std::shared_ptr<OEP> oep(new OEP(world,calc));
    const double energy=oep->value();

    real_function_3d slaterpot=oep->slater_potential(oep->get_calc()->amo);
    save(slaterpot,"slaterpotential");


//    oep->kinetic_energy_potential(oep->get_calc()->amo);
    oep->kinetic_energy_potential2(oep->get_calc()->amo);

    real_function_3d rhonemo=dot(world,oep->get_calc()->amo,oep->get_calc()->amo);
    oep->make_laplacian_density_oep(rhonemo);


    if (world.rank()==0) {
        printf("final energy   %12.8f\n", energy);
        printf("finished at time %.1f\n", wall_time());
    }

    finalize();
    return 0;
}
