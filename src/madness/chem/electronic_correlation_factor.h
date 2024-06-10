/*
 * electronic_correlation_factor.h
 *
 *  Created on: Jul 9, 2015
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_
#define SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_



#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
#include<madness/chem/molecule.h>
#include <iomanip>

namespace madness {
/// a class holding the electronic correlation factor for R12 theory
class CorrelationFactor {

    World& world;
    double _gamma;      ///< the correlation factor exp(-gamma r12)
    double dcut;        ///< the cutoff for the 1/r potential
    double lo;          ///< smallest length scale to be resolved

public:

    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor(World& world) : world(world), _gamma(-1.0), dcut(1.e-10),
        lo(1.e-10) {
    }

    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor(World& world, const double& gamma, const double dcut,
            const Molecule& molecule) : world(world), _gamma(gamma), dcut(dcut) {
        lo=1.e-6;//lo = molecule.smallest_length_scale();
//        if (world.rank()==0) {
//            if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
//            else if (gamma==0.0) print("constructed linear correlation factor");
//        }
    }
    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor(World& world, const double& gamma, const double dcut,
            const double lo) : world(world), _gamma(gamma), dcut(dcut), lo(lo) {
//        if (world.rank()==0) {
//            if (gamma>0.0) print("constructed correlation factor with gamma=",gamma);
//            else if (gamma==0.0) print("constructed linear correlation factor");
//        }
    }

    /// copy ctor
    CorrelationFactor(const CorrelationFactor& other) : world(other.world) {
        _gamma=other._gamma;
        dcut=other.dcut;
        lo=other.lo;
    }

    /// assignment; assume other's world is this world
    CorrelationFactor& operator=(const CorrelationFactor& other) {
        _gamma=other._gamma;
        dcut=other.dcut;
        lo=other.lo;
        return *this;
    }

    /// return the exponent of this correlation factor
    double gamma() const {return _gamma;}

    /// return the value of the correlation factor
    double operator()(const coord_6d& r) const {
        const double rr=r12(r);
        if (_gamma>0.0) return (1.0-exp(-_gamma*rr))/(2.0*_gamma);
        return 0.5*rr;
    }

    /// apply Kutzelnigg's regularized potential to an orbital product
    real_function_6d apply_U(const real_function_3d& phi_i, const real_function_3d& phi_j,
            const real_convolution_6d& op_mod, const bool symmetric=false) const {
	  if(not op_mod.modified()) MADNESS_EXCEPTION("ElectronicCorrelationFactor::apply_U, op_mod must be in modified_NS form",1);
	  const double thresh = FunctionDefaults<6>::get_thresh();
	  const bool debug = false;
	  if(symmetric) MADNESS_ASSERT((phi_i-phi_j).norm2() < FunctionDefaults<3>::get_thresh());

        real_function_6d result=real_factory_6d(world);

        for (int axis=0; axis<3; ++axis) {
            //if (world.rank()==0) print("working on axis",axis);
            real_derivative_3d D = free_space_derivative<double,3>(world, axis);
            const real_function_3d Di=(D(phi_i)).truncate();
            real_function_3d Dj;
            if(symmetric) Dj=madness::copy(Di);
            else Dj=(D(phi_j)).truncate();

            real_function_6d u=U1(axis);
            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                        .g12(u).particle1(copy(Di)).particle2(copy(phi_j)).thresh(thresh);
            tmp1.fill_cuspy_tree(op_mod).truncate();

            real_function_6d tmp2;
            if(symmetric) tmp2 = -1.0*swap_particles(tmp1);
            else{
            tmp2=CompositeFactory<double,6,3>(world)
                                    .g12(u).particle1(copy(phi_i)).particle2(copy(Dj)).thresh(thresh);
            tmp2.fill_cuspy_tree(op_mod).truncate();
            }

            result=result+(tmp1-tmp2).truncate();


            tmp1.clear();
            tmp2.clear();
            world.gop.fence();
            result.truncate().reduce_rank();
        }

        // include the purely local potential that (partially) cancels 1/r12
        if (_gamma>0.0) {
            fg_ func(_gamma,dcut);
            real_function_6d fg3=real_factory_6d(world).functor(func).is_on_demand();
            real_function_6d mul=CompositeFactory<double,6,3>(world)
                                .g12(fg3).particle1(copy(phi_i)).particle2(copy(phi_j)).thresh(thresh);;
            mul.fill_cuspy_tree(op_mod).truncate();
           // mul.print_size("mul");

            result=(result+mul).truncate().reduce_rank();
        }
        if(debug) result.print_size("Ue|ij>");
        return result;
    }

    /// return the U1 term of the correlation function
    real_function_6d U1(const int axis) const {
        U func(_gamma,axis,dcut);
        const real_function_6d u1=real_factory_6d(world)
                .functor(func).is_on_demand();
        return u1;
    }

    /// return the U1 term of the correlation function
    real_function_6d U2() const {
        if (world.rank()==0) print("U2 for the electronic correlation factor");
        if (world.rank()==0) print("is expensive -- do you really need it??");
        MADNESS_EXCEPTION("U2() not implemented, since it might be expensive",1);
        return real_factory_6d(world);
    }

    /// return the correlation factor as on-demand function
    real_function_6d f() const {
//        real_function_6d tmp=real_factory_6d(world).functor2(*this).is_on_demand();
        double thresh=FunctionDefaults<3>::get_thresh();
        real_function_6d tmp=TwoElectronFactory(world)
                .dcut(dcut).gamma(_gamma).f12().thresh(thresh);
        return tmp;
    }

    /// return f^2 as on-demand function
    real_function_6d f2() const {
        f2_ func(_gamma);
        real_function_6d tmp=real_factory_6d(world).functor(func).is_on_demand();
        return tmp;
    }

    /// return fg+sth as on-demand function
    real_function_6d fg() const {
        fg_ func(_gamma,dcut);
        real_function_6d tmp=real_factory_6d(world).functor(func).is_on_demand();
        return tmp;
    }

    /// return f/r as on-demand function
    real_function_6d f_over_r() const {
        f_over_r_ func(_gamma,dcut);
        real_function_6d tmp=real_factory_6d(world).functor(func).is_on_demand();
        return tmp;
    }

    /// return (\nabla f)^2 as on-demand functions
    real_function_6d nablaf2() const {
        nablaf2_ func(_gamma);
        real_function_6d tmp=real_factory_6d(world).functor(func).is_on_demand();
        return tmp;
    }

private:
    /// functor for the local potential (1-f12)/r12 + sth (doubly connected term of the commutator)

    /// TODO: turn this into coeffs directly
    struct fg_ : FunctionFunctorInterface<double,6> {
        double gamma;
        double dcut;
        fg_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
            MADNESS_ASSERT(gamma>0.0);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            return (1.0-e)*u(rr,dcut) + 0.5*gamma*e;
        }
    };

    /// functor for the local potential (1-f12)/r12
    struct f_over_r_ : FunctionFunctorInterface<double,6>  {
        double gamma;
        double dcut;
        f_over_r_(double gamma, double dcut) : gamma(gamma), dcut(dcut) {
            MADNESS_ASSERT(gamma>0.0);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            return (1.0-e)*u(rr,dcut)/(2.0*gamma);
        }
    };

    /// functor for the local part of the regularized potential: f12/r12*(r1-r2)(D1-D2)
    struct U : FunctionFunctorInterface<double,6>  {
        double gamma;
        int axis;
        double dcut;
        U(double gamma, int axis, double dcut) : gamma(gamma), axis(axis),
            dcut(dcut) {
            MADNESS_ASSERT(axis>=0 and axis<3);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const coord_3d vr12{r[0]-r[3],r[1]-r[4],r[2]-r[5]};
            const coord_3d N=unitvec(vr12);
            if (gamma>0.0) return -0.5*exp(-gamma*rr)*N[axis];
            MADNESS_EXCEPTION("no gamma in electronic corrfac::U1",1);
//          const double rr=r12(r);
//            const double g12=u(rr,dcut);
//            double a=0.5;
//            if (gamma>0.0) a=0.5*exp(-gamma*rr);
//            return -a*x12(r,axis) * g12;
        }
    };

    /// functor for the local potential (1-f12)^2
    struct f2_ : FunctionFunctorInterface<double,6>  {
        double gamma;
        f2_(double gamma) : gamma(gamma) {MADNESS_ASSERT(gamma>0.0);}
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double e=exp(-gamma*rr);
            const double f=(1.0-e)/(2.0*gamma);
            return f*f;
        }
    };

    /// functor for the local potential (\nabla f)^2
    struct nablaf2_ : FunctionFunctorInterface<double,6>  {
        double gamma;
        nablaf2_(double gamma) : gamma(gamma) {
            MADNESS_ASSERT(gamma>0.0);
            MADNESS_ASSERT(gamma==1.0);
        }
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const double f=exp(-2.0*gamma*rr)/(4.0*gamma*gamma);
            return f;
        }
    };

    /// Smoothed 1/r potential (c is the smoothing distance)
    static double u(double r, double c) {
        r = r/c;
        double r2 = r*r, pot;
        if (r > 6.5){
            pot = 1.0/r;
        } else if (r > 1e-2) {
            pot = erf(r)/r + exp(-r2)*0.56418958354775630;
        } else{
            pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
        }
        return pot/c;
    }

    static double r12(const coord_6d& r) {
        const double x12=r[0]-r[3];
        const double y12=r[1]-r[4];
        const double z12=r[2]-r[5];
        const double r12=sqrt(x12*x12 + y12*y12 + z12*z12);
        return r12;
    }
    static double x12(const coord_6d& r, const int axis) {
        return r[axis]-r[axis+3];
    }

	static coord_3d smoothed_unitvec(const coord_3d& xyz, double smoothing) {
//        if (smoothing==0.0) smoothing=molecule.get_eprec();
        // TODO:need to test this
        // reduce the smoothing for the unitvector
        //if (not (this->type()==None or this->type()==Two)) smoothing=sqrt(smoothing);
        const double r=xyz.normf();
        const double cutoff=smoothing;
        if (r>cutoff) {
            return 1.0/r*xyz;
        } else {
            const double xi=r/cutoff;
            const double xi2=xi*xi;
            const double xi3=xi*xi*xi;
//            const double nu21=0.5+1./32.*(45.*xi - 50.*xi3 + 21.*xi*xi*xi*xi*xi);
            const double nu22=0.5 + 1./64.*(105* xi - 175 *xi3 + 147* xi2*xi3 - 45* xi3*xi3*xi);
//            const double nu40=0.5 + 1./128.*(225 *xi - 350 *xi3 + 189*xi2*xi3);
            const double kk=2.*nu22-1.0;
            return kk/(r+1.e-15)*xyz;
        }
	}
};

/// a class holding the electronic correlation factor for R12 theory
/// CorrelationFactor2 = (1-0.5*exp(-gamma*r12), gamma=0.5
/// (CorrelationFactor + 1)*2.0 = CorrelationFactor2 (currently CorrelationFactor2 is only implemented for gamma=0.5 so use this gamma also on CorrelationFactor
class CorrelationFactor2 {

    World& world;
    double _gamma;      ///< the correlation factor exp(-gamma r12)
    typedef std::shared_ptr< FunctionFunctorInterface<double,6> > functorT;

public:

    double dcut;        ///< the cutoff for the 1/r potential
    double lo;          ///< smallest length scale to be resolved
    double vtol;        ///< initial projection threshold


    /// ctor, use negative gamma for linear correlation factor r12
    CorrelationFactor2(World& world) : world(world), _gamma(0.5), dcut(1.e-10),
        lo(1.e-10), vtol(FunctionDefaults<3>::get_thresh()*0.1) {
        MADNESS_ASSERT(_gamma==0.5);
    }

    /// return the exponent of this correlation factor
    double gamma() const {return _gamma;}

    real_function_6d function() const {
        functorT R=functorT(new R_functor(_gamma,1));
        return real_factory_6d(world).functor(R).is_on_demand();
    }

    real_function_6d square() const {
        functorT R2=functorT(new R_functor(_gamma,2));
        return real_factory_6d(world).functor(R2).is_on_demand();
    }

    real_function_6d inverse() const {
        functorT R=functorT(new R_functor(_gamma,-1));
        return real_factory_6d(world).functor(R).is_on_demand();
    }

    /// return the U1 term of the correlation function
    real_function_6d U1(const int axis) const {
        functorT U1f=functorT(new U1_functor(_gamma,axis));
        return real_factory_6d(world).functor(U1f).is_on_demand();
    }

    /// return the U2 term of the correlation function
    real_function_6d U2() const {
        functorT U2f=functorT(new U2_functor(_gamma));
        return real_factory_6d(world).functor(U2f).is_on_demand();
    }

    /// apply Kutzelnigg's regularized potential to an orbital product
    real_function_6d apply_U(const real_function_6d& psi, const double eps) const {
        const double bsh_thresh=1.e-7;

        real_function_6d result=real_factory_6d(world);

        real_convolution_6d op_mod = BSHOperator<6>(world, sqrt(-2*eps), lo,bsh_thresh);
        op_mod.modified()=true;

        for (int axis=0; axis<3; ++axis) {
            //if (world.rank()==0) print("working on axis",axis);
            real_derivative_6d D1 = free_space_derivative<double,6>(world, axis);
            real_derivative_6d D2 = free_space_derivative<double,6>(world, axis+3);
            const real_function_6d Drhs1=D1(psi).truncate();
            const real_function_6d Drhs2=D2(psi).truncate();

            const real_function_6d u1=U1(axis);

            real_function_6d tmp1=CompositeFactory<double,6,3>(world)
                                 .g12(u1).ket(copy(Drhs1));
            tmp1.fill_cuspy_tree(op_mod).truncate();

            real_function_6d tmp2=CompositeFactory<double,6,3>(world)
                                 .g12(u1).ket(copy(Drhs2));
            tmp2.fill_cuspy_tree(op_mod).truncate();
           // if (world.rank()==0) print("done with fill_tree");

            result=result+(tmp1-tmp2).truncate();
            tmp1.clear();
            tmp2.clear();
            world.gop.fence();
            result.truncate().reduce_rank();

           // if (world.rank()==0) printf("done with multiplication with U at ime %.1f\n",wall_time());
           // result.print_size("result");
        }

        real_function_6d u2=U2();
        real_function_6d r2=CompositeFactory<double,6,3>(world).ket(copy(psi))
                                    .g12(u2);
        r2.fill_tree(op_mod);
        result=(result+r2).truncate();
        return result;
    }


private:

    /// functor for the correlation factor R
    class R_functor : public FunctionFunctorInterface<double,6> {
        double gamma;
        int exponent;


    public:
        R_functor(double gamma, int e=1) : gamma(gamma), exponent(e) {
            MADNESS_ASSERT(gamma==0.5);
        }

        // only valid for gamma=1
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            double val=(1.0-0.5*exp(-gamma*rr));
            if (exponent==1) return val;
            else if (exponent==2) return val*val;
            else if (exponent==-1) return 1.0/val;
            else {
                MADNESS_EXCEPTION("fancy exponent in correlationfactor2",1);
            }
        }
    };

    /// functor for the U2 local potential
    class U2_functor : public FunctionFunctorInterface<double,6> {
        double gamma;

    public:
        U2_functor(double gamma) : gamma(gamma) {
            MADNESS_ASSERT(gamma==0.5);
        }

        // only valid for gamma=1
        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            // Taylor expansion for small r
            if (rr<1.e-4) { // valid for gamma==0.5, otherwise singular
                return (5./4.0 - rr + (35.0* rr*rr)/48.0 - (101.0*rr*rr*rr)/192.0);
            }
            const double egr=exp(-gamma*rr);
            return -(-8.*egr + 8.0 + rr*egr)/(4.0 *rr*egr - 8 *rr);
        }
    };

    /// functor for the U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}  potential

    /// the potential is given by
    /// U1 = -\frac{\vec\nabla_1 f_{12}}{f_{12}}
    ///    =  \frac{e^{-r12/2}{4-2e^{-r12/2}} \vec unitvec
    /// the derivative operators are not included
    class U1_functor : public FunctionFunctorInterface<double,6> {
        double gamma;
        int axis;

    public:
        U1_functor(double gamma, int axis) : gamma(gamma), axis(axis) {
            MADNESS_ASSERT(gamma==0.5);
            MADNESS_ASSERT(axis<3);
        }

        double operator()(const coord_6d& r) const {
            const double rr=r12(r);
            const coord_3d vr12{r[0]-r[3],r[1]-r[4],r[2]-r[5]};
            const coord_3d N=unitvec(vr12);
            // Taylor expansion for small r
            double val;
            if (rr<1.e-4) { // valid for gamma==0.5, otherwise singular
                val = 0.5 - 0.5*rr + 0.125*(3.*rr*rr) - (13.* rr*rr*rr)/48.0;
            } else {
                const double egr=exp(-gamma*rr);
                val=egr/(4.0-2.0*egr);
            }
            // NOTE the sign
            return -val*N[axis];
        }
    };

    /// helper function
    static double r12(const coord_6d& r) {
        const double x12=r[0]-r[3];
        const double y12=r[1]-r[4];
        const double z12=r[2]-r[5];
        const double r12=sqrt(x12*x12 + y12*y12 + z12*z12);
        return r12;
    }

};

}



#endif /* SRC_APPS_CHEM_ELECTRONIC_CORRELATION_FACTOR_H_ */
