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

/*!
  \file helium_mp2.cc
  \brief Solves the Hartree-Fock and MP2 equations for the helium atom
  \defgroup examplehehf Hartree-Fock and MP2 for the helium atom
  \ingroup examples
*/


//#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/lbdeux.h>

#include <iostream>



using namespace madness;

static const double dcut=1.e-5;
static const double r12cut=1.e-3;
static const double rcut=1.e-3;
static const double shift=0.0;

//static const double L = 32.0;   // box size
//static const long k = 6 ;        // wavelet order
//static const double thresh = 1e-3; // precision
//static const TensorType tt = TT_3D;
//static const long truncate_mode = 0;

//template<typename T>
//static std::string stringify(T arg) {
//	std::ostringstream o;
//	if (!(o << arg))
//		throw std::domain_error("stringify(double)");
//	return o.str();
//}

static double guess(const coord_3d& r) {
    const double x=r[0]+shift, y=r[1]+shift, z=r[2];
    return 6.0*exp(-sqrt(x*x+y*y+z*z+1e-8));
}

static double V(const coord_3d& r) {
    const double x=r[0]+shift, y=r[1]+shift, z=r[2];
    return -1.0/(sqrt(x*x+y*y+z*z+1e-10));
}

static double HO_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return x*x + y*y + z*z;
}

static double HO_6d(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double value=(x1*x1 + y1*y1 + z1*z1) + (x2*x2 + y2*y2 + z2*z2);
    return value;
}


static double gauss_3d(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double r2= x*x + y*y + z*z;
    const double norm=0.712705695388313;
    return norm*exp(-r2);
}

static double gauss_6d(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double r2=(x1*x1 + y1*y1 + z1*z1) + (x2*x2 + y2*y2 + z2*z2);
    const double norm=0.5;
    return norm*exp(-r2);
}


static double HO_vphi_3d(const coord_3d& r) {
    const double v=HO_3d(r);
    const double phi=gauss_3d(r);
    return v*phi;
}


static double Z2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -2.0/(sqrt(x*x+y*y+z*z+dcut*dcut));
}

static double Z1(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return -1.0/(sqrt(x*x+y*y+z*z+dcut*dcut));
}


static double V_1(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
	const double y1=r[1], y2=r[4];
	const double z1=r[2], z2=r[5];


	const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
	const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);

	const double value=-2.0/r1 - 2.0/r2;
	return value;


}

// Smoothed 1/r potential (c is the smoothing distance)
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

void distances(const coord_6d& r, double& r1, double& r2, double& r12) {
    const double x1=r[0], y1=r[1], z1=r[2];
    const double x2=r[3], y2=r[4], z2=r[5];
    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;
    r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
    r2 = sqrt(x2*x2 + y2*y2 + z2*z2);
    r12 = sqrt(xx*xx + yy*yy + zz*zz);
}

static double V(const coord_6d& r) {
    double r1, r2, r12;
    distances(r, r1, r2, r12);
    return -2.0*u(r1,rcut) - 2.0*u(r2,rcut) + 1.0*u(r12,r12cut);
}


static double coul(const coord_6d& r) {

//	// separation for 2-way decomposition (SVD; r1 -- r2)
//	const double x1=r[0], x2=r[3];
//	const double y1=r[1], y2=r[4];
//	const double z1=r[2], z2=r[5];
//
//	const double xx=x1-x2, yy=y1-y2, zz=z1-z2;
//
////	const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
////	const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
//	const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut*1.e6);

	double r1,r2,r12;
	distances(r, r1, r2, r12);
	const double value= 1.0*u(r12,r12cut);
//	const double value= + 1.0/r12;
	return value;
}


static double helium_pot(const coord_6d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
    const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    const double value=-2.0/r1 - 2.0/r2;// + 1.0/r12;
//    const double value= + 1.0/r12;
    return value;
}

// return the helium potential times the hylleraas function
static double V_times_phi(const coord_6d& r) {

	// separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];

    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    const double pot=-2.0/r1 - 2.0/r2;// + 1.0/r12;
    const double phi=exp(-1.8*(r1 + r2))*(1.0 + 0.5*r12);
    const double value=pot*phi;

    return value;
}


// according to McQuarrie
static double h_orbital(const coord_3d& r) {
    return exp(-sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6));
}


// according to McQuarrie
static double he_plus_orbital(const coord_3d& r) {
    return exp(-2.0*sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]+1e-6));
}

// according to Ed
static double he_orbital_3d(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);

    const double val=(0.00995312870402086*(31.05166416452748 - 7.405311261369526*r1 +
           r1*r1)*(4.335496673568937 + 0.24243181498262073*r1 +
           r1*r1)) * exp(-1.81607*r1);

    return val;
}

// according to McQuarrie
static double he_orbital_McQuarrie(const coord_3d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0];
    const double y1=r[1];
    const double z1=r[2];

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);

    const double val=exp(-(27.0/16.0)*r1);

    return val;
}


// according to Ed / McQuarrie
static double he_orbitals(const coord_6d& r) {

    // separation for 2-way decomposition (SVD; r1 -- r2)
    coord_3d r1;
    coord_3d r2;
    r1[0]=r[0];
    r1[1]=r[1];
    r1[2]=r[2];
    r2[0]=r[3];
    r2[1]=r[4];
    r2[2]=r[5];

//    const double val=he_orbital_3d(r1) * he_orbital_3d(r2);
    const double val=he_orbital_McQuarrie(r1) * he_orbital_McQuarrie(r2);

    return val;
}


static double hylleraas_3term(const coord_6d& r) {

    // E = -2.902432028988 E_h

    // separation for 2-way decomposition (SVD; r1 -- r2)
	const double x1=r[0], x2=r[3];
    const double y1=r[1], y2=r[4];
    const double z1=r[2], z2=r[5];
    const double xx=x1-x2, yy=y1-y2, zz=z1-z2;

    const double r1 = sqrt(x1*x1 + y1*y1 + z1*z1 + dcut*dcut);
    const double r2 = sqrt(x2*x2 + y2*y2 + z2*z2 + dcut*dcut);
    const double r12= sqrt(xx*xx + yy*yy + zz*zz + dcut*dcut);

    return -exp(-1.81607*(r1 + r2)) * (
    -1.33083943395992
    -0.388320016632985 * r12
    -0.174093511691879 *  ( r1*r1  + r2*r2  -2 * r1 * r2 )
    );
}


// Hylleraas 3-term minus He orbitals
static double he_correlation(const coord_6d& r) {
	return hylleraas_3term(r) - he_orbitals(r);
}



class YetAnotherWrapperClass {
    const real_function_6d& f;
    const real_function_3d f3;
    const long k;
    real_function_6d eri;
    const Tensor<double>& qx6;
    const Tensor<double>& qx3;
    Tensor<double> identity;

public:
    YetAnotherWrapperClass(const real_function_6d& f)
        : f(f)
        , f3(real_factory_3d(f.get_impl()->world))
        , k(f.k())
        , eri()
        , qx6(FunctionCommonData<double,6>::get(k).quad_x)
        , qx3(FunctionCommonData<double,3>::get(k).quad_x)
    {
        identity=Tensor<double>(k,k,k);
        identity=1.0;
        eri=ERIFactory<double,6>(f.get_impl()->world).dcut(1.e-8);
    }

    void operator()(const Key<6>& key, Tensor<double>& t) const {

#if 0
        // break key into particles
        const Vector<Translation, 6> l=key.translation();
        const Vector<Translation, 3> l1 {l[0],l[1],l[2]};
        const Vector<Translation, 3> l2 {l[3],l[4],l[5]};
        const Key<3> key1(key.level(),l1);
        const Key<3> key2(key.level(),l2);

//        Tensor<double> g12=eri.coeff(key);
        Tensor<double> g12=eri.get_impl()->coeffs2values(
		key,eri.get_impl()->get_functor()->coeff(key).full_tensor());
        Tensor<double> pot1(k,k,k), pot2(k,k,k);
        f3.get_impl()->fcube(key1,Z2,qx3,pot1);
        f3.get_impl()->fcube(key2,Z2,qx3,pot2);

        // direct product: V(1) * E(2) + E(1) * V(2)
//        g12+= outer(pot1,identity) + outer(identity,pot2);
        t.emul(g12);
#else
	Tensor<double> v(k,k,k,k,k,k);
       	f.get_impl()->fcube(key, helium_pot, qx6, v);
//       	f.get_impl()->fcube(key, coul, qx6, v);
//        t.emul(v);
        t=v;

#endif

    }

};

real_function_6d multiply_by_V(const real_function_6d& psi) {
    real_function_6d Vpsi = copy(psi);
    Vpsi.unaryop(YetAnotherWrapperClass(Vpsi));
    return Vpsi;
}

template<size_t NDIM>
void save_function(World& world, const Function<double,NDIM>& pair, const std::string name) {
    if (world.rank()==0) print("saving function",name);
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    ar & pair & FunctionDefaults<NDIM>::get_cell();
}

template<size_t NDIM>
void load_function(World& world, Function<double,NDIM>& pair, const std::string name) {
    archive::ParallelInputArchive ar(world, name.c_str());
//    archive::ParallelInputArchive ar(world, "restart");
    Tensor<double> cell;
    ar & pair & cell;
}


struct LBCost {
    double leaf_value;
    double parent_value;
    LBCost(double leaf_value=1.0, double parent_value=1.0)
        : leaf_value(leaf_value)
        , parent_value(parent_value)
    {}

    double operator()(const Key<6>& key, const FunctionNode<double,6>& node) const {
//        if (key.level() <= 1) {
//		return 100.0;
//        }
//        else {
//		return 1.0;
//	}
        if (node.is_leaf()) {
            return std::abs(node.coeff().rank());
        } else {
            return parent_value;
        }
    }
};


//struct true_op {
//    bool operator()(FunctionImpl<double,3>* impl, const Key<3>& key, const FunctionNode<double,3>& t) const {
//        return true;
//    }
//
//    template <typename Archive> void serialize (Archive& ar) {}
//};

struct true_op {
    bool operator()(FunctionImpl<double,6>* impl, const Key<6>& key, const FunctionNode<double,6>& t) const {
        return true;
    }

    template <typename Archive> void serialize (Archive& ar) {}
};


struct true_if_n_gt_op {
    long _l;
    true_if_n_gt_op() : _l(100) {}
    true_if_n_gt_op(const int level) : _l(level) {}
    bool operator()(FunctionImpl<double,6>* impl, const Key<6>& key, const FunctionNode<double,6>& t) const {
        return (key.level()>_l);
    }

    template <typename Archive> void serialize (Archive& ar) {}
};

/// Returns the box at level n that contains the given point in simulation coordinates
     Key<6> simpt2key(const Vector<double,6>& pt, Level n) {
         const std::size_t NDIM=6;
         Vector<Translation,NDIM> l;
         double twon = std::pow(2.0, double(n));
         for (std::size_t i=0; i<NDIM; ++i) {
             l[i] = Translation(twon*pt[i]);
         }
         return Key<NDIM>(n,l);
     }

//struct true_if_close_to_nucleus {
//    typedef Vector<double,6> coordT;
//    coordT nuc;
//    true_if_close_to_nucleus() {}
//    true_if_close_to_nucleus(coordT& a) : nuc(a) {}
//    bool operator()(FunctionImpl<double,6>* impl, const Key<6>& key, const FunctionNode<double,6>& t) const {
//        coordT simpt;
//        user_to_sim(nuc, simpt);
//        Key<6> specialkey = simpt2key(simpt, key.level());
//        return (specialkey.is_neighbor_of(key));
//    }
//
//    template <typename Archive> void serialize (Archive& ar) {
//        ar & nuc;
//    }
//
//};

/// iterate the helium wavefunction
/// @param[in]      world   world
/// @param[in]      Vpsi    -2.0*V*pair
/// @param[in/out]  psi     pair function
/// @param[in/out]  eps     energy
void iterate(World& world, real_function_6d& Vpsi, real_function_6d& psi, double& eps) {

    LoadBalanceDeux<6> lb(world);
    double ncoeff=std::pow(double(FunctionDefaults<6>::get_k()),double(6));
    lb.add_tree(Vpsi,LBCost(1.0,ncoeff));
    FunctionDefaults<6>::redistribute(world, lb.load_balance(2.0,false));
    if(world.rank() == 0) printf("redistributed at time   %.1fs\n", wall_time());

    MADNESS_ASSERT(eps<0.0);
    real_convolution_6d op = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
    op.modified()=false;

    if(world.rank() == 0) printf("starting convolution at time %.1fs\n", wall_time());
    real_function_6d tmp = op(Vpsi);
    if(world.rank() == 0) printf("ending convolution at time   %.1fs\n", wall_time());

    tmp.print_size("tmp before truncation");
    tmp.truncate();
    tmp.print_size("tmp after truncation");
    if(world.rank() == 0) printf("truncated at time   %.1fs\n", wall_time());

    double norm = tmp.norm2();
    real_function_6d r = tmp-psi;
    double rnorm=r.norm2();
    double eps_new = eps -0.5 * inner(r,Vpsi)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("eps_new",eps_new);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

void iterate(World& world, const real_function_3d& V, real_function_3d& psi, double& eps) {

	real_function_3d Vpsi = (V*psi);
    Vpsi.scale(-2.0).truncate();

    real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.0001, 1e-6);
    real_function_3d tmp;
//
//    // direct Vpsi
//    if (0) {
//    	// set the impl of arg to get the structure of the target tree
//    	real_function_3d arg=CompositeFactory<double,3,3>(world)
//    	        .ket(copy(Vpsi).get_impl())
//    			.muster(copy(Vpsi).get_impl());
//
//    	tmp = op(arg).truncate();
//    }
//
//    // composite Vpsi
//    else {
////        psi.refine_general(true_op());
//        print("psi.size()",psi.size(),psi.tree_size());
//    	real_function_3d arg=CompositeFactory<double,3,3>(world)
////    	        .V_for_particle1(copy(V).get_impl())
//    	        .ket(copy(Vpsi).get_impl())
//				.muster(copy(Vpsi).get_impl());
//
//    	tmp = op(arg);
//    	print("tmp.size()",tmp.size(),tmp.norm2());
////    	tmp.scale(-2.0);
//    	tmp.truncate();
//        print("tmp.truncated.size()",tmp.size());
//
//    }

    // conventional
//    tmp=op(Vpsi).truncate();

    double norm = tmp.norm2();
    real_function_3d r = tmp-psi;
    double rnorm = r.norm2();
    double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
    if (world.rank() == 0) {
        print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        print("eps_new",eps_new);
    }
    psi = tmp.scale(1.0/norm);
    eps = eps_new;
}

void compute_energy(World& world, const real_function_3d& psi, const real_function_3d& pot, double& ke, double& pe) {


    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
    	real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    	real_function_3d dpsi = D(psi);
    	kinetic_energy += 0.5*inner(dpsi,dpsi);
    }
    ke=kinetic_energy;

    pe=inner(psi,pot*psi);
    if(world.rank() == 0) {
        printf("compute the energy at time   %.1fs\n", wall_time());
        printf("kinetic energy      %12.8f\n", ke);
        printf("potential energy    %12.8f\n", pe);
        printf("total energy        %12.8f\n", pe+ke);
    }
}



void compute_energy(World& world, const real_function_6d& pair,
		const real_function_3d& pot1, const real_function_3d& pot2, double& ke, double& pe) {

	// compute kinetic energy
	ke=0.0;
	if (1) {
		for (int axis=0; axis<6; axis++) {
			real_derivative_6d D = free_space_derivative<double,6>(world, axis);
			real_function_6d dpsi = D(pair);
//            double a=0.5*inner(dpsi,dpsi);
			double aa=dpsi.norm2();
			double a=0.5*aa*aa;
			ke += a;
			if (world.rank()==0) print("done with axis",axis, a);
		}
	}
	if (world.rank()==0) {
	    print("kinetic energy:", ke);
	    printf("\nkinetic at time %.1fs\n\n", wall_time());
	}

	// compute potential energy
	pe=0.0;
	if (1) {

		// two-electron interaction potential
		real_function_6d eri=ERIFactory<double,6>(world).dcut(1.e-8);


		real_function_6d v11=CompositeFactory<double,6,3>(world)
				.ket(copy(pair))
				.g12(eri)
				.V_for_particle1(copy(pot1))
				.V_for_particle2(copy(pot2))
				;

		// make the tree (optional!)
        if (1) {
		    // this is a dummy convolution
            real_convolution_6d op = BSHOperator<6>(world, -0.1, 0.00001, 1e-6);
            op.modified()=true;
            v11.fill_tree(op);
        }

		double a=inner(pair,v11);
        if (world.rank()==0) print("<phi|V_tot|phi> ", a);
		pe=a;
	} else {
		pe=-2.0*ke;
	}
	if (world.rank()==0) print("total energy",ke+pe);

	if(world.rank() == 0) printf("\npotential at time %.1fs\n\n", wall_time());

}


void solve(World& world, real_function_6d& pair, double& energy, long maxiter, double dcut) {

    const long k=FunctionDefaults<6>::get_k();
    const double thresh=FunctionDefaults<6>::get_thresh();

	if (world.rank()==0) {
		print("solving the helium atom with parameters");
		print("energy   ",energy);
		print("dcut     ",dcut);
		print("k        ",k);
		print("thresh   ",thresh);
	}

	// one-electron potential
	real_function_3d pot1=real_factory_3d(world).f(Z2);
	real_function_3d pot2=real_factory_3d(world).f(Z2);

	if(world.rank() == 0) printf("\nproject at time %.1fs\n\n", wall_time());

	for (long i=0; i<maxiter; i++) {

		// two-electron interaction potential
		real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

		real_function_6d vphi=CompositeFactory<double,6,3>(world)
							.ket(copy(pair))
							.g12(eri)
							.V_for_particle1(copy(pot1))
							.V_for_particle2(copy(pot2))
							;
//		vphi.get_impl()->make_Vphi();
		MADNESS_EXCEPTION("fix solve in helium_mp2",1);

		long tree_size1=pair.tree_size();
        long tree_size2=vphi.tree_size();
        long size2=vphi.size();
        if (world.rank()==0) print("refined pair",tree_size1,tree_size2,size2);


        // plot xy plane containing the origin
        for (int ii=0; ii<20; ii++) {
            coord_6d fix_coord(0.0);
            // electron 2:
            fix_coord[3]=0.0;
            fix_coord[4]=0.1+0.1*ii;
            fix_coord[5]=0.0;

            Tensor<double> cell(6,2);
            cell(Slice(_),0)=-2.0;
            cell(Slice(_),1)= 2.0;
            std::string filename="plot_plane_d"+stringify(fix_coord[4])+"_k"+stringify(k)+"_eps"+stringify(thresh);
            vphi.get_impl()->print_plane(filename,"xy",fix_coord);

        }
        std::string name="vphi_k"+stringify(k)+"_e"+stringify(thresh)+"_it"+stringify(i);
		double L=FunctionDefaults<6>::get_cell_width()[0];
		coord_6d lo(0.0), hi(0.0);
		lo[0]=-L/2;
		hi[0]=L/2;
		for (int ii=-5; ii<6; ii++) {
		    lo[3]=hi[3]=double(ii);
		    trajectory<6> line=trajectory<6>::line2(lo,hi,601);
		    plot_along<6>(world,line,vphi,(name+"lineplot"+stringify(ii)));
            plot_along<6>(world,line,helium_pot,(name+"lineplot_coul"+stringify(ii)));
		}


		iterate(world,vphi,pair,energy);
		pair.reconstruct();

		long tree_size=pair.tree_size();
		long size=pair.size();
		if(world.rank() == 0) print("pair.tree_size() in iteration",i,":",tree_size);
		if(world.rank() == 0) print("pair.size() in iteration",i,     ":",size);
		name="restart_k"+stringify(k)+"_e"+stringify(thresh)+"_it"+stringify(i);
		save_function(world,pair,name);


		if (i%3==0) {
			double ke, pe;
			compute_energy(world,pair,pot1,pot2,ke,pe);
			if (world.rank()==0) print("virial ratio in iteration  :",i, pe/ke, ke+pe);
		}
	}


    // plot xy plane containing the origin
    for (int i=0; i<20; i++) {
        coord_6d fix_coord(0.0);
        // electron 2:
        fix_coord[3]=0.0;
        fix_coord[4]=0.1+0.1*i;
        fix_coord[5]=0.0;

        pair.reconstruct();
        std::shared_ptr<FunctionImpl<double,6> > pair_impl=pair.get_impl();
        Tensor<double> cell(6,2);
        cell(Slice(_),0)=-2.0;
        cell(Slice(_),1)= 2.0;
        std::string filename="plot_plane_d"+stringify(fix_coord[4])+"_k"+stringify(k)+"_eps"+stringify(thresh);
        pair_impl->print_plane(filename,"xy",fix_coord);

    }

	if (world.rank()==0) print("finished",maxiter,"iterations");
	double ke, pe;
	compute_energy(world,pair,pot1,pot2,ke,pe);
	if (world.rank()==0) print("virial ratio   :", pe/ke, ke+pe);


}


// test the modified NS form
void test_modified(World& world) {

    real_function_3d V=real_factory_3d(world).f(Z1);
    real_function_3d psi=real_factory_3d(world).f(gauss_3d);
    psi.scale(1.0/inner(psi,psi));

    double energy=-0.5;
    double eps;

    coord_3d lo, hi;
    const double L=FunctionDefaults<3>::get_cell_width()[0];
    lo[0]=-L/2; hi[0]=L/2;
    trajectory<3> traj=trajectory<3>::line2(lo,hi,201);
    plot_along(world,traj,psi,"name");


    if(world.rank() == 0) printf("starting at time %.1fs\n\n", wall_time());
    for (int i=0; i<20; i++) {

        print("");
        print("entering iteration",i);
        eps=energy;
        real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.0001, 1e-6);

        op.modified()=true;
        op.doleaves=false;
        print("operator modified",op.modified());

        real_function_3d Vpsi = (V*psi);
        Vpsi.scale(-2.0).truncate();
//        real_function_3d tmp=op(Vpsi).truncate();
        real_function_3d tmp;
        MADNESS_EXCEPTION("op() only with NDIM==6",1);
        psi=tmp;
        psi.scale(1.0/psi.norm2());
        plot_along(world,traj,psi,"name"+stringify(i));

        double ke,pe;
        compute_energy(world,psi,V,ke,pe);

    }
    if(world.rank() == 0) printf("finished at time %.1fs\n\n", wall_time());


}



// test the modified NS form
void test_recursive_application(World& world) {

    if (world.rank()==0) print("in test_recursive_application");

    // one orbital at a time
//    real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
    real_function_3d orbital=real_factory_3d(world).f(h_orbital);

    double norm=inner(orbital,orbital);
    orbital.scale(1.0/sqrt(norm));

    real_function_6d pair=hartree_product(orbital,orbital);
    double norm2=inner(pair,pair);
    if (world.rank()==0) print("<phi | phi>",norm2);
    double ke,pe,eps;
    const real_function_3d pot1=real_factory_3d(world).f(Z1);
    const real_function_3d pot2=real_factory_3d(world).f(Z1);


    {
        // 3d
        compute_energy(world,orbital,pot1,ke,pe);
        eps=-0.5;

        if (1) {
            real_function_3d vphi=-2.0*(orbital*pot1);
            real_convolution_3d op2 = BSHOperator<3>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op2.modified()=false;
            real_function_3d tmp;//=op2(vphi);
            MADNESS_EXCEPTION("op() only with NDIM==6",1);
            tmp.scale(1.0/sqrt(tmp.norm2()));
            double ke, pe;
            compute_energy(world,tmp,pot1,ke,pe);

            coord_3d lo, hi;
            const double L=FunctionDefaults<3>::get_cell_width()[0];
            lo[0]=-L/2; hi[0]=L/2;
            trajectory<3> traj=trajectory<3>::line2(lo,hi,201);
            plot_along(world,traj,tmp,"name0");

        }

        {
            real_function_3d two_phi=orbital.scale(-2.0);

            real_function_3d phi2=CompositeFactory<double,3,3>(world)
                                .ket(copy(two_phi))
                                .V_for_particle1(copy(pot1))
                                ;
            {
                real_convolution_3d op = BSHOperator<3>(world, sqrt(-2*eps), 0.00001, 1e-6);
                op.modified()=true;
                phi2.fill_tree(op);
            }

            real_convolution_3d op = BSHOperator<3>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op.modified()=false;
//            phi2=op(phi2);
            phi2.scale(1.0/sqrt(phi2.norm2()));

            compute_energy(world,phi2,pot1,ke,pe);


            coord_3d lo, hi;
            const double L=FunctionDefaults<3>::get_cell_width()[0];
            lo[0]=-L/2; hi[0]=L/2;
            trajectory<3> traj=trajectory<3>::line2(lo,hi,201);
            plot_along(world,traj,phi2,"name1");

        }



    }

//    double v=inner(pair,vphi);
//    if (world.rank()==0) print("inner(pair,vphi",v);


}


void test_adaptive_tree(World& world, const bool restart, const std::string restartname) {

    if (world.rank()==0) print("in test_adaptive_tree");

    double thresh=FunctionDefaults<6>::get_thresh();
//    const double L=FunctionDefaults<6>::get_cell_width()[0];
    long k=FunctionDefaults<6>::get_k();

    // one orbital at a time
    real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
//    real_function_3d orbital=real_factory_3d(world).f(he_plus_orbital);
    orbital.scale(1.0/orbital.norm2());

    const real_function_3d pot1=real_factory_3d(world).f(Z2);
    const real_function_3d pot2=real_factory_3d(world).f(Z2);

    // get he+ energy
//    for (int i=0; i<20; ++i) {
//        double ke,pe;
//        compute_energy(world,orbital,pot1,ke,pe);
//        double eps=ke+pe;
//        std::string name="orbital_it"+stringify(i);
//        save_function(world,orbital,name);
//        iterate(world,pot1,orbital,eps);
//    }
//    orbital=real_factory_3d(world).f(he_orbital_McQuarrie);

    double en1=inner(orbital,pot1*orbital);
    if (world.rank()==0) printf("1e pot energy computed as < phi | V phi>: %12.8f\n",en1);

    double ke,pe,eps;
    compute_energy(world,orbital,pot1,ke,pe);

    real_function_6d pair;
    if (restart) {
        load_function(world,pair,restartname);
    } else {
        pair=hartree_product(orbital,orbital);
    }
    double norm2=inner(pair,pair);
    if (world.rank()==0) print("<phi | phi>",norm2);

    pair=pair=symmetrize(pair,"sy_particle");

    compute_energy(world,pair,pot1,pot2,ke,pe);
    eps=ke+pe;

//    real_function_6d vpair=hartree_product(orbital*pot1,orbital);
//    double en=inner(pair,vpair);
//    if (world.rank()==0) printf("2e pot energy computed as < phi | V phi>: %12.8f\n",en);
//
//    long size,tsize;
//    size=vpair.size();
//    tsize=vpair.tree_size();
//    if (world.rank()==0) print("done computing pot size",size,tsize);

    for (int ii=0; ii<10; ++ii) {



        // two-electron interaction potential
        real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

        real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                            .ket(copy(pair))
                                            .g12(eri)
                                            .V_for_particle1(copy(pot1))
                                            .V_for_particle2(copy(pot2))
                                            ;

        // make the tree
        if (1) {
            real_convolution_6d op = BSHOperator<6>(world, sqrt(-2*eps), 0.00001, 1e-6);
            op.modified()=true;
            vphi.fill_tree(op);
            vphi.scale(-2.0);
            vphi.print_size("vphi");

        }

        // convolute
        {
            iterate(world,vphi,pair,eps);

        }

        {
            const std::string appendix="_k"+stringify(k)+"_eps"+stringify(thresh)+"_it"+stringify(ii);
            std::string filename="pair"+appendix;
            save_function(world,pair,filename);
            pair.print_size("pair");
        }

        compute_energy(world,pair,pot1,pot2,ke,pe);
//        eps=ke+pe;


        double norm=pair.norm2();
        if (world.rank()==0) printf("norm(pair) : %12.8f",norm);
        if(world.rank() == 0) printf("\nfinished iteration %d at time %.1fs\n\n", ii, wall_time());

    }
}


void test_truncation(World& world) {

    const real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
    const real_function_3d pot1=real_factory_3d(world).f(Z2);
    const real_function_3d prod=orbital*pot1;
    prod.print_size("prod");

    real_function_3d a1,a2;
    {
        a1=copy(prod);
        a1.reconstruct();
        a1.truncate();
        a1.reconstruct();
        a1.print_size("a1");
    }
    {
        a2=copy(prod);
        a2.compress();
        a2.truncate();
        a2.reconstruct();
        a2.print_size("a2");
    }
    real_function_3d a=a1-a2;
    print("norm(a)",a.norm2());


}


void test_compress(World& world) {

    const real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);
    const real_function_3d pot1=real_factory_3d(world).f(Z2);
    const real_function_3d prod=orbital*pot1;

    real_function_6d pair=hartree_product(orbital,orbital);
    real_function_6d vpair=hartree_product(prod,prod);
    pair.print_size("pair");
    double ke,pe;

    real_function_6d a1,a2,a3;
//    {
//        a1=copy(pair);
//        a1.compress();
//        a1.reconstruct();
//        a1.print_size("a1");
//    }
    {
        a2=copy(vpair);
        a2.print_size("a2");
    }
    {
        a3=copy(pair);
        real_function_6d vphi=CompositeFactory<double,6,3>(world)
                                            .ket(copy(a3))
//                                            .g12(eri.get_impl())
                                            .V_for_particle1(copy(pot1))
                                            .V_for_particle2(copy(pot1))
                                            ;

        // make the tree
        real_convolution_6d op = BSHOperator<6>(world, sqrt(1.0), 0.00001, 1e-6);
        op.modified()=true;
        vphi.fill_tree(op);
        a3=vphi;
        a3.print_size("a3");
    }
    real_function_6d a=a3-a2;
    print("norm(a)",a.norm2());


}


int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    std::cout.precision(6);

    // defaults
    double L = 16;   // box size
    long k = 4 ;        // wavelet order
    double thresh = 1.e-2; // precision
    TensorType tt=TT_2D;

    // set the parameters
    bool restart=false;
    std::string restart_name;

    for(int i = 1; i < argc; i++) {
        const std::string arg=argv[i];

        // break parameters into key and val
        size_t pos=arg.find("=");
        std::string key=arg.substr(0,pos);
        std::string val=arg.substr(pos+1);

        if (key=="size") L=atof(val.c_str());               // usage: size=10
        if (key=="k") k=atoi(val.c_str());                  // usage: k=5
        if (key=="thresh") thresh=atof(val.c_str());        // usage: thresh=1.e-3
        if (key=="TT") {
            if (val=="TT_2D") tt=TT_2D;
            else if (val=="TT_3D") tt=TT_3D;
            else if (val=="TT_FULL") tt=TT_FULL;
            else {
                print("arg",arg, "key",key,"val",val);
                MADNESS_EXCEPTION("confused tensor type",0);
            }

        }
        if (key=="restart") {                               // usage: restart=path/to/mo_file
            restart_name=stringify(val);
            restart=true;
        }
    }



    if (world.rank()==0) {
        print("restart mode",restart," restart_name=", restart_name);
        printf("\nstarting at time %.1fs\n\n", wall_time());
    }


    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
    FunctionDefaults<6>::set_tensor_type(tt);
    FunctionDefaults<6>::set_apply_randomize(true);



    if (world.rank()==0) {
        print("polynomial order:  ", FunctionDefaults<6>::get_k());
        print("threshold:         ", FunctionDefaults<6>::get_thresh());
        print("cell size:         ", L);
        print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
        print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());
        print("");
        print("orthogonalization  ", OrthoMethod());
        print("facReduce          ", GenTensor<double>::fac_reduce());
        print("max displacement   ", Displacements<6>::bmax_default());
        print("apply randomize    ", FunctionDefaults<6>::get_apply_randomize());
        print("world.size()       ", world.size());
        print("");
    }

    if (world.rank()==0) {
        print("size consistency of the 6d green's function?");
        print("");
    }

    if (0) test_compress(world);
    if (0) test_truncation(world);
    if (0) test_modified(world);
    if (0) test_recursive_application(world);
    if (1) test_adaptive_tree(world,restart,restart_name);

    if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    finalize();

    return 0;


    // helium
#if 0

    if (0) {
        // compute the energy of Hylleraas
        real_function_6d pair=real_factory_6d(world).f(hylleraas_3term);
        long tree_size=pair.tree_size();
        long size=pair.size();
        if (world.rank()==0) print("Hylleraas-3term;  tree size",tree_size);
        if (world.rank()==0) print("Hylleraas-3term;  size     ",size);

        // normalize pair function
        double norm=inner(pair,pair);
        pair.scale(1.0/sqrt(norm));

        tree_size=pair.tree_size();
        size=pair.size();
        if (world.rank()==0) print("Hylleraas-3term;  tree size",tree_size);
        if (world.rank()==0) print("Hylleraas-3term;  size     ",size);

        // initial energy
        double ke,pe;
        real_function_3d pot1=real_factory_3d(world).f(Z2);
        real_function_3d pot2=real_factory_3d(world).f(Z2);

        compute_energy(world,pair,pot1,pot2,ke,pe);

        return 0;
    }


    // one orbital at a time
	real_function_3d orbital=real_factory_3d(world).f(he_orbital_McQuarrie);

	{
		long tree_size=orbital.tree_size();
		long size=orbital.size();
		if (world.rank()==0) {
		    print("orbital.tree_size()",tree_size);
		    print("orbital.size()     ",size);
		}
	}


    {
    	double norm=inner(orbital,orbital);
    	orbital.scale(1.0/sqrt(norm));
//	orbital.print_tree();

    	// compute kinetic energy
    	double kinetic_energy = 0.0;
    	for (int axis=0; axis<3; axis++) {
    		real_derivative_3d D = free_space_derivative<double,3>(world, axis);
    		real_function_3d dpsi = D(orbital);
    		kinetic_energy += 0.5*inner(dpsi,dpsi);
    	}
    	if (world.rank()==0) print("kinetic energy/electron  :", kinetic_energy);

    	// compute potential energy
    	real_function_3d one_el_pot=real_factory_3d(world).f(Z2);
    	double pe=inner(orbital,one_el_pot*orbital);
    	if (world.rank()==0) print("potential energy/electron:",pe);
    }


    real_function_6d pair;
    double energy;
    // where we start and where we end
    double current_thresh=1.e-2;
    const double max_thresh=FunctionDefaults<6>::get_thresh();

    if (not restart) {
        pair=hartree_product(orbital,orbital);

        LoadBalanceDeux<6> lb(world);
        double ncoeff=std::pow(FunctionDefaults<6>::get_k(),6);
        lb.add_tree(pair,LBCost(1.0,ncoeff));
        FunctionDefaults<6>::redistribute(world, lb.load_balance(2.0,false));

        // normalize pair function
        double norm=inner(pair,pair);
        pair.scale(1.0/sqrt(norm));
    } else {
        load_function(world,pair,restart_name);
    	pair=project(pair,FunctionDefaults<6>::get_k(),pair.thresh());	// change to appropriate k
        if (world.rank()==0) {
            print("restart from file ",restart_name);
            print("k                 ",pair.k());
            print("thresh            ",pair.thresh());
//            print("read energy       ",energy);
        }
        current_thresh=pair.thresh();
    }
    if(world.rank() == 0) printf("\npair function at time %.1fs\n\n", wall_time());

    {
    	long tree_size=pair.tree_size();
    	long size=pair.size();
    	if (world.rank()==0) {
    	    print("pair.tree_size()",tree_size);
    	    print("pair.size()     ",size);
    	}
    }

    // initial energy
//    if (not restart) {
        double ke,pe;
        real_function_3d pot1=real_factory_3d(world).f(Z2);
        real_function_3d pot2=real_factory_3d(world).f(Z2);
//        // Coulomb potential
//        real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
////        real_function_3d orbital2=real_factory_3d(world).f(he_orbital_McQuarrie);
//        real_function_3d rho = 0.5*square(orbital).truncate();
//        real_function_3d coulpot = op(rho).truncate();
//        pot1+=coulpot;
//        pot2+=coulpot;


        compute_energy(world,pair,pot1,pot2,ke,pe);
        energy=ke+pe;
        if (world.rank()==0) print("computed energy   ",energy);
//    }

    // solve for thresh=1.e-3
    if (current_thresh>1.e-3) {
        FunctionDefaults<6>::set_thresh(1.e-3);
    	real_function_6d pair1=project(pair,FunctionDefaults<6>::get_k(),1.e-3);
    	solve(world,pair1,energy,4,1.e-4);
    	pair=pair1;
    	current_thresh=1.e-3;
    }

    // solve for thresh=1.e-4
    if ((current_thresh>1.e-4) and (max_thresh<9.e-4)) {
        FunctionDefaults<6>::set_thresh(1.e-4);
    	real_function_6d pair1=project(pair,FunctionDefaults<6>::get_k(),1.e-4);
    	solve(world,pair1,energy,8,1.e-4);
    	pair=pair1;
    	current_thresh=1.e-4;
    }

    // solve for thresh=1.e-5
    if ((current_thresh>1.e-5) and (max_thresh<9.e-5)) {
        FunctionDefaults<6>::set_thresh(1.e-5);
    	real_function_6d pair1=project(pair,FunctionDefaults<6>::get_k(),1.e-5);
    	solve(world,pair1,energy,8,1.e-8);
    	pair=pair1;
        current_thresh=1.e-5;
    }

    // solve for thresh=1.e-6
    if (max_thresh<9.e-6) {
        FunctionDefaults<6>::set_thresh(1.e-6);
    	real_function_6d pair1=project(pair,FunctionDefaults<6>::get_k(),1.e-6);
    	solve(world,pair1,energy,8,1.e-8);
    	pair=pair1;
    }

    // finally solve
    FunctionDefaults<6>::set_thresh(max_thresh);
    solve(world,pair,energy,10,1.e-8);

#endif


     // hylleraas for ACS
#if 0

     FunctionDefaults<3>::set_k(k);
     FunctionDefaults<3>::set_thresh(thresh);
     FunctionDefaults<3>::set_cubic_cell(-L/2,L/2);

     FunctionDefaults<6>::set_k(k);
     FunctionDefaults<6>::set_thresh(thresh);
     FunctionDefaults<6>::set_cubic_cell(-L/2,L/2);
     FunctionDefaults<6>::set_tensor_type(TT_2D);


     print("polynomial order:  ", FunctionDefaults<6>::get_k());
     print("threshold:         ", FunctionDefaults<6>::get_thresh());
     print("cell size:         ", L);
     print("truncation mode:   ", FunctionDefaults<6>::get_truncate_mode());
     print("tensor type:       ", FunctionDefaults<6>::get_tensor_type());

	 // plot circle with radius 0.5, with electron 2 at (0.0, 0.5, 0.0)
     if (0) {
	     coord_3d el2(0.0);
	     const double phi=1.0;
	     const double radius=0.5;
	     el2[0]=radius * sin(phi);
	     el2[1]=radius * cos(phi);
 		 std::string filename="hylleraas_plotfile_k"+stringify(k)+"_eps"+stringify(thresh);
		 trajectory<6> traj(0.5,el2,201);
		 plot_along<6>(world,traj,hylleraas_3term,filename);
		 print("plotting done");
     }

     // pair function
//     real_function_6d pair=real_factory_6d(world).f(he_correlation);
     real_function_6d pair=real_factory_6d(world).f(hylleraas_3term);
//     print("Hylleraas correlation part: 3-term minus Ed's HF");
     print("Hylleraas-3term;  tree size",pair.tree_size());
     print("Hylleraas-3term;  size     ",pair.size());

     // normalize pair function
     double norm=inner(pair,pair);
     print("norm(ij_pair)",norm);
     pair.scale(1.0/sqrt(norm));

     print("Hylleraas-3term;  tree size",pair.tree_size());
     print("Hylleraas-3term;  size     ",pair.size());

	 // plot circle with radius 0.5, with electron 2 at (0.0, 0.5, 0.0)
     {
		 coord_3d el2(0.0);
		 const double phi=1.0;
		 const double radius=0.5;
		 el2[0]=radius * sin(phi);
		 el2[1]=radius * cos(phi);
		 std::string filename="plotfile_k"+stringify(k)+"_eps"+stringify(thresh);
		 trajectory<6> traj(0.5,el2,201);
		 plot_along<6>(world,traj,pair,filename);
     }

     // plot xy plane containing the origin
     for (int i=0; i<20; i++) {
    	 coord_6d fix_coord(0.0);
    	 // electron 2:
    	 fix_coord[3]=0.0;
    	 fix_coord[4]=0.1+0.1*i;
    	 fix_coord[5]=0.0;

    	 std::shared_ptr<FunctionImpl<double,6> > pair_impl=pair.get_impl();
    	 Tensor<double> cell(6,2);
    	 cell(Slice(_),0)=-2.0;
    	 cell(Slice(_),1)= 2.0;
    	 std::string filename="plot_plane_d"+stringify(fix_coord[4])+"_k"+stringify(k)+"_eps"+stringify(thresh);
    	 pair_impl->print_plane(filename,"xy",fix_coord,cell);

     }
     if (1) {
		 // plot function thru xy plane containing the origin
		 std::shared_ptr<FunctionImpl<double,6> > pair_impl=pair.get_impl();
		 Tensor<double> cell(6,2);
		 cell(Slice(_),0)=-2.0;
		 cell(Slice(_),1)= 2.0;
		 std::string filename="plot_function_k"+stringify(k)+"_eps"+stringify(thresh);
		 FILE * f = fopen(filename.c_str(), "w");
		 if(!f) MADNESS_EXCEPTION("plotvtk: failed to open the plot file", 0);

         fprintf(f,"\\psset{xunit=4cm}\n");
         fprintf(f,"\\psset{yunit=4cm}\n");
         fprintf(f,"\\begin{pspicture}(0,0)(0,4)\n");
         fprintf(f,"\\pslinewidth=0.05pt\n");


		 for (int i=-100; i<100; i++) {
			 for (int j=-100; j<100; j++) {
				 coord_6d fix_coord(0.0);
				 // electron 2:
				 fix_coord[0]=0.0+i*0.05;
				 fix_coord[1]=0.0+j*0.05;
				 fix_coord[2]=0.0;
				 fix_coord[3]=0.0;
				 fix_coord[4]=0.5;
				 fix_coord[5]=0.0;
				 print("coord",fix_coord);
				 fprintf(f,"%12.8f   %12.8f   %12.8f\n",fix_coord[0],fix_coord[1],pair(fix_coord));
			 }
		 }

         fprintf(f,"\\end{pspicture}\n");
		 fclose(f);

     }

#if 0

     // compute kinetic energy
     double kinetic_energy = 0.0;
     for (int axis=0; axis<6; axis++) {
     	real_derivative_6d D = free_space_derivative<double,6>(world, axis);
     	real_function_6d dpsi = D(pair);
     	kinetic_energy += 0.5*inner(dpsi,dpsi);
     	print("done with axis",axis);
     }
     print("kinetic energy:", kinetic_energy);

     if(world.rank() == 0) printf("\nkinetic at time %.1fs\n\n", wall_time());


     // compute potential energy
     double potential_energy=0.0;
     {
 		// doomed copy of pair, to save pair
 		real_function_6d copy_of_pair=copy(pair);

 		// two-electron interaction potential
 		real_function_6d eri=ERIFactory<double,6>(world).dcut(dcut);

 	    // one-electron potential
 	    real_function_3d pot1=real_factory_3d(world).f(Z2);
 	    real_function_3d pot2=real_factory_3d(world).f(Z2);


 		real_function_6d v11=CompositeFactory<double,6,3>(world)
 				.ket(copy_of_pair.get_impl())
 				.g12(eri.get_impl())
 				.V_for_particle1(pot1.get_impl())
 				.V_for_particle2(pot2.get_impl())
 				;


 		double a=inner(pair,v11);
 		print("<phi|V_tot|phi> ", a);
 		potential_energy=a;
     }


     print("for Hylleraas 3-term");
     print("total energy   :",kinetic_energy+potential_energy);
     print("expected energy:",-2.902432);
#endif



#endif


     if(world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
     world.gop.fence();

     finalize();
     return 0;
}


