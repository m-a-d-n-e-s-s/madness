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

/// \file SCFOperators.h
/// \brief Operators for the molecular HF and DFT code
/// \defgroup chem The molecular density functional and Hartree-Fock code


#ifndef MADNESS_CHEM_SCFOPERATORS_H_
#define MADNESS_CHEM_SCFOPERATORS_H_

#include <madness.h>
#include <madness/mra/macrotaskq.h> // otherwise issues with install

namespace madness {

// forward declaration
class SCF;
class Nemo;
class NemoBase;
class OEP;
class NuclearCorrelationFactor;
class XCfunctional;
struct nemo_u1_functors;
class MacroTaskQ;
class Molecule;

typedef std::vector<real_function_3d> vecfuncT;

template<typename T, std::size_t NDIM>
class SCFOperatorBase {

public:
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;
    mutable nlohmann::json statistics;

    SCFOperatorBase() = default;
    SCFOperatorBase(std::shared_ptr<MacroTaskQ> taskq) : taskq(taskq) {}

    virtual ~SCFOperatorBase() {}

    std::shared_ptr<MacroTaskQ> taskq=0;

    /// print some information about this operator
    virtual std::string info() const = 0;

    /// apply this operator on the argument function
    ///
    /// \param  ket the argument function
    /// \return op(ket)
    virtual functionT operator()(const functionT& ket) const = 0;

    /// apply this operator on the argument vector of functions

    /// \param vket     argument vector
    /// \return         op(vket)
    virtual vecfuncT operator()(const vecfuncT& vket) const = 0;

    /// compute the matrix element <bra | op | ket>

    /// \param bra  bra state
    /// \param ket  ket state
    /// \return     the matrix element <bra | op | ket>
    virtual T operator()(const functionT& bra, const functionT& ket) const = 0;

    /// compute the matrix <vbra | op | vket>

    /// \param vbra  vector of bra states
    /// \param vket  vector of ket states
    /// \return     the matrix <vbra | op | vket>
    virtual tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const = 0;

};

template<typename T, std::size_t NDIM>
class Exchange : public SCFOperatorBase<T,NDIM> {
public:

    class ExchangeImpl;
    using implT = std::shared_ptr<ExchangeImpl>;
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;
private:
    implT impl;

public:
    enum ExchangeAlgorithm {
        small_memory, large_memory, multiworld_efficient, multiworld_efficient_row, fetch_compute
    };
    // print out algorithm
    friend std::ostream& operator<<(std::ostream& os, const ExchangeAlgorithm& alg) {
        switch (alg) {
            case small_memory:
                os << "smallmem";
                break;
            case large_memory:
                os << "largemem";
                break;
            case multiworld_efficient:
                os << "multiworld";
                break;
            case multiworld_efficient_row:
                os << "multiworld_row";
                break;
            case fetch_compute:
                os << "fetch_compute";
                break;
            default:
                os << "unknown algorithm";
        }
        return os;
    }

    static std::string to_string(const ExchangeAlgorithm alg) {
        std::stringstream ss;
        ss << alg;
        return ss.str();
    }

    static ExchangeAlgorithm string2algorithm(const std::string& alg_string) {
        ExchangeAlgorithm alg;
        std::string alg_lc=commandlineparser::tolower(alg_string);
        if (alg_lc=="smallmem") alg=small_memory;
        else if (alg_lc=="largemem") alg=large_memory;
        else if (alg_lc=="multiworld") alg=multiworld_efficient;
        else if (alg_lc=="multiworld_row") alg=multiworld_efficient_row;
        else if (alg_lc=="fetch_compute") alg=fetch_compute;
        else {
            std::string msg="unknown Exchange algorithm: "+alg_string;
            MADNESS_EXCEPTION(msg.c_str(),1);
        }
        return alg;
    }

    Exchange(World& world, const double lo, const double thresh=FunctionDefaults<NDIM>::get_thresh());

    /// ctor with a conventional calculation
    Exchange(World& world, const SCF *calc, const int ispin);

    /// ctor with a nemo calculation
    Exchange(World& world, const Nemo *nemo, const int ispin);

    std::string info() const {return "K";}

    bool is_symmetric() const;

    Exchange& set_symmetric(const bool flag);

    Exchange& set_algorithm(const ExchangeAlgorithm& alg);

    /// how the cloud will handle the data
    Exchange& set_macro_task_info(const MacroTaskInfo& info);
    Exchange& set_macro_task_info(const std::vector<std::string>& info) {
        impl->set_macro_task_info(info);
        return *this;
    }


    Exchange& set_printlevel(const long& level);

    /// batches per rank in the owner-pinned symmetric partition (>= 1)
    Exchange& set_batch_granularity(const long level);

    /// 1 = gather tile results per subworld, 2 = also gather per node first
    Exchange& set_accumulation_mode(const int mode);

    /// place tasks by their measured cost rather than by counting them
    Exchange& set_cost_aware_assignment(const bool flag);

    Exchange& set_taskq(std::shared_ptr<MacroTaskQ> taskq1) {
        this->taskq=taskq1;
        return *this;
    }

    Exchange& set_bra_and_ket(const vecfuncT& bra, const vecfuncT& ket);

    Function<T, NDIM> operator()(const Function<T, NDIM>& ket) const {
        vecfuncT vket(1, ket);
        vecfuncT vKket = this->operator()(vket);
        return vKket[0];
    }

    /// apply the exchange operator on a vector of functions

    /// note that only one spin is used (either alpha or beta orbitals)
    /// @param[in]  vket       the orbitals |i> that the operator is applied on
    /// @return     a vector of orbitals  K| i>
    vecfuncT operator()(const vecfuncT& vket) const;

    /// compute the matrix element <bra | K | ket>

    /// @param[in]  bra    real_function_3d, the bra state
    /// @param[in]  ket    real_function_3d, the ket state
    T operator()(const Function<T, NDIM>& bra, const Function<T, NDIM>& ket) const {
        return inner(bra, this->operator()(ket));
    }

    /// compute the matrix < vbra | K | vket >

    /// @param[in]  vbra    vector of real_function_3d, the set of bra states
    /// @param[in]  vket    vector of real_function_3d, the set of ket states
    /// @return K_ij
    Tensor<T> operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        vecfuncT vKket = this->operator()(vket);
        World& world=vket[0].world();
        auto result = matrix_inner(world, vbra, vKket);
        return result;
    }


};


template<typename T, std::size_t NDIM>
class Kinetic : public SCFOperatorBase<T,NDIM> {
    typedef DistributedMatrix<T> distmatT;
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;

public:
    Kinetic(World& world) : world(world) {
        gradop = gradient_operator<T,NDIM>(world);
    }

    std::string info() const {return "T";}

    functionT operator()(const functionT& ket) const {
        MADNESS_EXCEPTION("do not apply the kinetic energy operator on a function!",1);
        return ket;
    }

    vecfuncT operator()(const vecfuncT& vket) const {
        MADNESS_EXCEPTION("do not apply the kinetic energy operator on a function!",1);
        return vket;
    }

    T operator()(const functionT& bra, const functionT& ket) const {
        vecfuncT vbra(1,bra), vket(1,ket);
        Tensor<T> tmat=this->operator()(vbra,vket);
        return tmat(0l,0l);
    }

    tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        distmatT dkinetic;
        if (&vbra==&vket) {
            dkinetic = kinetic_energy_matrix(world,vbra);
        } else {
            dkinetic = kinetic_energy_matrix(world,vbra,vket);
        }
        tensorT kinetic(vbra.size(),vket.size());
        dkinetic.copy_to_replicated(kinetic);
        return kinetic;
    }

private:
    World& world;
    std::vector< std::shared_ptr<Derivative<T,NDIM> > > gradop;

    distmatT kinetic_energy_matrix(World & world, const vecfuncT & v) const;
    distmatT kinetic_energy_matrix(World & world, const vecfuncT & vbra,
            const vecfuncT & vket) const;

};


template<typename T, std::size_t NDIM>
class DerivativeOperator : public SCFOperatorBase<T,NDIM> {
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;

public:

    DerivativeOperator(World& world, const int axis1) : world(world), axis(axis1) {
        gradop = free_space_derivative<T,NDIM>(world, axis);
    }

    std::string info() const {return "D";}

    functionT operator()(const functionT& ket) const {
        vecfuncT vket(1,ket);
        return this->operator()(vket)[0];
    }

    vecfuncT operator()(const vecfuncT& vket) const {
        vecfuncT dvket=apply(world, gradop, vket, false);
        world.gop.fence();
        return dvket;
    }

    T operator()(const functionT& bra, const functionT& ket) const {
        vecfuncT vbra(1,bra), vket(1,ket);
        Tensor<T> tmat=this->operator()(vbra,vket);
        return tmat(0l,0l);
    }

    tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        vecfuncT dvket=this->operator()(vket);
        return matrix_inner(world,vbra,dvket, bra_equiv_ket);
    }

private:
    World& world;
    int axis;
    Derivative<T,NDIM> gradop;

};


/// the Laplacian operator: \sum_i \nabla^2_i

/// note that the application of the Laplacian operator is in general
/// unstable and very sensitive to noise and cusps in the argument.
///
/// !!! BE SURE YOU KNOW WHAT YOU ARE DOING !!!
///
/// For computing matrix elements, which is reasonably stable, we refer
template<typename T, std::size_t NDIM>
class Laplacian : public SCFOperatorBase<T,NDIM> {
    typedef Function<T,NDIM> functionT;
    typedef std::vector<functionT> vecfuncT;
    typedef Tensor<T> tensorT;

public:

    Laplacian(World& world, const double e=0.0) : world(world), eps(e) {
        gradop = gradient_operator<T,NDIM>(world);
    }

    std::string info() const {return "D^2";}

    functionT operator()(const functionT& ket) const {
        vecfuncT vket(1,ket);
        return this->operator()(vket)[0];
    }

    vecfuncT operator()(const vecfuncT& vket) const;

    T operator()(const functionT& bra, const functionT& ket) const {
        vecfuncT vbra(1,bra), vket(1,ket);
        Tensor<T> tmat=this->operator()(vbra,vket);
        return tmat(0l,0l);
    }

    tensorT operator()(const vecfuncT& vbra, const vecfuncT& vket) const {
        Kinetic<T,NDIM> t(world);
        return -2.0*t(vbra,vket);
    }

private:
    World& world;
    std::vector< std::shared_ptr< Derivative<T,NDIM> > > gradop;
    double eps;
};



template<typename T, std::size_t NDIM>
class Coulomb : public SCFOperatorBase<T,NDIM> {
public:

    class MacroTaskCoulomb : public MacroTaskOperationBase {
    public:
        // you need to define the exact argument(s) of operator() as tuple
        typedef std::tuple<const Function<double,NDIM>&, const std::vector<Function<T,NDIM>> &> argtupleT;

        using resultT = std::vector<Function<T,NDIM>>;

        class MacroTaskPartitionerCoulomb : public MacroTaskPartitioner {
        public:
            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {
                partitionT p={std::pair(Batch(_,_),1.0)};
                return p;
            }
        };

        MacroTaskCoulomb() {
            partitioner.reset(new MacroTaskPartitionerCoulomb());
        }

        // you need to define an empty constructor for the result
        // resultT must implement operator+=(const resultT&)
        resultT allocator(World &world, const argtupleT &argtuple) const {
            std::size_t n = std::get<1>(argtuple).size();
            resultT result = zero_functions_compressed<T,NDIM>(world, n);
            return result;
        }

        resultT operator()(const Function<double,NDIM>& vcoul, const std::vector<Function<T,NDIM>> &arg) const {
            return truncate(vcoul * arg);
        }
    };

    /// default empty ctor
    Coulomb(World& world) : world(world) {};

    /// default empty ctor
    Coulomb(World& world, const double lo, const double thresh=FunctionDefaults<3>::get_thresh()) : world(world) {
        reset_poisson_operator_ptr(lo,thresh);
    };

    /// ctor with an SCF calculation providing the MOs and density
    Coulomb(World& world, const SCF* calc);

    /// ctor with a Nemo calculation providing the MOs and density
    Coulomb(World& world, const Nemo* nemo);

    std::string info() const {return "J";}

    Coulomb& set_taskq(std::shared_ptr<MacroTaskQ> taskq1) {
        this->taskq=taskq1;
        return *this;
    }

    void reset_poisson_operator_ptr(const double lo, const double econv);

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
        std::vector<Function<T,NDIM> > vket(1,ket);
        return this->operator()(vket)[0];
    }

    std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vket) const {
        MacroTaskCoulomb t;
        World& world=vket.front().world();
        MacroTask task(world, t, this->taskq);
        auto result=task(vcoul,vket);
        return result;
    }

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,vcoul*ket);
    }

    Tensor<T> operator()(const std::vector<Function<T,NDIM> >& vbra,
    		const std::vector<Function<T,NDIM> >& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        std::vector<Function<T,NDIM> > vJket;
        for (std::size_t i=0; i<vket.size(); ++i) {
            vJket.push_back(this->operator()(vket[i]));
        }
        return matrix_inner(world,vbra,vJket,bra_equiv_ket);
    }

    /// getter for the Coulomb potential
    const real_function_3d& potential() const {return vcoul;}

    /// setter for the Coulomb potential
    real_function_3d& potential() {return vcoul;}

    real_function_3d compute_density(const SCF* calc) const;

    /// given a density compute the Coulomb potential

    /// this function uses a newly constructed Poisson operator. Note that
    /// the accuracy parameters must be consistent with the exchange operator.
    Function<T,NDIM> compute_potential(const Function<T,NDIM>& density) const {
    	return (*poisson)(density).truncate();
    }

    /// given a set of MOs in an SCF calculation, compute the Coulomb potential

    /// this function uses the Poisson operator of the SCF calculation
    real_function_3d compute_potential(const SCF* calc) const;

    /// given a set of MOs in an SCF calculation, compute the Coulomb potential

    /// this function uses the Poisson operator of the SCF calculation
    real_function_3d compute_potential(const Nemo* nemo) const;

private:
    World& world;
    std::shared_ptr<real_convolution_3d> poisson;
    double lo=1.e-4;
    real_function_3d vcoul; ///< the coulomb potential
};


template<typename T, std::size_t NDIM>
class Nuclear : public SCFOperatorBase<T,NDIM> {
public:

    Nuclear(World& world, const SCF* calc);

    Nuclear(World& world, const NemoBase* nemo);

    /// simple constructor takes a molecule, no nuclear correlation factor or core potentials
    Nuclear(World& world, const Molecule& molecule);

    Nuclear(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf)
        : world(world), ncf(ncf) {}

    std::string info() const {return "Vnuc";}

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
        std::vector<Function<T,NDIM> > vket(1,ket);
        return this->operator()(vket)[0];
    }

    std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vket) const;

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<T> operator()(const  std::vector<Function<T,NDIM> >& vbra,
    		const  std::vector<Function<T,NDIM> >& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        std::vector<Function<T,NDIM> > vVket=this->operator()(vket);
        return matrix_inner(world,vbra,vVket,bra_equiv_ket);
    }

private:
    World& world;
    std::shared_ptr<NuclearCorrelationFactor> ncf;

};


/// the z component of the angular momentum

/// takes real and complex functions as input, will return complex functions
template<typename T, std::size_t NDIM>
class Lz : public SCFOperatorBase<T,NDIM> {
private:
    World& world;
public:

	bool use_bsplines=true;

	Lz(World& world, bool use_bspline_derivative=true) : world(world), use_bsplines(use_bspline_derivative) {};

    std::string info() const {return "Lz";}


    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
		std::vector<Function<T,NDIM> > vket(1,ket);
        return this->operator()(vket)[0];
    }

    std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vket) const {

		// the operator in cartesian components as
		// L_z =  - i (x del_y - y del_x)

		if (vket.size()==0) return std::vector<complex_function_3d>(0);

	    real_function_3d x=real_factory_3d(world).functor([] (const coord_3d& r) {return r[0];});
	    real_function_3d y=real_factory_3d(world).functor([] (const coord_3d& r) {return r[1];});

	    Derivative<T,NDIM> Dx = free_space_derivative<T,NDIM>(world, 0);
		Derivative<T,NDIM> Dy = free_space_derivative<T,NDIM>(world, 1);
		if (use_bsplines) {
			Dx.set_bspline1();
			Dy.set_bspline1();
		}

		reconstruct(world,vket,true);
	    std::vector<Function<T,NDIM> > delx=apply(world,Dx,vket,false);
	    std::vector<Function<T,NDIM> > dely=apply(world,Dy,vket,true);

	    std::vector<Function<T,NDIM> > result1=x*dely - y*delx;
	    std::vector<complex_function_3d> cresult1=convert<T,double_complex,NDIM>(world,result1);
	    std::vector<complex_function_3d> result=double_complex(0.0,-1.0)*cresult1;
		return result;
	}

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<T> operator()(const std::vector<Function<T,NDIM> >& vbra,
    		const std::vector<Function<T,NDIM> >& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        std::vector<complex_function_3d> vVket=this->operator()(vket);
        return matrix_inner(world,vbra,vVket,bra_equiv_ket);
    }

};



/// derivative of the (regularized) nuclear potential wrt nuclear displacements
template<typename T, std::size_t NDIM>
class DNuclear : public SCFOperatorBase<T,NDIM> {
public:

    DNuclear(World& world, const SCF* calc, const int iatom, const int iaxis);

    DNuclear(World& world, const Nemo* nemo, const int iatom, const int iaxis);

    DNuclear(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf,
            const int iatom, const int iaxis)
        : world(world), ncf(ncf), iatom(iatom), iaxis(iaxis) {}

    std::string info() const {return "DVnuc";}

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
        std::vector<Function<T,NDIM>> vket(1,ket);
        return this->operator()(vket)[0];
    }

    std::vector<Function<T,NDIM>> operator()(const std::vector<Function<T,NDIM>>& vket) const;

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,this->operator()(ket));
    }

    Tensor<T> operator()(const std::vector<Function<T,NDIM>>& vbra, const std::vector<Function<T,NDIM>>& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        std::vector<Function<T,NDIM>> vVket=this->operator()(vket);
        return matrix_inner(world,vbra,vVket,bra_equiv_ket);
    }

private:
    World& world;
    std::shared_ptr<NuclearCorrelationFactor> ncf;
    int iatom;  ///< index of the atom which is displaced
    int iaxis;  ///< x,y,z component of the atom

};

template<typename T, std::size_t NDIM>
class LocalPotentialOperator : public SCFOperatorBase<T,NDIM> {
public:
    LocalPotentialOperator(World& world) : world(world) {};
    LocalPotentialOperator(World& world, const std::string info, const Function<T,NDIM> potential)
            : world(world), info_str(info), potential(potential) {};

    std::string info() const {return info_str;}

    void set_info(const std::string new_info) {
        info_str=new_info;
    }

    void set_potential(const Function<T,NDIM>& new_potential) {
        potential=copy(new_potential);
    }

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
        return (potential*ket).truncate();
    }

    std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vket) const {
        return truncate(potential*vket);
    }

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        return inner(bra,potential*ket);
    }

    Tensor<T> operator()(const std::vector<Function<T,NDIM> >& vbra,
                         const std::vector<Function<T,NDIM> >& vket) const {
        const auto bra_equiv_ket = &vbra == &vket;
        return matrix_inner(world,vbra,potential*vket,bra_equiv_ket);
    }

private:
    World& world;
    std::string info_str="Vlocal";
    Function<T,NDIM> potential;
};

/// operator class for the handling of DFT exchange-correlation functionals
/// divergence of a real vector field, honouring the `dft_deriv` setting

/// vmra.h's div() hardcodes a plain ABGV gradient, so any caller that wants
/// `dft_deriv bspline`/`ble` to reach a divergence must come through here.
/// Shared by XCOperator::div_dft_deriv and by nemo's flux_bsh_term, which
/// otherwise disagreed on the derivative operator and made the two weak-form
/// routes incomparable -- see XC-NEMO-MEMORY-INVESTIGATION.md 9.4.
/// `do_refine` mirrors vmra.h's div(v, do_refine): refining before
/// differentiating helps a flux built from projected densities, and hurts one
/// that is already the smooth output of a Green's-function apply -- measured 4x
/// wall and 2x memory on He/PBE when switched on at the flux_bsh_term call site.
real_function_3d div_dft_deriv(World& world, const std::vector<real_function_3d>& v,
                               const std::string& dft_deriv, bool do_refine = true);

/// how XCOperator::set_tau() evaluates the two U1 terms of the product rule

/// With a nuclear correlation factor psi = R F and
/// |grad psi|^2 = R^2 (|grad F|^2 - 2 F U1.grad F + |U1|^2 F^2). U1 = -grad(R)/R is
/// analytic but componentwise non-smooth at each nucleus (U1_x ~ x/r), so carrying
/// it as an MRA Function costs depth ~18 and every product with it inherits that
/// depth -- which refine_to_common_level then imposes on every xc intermediate.
/// See XC-NEMO-MEMORY-INVESTIGATION.md 9.2 and 9.8.
enum class TauU1 {
    from_env,    ///< MAD_XC_TAU_NO_U1 / MAD_XC_TAU_U1_POINTWISE decide; the default
    mra,         ///< U1 and |U1|^2 projected into Functions and multiplied (shipped route)
    pointwise,   ///< U1 evaluated from its functor at the quadrature points of the orbital tree
    none         ///< drop both U1 terms -- DIAGNOSTIC ONLY, gives a wrong tau
};

template<typename T, std::size_t NDIM>
class XCOperator : public SCFOperatorBase<T,NDIM> {
public:

    /// default ctor without information about the XC functional
    XCOperator(World& world) : world(world), nbeta(0), ispin(0),
        extra_truncation(FunctionDefaults<3>::get_thresh()*0.01) {}

    /// custom ctor with information about the XC functional
    XCOperator(World& world, std::string xc_data, const bool spin_polarized,
            const real_function_3d& arho, const real_function_3d& brho,
            std::string deriv="abgv");

    /// custom ctor with the XC functional
    XCOperator(World& world, std::shared_ptr<XCfunctional> xc,
               const bool spin_polarized,
               int ispin,
               int nbeta,
               const real_function_3d& arho, const real_function_3d& brho,
               std::string deriv="abgv");

    /// ctor with an SCF calculation, will initialize the necessary intermediates
    XCOperator(World& world, const SCF* scf, int ispin=0, std::string deriv="abgv");

    /// ctor with a Nemo calculation, will initialize the necessary intermediates
    XCOperator(World& world, const Nemo* nemo, int ispin=0);

    /// ctor with an SCF calculation, will initialize the necessary intermediates
    XCOperator(World& world, const SCF* scf, const real_function_3d& arho,
            const real_function_3d& brho, int ispin=0, std::string deriv="abgv");

    /// ctor with an Nemo calculation, will initialize the necessary intermediates
    XCOperator(World& world, const Nemo* scf, const real_function_3d& arho,
            const real_function_3d& brho, int ispin=0);

    std::string info() const {return "Vxc";}

    XCOperator& set_extra_truncation(const double& fac) {
        extra_truncation=fac;
        if (world.rank()==0)
            print("set extra truncation in XCOperator to", extra_truncation);
        return *this;
    }

    /// set the spin state this operator is acting on
    void set_ispin(const int i) const {ispin=i;}

    /// print the meta-gga de/dtau range each time the operator is applied
    XCOperator& set_print_level(const int p) {print_level=p; return *this;}

    /// apply the xc potential on a set of orbitals
    std::vector<Function<T,NDIM> > operator()(const std::vector<Function<T,NDIM> >& vket) const;

    /// apply the xc potential on an orbitals
    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
    	std::vector<Function<T,3> > vket(1,ket);
    	std::vector<Function<T,3> > vKket=this->operator()(vket);
        return vKket[0];
    }

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
         MADNESS_EXCEPTION("no implementation of matrix elements of the xc operator",1);
    };

    Tensor<T> operator()(const std::vector<Function<T,NDIM>>& vbra, const std::vector<Function<T,NDIM>>& vket) const {
        MADNESS_EXCEPTION("no implementation of matrix elements of the xc operator", 1);
    }

    /// compute the xc energy using the precomputed intermediates vf and delrho
    double compute_xc_energy() const;

    /// return the local xc potential
    real_function_3d make_xc_potential() const;

    /// true if the functional contributes a non-multiplicative (meta-gga) term
    bool has_tau_term() const;

    /// compute the kinetic energy density and add it to the intermediates

    /// tau is orbital-dependent, so unlike the density it cannot be recovered
    /// from what the ctors are given -- it has to be supplied separately. Call
    /// this after construction and before make_xc_potential() whenever
    /// has_tau_term() is true; make_xc_potential() throws otherwise.
    /// The occupation numbers are required, not optional: amo/bmo may carry
    /// virtual orbitals (occupation zero), which contribute nothing to the
    /// density but would inflate an unweighted sum of |grad psi|^2, and
    /// occupations may be fractional.
    /// With a nuclear correlation factor the vectors are the nemos F, psi = R F,
    /// and the product rule grad psi = R (grad F - U1 F) is applied here -- the
    /// cusp stays in the analytic U1 and only the cusp-free F is differentiated.
    /// @param[in]  amo  alpha orbitals (nemos if a nuclear correlation factor is set)
    /// @param[in]  aocc occupation numbers of amo
    /// @param[in]  bmo  beta orbitals, ignored if the calculation is spin-restricted
    /// @param[in]  bocc occupation numbers of bmo
    /// @param[in]  u1mode how the two U1 terms are evaluated; see TauU1
    void set_tau(const vecfuncT& amo, const Tensor<double>& aocc,
                 const vecfuncT& bmo=vecfuncT(),
                 const Tensor<double>& bocc=Tensor<double>(),
                 TauU1 u1mode=TauU1::from_env) const;

    /// inject a kinetic energy density computed elsewhere

    /// For diagnostics: it lets one potential be assembled from several taus on a
    /// single density, so a difference between two potentials is the tau and not
    /// the orbitals. Note make_libxc_args floors tau from below at the von
    /// Weizsaecker bound rho chi/8, so tau = 0 does not reach libxc as zero.
    /// @param[in]  taua alpha kinetic energy density
    /// @param[in]  taub beta kinetic energy density; required iff spin-polarized
    ///                  with beta electrons
    void set_tau(const real_function_3d& taua,
                 const real_function_3d& taub=real_function_3d()) const;

    /// true once set_tau() has supplied tau, by either route

    /// The pointwise route never builds a tau Function, so the presence of
    /// enum_taua is not the right test any more.
    bool has_tau_args() const;

    /// the four analytic U1 quantities the xc ops evaluate at their quadrature points

    /// Empty unless the pointwise route is in use. Handed to xc_functional /
    /// xc_potential rather than projected into xc_args, because a projected product
    /// with U1 is exactly what rings.
    nemo_u1_functors make_u1_functors() const;

    /// the kinetic energy density of one spin channel, as set by set_tau()

    /// exposed for diagnostics and for the exact check int(tau) == T.
    /// On the pointwise route tau is not stored as a Function: this assembles one,
    /// which reintroduces the projection the route avoids, so the result is a
    /// diagnostic and not what the functional saw near a nucleus.
    real_function_3d get_tau(const int spin=0) const;

    /// de/dtau, as computed by make_xc_potential()
    real_function_3d get_vtau() const {return vtau;}

    /// 2 de/dsigma (same spin), if MAD_XC_EXPORT_DFDS asked for it

    /// Diagnostic. Uninitialized unless the switch is on -- it is the coefficient of
    /// the semilocal flux X = u grad(rho), i.e. where div(X)'s 1/r originates, and
    /// nothing in the potential reads it back. Requires make_xc_potential() first.
    real_function_3d get_dfdsigma() const {return dfdsigma;}

    /// opt in to the weak form, if MAD_XC_WEAK_GGA asks for it

    /// Load-bearing: in weak form make_xc_potential() returns only de/drho, so a
    /// caller that does not also apply weak_xc_terms()/weak_xc_matrix() would
    /// silently drop the whole semilocal contribution and return a plausible but
    /// wrong energy. Only Nemo::compute_nemo_potentials implements the split, so
    /// only it opts in; SCF, OEP, TDHF and the response kernels keep the
    /// multiplicative potential whatever the environment says.
    XCOperator& allow_weak_form() {weak_form_ok=true; return *this;}

    /// true if this operator is running in weak form, i.e. make_xc_potential()
    /// returns only de/drho and the flux is carried separately
    bool is_weak_form() const;

    /// the semilocal flux X_sigma = 2 (de/dsigma_ss) grad(rho_s) + (de/dsigma_ab) grad(rho_s')

    /// only assigned in weak form, by make_xc_potential(). The same-spin and
    /// cross-spin contributions are summed: they enter the potential as a single
    /// divergence, so nothing needs them apart.
    const vecfuncT& get_semilocal_flux() const {return semilocal_flux;}

    /// weak-form split of the non-multiplicative xc terms

    /// Writes the decomposition
    /// \f[
    ///   \hat v_{xc}^{semilocal+\tau}\psi_i = \mathrm{mult}_i - \nabla\cdot\mathbf Y_i
    /// \f]
    /// with (nemo kets \f$ F_i \f$, \f$ W_i = v_\tau(\nabla F_i - U_1 F_i) \f$)
    /// \f[
    ///   \mathrm{mult}_i = \mathbf X\cdot\nabla F_i + \tfrac12 U_1\cdot\mathbf W_i,
    ///   \qquad \mathbf Y_i = \mathbf X F_i + \tfrac12 \mathbf W_i .
    /// \f]
    /// `mult` is an ordinary multiplicative-style term and goes into \f$ V\psi \f$;
    /// `flux` is what the caller pushes through the Green's function, as
    /// \f$ 2\nabla\cdot(G*\mathbf Y_i) \f$, so that neither \f$\mathbf X\f$ nor
    /// \f$ v_\tau \f$ is ever differentiated. Requires make_xc_potential() first.
    void weak_xc_terms(const std::vector<Function<T,NDIM> >& vket,
                       std::vector<Function<T,NDIM> >& mult,
                       std::vector<std::vector<Function<T,NDIM> > >& flux) const;

    /// xc contribution to the Fock matrix in weak form

    /// \f[
    ///   F_{ij} = \langle\psi_i|\tfrac{\partial f}{\partial\rho}|\psi_j\rangle
    ///          + \int \mathbf X\cdot\nabla(\psi_i\psi_j)
    ///          + \tfrac12\int v_\tau \nabla\psi_i\cdot\nabla\psi_j
    /// \f]
    /// symmetric by construction, and no derivative touches \f$\mathbf X\f$ or
    /// \f$ v_\tau \f$. With a nuclear correlation factor the kets are the nemos and
    /// \f$ \nabla(\psi_i\psi_j) = R^2[F_j\nabla F_i + F_i\nabla F_j - 2U_1F_iF_j] \f$,
    /// so the cusp stays in the analytic \f$ U_1 \f$.
    /// @param[in] vket    the orbitals (nemos if an ncf is set)
    /// @param[in] v_local de/drho, i.e. what make_xc_potential() returned
    Tensor<T> weak_xc_matrix(const std::vector<Function<T,NDIM> >& vket,
                             const real_function_3d& v_local) const;

    /// apply the non-multiplicative meta-gga term on a set of orbitals

    /// \f[
    ///   \hat v_\tau \psi_i = -\frac{1}{2}\nabla\cdot
    ///        \left(\frac{\partial e_{xc}}{\partial\tau_\sigma}\nabla\psi_i\right)
    /// \f]
    /// evaluated as \f$ -\frac{1}{2}\sum_x D_x(v_\tau D_x\psi_i) \f$, so that
    /// \f$ v_\tau \f$ is only ever multiplied and never differentiated. That
    /// nested form is also self-adjoint by construction, while the expanded
    /// \f$ -\frac{1}{2}(v_\tau\nabla^2\psi + \nabla v_\tau\cdot\nabla\psi) \f$
    /// is symmetric only up to discretization error and needs \f$\nabla^2\psi\f$.
    /// Requires make_xc_potential() to have been called first, which is where
    /// \f$ v_\tau \f$ is computed.
    ///
    /// With a nuclear correlation factor the kets are the nemos \f$ F \f$ and the
    /// result is the term already divided by \f$ R \f$, which is what nemo's
    /// equations for \f$ F \f$ need. The \f$ R \f$ factors cancel analytically
    /// (see the implementation), so nothing is ever divided by \f$ R \f$.
    std::vector<Function<T,NDIM> > apply_tau_term(const std::vector<Function<T,NDIM> >& vket) const;

    /// construct the xc kernel and apply it directly on the (response) density

    /// the xc kernel is the second derivative of the xc functions wrt the density
    /// @param[in]  density the (response) density on which the kernel is applied
    /// @return     kernel * density
    real_function_3d apply_xc_kernel(const real_function_3d& density,
            const vecfuncT grad_dens_pt=vecfuncT()) const;

private:

    /// the world
    World& world;

    /// which derivative operator to use
    std::string dft_deriv;

    /// print level; >=2 logs the meta-gga de/dtau range
    int print_level=0;

public:
    /// interface to the actual XC functionals
    std::shared_ptr<XCfunctional> xc;

private:
    /// number of beta orbitals
    int nbeta;

    /// the XC functionals depend on the spin of the orbitals they act on
    mutable int ispin;

    /// additional truncation for the densities in the XC kernel

    /// the densities in the DFT kernal are processed as their inverses,
    /// so noise in the small density regions might amplify and lead to inaccurate
    /// results. Extra truncation will tighten the truncation threshold by a
    /// specified factor, default is 0.01.
    double extra_truncation;

    /// the nuclear correlation factor, if it exists, for computing derivatives for GGA
    std::shared_ptr<NuclearCorrelationFactor> ncf;

    /// functions that are need for the computation of the XC operator

    /// the ordering of the intermediates is fixed, but the code can handle
    /// non-initialized functions, so if e.g. no GGA is requested, all the
    /// corresponding vector components may be left empty.
    /// For the ordering of the intermediates see xcfunctional::xc_arg
    mutable vecfuncT xc_args;

    /// caller has opted in to the weak form, see allow_weak_form()
    bool weak_form_ok=false;

    /// the regularized (nemo) alpha density rho_a/R^2, if the ctor was given it

    /// needed by algorithm A'' to build grad(rho) = R^2 (grad n - 2 U1 n), i.e.
    /// make_ddensity's route, whose divergence conserves; the alternative
    /// rho*zeta carries zeta's tail garbage to the box wall.
    mutable real_function_3d arho_reg;

    /// the semilocal flux, assigned by make_xc_potential() in weak form only
    mutable vecfuncT semilocal_flux;

    /// de/dtau, the prefactor of the non-multiplicative meta-gga term

    /// falls out of the same pointwise pass as the multiplicative potential, so
    /// it is stashed by make_xc_potential() rather than recomputed
    mutable real_function_3d vtau;

    /// 2 de/dsigma, same spin; only assigned when MAD_XC_EXPORT_DFDS is set
    mutable real_function_3d dfdsigma;

    /// gradient operator honouring dft_deriv, for the meta-gga term
    std::shared_ptr<Derivative<T,NDIM> > make_derivative(const int axis) const;

    /// divergence of a vector field, honouring dft_deriv

    /// vmra.h's div() hardcodes a plain ABGV gradient, so the `dft_deriv` setting
    /// reached grad(log rho) -- the well-conditioned derivative -- and never the
    /// flux divergence, which is the badly conditioned one.
    real_function_3d div_dft_deriv(const vecfuncT& v) const;


    /// compute the intermediates for the XC functionals

    /// @param[in]  arho    density of the alpha orbitals
    /// @param[in]  brho    density of the beta orbitals (necessary only if spin-polarized)
    /// @param[in]  arho_reg  regularized alpha density rho_a/R^2, nemo only, optional
    /// @param[in]  brho_reg  regularized beta density rho_b/R^2, nemo only, optional
    /// @return xc_args vector of intermediates as described above

    /// If the regularized densities are supplied, zeta = grad log(rho) is built as
    /// grad log(rho_reg) - 2 U1 -- exact, and it keeps the nuclear cusp of
    /// rho = R^2 rho_reg out from under the numerical derivative.
    vecfuncT prep_xc_args(const real_function_3d& arho, const real_function_3d& brho,
                          const real_function_3d& arho_reg = real_function_3d(),
                          const real_function_3d& brho_reg = real_function_3d()) const;

    /// compute the intermediates for the XC functionals

    /// @param[in]  dens_pt     perturbed densities from CPHF or TDDFT equations
    /// @param[in,out] xc_args   vector of intermediates as described above
    /// @param[out] ddens_pt    xyz-derivatives of dens_pt
    void prep_xc_args_response(const real_function_3d& dens_pt,
            vecfuncT& xc_args, vecfuncT& ddens_pt) const;

    /// check if the intermediates are initialized
    bool is_initialized() const {
        return (xc_args.size()>0);
    }

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

    /// simple structure to take the pointwise exponential of a function, shifted by +14
    struct expme{
        typedef double resultT;
        struct expme1 {
            double operator()(const double& val) {return exp(val-14.0);}
        };
        Tensor<double> operator()(const Key<3>& key, const Tensor<double>& val) const {
            Tensor<double> result=copy(val);
            expme1 op;
            return result.unaryop(op);
        }

        template <typename Archive>
        void serialize(Archive& ar) {}

    };
};

/// Computes matrix representation of the Fock operator
template<typename T, std::size_t NDIM>
class Fock : public SCFOperatorBase<T,NDIM> {
public:
    Fock(World& world) : world(world) {}

    Fock(World& world, const Nemo* nemo);
    Fock(World& world, const OEP* nemo);
    Fock(World& world, const NemoBase* nemo);

    /// pretty print what this is actually computing
    std::string info() const {
        std::string s;
        for (auto& op : operators) {
            double number=std::get<0>(op.second);
            if (number==-1.0) {
                s+=" - ";
            } else if (number!=1.0) {
                std::stringstream snumber;
                snumber << std::fixed << std::setw(2) << number;
                s+=" "+snumber.str()+ " ";
            } else {
                MADNESS_CHECK(number==1.0);
                s+=" + ";
            }
            s+=op.first;
        }
        return s;
    }

    /// add an operator with default prefactor 1.0
    void add_operator(std::string name, std::shared_ptr<SCFOperatorBase<T,NDIM>> new_op) {
        operators.insert({name,valueT(1.0,new_op)});
    }

    /// add an operator with custom prefactor (e.g. -1.0 for the exchange, supposedly)
    void add_operator(std::string name, std::tuple<double,std::shared_ptr<SCFOperatorBase<T,NDIM>>> new_op) {
        operators.insert({name,new_op});
    }

    /// remove operator, returns 0 if no operator was found
    int remove_operator(std::string name) {
        return operators.erase(name);
    }

    Function<T,NDIM> operator()(const Function<T,NDIM>& ket) const {
      MADNESS_EXCEPTION("Fock(ket) not yet implemented",1);
      Function<T,NDIM> result;
      return result;
    }

    std::vector<Function<T,NDIM>> operator()(const std::vector<Function<T,NDIM>>& vket) const {
        // make sure T is not part of the Fock operator, it's numerically unstable!
        MADNESS_CHECK(operators.count("T")==0);
        std::vector<Function<T,NDIM>> result = zero_functions_compressed<T, NDIM>(world, vket.size());
        for (const auto& op : operators) {
            result+=std::get<0>(op.second) * (*std::get<1>(op.second))(vket);
        }
        return result;
    }

    T operator()(const Function<T,NDIM>& bra, const Function<T,NDIM>& ket) const {
        std::vector<Function<T,NDIM>> vbra(1,bra), vket(1,ket);
        return (*this)(vbra,vket)(0,0);
    }

    /// compute the Fock matrix by summing up all contributions
    Tensor<T> operator()(const std::vector<Function<T,NDIM>>& vbra, const std::vector<Function<T,NDIM>>& vket) const {
        return this->operator()(vbra,vket,false);
    }

    /// compute the Fock matrix by summing up all contributions
    Tensor<T> operator()(const std::vector<Function<T,NDIM>>& vbra, const std::vector<Function<T,NDIM>>& vket,
            const bool symmetric) const {
        Tensor<T> fock(vbra.size(),vket.size());
        for (const auto& op : operators) {
            Tensor<T> tmp=std::get<0>(op.second) * (*std::get<1>(op.second))(vbra,vket);
//            print("Operator",std::get<1>(op.second)->info());
//            print(tmp);
            fock+=tmp;
        }
        return fock;
    }


private:
    /// the world
    World& world;

    /// type defining Fock operator contribution including prefactor
    typedef std::tuple<double,std::shared_ptr<SCFOperatorBase<T,NDIM> > > valueT;

    /// all the Fock operator contribution
    std::map<std::string,valueT> operators;
};

}
#endif /* MADNESS_CHEM_SCFOPERATORS_H_ */
