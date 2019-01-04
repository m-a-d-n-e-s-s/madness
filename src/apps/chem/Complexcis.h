/*
 * Complexcis.h
 *
 *  Created on: 21 Nov 2018
 *      Author: fbischoff
 */

#ifndef SRC_APPS_CHEM_COMPLEXCIS_H_
#define SRC_APPS_CHEM_COMPLEXCIS_H_

#include <chem/Nemocomplex.h>
#include <chem/projector.h>
#include<tuple>


namespace madness {



class Complex_CIS_Parameters : public CalculationParametersBase {
public:
	enum parameterenum {guess_excitation_operators_,exops_,freeze_,guess_excitations_,thresh_,maxiter_,omega_,
						swap_ab_};

	/// the parameters with the enum key, the constructor taking the input file key and a default value
	ParameterMap params={
        		init<std::string>(guess_excitation_operators_,{"guess_excitation_operators",{"dipole+"}}),
        		init<std::vector<std::string> >(exops_,{"exops",{"x 1.0","y 1.0","z 1.0","x 2.0 , y 2.0 , z 2.0"}}),
        		init<int>(freeze_,{"freeze",0}),
        		init<int>(guess_excitations_,{"guess_excitations",4}),
        		init<double>(thresh_,{"thresh",FunctionDefaults<3>::get_thresh()}),
        		init<double>(omega_,{"omega",0.0}),
        		init<int>(maxiter_,{"maxiter",10}),
				init<bool>(swap_ab_,{"swap_ab",false})
    };

	/// ctor reading out the input file
	Complex_CIS_Parameters(World& world) {

		// read input file
		read(world,"input","response",params);

		// set derived values
//		params[param2_].set_derived_value(this->get<int>(param1_)*10.0);

		// print final parameters
//		if (world.rank()==0) print(params,"Our parameters");
	}

	std::string guess_excitation_operators() const {return get<std::string>(guess_excitation_operators_);};
	std::vector<std::string> exops() const {return get<std::vector<std::string> >(exops_);};
	int freeze() const {return get<int>(freeze_);};
	int guess_excitations() const {return get<int>(guess_excitations_);};
	double thresh() const {return get<double>(thresh_);};
	double omega() const {return get<double>(omega_);};
	int maxiter() const {return get<int>(maxiter_);};
	bool swap_ab() const {return get<bool>(swap_ab_);};


	/// return the value of the parameter
	template<typename T>
	T get(parameterenum k) const {
		if (params.find(int(k))!=params.end()) {
			return params.find(int(k))->second.get_parameter<T>().get();
		} else {
			MADNESS_EXCEPTION("could not fine parameter ",1);
		}
	}
};


class Complex_cis {
public:

	struct root {
		std::vector<complex_function_3d> afunction;
		std::vector<complex_function_3d> bfunction;
		std::vector<complex_function_3d> apot;
		std::vector<complex_function_3d> bpot;
		double omega=0.0;
		double delta=0.0;				// last wave function error
		double energy_change=0.0;		// last energy_change
	};

	static std::vector<root> transform(World& world,
			  const std::vector<root>& v,
			  const Tensor<double_complex>& c,
			  bool fence=true) {

		int n = v.size();  // n is the old dimension
		int m = c.dim(1);  // m is the new dimension

		int nocca=v.front().afunction.size();
		int noccb=v.front().bfunction.size();

		std::vector<root> result(m);

		for (auto&& r : result) {
			r.afunction= zero_functions_compressed<double_complex,3>(world, nocca, false);
			r.bfunction= zero_functions_compressed<double_complex,3>(world, noccb, false);
			r.apot= zero_functions_compressed<double_complex,3>(world, nocca, false);
			r.bpot= zero_functions_compressed<double_complex,3>(world, noccb, false);
		}
		for (auto&& vv : v) {
			compress(world, vv.afunction,false);
			compress(world, vv.bfunction,false);
			compress(world, vv.apot,false);
			compress(world, vv.bpot,false);
		}
		world.gop.fence();

		for (int i=0; i<m; ++i) {
			for (int j=0; j<n; ++j) {
				gaxpy(world,double_complex(1.0),result[i].afunction,c(j,i),v[j].afunction,false);
				gaxpy(world,double_complex(1.0),result[i].bfunction,c(j,i),v[j].bfunction,false);
				if (v[j].apot.size()>0) gaxpy(world,double_complex(1.0),result[i].apot,c(j,i),v[j].apot,false);
				if (v[j].bpot.size()>0) gaxpy(world,double_complex(1.0),result[i].bpot,c(j,i),v[j].bpot,false);
			}
		}

		if (fence) world.gop.fence();
		return result;
	}



	Complex_cis(World& w, Nemo_complex& n) : world(w), cis_param(world), nemo(n),
		Qa(world,conj(world,nemo.amo),nemo.amo), Qb(world,conj(world,nemo.bmo),nemo.bmo) {
		print("Qa projector",Qa.get_ket_vector().size());
		print("Qb projector",Qb.get_ket_vector().size());

	}

	virtual ~Complex_cis() {};

	double value();

	void iterate(std::vector<root>& roots) const;

	void compute_potentials(std::vector<root>& roots, const real_function_3d& totdens) const;

	std::vector<complex_function_3d> compute_residuals(root& root) const;

	void update_roots(std::vector<root>& aroot, std::vector<root>& broot, std::vector<root>& troot) const;

	std::vector<root> read_guess(const std::string spin) const;

	std::vector<root> make_guess() const;

	void canonicalize(const std::vector<complex_function_3d>& mo, const real_function_3d& density,
			std::vector<complex_function_3d>& virtuals, Tensor<double>& veps) const;

	Tensor<double_complex> make_CIS_matrix(const Tensor<double>& veps, const Tensor<double>& oeps) const;

	Tensor<double_complex> compute_fock_pt(const std::vector<root>& roots) const;


	/// return the active orbitals only
	Tensor<double> noct(const Tensor<double>& eps) const {
		if (eps.size()<=cis_param.freeze()) return Tensor<double>();
		return eps(Slice(cis_param.freeze(),-1,1));
	}

	/// return the active orbitals only
	std::vector<complex_function_3d> active_mo(const std::vector<complex_function_3d>& mo) const {
		if (mo.size()<=cis_param.freeze()) return std::vector<complex_function_3d> ();
		std::vector<complex_function_3d> result;
		result.insert(result.end(),mo.begin()+cis_param.freeze(),mo.end());
		return result;
	}

	/// little helper function
	template<typename T>
	static Tensor<T> concatenate(const Tensor<T>& t1, const Tensor<T>& t2) {
		MADNESS_ASSERT(t1.ndim()==1 or t1.size()==0);
		MADNESS_ASSERT(t2.ndim()==1 or t2.size()==0);
		Tensor<T> result(t1.size()+t2.size());
		if (t1.size()>0) result(Slice(0,t1.size()-1,1))=t1;	// slices count inclusive
		if (t2.size()>0) result(Slice(t1.size(),-1,1))=t2;
		return result;
	}

	template<typename T>
	static std::tuple<Tensor<T>, Tensor<T> > split(const Tensor<T>& rhs, int dim1) {

		MADNESS_ASSERT(dim1>=0 and dim1<=rhs.size());
		if (dim1==0) return std::make_tuple(Tensor<T>(),rhs);
		if (dim1==rhs.size()) return std::make_tuple(rhs,Tensor<T>());

		Tensor<T> t1=rhs(Slice(0,dim1-1,1));
		Tensor<T> t2=rhs(Slice(dim1,-1,1));
		return std::make_tuple(t1,t2);
	}

	void orthonormalize(std::vector<root>& roots, const Tensor<double_complex>& fock_pt_a) const;

	void orthonormalize(std::vector<root>& roots) const;

	void normalize(std::vector<root>& roots) const;

	void compare_to_file(const std::vector<complex_function_3d>& rhs, const std::string name) const {
		if (nemo.cparam.spin_restricted) {
			save_function(rhs,name);

		} else {
			std::vector<complex_function_3d> rhs_file=zero_functions_compressed<double_complex,3>(world,rhs.size());
			std::vector<complex_function_3d> rhs_file1;
			load_function(world,rhs_file1,name);
			for (int i=0; i<rhs_file1.size(); ++i) rhs_file[i]=rhs_file1[i];
			std::vector<double> dnorm=norm2s(world,rhs-rhs_file);
			print(name,"diffnorm",dnorm);
		}
	}


	static std::tuple<std::vector<complex_function_3d>, std::vector<complex_function_3d> >
	split(const std::vector<complex_function_3d>& rhs, int dim1) {

		MADNESS_ASSERT(dim1>=0 and dim1<=rhs.size());
		if (dim1==0) return std::make_tuple(std::vector<complex_function_3d>(),rhs);
		if (dim1==rhs.size()) return std::make_tuple(rhs,std::vector<complex_function_3d>());

		std::vector<complex_function_3d> t1,t2;
		copy(rhs.begin(),rhs.begin()+dim1,back_inserter(t1));
		copy(rhs.begin()+dim1,rhs.end(),back_inserter(t2));

		return std::make_tuple(t1,t2);
	}

	/// the world
	World& world;

	/// the parameters
	Complex_CIS_Parameters cis_param;

	/// the reference
	Nemo_complex nemo;

	/// orthogonality projector
	QProjector<double_complex,3> Qa, Qb;

	/// the x vectors
	std::vector<root> roots;

};

} /* namespace madness */

#endif /* SRC_APPS_CHEM_COMPLEXCIS_H_ */
