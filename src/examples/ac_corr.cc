/*
 * ac_corr.cpp
 *
 *  Created on: Nov 17, 2016
 *      Author: msahre
 */

#include <madness/world/MADworld.h>
#include <madness/mra/mra.h>
#include <madness/mra/funcplot.h>
#include <exception>
#include <iterator>
#include <list>
//#include <chem/molecule.h>
#include <madness/world/info.h>
#include <chem/AC.h>



using namespace madness;

/// Asymptotic correction for DFT. In the correction the xc-potential is replaced by an 1/r term far away from the nuclei
/// to give the correct asymptotic behavior. Close to the nuclei the standard xc-potential is used. The transition between
/// the different potentials is achieved via a linear interpolation.
/// This is a test code to compute the corrected potential in 1D or 2D.
/// The density functional is computed using the slater potential. The orbitals have to be inserted manually in the xc_functor
/// class in the density() function.
/// The molecule/atom has to be initialized manually in the main function and has to be stored as a vector of atom_information.


/// Functor for the exchange correlation potential
template<unsigned long int NDIM>
class xc_functor : public FunctionFunctorInterface<double,NDIM>{

public:
	xc_functor(){}
	xc_functor(std::vector<atom_information<NDIM> > atoms): atoms(atoms){}

	/// returns the slater potential
	double operator ()(const Vector<double, NDIM> &r)const{
		return xc_potential(r);
	}

private:
	/// computes distance between coordinate of electron and coordinate of nucleus
	double get_distance(Vector<double,NDIM> elec, Vector<double,NDIM> nuc) const{
		double distance = 0.0;
		for(unsigned i = 0; i < NDIM; i++){
			distance += (elec[i]-nuc[i])*(elec[i]-nuc[i]);
		}
		distance = sqrt(distance);

		return distance;
	}

	/// forms the density
	double density (const Vector<double, NDIM> &r) const{
		double dist = get_distance(r, atoms[0].coord)*get_distance(r, atoms[0].coord);

		//double sto2g_h = 0.430128498*exp(-1.309756377*get_distance(r, atoms[0].coord)*get_distance(r, atoms[0].coord))+ 0.678913531*exp(-0.233135974*get_distance(r, atoms[0].coord)*get_distance(r, atoms[0].coord));
		//double density = sto2g_h*sto2g_h;
		//double density = sto2g_h*sto2g_h;
		double sto2g_he = 0.4301280*exp(-2.4328790*get_distance(r, atoms[0].coord)*get_distance(r, atoms[0].coord))+ 0.6789140*exp(-0.4330510*get_distance(r, atoms[0].coord)*get_distance(r, atoms[0].coord));
		double density = 1.0*sto2g_he*sto2g_he;
		//double be_1s = 0.4301280*exp(-11.5356690*dist)+ 0.6789140*exp(-2.0533430*dist);
		//double be_2s = 0.0494720*exp(-0.5081630*dist)+ 0.9637820s*exp(-0.1288840*dist);

		return 2.0*density;
	}

	/// equation for the slater potential
	double xc_potential(const Vector<double, NDIM> &r) const{
		return -pow((3.0/M_PI)*density(r), (1.0/3.0));
	}

	/// Needed information about the molecule to apply asymtotic correction
	std::vector<atom_information<NDIM> > atoms;

};



int main(int argc, char** argv){

	const long k = 7;
	const double thresh = 1.e-5;
	const double L = 50.0;

	/// interval limits for interpolation region
    double lim1 = 3.0;
    double lim2 = 4.0;

    madness::initialize(argc,argv);
    madness::World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);
    madness::print("Hello from processor",world.rank());

    ///////////////////////////////////////////////////////////////////
    //																 //
    //                      Test for 1D functionals					 //
    //																 //
    ///////////////////////////////////////////////////////////////////

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_truncate_mode(3);

    /// Coordinates of atoms
	// atom 1
    Vector<double, 1> coord1;
	coord1=0.0;
//	// atom 2
//	Vector<double, 1> coord2;
//	coord2=1.0;
//	//atom 3
//	Vector<double, 1> coord3;
//	coord3=2.0;

	/// Vector for the atom information for all atoms of the molecule
	std::vector<atom_information<1> > atomvec;

	/// atom information atom 1
	atom_information<1> atom1;
	atom1.coord = coord1;
	atom1.charge = 2;
	atom1.R1 = slater_radius(atom1.charge)*lim1;
	atom1.R2 = slater_radius(atom1.charge)*lim2;

//	/// atom information atom 2
//	atom_information<1> atom2;
//	atom2.coord = coord2;
//	atom2.charge = 2;
//	atom2.R1 = slater_radius(atom2.charge)*lim1;
//	atom2.R2 = slater_radius(atom2.charge)*lim2;

//	/// atom information atom 3
//	atom_information<1> atom3;
//	atom3.coord = coord3;
//	atom3.charge = 2;
//	atom3.R1 = slater_radius(atom3.charge)*lim1;
//	atom3.R2 = slater_radius(atom3.charge)*lim2;

	/// put atoms in atom_information vector
	atomvec.push_back(atom1);
	//atomvec.push_back(atom2);
	//atomvec.push_back(atom3);

	/// Create ACParameters object (necessary to use AC class)
	ACParameters<1> param;
	param.atoms_ = atomvec;
	param.R1_ = lim1;
	param.R2_ = lim2;
	param.dft_coefficient_ = 1.0;
	param.e_ion_ = 0.0;
	param.eh_ = 0.0;
	param.interpolation_scheme_ = "linear";
	param.num_elec_ = 2.0;
	param.use_mult_ = true;

	/// create an object of the AC class to calculate the correction
	AC<1> ac(param);

	/// make xc_potential 1D
	/// WARNING: You have to change the code of the density function in the XC functor class to get the right density for your molecule
	std::cout << "Computing density...\n";
	std::shared_ptr<FunctionFunctorInterface<double, 1> > xc_ptr(new xc_functor<1>(atomvec));
	real_function_1d xc_pot = real_factory_1d(world).functor(xc_ptr);
	real_function_1d xc_ref = real_factory_1d(world).functor(xc_ptr);
	plot_plane(world, xc_pot, "test_xc_potential");

	/// apply correction
	xc_pot = ac.apply(xc_pot);
	std::cout << "Plotting corrected potential...\n";
	plot_plane(world, xc_pot, "test_total");

	/// difference between corrected and uncorrected standard potential
	real_function_1d int_total = xc_pot - xc_ref;
	plot_plane(world, int_total, "test_diff");


    ///////////////////////////////////////////////////////////////////
    //																 //
    //                      Test for 2D functionals					 //
    //																 //
    ///////////////////////////////////////////////////////////////////
//    FunctionDefaults<2>::set_k(k);
//    FunctionDefaults<2>::set_thresh(thresh);
//    FunctionDefaults<2>::set_cubic_cell(-L,L);
//    FunctionDefaults<2>::set_truncate_mode(3);

    madness::finalize();
	return 0;
}




