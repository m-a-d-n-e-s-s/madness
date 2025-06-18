/*
 *
 *   Written by: bsundahl and jscanderson
 *   Date: A long time ago...
 *
 * TODO:
 *  - Reading input from NWChem needs to be updated for the spin-restricted version of the code
 *  - Try setting truncate mode to 0 just for the final calculation of the energy 
 *  - 
 *
 */ 

#include "DF.h"
//#include "Plot_VTK.h"
#include "fcwf.h"
#include <madness/chem/potentialmanager.h>

using namespace madness;

//Just a function for f(x,y,z) = x
double myxfunc(const madness::coord_3d& r){
     return r[0];
}

//Just a function for f(x,y,z) = y
double myyfunc(const madness::coord_3d& r){
     return r[1];
}

// A function that constructs a cost tree, which is a heuristic used for load balancing
template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() < 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};

// Pulled from SCF.cc, starts a timer
void DF::start_timer(World& world)
{
   world.gop.fence();
   ttt.push_back(wall_time());
   sss.push_back(cpu_time());
}

// Used by end_timer
double DF::pop(std::vector<double>& v)
{
   double x = v.back();
   v.pop_back();
   return x;
}

// Stops a timer
Tensor<double> DF::end_timer(World& world)
{
   Tensor<double> times(2);
   times[0] = wall_time() - pop(ttt);
   times[1] = cpu_time() - pop(sss);
   return times;
}

//reports back with time spent in calculation so far (current time - the first time in ttt,sss)
Tensor<double> DF::get_times(World& world){
     Tensor<double> times(2);
     times[0] = wall_time() - ttt[0];
     times[1] = cpu_time() - sss[0];
     return times;
}

//Changes a real function to a complex function
template <typename Q, int NDIM>
struct function_real2complex_op
{
  typedef std::complex<Q> resultT;
  Tensor<resultT> operator()(const Key<NDIM>& key, const Tensor<Q>& t) const
  {
    Tensor<resultT> result(t.ndim(), t.dims());
    BINARY_OPTIMIZED_ITERATOR(const Q, t, resultT, result, *_p1 = resultT(*_p0,0.0););
    return result;
  }
  template <typename Archive>
  void serialize(Archive& ar) {}
};

Function<std::complex<double>,3> function_real2complex(const Function<double,3>& r)
{
  return unary_op_coeffs(r, function_real2complex_op<double,3>());
}

/// Factory function generating operator for convolution with grad(bsh) in 3D
/// Returns a 3-vector containing the convolution operator for the
/// x, y, and z components of grad(bsh)
//Mostly copied from a function written by W. Scott Thornton
static
inline
std::vector< std::shared_ptr< SeparatedConvolution<double,3> > >
GradBSHOperator_Joel(World& world,
                    double mu,
                    double lo,
                    double eps,
                    const BoundaryConditions<3>& bc=FunctionDefaults<3>::get_bc(),
                    int k=FunctionDefaults<3>::get_k())
{
    typedef SeparatedConvolution<double,3> real_convolution_3d;
    typedef std::shared_ptr<real_convolution_3d> real_convolution_3d_ptr;
    const double pi = constants::pi;
    const Tensor<double> width = FunctionDefaults<3>::get_cell_width();
    double hi = width.normf(); // Diagonal width of cell
    const bool isperiodicsum = (bc(0,0)==BC_PERIODIC);
    if (isperiodicsum) hi *= 100; // Extend range for periodic summation

    GFit<double,3> fit=GFit<double,3>::BSHFit(mu,lo,hi,eps,false);
    Tensor<double> coeff=fit.coeffs();
    Tensor<double> expnt=fit.exponents();

     // Stuff that Joel added
     // Go through coeff and expnt and lump together terms with higher exponents than we want to use
     // --------------------------------------------------------------------------------------------

     // First, select how large of an exponent we're going to keep
     double max_expnt = 10000.0; //Just taking a stab at this for now

     // Then need to truncate coeff and expnt into new Tensors
    int rank = coeff.dim(0);
    double Cdelta = 0.0;
    //double max_kept;
    int max_j=0;
    for(int j = 0; j < rank; j++){
         if(expnt[j] > max_expnt){
               Cdelta += coeff[j]*std::pow(constants::pi/expnt[j],1.5);
         }
         else{
	   //max_kept = expnt[j];
               max_j = j;
               break;
         } 
    }
     coeff = coeff(Slice(max_j,-1));
     expnt = expnt(Slice(max_j,-1));

     // Then calculate what the new coefficient needs to be out front
     coeff[0] = coeff[0] + Cdelta * std::pow(expnt[0]/constants::pi,1.5); 

     //reset rank because we use it below
     rank = coeff.dim(0);
     //----------------------------------------------------------------------------------------------

    if (bc(0,0) == BC_PERIODIC) {
        fit.truncate_periodic_expansion(coeff, expnt, width.max(), true);
    }

    std::vector<real_convolution_3d_ptr> gradG(3);

    for (int dir=0; dir<3; dir++) {
        std::vector< ConvolutionND<double,3> > ops(rank);
        for (int mu=0; mu<rank; mu++) {
            // We cache the normalized operator so the factor is the value we must multiply
            // by to recover the coeff we want.
            double c = std::pow(sqrt(expnt(mu)/pi),3); // Normalization coeff
            ops[mu].setfac(coeff(mu)/c/width[dir]);

            for (int d=0; d<3; d++) {
                if (d != dir)
                    ops[mu].setop(d,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[d]*width[d], 0, isperiodicsum));
            }
            ops[mu].setop(dir,GaussianConvolution1DCache<double>::get(k, expnt(mu)*width[dir]*width[dir], 1, isperiodicsum));
        }
        gradG[dir] = real_convolution_3d_ptr(new SeparatedConvolution<double,3>(world, ops));
    }

    return gradG;
}


//Stolen from SCF.cc to aid in orthonormalization. Used in orthogonalize_inplace.
Tensor<std::complex<double>> Q2(const Tensor<std::complex<double>>& s) {
    Tensor<std::complex<double>> Q = -0.5*s;
    for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
    return Q;
}

//generic f(r)=||r|| function for calculation of the radial expectation value
double myr(const coord_3d& r){
     return std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

//Creates the fermi nuclear potential from the charge distribution. Also calculates the nuclear repulsion energy
void DF::make_fermi_potential(World& world, real_convolution_3d& op, real_function_3d& potential, double& nuclear_repulsion_energy){
     if(world.rank()==0) print("\n***Making a Fermi Potential***");
     
     //Get list of atom coordinates
     auto molecule = Init_params.molecule;

     //variables for upcoming loop
     real_function_3d temp;
     double tempnorm;
     potential = madness::real_factory_3d(world);

     //Go through the atoms in the molecule and construct the total charge distribution due to all nuclei
     const double safety = 0.1;
     double vtol = madness::FunctionDefaults<3>::get_thresh() * safety;
     for(auto&& atom : molecule.get_atoms()){
          madness::FermiPotentialFunctor rho(atom);
          temp = real_factory_3d(world).functor(rho).thresh(vtol);
          tempnorm = temp.trace();
          temp.scale(-1. * atom.atomic_number / tempnorm);
          potential += temp;
     }

     //Potential is found by application of the coulomb operator to the charge distribution
     potential = apply(op,potential);

     //Calculate the nuclear repulsion energy
     //It doesn't change iteration to iteration, so we want to calculate it once and store the result
     //We calculate it inside this function because here we already have access to the nuclear charges and coordinates
     nuclear_repulsion_energy = 0.0;
     double rr;
     for(unsigned int m = 0; m < molecule.natom(); m++){
          auto& atom_m = molecule.get_atom(m);
          for(unsigned int n = m+1; n < molecule.natom(); n++){
               auto& atom_n = molecule.get_atom(n);
               coord_3d dist = atom_m.get_coords() - atom_n.get_coords();
               rr = std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
               nuclear_repulsion_energy += atom_m.atomic_number * atom_n.atomic_number / rr;
          }
     }

}

// Collective constructor
// Design of this was mostly taken from Bryan Sundahl's response code
DF::DF(World & world, const char* filename) : DF(world, (world.rank() == 0 ? std::make_shared<std::ifstream>(filename) : nullptr))
{}

// Constructor that actually does stuff
// Design of this was mostly taken from Bryan Sundahl's response code
DF::DF(World & world,std::shared_ptr<std::istream> input) {
     // Start a timer
     start_timer(world);

     // Try and open input file
     if(world.rank() == 0){
          if (input->fail()) MADNESS_EXCEPTION("Dirac Fock failed to open input stream", 0);
   
          // Welcome user 
          print("\n   Preparing to solve the Dirac-Hartree-Fock equations.\n"); 

          // Read input files
          DFparams.read(*input);

          // Print out what was read in
          DFparams.print_params();
     }

     // Broadcast to all other nodes
     world.gop.broadcast_serializable(DFparams, 0);

     // Read in archive, but first find out if we're reading an nwchem file or other archive
     if(DFparams.nwchem){
          Init_params.readnw(world, DFparams.archive, DFparams.Krestricted);
     }
     else{
          Init_params.read(world, DFparams.archive, DFparams.restart, DFparams.Krestricted);
     }

     //print initialization parameters and molecule geometry
     if(world.rank() == 0){
          Init_params.print_params();
          print_molecule(world);
     }   

     // Set some function defaults   
     FunctionDefaults<3>::set_thresh(DFparams.thresh); //Always use user-specified thresh

     //Truncate mode 0 is preferrable, but currently way too expensive for relativistic calculations, so for now use 1
     FunctionDefaults<3>::set_truncate_mode(0);   

     //If user requests different k, then project functions
     if(DFparams.k != Init_params.order){

          //set function default
          FunctionDefaults<3>::set_k(DFparams.k);

          //Loop over orbitals
          for(unsigned int i = 0; i < Init_params.num_occupied; i++){

               //Look over 4 components of orbitals
               for(unsigned int j = 0; j < 4; j++){

                    //project
                    Init_params.orbitals[i][j] = project(Init_params.orbitals[i][j], FunctionDefaults<3>::get_k(), DFparams.thresh, false);
               }

               //fence and truncate
               world.gop.fence();
               Init_params.orbitals[i].truncate();
          }
     }

     //Set local orbitals and energies to those from the archive
     energies = Init_params.energies;
     occupieds = Init_params.orbitals;
     total_energy = Init_params.Init_total_energy;
     //If nonrelativistic calculation was spinrestricted then we're doing a closed shell calculation 
     //This is a little incorrect, as we're equating two separate concepts, but it works for now.
     closed_shell = Init_params.closed_shell;

     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("Preparation complete: ", times[0]);

}

//returns a new Fcwf that is the result of applying the Dirac free-particle hamiltonian on psi
Fcwf apply_T(World& world, Fcwf& psi){
     double myc = 137.0359895; //speed of light in atomic units
     std::complex<double> myi(0,1);
     complex_derivative_3d Dx(world,0);
     complex_derivative_3d Dy(world,1);
     complex_derivative_3d Dz(world,2);
     Fcwf Tpsi(world);
     
     //reconstruct psi
     psi.reconstruct();

     //take derivatives
     Fcwf psix = apply(world,Dx,psi); 
     Fcwf psiy = apply(world,Dy,psi);
     Fcwf psiz = apply(world,Dz,psi); 

     //compress
     psix.compress();
     psiy.compress();
     psiz.compress();
     psi.compress();

     //combine to calculate application of T
     Tpsi[0] = psiz[2] + psix[3] - myi*psiy[3];
     Tpsi[1] = psix[2] + myi*psiy[2] - psiz[3];
     Tpsi[2] = (myc*myc)*(psiz[0] + psix[1] - myi*psiy[1] - (2.0*myi)*psi[2]);
     Tpsi[3] = (myc*myc)*(psix[0] + myi*psiy[0] - psiz[1] - (2.0*myi)*psi[3]);

     return Tpsi * (-myi);
}

//function to calculate the kinetic + rest energy expectation value using Dirac Hamiltonian c*\alpha*p+\Beta*m*c*c
double DF::rele(World& world, Fcwf& psi){
     Fcwf Tpsi = apply_T(world, psi);
     std::complex<double> energy  = inner(psi, Tpsi);
     return energy.real();
}

//Calculates K*psi for each psi in the orbitals vector and stores them in result
void DF::exchange(World& world, real_convolution_3d& op, std::vector<Fcwf>& Kpsis){

     //start timer
     start_timer(world);

     //zero out Kpsis
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          Kpsis[i] = Fcwf(world);
     }

     //reconstruct
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          occupieds[i].reconstruct();
     }

     //Calculate and accumulate exchange contributions
     unsigned int n = Init_params.num_occupied;
     double myc = 137.0359895; //speed of light in atomic units

     //Calculates exchange contributions from the orbitals that we have stored
     //Loop through orbitals phi_i, computing K(phi_i), and while we're at it, use symmetry to start calculating contributions to later orbitals
     for(unsigned int i = 0; i < n; i++){

          //temp will hold the contributions (results of the coulomb operator application) that we need to finish calculation of K(phi_i)
          //For the first iteration, we need all n "contributions," so temp has length n, but this will decrease by 1 each iteration
          std::vector<complex_function_3d> temp(n-i);
          for(unsigned int j = 0; j < n-i; j++){
               temp[j] = complex_factory_3d(world);    
          }
          compress(world, temp);

          //break up ith through nth orbitals into their components to facilitate use of vmra functions
          std::vector<complex_function_3d> temp0(n-i);
          std::vector<complex_function_3d> temp1(n-i);
          std::vector<complex_function_3d> temp2(n-i);
          std::vector<complex_function_3d> temp3(n-i);
          for(unsigned int j = i; j < n; j++){
               temp0[j-i] = occupieds[j][0];
               temp1[j-i] = occupieds[j][1];
               temp2[j-i] = occupieds[j][2];
               temp3[j-i] = occupieds[j][3];
          }

          //These gaxpy calls accomplish the (\phi_j^\dagger)(\phi_i) in the numerator 
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][0]*conj(world,temp0));
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][1]*conj(world,temp1));
          gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][2]*conj(world,temp2));
          gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][3]*conj(world,temp3));

          //truncate before apply phase
          truncate(world, temp);

          //apply coulomb operator
          temp = apply(world, op, temp);
          
          //truncate again
          truncate(world,temp);

          //Now multiply by phi_j's and accumulate to K(phi_i)
          Kpsis[i][0] += sum(world, mul(world, temp, temp0));
          Kpsis[i][1] += sum(world, mul(world, temp, temp1));
          Kpsis[i][2] += sum(world, mul(world, temp, temp2));
          Kpsis[i][3] += sum(world, mul(world, temp, temp3));
          

          //Everything in temp can be used to accumulate a small part of K(phi_k) for k in [i+1,n] so we avoid calculating the same quantity on future iterations
          
          //First, take the complex conjugate of the contributions we have
          temp = conj(world, temp);

          //multiply by phi_i
          temp0 = occupieds[i][0]*temp;
          temp1 = occupieds[i][1]*temp;
          temp2 = occupieds[i][2]*temp;
          temp3 = occupieds[i][3]*temp;

          //acummulate
          for(unsigned int j = i+1; j < n; j++){
               Kpsis[j][0] += temp0[j-i];
               Kpsis[j][1] += temp1[j-i];
               Kpsis[j][2] += temp2[j-i];
               Kpsis[j][3] += temp3[j-i];
          }
         
     }

     //If our calculation is Kramers-restricted we need exchange contributions from the time-reversed orbitals that we don't explicitly store.
     //This will look very similar to the above loop (over i), but with rearrangements to the temporary vectors
     //Later, these can probably be combined into a single loop
     
     if(DFparams.Krestricted){
          int num_contrib=n;
          if(!closed_shell){
               num_contrib = n-1;
          }

          //Loop through orbitals phi_i, computing K(phi_i), and while we're at it, use symmetry to start calculating contributions to later orbitals
          for(unsigned int i = 0; i < num_contrib; i++){

               //temp will hold the contributions (results of the coulomb operator application) that we need to finish calculation of K(phi_i)
               //For the first iteration, we need all n "contributions," so temp has length n, but this will decrease by 1 each iteration
               std::vector<complex_function_3d> temp(num_contrib-i);
               for(unsigned int j = 0; j < num_contrib-i; j++){
                    temp[j] = complex_factory_3d(world);    
               }
               compress(world, temp);

               //break up the time-reversals of the ith through nth orbitals into their components to facilitate use of vmra functions
               std::vector<complex_function_3d> temp0(num_contrib-i);
               std::vector<complex_function_3d> temp1(num_contrib-i);
               std::vector<complex_function_3d> temp2(num_contrib-i);
               std::vector<complex_function_3d> temp3(num_contrib-i);
               for(unsigned int j = i; j < num_contrib; j++){
                    //the next four lines accomplish the rearrangement needed to get the time reversal rather than the orbital itself, but skip the complex conjugation, which will come later
                    temp0[j-i] = -1.0*occupieds[j][1];
                    temp1[j-i] = occupieds[j][0];
                    temp2[j-i] = -1.0*occupieds[j][3];
                    temp3[j-i] = occupieds[j][2];
               }

               //These gaxpy calls accomplish the (\phi_j^T)(\phi_i) in the numerator. Conjugation of phi_j is left to later
               gaxpy(world, 1.0, temp, 1.0, occupieds[i][0]*temp0);
               gaxpy(world, 1.0, temp, 1.0, occupieds[i][1]*temp1);
               gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][2]*temp2);
               gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][3]*temp3);

               //truncate before apply phase
               truncate(world, temp);

               //apply coulomb operator
               temp = apply(world, op, temp);
               
               //truncate again
               truncate(world,temp);

               //Now multiply by phi_j's and accumulate to K(phi_i)
               //Here is where we put the complex conjugation we've left out. This allows for one conjugation here instead of two conjugations earlier.
               Kpsis[i][0] += sum(world, mul(world, temp, conj(world,temp0)));
               Kpsis[i][1] += sum(world, mul(world, temp, conj(world,temp1)));
               Kpsis[i][2] += sum(world, mul(world, temp, conj(world,temp2)));
               Kpsis[i][3] += sum(world, mul(world, temp, conj(world,temp3)));
               
               //Now for the next part (accumulating the n-i "symmetric" contributions), we already have the complex conjugate of temp, so we can go straight to multiplication by the time-reversal of phi_i
               temp0 = conj(occupieds[i][1])*temp;
               temp1 = -1.0*conj(occupieds[i][0])*temp;
               temp2 = conj(occupieds[i][3])*temp;
               temp3 = -1.0*conj(occupieds[i][2])*temp;

               //accumulate
               for(unsigned int j = i+1; j < num_contrib; j++){
                    Kpsis[j][0] += temp0[j-i];
                    Kpsis[j][1] += temp1[j-i];
                    Kpsis[j][2] += temp2[j-i];
                    Kpsis[j][3] += temp3[j-i];
               }
               
          }

          //For open shell calculations need to calculate final contributions for the singly-occupied orbital
          //This will look very similar to the two above loops
          //Could put this as an addition to the first iteration of the loop above for more concise, though less readable code
          if(!closed_shell){
               std::vector<complex_function_3d> temp(n-1);
               for(unsigned int j = 0; j < n-1; j++){
                    temp[j] = complex_factory_3d(world);    
               }
               compress(world, temp);

               std::vector<complex_function_3d> temp0(n-1);
               std::vector<complex_function_3d> temp1(n-1);
               std::vector<complex_function_3d> temp2(n-1);
               std::vector<complex_function_3d> temp3(n-1);
               for(unsigned int j = 0; j < n-1; j++){
                    temp0[j] = -1.0*occupieds[j][1];
                    temp1[j] = occupieds[j][0];
                    temp2[j] = -1.0*occupieds[j][3];
                    temp3[j] = occupieds[j][2];
               }

               gaxpy(world, 1.0, temp, 1.0, occupieds[n-1][0]*temp0);
               gaxpy(world, 1.0, temp, 1.0, occupieds[n-1][1]*temp1);
               gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[n-1][2]*temp2);
               gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[n-1][3]*temp3);

               truncate(world, temp);

               temp = apply(world, op, temp);
               
               truncate(world,temp);

               Kpsis[n-1][0] += sum(world, mul(world, temp, conj(world,temp0)));
               Kpsis[n-1][1] += sum(world, mul(world, temp, conj(world,temp1)));
               Kpsis[n-1][2] += sum(world, mul(world, temp, conj(world,temp2)));
               Kpsis[n-1][3] += sum(world, mul(world, temp, conj(world,temp3)));
               
          }
     }

     //Truncate
     for(unsigned int i=0; i < n; i++) Kpsis[i].truncate();

     //Report time
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
}

//Diagonalize occupieds in the Fock space of occupieds. 
//occupieds is transformed in place.
//This requires Kpsis to be precomputed.
//Kpsis are transformed in place as well.
void DF::diagonalize(World& world, real_function_3d& myV, real_convolution_3d& op, std::vector<Fcwf>& Kpsis){


     if(world.rank()==0) print("\n***Diagonalizing***");
     start_timer(world);

     //Create a few integers to help keep track of loops
     unsigned int n = Init_params.num_occupied; //number of orbitals we have stored
     unsigned int np = closed_shell ? n : n-1; //number of PAIRS we want in computations
     unsigned int m = DFparams.Krestricted ? n+np : n; //Size of matrix to use

     //Initialize tensors and vectors to store temporary wavefunctions needed for computation
     Tensor<std::complex<double>> fock(m, m);
     Tensor<std::complex<double>> overlap(m, m);
     Tensor<std::complex<double>> U(m,m);
     Tensor<double> evals(m);
     std::vector<Fcwf> temp_orbitals;
     std::vector<Fcwf> kramers_pairs;

     //Make a permutation matrix for use later (in the Kramers-restricted case)
     //Fock matrix in the Kramers-restricted case is built with orbitals in the following order: 
     //
     //stored (doubly occupied) orbitals, 
     //singly-occupied orbital (if open shell), 
     //then kramers pairs of doubly-occupied orbitals
     //
     //However, before diagonalization we want the fock matrix as if our orbitals were in order of increasing energy.
     //This permutation matrix will accomplish the column and row swapping necessary to reorder
     Tensor<double> P(m,m);
     if(DFparams.Krestricted){
          for(unsigned int j=0; j < np; j++){
               P(j,2*j) = 1;
               P(n+j,2*j+1) = 1;
          }
          if(!closed_shell) P(n-1,m-1) = 1;
     }

     //Also make the vector of Kramers Pairs for use later (if Kramers-restricted)
     if(DFparams.Krestricted){
         for(unsigned int j = 0; j < np; j++){
               kramers_pairs.push_back(occupieds[j].KramersPair());
          }
     }
     if(world.rank()==0) print("     Forming Matrices");
     start_timer(world);
     
     ////Form the Fock Matrix
     //
     //In the Kramers-unrestricted case, the Fock matrix is simply nxn (or mxm, as n=m)
     //
     //In the Kramers restricted case, construction of the Fock matrix is a bit more complicated:
     //
     //With the orbitals in the order described above, the Fock matrix can be thought of as a 2x2 block matrix.
     //
     //The dimensions are (using F to indicate blocks of the Fock matrix):
     //F_11: n x n
     //F_12: n x np
     //F_21: np x n
     //F_22: np x np
     //
     //If the system is closed shell, then n=np.
     //
     //To save on computation we explicity calculate F_11 and F_21, and then use these to calculate the remaining blocks
     //
     //We do this by first making a vector of F*psi for all psi in the orbitals we have stored

     //calculate potential due to nuclei and mean field
     if(world.rank() == 0) print("          Adding (V+J)psi");
     real_function_3d rho = real_factory_3d(world);
     double fac = (DFparams.Krestricted ? 2.0 : 1.0);
     for(unsigned int j = 0; j < np; j++){
          rho += fac*squaremod(occupieds[j]);
     }
     if(!closed_shell) rho += squaremod(occupieds[n-1]);
     real_function_3d potential = myV + apply(op,rho);
     potential.truncate();

     //add in coulomb parts to neworbitals
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals.push_back(occupieds[j]*potential); //add in coulomb term
     }

     if(world.rank() == 0) print("          Subtracting K*psi");
     //Move Kpsis to new orbitals, as they are part of the fock operator
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] -= Kpsis[j]; //Must be subtracted because exchange function doesn't include the negative.
     }

     //add in T_psi
     if(world.rank()==0) print("          Adding T*psi");
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] += apply_T(world, occupieds[j]);  //add in "kinetic" term
     }

     //Now that we have F*psi (temp_orbitals), we can get on with integration
     if(world.rank()==0) print("          Integrating to form Fock Matrix");
     start_timer(world);

     //Now compute the fock matrix
     Tensor<std::complex<double>> tempmatrix = matrix_inner(world, occupieds, temp_orbitals);
     if(DFparams.Krestricted){
          fock(Slice(0,n-1),Slice(0,n-1)) = copy(tempmatrix);
          fock(Slice(n,m-1),Slice(n,m-1)) = conj(tempmatrix(Slice(0,np-1),Slice(0,np-1)));
          tempmatrix = matrix_inner(world,kramers_pairs,temp_orbitals);
          fock(Slice(n,m-1),Slice(0,n-1)) = copy(tempmatrix);
          fock(Slice(0,n-1),Slice(n,m-1)) = conj_transpose(tempmatrix);
     }
     else{
          fock = tempmatrix;
     }

     //permute and symmetrize
     if(DFparams.Krestricted) fock = inner(transpose(P),inner(fock,P));
     fock = (1.0/2.0)*(fock + conj_transpose(fock));
     
     //End timer for Fock matrix calculation
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("               ", times[0]); 


     ////and the overlap matrix
     if(world.rank()==0) print("          Integrating to form Overlap Matrix");
     start_timer(world);
     if(DFparams.Krestricted){
          overlap(Slice(0,n-1),Slice(0,n-1)) = matrix_inner(world,occupieds,occupieds);
          overlap(Slice(0,n-1),Slice(n,m-1)) = matrix_inner(world,occupieds,kramers_pairs);
          overlap(Slice(n,m-1),Slice(0,n-1)) = matrix_inner(world,kramers_pairs,occupieds);
          overlap(Slice(n,m-1),Slice(n,m-1)) = matrix_inner(world,kramers_pairs,kramers_pairs);
     }
     else{
          overlap = matrix_inner(world,occupieds,occupieds);
     }

     //permute and symmetrize
     if(DFparams.Krestricted) overlap = inner(transpose(P),inner(overlap,P));
     overlap = (1.0/2.0)*(overlap + conj_transpose(overlap));

     //End timers for overlap calculation and total matrix formation
     times = end_timer(world);
     if(world.rank()==0) print("               ", times[0]);
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     //Now call eigensolver
     if(world.rank()==0) print("     Eigensolver");
     start_timer(world);
     sygv(fock, overlap, 1, U, evals);
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     //Before applying the transformation, fix arbitrary rotations introduced by the eigensolver. 
     //This is done in three steps:
     // 1) Column-swapping in order to obtain a diagonally dominant matrix
     // 2) Removal of complex phases, so the diagonal is comprised of positive real numbers
     // 3) Identification of clusters with the same eigenvalue (within some tolerance) and removal of arbitrary
     //    rotations within the eigenspaces formed by those clusters
     //
     //This is largely stolen from SCF.cc, but heavily changed for FCWFs
     if(world.rank()==0) print("     Removing Rotations");
     start_timer(world);

     double thresh_degenerate = DFparams.thresh*10.0; //threshold for determining if eigenvalues are equal

     //swap columns for a diagonally dominant matrix
     bool switched = true;
     while (switched) {
          switched = false;
          for (unsigned int kk = 0; kk < m; kk++) {
               for (unsigned int j = kk + 1; j < m; j++) {
                    double sold = std::real(U(kk,kk)*std::conj(U(kk,kk))) + std::real(U(j,j)*std::conj(U(j,j)));
                    double snew = std::real(U(kk,j)*std::conj(U(kk,j))) + std::real(U(j,kk)*std::conj(U(j,kk)));
                    if (snew > sold and not ((evals[j] - evals[kk]) > thresh_degenerate * std::fabs(evals[kk])) ) {
                         if(world.rank()==0){
                              print("          swapping columns ", kk+1, " and ", j+1);
                         }
                         Tensor<std::complex<double>> tmp = copy(U(_, kk));
                         U(_, kk) = U(_, j);
                         U(_, j) = tmp;
                         std::swap(evals[kk], evals[j]);
                         switched = true;
                    }
               }
          }
     }

     // Fix phases.
     for (unsigned int kk = 0; kk < m; ++kk)
          U(_, kk).scale(std::conj(U(kk,kk))/std::abs(U(kk,kk)));

     //Find clusters of degenerate eigenvalues and rotate eigenvectors to maximize overlap with previous ones
     unsigned int ilo = 0; // first element of cluster
     if(world.rank()==0) print("          Degeneracy threshold: ",thresh_degenerate);

     while (ilo < m - 1) {
         unsigned int ihi = ilo;
         while (fabs(evals[ilo] - evals[ihi + 1]) < thresh_degenerate * std::fabs(evals[ilo])){
             ++ihi;
             if (ihi == m - 1) break;
         }
         unsigned int nclus = ihi - ilo + 1; //size of the cluster
         if (nclus > 1) {
              if(world.rank()==0){
                    //some printing to tell the user about the clusters found
                    print("          found cluster from ", ilo + 1, " to " , ihi + 1);
                    for(unsigned int kk = ilo; kk <= ihi; kk++){
                         print("               ",evals[kk]);
                    }
              }

              //Use the polar decomposition to undo rotations:
              //For a description of how this works see Matrix Computations by Golub and Van Loan, 4th Ed., Section 6.4.1
              Tensor<std::complex<double>> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));
              Tensor<std::complex<double>> VH(nclus,nclus);
              Tensor<std::complex<double>> W(nclus,nclus);
              Tensor<double> sigma(nclus);

              svd(q, W, sigma, VH);

              //W*VH is the rotation part of q. Undo it by taking the adjoint and right-multiplying
              q = conj_transpose(inner(W,VH));
              U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
              
         }
         ilo = ihi + 1;
     }
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     //Now undo permutations before the transformation phase
     if(DFparams.Krestricted){
          U = inner(P,inner(U,transpose(P)));
          evals = inner(evals,transpose(P));
     }

     if(world.rank()==0) print("     Applying Transformation");
     start_timer(world);

     //Apply the transformation to the orbitals by right-multiplying
     //(if Krestricted) Because we have n orbitals but the matrix has n+np rows, this happens in two stages
     if(DFparams.Krestricted){
          tempmatrix = U(Slice(0,n-1),Slice(0,n-1));
          occupieds = transform(world, occupieds, tempmatrix);
          tempmatrix = U(Slice(n,m-1),Slice(0,n-1));
          occupieds += transform(world, kramers_pairs, tempmatrix);
     }
     else{
          occupieds = transform(world,occupieds,U);
     }

     ////Apply the transformation to the Exchange-applied functions as well
     if(DFparams.Krestricted){
          for(unsigned int j = 0; j < np; j++){
               kramers_pairs[j] = Kpsis[j].KramersPair();
          }
          tempmatrix = U(Slice(0,n-1),Slice(0,n-1));
          Kpsis = transform(world, Kpsis, tempmatrix);
          tempmatrix = U(Slice(n,m-1),Slice(0,n-1));
          Kpsis += transform(world, kramers_pairs, tempmatrix);
     }
     else{
          Kpsis = transform(world,Kpsis,U);
     }

     //truncate
     for(int kk = 0; kk < n; kk++){
           Kpsis[kk].truncate();
           occupieds[kk].truncate();
     }

     //Set energies = evals
     energies = evals(Slice(0,n-1));

     //End timers for this function
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

}

//orthogonalize occupieds in place. This function mimics a similar one from SCF.cc
void DF::orthogonalize_inplace(World& world){

     unsigned int n = occupieds.size();
     double maxq;

     //normalize beforehand
     for(unsigned int i = 0; i < n; i++){
          occupieds[i].normalize();
     }

     //Basically stolen from SCF.cc. Orthogonalization based on Taylor expansion of (overlap)^(1/2)
     do{
          Tensor<std::complex<double>> Q = Q2(matrix_inner(world,occupieds,occupieds));
          maxq = 0.0;
          for(unsigned int j=1; j<n; j++){
               for(unsigned int i=0; i<j; i++){
                    maxq = std::max(std::abs(Q(j,i)),maxq);
               }
          }
          occupieds = transform(world, occupieds, Q);
     } while (maxq>0.01);

     //normalize afterward
     for(unsigned int i = 0; i < n; i++){
          occupieds[i].normalize();
     }

}

//Apply's Green's function to Vpsi (a Fcwf). Overwrites Vpsi with new Fcwf
//Use of this function has largely been replaced by apply_BSH_new, but
//this one is kept in case one wants to avoid use of the derivative operator.
void apply_BSH(World& world, Fcwf& Vpsi, double& eps, double& small, double& thresh){

     //necessary constants
     double myc = 137.0359895; //speed of light
     double c2 = myc*myc; //speed of light squared
     std::complex<double> myi(0,1); //imaginary number i
     std::complex<double> ic = myi*myc; //i*c
    
     //calculate exponent for equivalent BSH operator
     double mu = std::sqrt(-(2*eps*c2+eps*eps)/c2);

     world.gop.fence();
     
     //pointer to BSH operator
     std::shared_ptr<real_convolution_3d> op = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world,mu,small,thresh)); 

     //vector of pointers to the gradient of the BSH operator
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator(world, mu, small, thresh); 

     //Subsitute the below line for the above one if you want to screen operator coefficients
     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator_Joel(world, mu, 1e-8, thresh); 

     //operators are copied and organized into a vector for use of vmra's more efficient functions
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> allops(16);
     for(unsigned int i = 0; i < 4; i++){
          allops[i] = op;
          allops[4+i] = op3[0];
          allops[8+i] = op3[1];
          allops[12+i] = op3[2];
     }

     //create intermediate functions necessary to compute new components
     std::vector<complex_function_3d> temp(16);
     for(unsigned int i = 0; i < 4; i++){
          temp[i] = Vpsi[i];
          temp[4+i] = Vpsi[i];
          temp[8+i] = Vpsi[i];
          temp[12+i] = Vpsi[i];
     }

     //vector apply accomplishes all 16 necessary operator applications
     temp = apply(world, allops, temp);

     //All four components of the result fcwf can be written as a linear combination of the functions in temp now
     //Easy way to do this is to write out a transformation tensor and then call transform
     Tensor<std::complex<double>> U(16,4);
     U(0,0) = 2*c2+eps; U(14,0) = -ic; U(7,0) = -ic; U(11,0) = -myc;
     U(1,1) = 2*c2+eps; U(6,1) = -ic; U(10,1) = myc; U(15,1) = ic;
     U(12,2) = -ic; U(5,2) = -ic; U(9,2) = -myc; U(2,2) = eps;
     U(4,3) = -ic; U(8,3) = myc; U(13,3) = ic; U(3,3) = eps;
     U *= (1.0/c2);

     //Apply transformation tensor to temp (a vector of functions) and store back in Vpsi (a FCWF)
     temp = transform(world, temp, U);
     Vpsi[0] = temp[0];
     Vpsi[1] = temp[1];
     Vpsi[2] = temp[2];
     Vpsi[3] = temp[3];

}

//This function applies the Dirac Green's function to Vpsi (a Fcwf), overwriting Vpsi with a new Fcwf.
//
//Instead of using the derivative of the Green's function, this version first applies the
//nonrelativistic Green's function, then applies (H_D + eps) to the result.
//
//Analytically this is equivalent to apply_BSH, but faster because the application of
//the derivative operator is faster than application of an integral operator.
//
//Empirically this has resulted in no decrease in accuracy, despite reliance on the "noisier" derivative operator
void apply_BSH_new(World& world, Fcwf& Vpsi, double& eps, double& small, double& thresh){

     //necessary constants
     double myc = 137.0359895; //speed of light
     double c2 = myc*myc; //speed of light squared
     std::complex<double> myi(0,1); //imaginary number i
     //std::complex<double> ic = myi*myc; //i*c
    
     //calculate exponent for equivalent BSH operator
     double mu = std::sqrt(-(2*eps*c2+eps*eps)/c2);

     world.gop.fence();

     //create BSH operator
     real_convolution_3d op = BSHOperator3D(world, mu,small,thresh); 

     //Apply BSH operator to Vpsi
     Vpsi = apply(world, op, Vpsi);

     //Apply (1/c^2)(H_D + eps) to Vpsi. Using apply_T for convenience, but this requires adding 2c^2Vpsi
     Vpsi = apply_T(world, Vpsi)*(1.0/c2) + Vpsi * ((eps+2*c2)/c2);

}

// Small function to print geometry of a molecule nicely
// Straight up stolen from Bryan
void DF::print_molecule(World &world)
{
   if(world.rank() == 0)
   {
      // Precision is set to 10 coming in, drop it to 5
      std::cout.precision(5);
      std::cout << std::fixed;

      // First get atoms
      const std::vector<Atom> atoms = Init_params.molecule.get_atoms();
      int num_atoms = atoms.size();

      // Now print
      print("\n   Geometry Information");
      print("   --------------------\n");
      print("   Units: a.u.\n");
      print(" Atom            x                 y                 z");
      print("----------------------------------------------------------------");
      for(int j = 0; j < num_atoms; j++)
      {
           Vector<double,3> coords = atoms[j].get_coords();
           std::cout << std::setw(3) << atomic_number_to_symbol(atoms[j].get_atomic_number());
           std::cout << std::setw(18) << std::right << coords[0] << std::setw(18) << coords[1] << std::setw(18) << coords[2] << endl;
      }
      print("");

      // Reset precision
      std::cout.precision(10);
      std::cout << std::scientific;
   }
}

//Saves everything necessary to restart a DFdriver job.
void DF::saveDF(World& world){

     //get current time
     Tensor<double> times = get_times(world);

     //start timer
     start_timer(world);

     //print time of save
     if(world.rank()==0) print("\n***Saving at time: ",times[0]," ***");

     //Create archive and save the following:
     // 1) Total energy (double)
     // 2) Krestricted (boolean)
     // 3) closed_shell (boolean)
     // 3) number of occupied orbitals (int)
     // 4) orbital energies (vector of doubles)
     // 5) box size (double)
     // 6) wavelet order (int)
     // 7) molecule (molecule)
     // 8) occupied orbitals as complex functions
     try{
          //create archive
          archive::ParallelOutputArchive output(world, DFparams.savefile.c_str(), 1);

          //save simulation parameters and calculated properties
          output & total_energy & DFparams.Krestricted & closed_shell & Init_params.num_occupied & energies & Init_params.L & Init_params.order & Init_params.molecule;

          //Save orbitals
          //Loop over all occupied orbitals
          for(unsigned int i = 0; i < Init_params.num_occupied; i++){
               //Loop over four components
               for(int j = 0; j < 4; j++){
                    output & occupieds[i][j];
               }
          }
     }
     catch(const char* s){
          if(world.rank()==0) print("Failed to save DF restart data with error message:", s);
     }

     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
}

//Creates the (Gaussian) nuclear potential from the molecule object
void DF::make_gaussian_potential(World& world, real_function_3d& potential){
     if(world.rank()==0) print("\n***Making a Gaussian Potential***");
     GaussianPotentialFunctor Vfunctor(Init_params.molecule);
     potential = real_factory_3d(world).functor(Vfunctor).truncate_mode(0).truncate_on_project();
}

//Creates the (Gaussian) nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
void DF::make_gaussian_potential(World& world, real_function_3d& potential, double& nuclear_repulsion_energy){
     if(world.rank()==0) print("\n***Making a Gaussian Potential***");
     auto molecule = Init_params.molecule;

     GaussianPotentialFunctor Vfunctor(molecule);
     potential = real_factory_3d(world).functor(Vfunctor).truncate_mode(0).truncate_on_project();

     nuclear_repulsion_energy = 0.0;
     double rr;
     for(int m = 0; m < molecule.natom(); m++){
          auto& atom_m = molecule.get_atom(m);
          for(int n = m+1; n < molecule.natom(); n++){
               auto& atom_n = molecule.get_atom(n);
               coord_3d dist = atom_m.get_coords() - atom_n.get_coords();
               rr = std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
               nuclear_repulsion_energy += atom_m.atomic_number * atom_n.atomic_number / rr;
          }
     }

}

//Own version of load balancing for DF. Load balance on the functions as well as the nuclear potential
void DF::DF_load_balance(World& world, real_function_3d& Vnuc){
     if(world.rank()==0) print("\n***Load Balancing***");
     start_timer(world);
     
     //create Load balance object
     LoadBalanceDeux<3> lb(world);

     //Add functions that we want to load balance based on
     lb.add_tree(Vnuc, lbcost<double,3>(12.0,96.0),true);
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          for(int kk = 0; kk < 4; kk++){
               lb.add_tree(occupieds[j][kk], lbcost<std::complex<double>,3>(24.0,192.0),true);
          }
     }

     //Redistribute
     FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

     //End timers
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
     
}

//Function to output 2*n lineplots. For each orbital we make one "density" from the large component
//and one from the small component, evaluated on a grid with npt points from 0 to endpnt on the x axis
void DF::make_component_lineplots(World& world, const char* filename1, const char* filename2, int npt, double endpnt){

     //vectors to store densities
     std::vector<real_function_3d> large_densities;
     std::vector<real_function_3d> small_densities;

     //Push back densities
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities.push_back(squaremod_large(occupieds[i]));
          small_densities.push_back(squaremod_small(occupieds[i]));
     }

     //reconstruct
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities[i].reconstruct();
          small_densities[i].reconstruct();
     }

     //create lineplots
     double h = endpnt*(1.0/(npt-1));
     if(world.rank()==0){
          //open files
          FILE* file1 = fopen(filename1,"w");
          FILE* file2 = fopen(filename2,"w");
          if(!file1) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          if(!file2) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);

          //Loop through occupied orbitals
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               //Loop through desired grid
               for(int j=0; j<npt; ++j){
                    //calculate next grid point
                    coord_3d r({j*h,0.0,0.0});

                    //write grid point x value
                    fprintf(file1,"%.14e ", r[0]);
                    fprintf(file2,"%.14e ", r[0]);

                    //write function values (of the densities)
                    plot_line_print_value(file1, large_densities[i].eval(r));
                    plot_line_print_value(file2, small_densities[i].eval(r));

                    //write newlines
                    fprintf(file1,"\n");
                    fprintf(file2,"\n");
               }
               //newlines between orbitals
               fprintf(file1,"\n");
               fprintf(file2,"\n");
          }
          //close files
          fclose(file1);
          fclose(file2);
     }
     world.gop.fence();

}

//Same as above, but uses a log scale for the horizontal axis
void DF::make_component_logplots(World& world, const char* filename1, const char* filename2, int npt, int startpnt, int endpnt){
     
     //vectors to store densities
     std::vector<real_function_3d> large_densities;
     std::vector<real_function_3d> small_densities;

     //Push back densities
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities.push_back(squaremod_large(occupieds[i]));
          small_densities.push_back(squaremod_small(occupieds[i]));
     }

     //reconstruct
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities[i].reconstruct();
          small_densities[i].reconstruct();
     }

     //create lineplots
     double h = (endpnt-startpnt)*(1.0/(npt-1));
     if(world.rank()==0){
          //open files
          FILE* file1 = fopen(filename1,"w");
          FILE* file2 = fopen(filename2,"w");
          if(!file1) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          if(!file2) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);

          //Loop through occupied orbitals
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               //Loop through desired grid
               for(int j=0; j<npt; ++j){
                    //calculate next grid point
                    coord_3d r({std::pow(10,startpnt + j*h),0.0,0.0});

                    //write grid point x value
                    fprintf(file1,"%.14e ", r[0]);
                    fprintf(file2,"%.14e ", r[0]);

                    //write function values (of the densities)
                    plot_line_print_value(file1, large_densities[i].eval(r));
                    plot_line_print_value(file2, small_densities[i].eval(r));

                    //write newlines
                    fprintf(file1,"\n");
                    fprintf(file2,"\n");
               }
               //newlines between orbitals
               fprintf(file1,"\n");
               fprintf(file2,"\n");
          }
          //close files
          fclose(file1);
          fclose(file2);
     }
     world.gop.fence();

}

//Another lineplotting function, but write the entire density, instead of one divided into large and small components
void DF::make_density_lineplots(World& world, const char* filename, int npt, double endpnt){

     //vector to store densities
     std::vector<real_function_3d> densities;

     //Push back densities
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          densities.push_back(squaremod(occupieds[i]));
     }

     //reconstruct
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          densities[i].reconstruct();
     }

     //create lineplots
     double h = endpnt*(1.0/(npt-1));
     if(world.rank()==0){
          //open file
          FILE* file = fopen(filename,"w");
          if(!file) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          
          //Loop through occupied orbitals
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               //Loop through desired grid
               for(int j=0; j<npt; ++j){
                    //calculate next grid point
                    coord_3d r({j*h,0.0,0.0});

                    //write grid point x value
                    fprintf(file,"%.14e ", r[0]);

                    //write density function value
                    plot_line_print_value(file, densities[i].eval(r));

                    //newline
                    fprintf(file,"\n");
               }
               //newline between orbitals
               fprintf(file,"\n");
          }
          //close file
          fclose(file);
     }
     world.gop.fence();

}

//One complete iteration of the Dirac-Hartree-Fock solver
bool DF::iterate(World& world, real_function_3d& V, real_convolution_3d& op, real_function_3d& JandV, std::vector<Fcwf>& Kpsis, XNonlinearSolver<std::vector<Fcwf>, std::complex<double>, Fcwf_vector_allocator>& kainsolver, double& tolerance, int& iteration_number, double& nuclear_repulsion_energy){

     //Get and print the time of this iteration's start, and start a timer
     Tensor<double> times = get_times(world);
     if(world.rank()==0) print("\n\n\nIteration: ", iteration_number, " at ",times[0]);
     if(world.rank()==0) print("--------------");
     start_timer(world);

     //A vector to hold residuals after BSH application
     std::vector<Fcwf> Residuals;

     //Norm of a residual. We use these norms one-at-a-time, so no vector is needed
     double residualnorm;

     //A working FCWF that will have multiple uses
     Fcwf temp_function(world);

     //Boolean used in while loop.
     bool iterate_again = false; //Initialize to false = assume iterations will stop

     //First diagonalize the occupied orbitals in the Fock space (of occupied orbitals). Also transforms Kpsis.
     diagonalize(world, V, op, Kpsis);

     //Diagonalization forces us to recompute density
     real_function_3d rho = real_factory_3d(world);
     double fac = (DFparams.Krestricted ? 2 : 1);
     if(closed_shell){
          for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
     }
     else{
          for(unsigned int kk = 0; kk < Init_params.num_occupied-1; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
          rho += squaremod(occupieds[Init_params.num_occupied-1]);
     }

     //Apply coulomb operator to density and combine with V to get total potential
     JandV = V + apply(op,rho); 
     JandV.truncate();


     //Apply Dirac BSH to each psi
     if(world.rank()==0) print("\n***Applying BSH operator***");
     start_timer(world);
     double maxresidual = -1.0;
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){

          //construct the function to which we will apply the BSH
          temp_function = occupieds[j]*JandV;
          temp_function.scale(-1.0);
          temp_function += Kpsis[j];
          temp_function.truncate();

          //temp_function now holds (K-V-J)psi, so apply the BSH
          apply_BSH_new(world,  temp_function, energies[j], DFparams.small, DFparams.thresh);

          //truncate
          temp_function.truncate();

          //Now calculate the residual
          residualnorm = (occupieds[j] - temp_function).norm2();

          //Print residual norm for user to keep track
          if(world.rank()==0) printf("     Orbital: %3i,  Resid: %.10e\n",j+1, residualnorm);

          //Keep track of the maximum residual
          maxresidual = std::max(maxresidual, residualnorm);

          //Store residual function if we're using KAIN. Not necessary if we're not using kain
          //We don't use KAIN on the first iteration.
          if(iteration_number != 1 and DFparams.kain){
               Residuals.push_back(occupieds[j] - temp_function);
          }
          else{
               //if not using KAIN, then the result is the new orbital, with some step restriction
               residualnorm = (occupieds[j] - temp_function).norm2(); 
               if(residualnorm > DFparams.maxrotn){
                    double s = DFparams.maxrotn / residualnorm;
                    if(world.rank()==0) print("     restricting step for orbital: ", j+1);
                    occupieds[j] = temp_function*s + occupieds[j]*(1.0-s);
               }
               else{
                    occupieds[j] = temp_function; 
               }
          }
     }

     //Print max residual and the tolerance that we're using
     if(world.rank()==0) printf("                max Resid: %.10e\n",maxresidual);
     if(world.rank()==0) printf("                tolerance: %.10e\n",tolerance);

     //End timers
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //If any residual is still larger than the tolerance then we need to iterate again.
     //Can just enforce this on the max residual
     if(maxresidual > tolerance) iterate_again = true;

     //Apply the kain solver, if called for
     if(iteration_number != 1 and DFparams.kain){
          if(world.rank()==0) print("\n***Applying KAIN Solver***");
          start_timer(world);

          //Implement KAIN, but enforce a step restriction to KAIN doesn't take too large of a step
          Residuals = kainsolver.update(occupieds, Residuals); //Using Residuals for updated orbitals to save storage
          for(unsigned int i=0; i < Init_params.num_occupied; i++){

               //see how big of a step KAIN took for each orbital
               residualnorm = (occupieds[i]-Residuals[i]).norm2();

               //Restrict the step taken by KAIN if it's too big
               //This code is basically stolen from SCF.cc
               if(residualnorm > DFparams.maxrotn){
                    double s = DFparams.maxrotn / residualnorm;
                    if(world.rank()==0) print("     restricting step for orbital: ", i+1);
                    occupieds[i] = Residuals[i]*s + occupieds[i]*(1.0-s);
               }
               else{
                    occupieds[i] = Residuals[i];
               }
          }

          //End timers
          times = end_timer(world);
          if(world.rank()==0) print("     ", times[0]);
     }

     //truncate after BSH(+KAIN) application
     for(unsigned int i = 0; i < Init_params.num_occupied; i++) occupieds[i].truncate();

     //orthonormalize
     if(world.rank()==0) print("\n***Orthonormalizing***");
     start_timer(world);

     orthogonalize_inplace(world);

     //truncate here and normalize again
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
           occupieds[i].truncate();
           occupieds[i].normalize();
     }

     //End orthonormalization timer
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Now calculate new exchange functions. Has timer built in.
     if(world.rank()==0) print("\n***Recalculating Exchange***");
     exchange(world, op, Kpsis);

     //Calculate new J+V term
     if(world.rank()==0) print("\n***Recalculating Coulomb***");
     start_timer(world);
     
     //Make density
     rho = real_factory_3d(world);
     if(closed_shell){
          for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
     }
     else{
          for(unsigned int kk = 0; kk < Init_params.num_occupied-1; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
          rho += squaremod(occupieds[Init_params.num_occupied-1]);
     }
     
     //Apply coulomb operator to density and combine with V to get total potential
     JandV = V + apply(op,rho); 
     JandV.truncate();
     
     //End timer for coulombic potential calculation
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Calculate and print total energy each iteration, as well as a breakdown of different contributions
     if(world.rank()==0){
          print("\n***Printing Current Energies***");
     }

     start_timer(world);

     double kinetic_energy = 0.0;
     double coulomb_energy = 0.0;
     double exchange_energy = 0.0;
     double nuclear_attraction_energy = 0.0;
     double old_total_energy = total_energy;
     double myc = 137.0359895; //speed of light
     Tensor<double> nuclear_attraction_tensor;
     Tensor<double> coulomb_tensor;
     Tensor<double> exchange_tensor;

     //Potential due to the density only (no nuclear potential)
     real_function_3d Jop = apply(op,rho);

     //Compute kinetic energy contributions
     if(closed_shell){
          for(unsigned int j = 0; j < Init_params.num_occupied; j++){
               kinetic_energy += fac*rele(world, occupieds[j]);
          }
     }
     else{
          for(unsigned int j = 0; j < Init_params.num_occupied-1; j++){
               kinetic_energy += fac*rele(world, occupieds[j]);
          }
          kinetic_energy += rele(world, occupieds[Init_params.num_occupied-1]);
     }
          
     //Compute electron-nuclear attraction energy contributions, taking advantage of vmra's inner
     std::vector<complex_function_3d> occupieds1(Init_params.num_occupied);
     std::vector<complex_function_3d> occupieds2(Init_params.num_occupied);
     std::vector<complex_function_3d> occupieds3(Init_params.num_occupied);
     std::vector<complex_function_3d> occupieds4(Init_params.num_occupied);
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          occupieds1[i] = occupieds[i][0];
          occupieds2[i] = occupieds[i][1];
          occupieds3[i] = occupieds[i][2];
          occupieds4[i] = occupieds[i][3];
     }
     nuclear_attraction_tensor = real(inner(world,occupieds1,mul(world,V,occupieds1)));
     nuclear_attraction_tensor += real(inner(world,occupieds2,mul(world,V,occupieds2)));
     nuclear_attraction_tensor += real(inner(world,occupieds3,mul(world,V,occupieds3)))*(1.0/(myc*myc));
     nuclear_attraction_tensor += real(inner(world,occupieds4,mul(world,V,occupieds4)))*(1.0/(myc*myc));
     nuclear_attraction_energy = fac*nuclear_attraction_tensor.sum();
     if(!closed_shell and DFparams.Krestricted) nuclear_attraction_energy -= nuclear_attraction_tensor(Init_params.num_occupied-1);
     
     //Compute electron-electron repulsion energy contribution, again using vmra
     //Regarding use of fac here:
     //   Normally, this sum runs over all *unique* pairs of orbitals
     //   Our sum runs over all pairs of orbitals, so we divide by 2.
     //   In the Kramers-restricted case, we multiply by 2 and this cancels.
     //   In the Kramers-restricted open-shell case, the doubling doesn't apply to the last orbital, so we subtract the excess
     coulomb_tensor = real(inner(world,occupieds1,mul(world,Jop,occupieds1)));
     coulomb_tensor += real(inner(world,occupieds2,mul(world,Jop,occupieds2)));
     coulomb_tensor += real(inner(world,occupieds3,mul(world,Jop,occupieds3)))*(1.0/(myc*myc));
     coulomb_tensor += real(inner(world,occupieds4,mul(world,Jop,occupieds4)))*(1.0/(myc*myc));
     coulomb_energy = coulomb_tensor.sum()*(fac/2.0);
     if(!closed_shell and DFparams.Krestricted) coulomb_energy -= 0.5*coulomb_tensor(Init_params.num_occupied-1);
     
     //Calculate Exchange energy contribution, with some similar logic to the above
     std::vector<complex_function_3d> Kpsis1(Init_params.num_occupied);
     std::vector<complex_function_3d> Kpsis2(Init_params.num_occupied);
     std::vector<complex_function_3d> Kpsis3(Init_params.num_occupied);
     std::vector<complex_function_3d> Kpsis4(Init_params.num_occupied);
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          Kpsis1[i] = Kpsis[i][0];
          Kpsis2[i] = Kpsis[i][1];
          Kpsis3[i] = Kpsis[i][2];
          Kpsis4[i] = Kpsis[i][3];
     }
     exchange_tensor = real(inner(world,occupieds1,Kpsis1));
     exchange_tensor += real(inner(world,occupieds2,Kpsis2));
     exchange_tensor += real(inner(world,occupieds3,Kpsis3))*(1.0/(myc*myc));
     exchange_tensor += real(inner(world,occupieds4,Kpsis4))*(1.0/(myc*myc));
     exchange_energy = exchange_tensor.sum()*(fac/2.0);
     if(!closed_shell and DFparams.Krestricted) exchange_energy -= 0.5*exchange_tensor(Init_params.num_occupied-1);

     //compute total energy using the above computed contributions
     total_energy = kinetic_energy + coulomb_energy - exchange_energy + nuclear_attraction_energy + nuclear_repulsion_energy;

     //End timers for total energy calculation
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Now we calculate other quantities for later printing:
     // 1) expectation values for r
     // 2) number of coefficients used to represent each orbital
     // 3) maximum depth used in representing each orbital
     // 4) norms of each component of each orbital

     //Need function f(r)=r for the expectation values
     real_function_3d rfunc = real_factory_3d(world).f(myr);

     //Now loop through occupied orbitals, computing desired quantities (above)
     std::vector<double> r_expec_vec(Init_params.num_occupied);
     std::vector<int> numcoeffs_vec(Init_params.num_occupied);
     std::vector<int> maxdepth_vec(Init_params.num_occupied);
     std::vector<double> comp1norm(Init_params.num_occupied);
     std::vector<double> comp2norm(Init_params.num_occupied);
     std::vector<double> comp3norm(Init_params.num_occupied);
     std::vector<double> comp4norm(Init_params.num_occupied);
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){

          //calculate <r>
          double r_expec = std::real(inner(occupieds[j], occupieds[j]*rfunc));

          //find number of coefficients (sum over all components)
          int numcoeffs = occupieds[j][0].size() + occupieds[j][1].size() + occupieds[j][2].size() + occupieds[j][3].size();

          //find maximum depth (max over all components)
          int maxdepth = std::max(occupieds[j][0].max_depth(), occupieds[j][1].max_depth());
          maxdepth = std::max(int(occupieds[j][2].max_depth()), maxdepth);
          maxdepth = std::max(int(occupieds[j][3].max_depth()), maxdepth);

          //Store in vectors
          r_expec_vec[j] = r_expec;
          numcoeffs_vec[j] = numcoeffs;
          maxdepth_vec[j] = maxdepth;

          //Compute and store norms of each component. Remember to scale small component by c
          comp1norm[j] = occupieds[j][0].norm2();
          comp2norm[j] = occupieds[j][1].norm2();
          comp3norm[j] = occupieds[j][2].norm2()/myc;
          comp4norm[j] = occupieds[j][3].norm2()/myc;

     }

     //Print everything
     if(world.rank()==0){
          print("Orbital                 Energy            <r>     No. coeffs     Max depth        ||1||        ||2||        ||3||        ||4||");
          print("------------------------------------------------------------------------------------------------------------------------------");
          for(unsigned int j = 0; j < Init_params.num_occupied; j++){
               printf("%7i %22.10e %14.5e %14i %13i %12.5e %12.5e %12.5e %12.5e\n",j+1,energies[j],r_expec_vec[j],numcoeffs_vec[j],maxdepth_vec[j],comp1norm[j],comp2norm[j],comp3norm[j],comp4norm[j]);
          }
     }

     //Print aggregate information
     if(world.rank()==0){
          print("\n              Kinetic Energy: ",kinetic_energy);
          print("             Coulomb  Energy: ",coulomb_energy);
          print("             Exchange Energy: ",exchange_energy);
          print("   Nuclear Attraction Energy: ",nuclear_attraction_energy);
          print("    Nuclear Repulsion Energy: ",nuclear_repulsion_energy);
          print("                Total Energy: ",total_energy);
          print("       Total Energy Residual: ", std::fabs(total_energy - old_total_energy));
     }

     //final truncatation of orbitals now that we've computed properties
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          occupieds[j].truncate();
     }

     times = end_timer(world);
     if(world.rank()==0) print("     Iteration time:", times[0]);
 
     return iterate_again;
}

// Solves for the ground state Dirac Hartree Fock orbitals
void DF::solve_occupied(World & world)
{

     //State what we're doing here
     if(world.rank()==0){
          if(DFparams.Krestricted){
               if(closed_shell) print("\nSolving for ", Init_params.num_occupied, " doubly-occupied orbitals\n------------------------------\n");
               else print("\nSolving for ", Init_params.num_occupied-1, " doubly-occupied, 1 singly-occupied orbitals\n------------------------------\n");
          }
          else{
               print("\nSolving for ", Init_params.num_occupied, " single-occupied orbitals\n------------------------------\n");
          }
     }

     //Will need a coulomb operator
     real_convolution_3d op = CoulombOperator(world,DFparams.small,DFparams.thresh);

     //allocator is useful to have, but also required for use of KAIN
     Fcwf_vector_allocator allocator(world,Init_params.num_occupied);

     //initialize kain solver
     XNonlinearSolver<std::vector<Fcwf>, std::complex<double>, Fcwf_vector_allocator> kainsolver(allocator);
     kainsolver.set_maxsub(DFparams.maxsub);

     //normalize initial guesses
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          occupieds[i].normalize();
     }

     //Form nuclear potential
     real_function_3d Vnuc;
     double nuclear_repulsion_energy;
     if(DFparams.nucleus == 1){
          make_fermi_potential(world, op, Vnuc, nuclear_repulsion_energy);
     }
     else{
          make_gaussian_potential(world, Vnuc, nuclear_repulsion_energy);
     }

     //Initial load balance
     DF_load_balance(world, Vnuc);

     //Initialize vector of Fcwfs to hold exchange applied to the orbitals
     std::vector<Fcwf> Kpsis = allocator();

     //Calculate initial exchange
     if(world.rank()==0) print("\n***Calculating Initial Exchange***");
     exchange(world, op, Kpsis);

     //Calculate initial J+V term
     if(world.rank()==0) print("\n***Calculating Initial Coulomb***");
     start_timer(world);
     real_function_3d rho = real_factory_3d(world);
     double fac = (DFparams.Krestricted ? 2 : 1);
     if(closed_shell){
          for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
     }
     else{
          for(unsigned int kk = 0; kk < Init_params.num_occupied-1; kk++){
               rho += fac*squaremod(occupieds[kk]);
          }
          rho += squaremod(occupieds[Init_params.num_occupied-1]);
     }
     real_function_3d JandV = Vnuc + apply(op,rho); 
     JandV.truncate();
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Set tolerance for residuals
     double tol = 50.0*DFparams.thresh; 

     //Now time to start iterating
     bool keep_going = true;
     int iteration_number = 1;
     while((keep_going and iteration_number <= DFparams.max_iter) or iteration_number <= DFparams.min_iter){
          keep_going = iterate(world, Vnuc, op, JandV, Kpsis, kainsolver, tol, iteration_number, nuclear_repulsion_energy);
          
          //Load balance and save between iterations
          if(keep_going and iteration_number <= DFparams.lb_iter) DF_load_balance(world, Vnuc);
          if(DFparams.do_save) saveDF(world);

          //Increment iteration counter
          iteration_number++;
     }

     ////Calculation of Effective Electric Field:
     //if(world.rank()==0) print("Effective Electric Field calculation");
     //std::complex<double> myi(0,1);
     //std::complex<double> one(1,0);
     //real_derivative_3d Dx(world,0);
     //real_derivative_3d Dy(world,1);
     //real_derivative_3d Dz(world,2);
     //double Eeff(0.0);
     //for(unsigned int j; j < Init_params.num_occupied; j++){
     //     real_function_3d LL(world);
     //     for(unsigned int kk; kk < Init_params.num_occupied; kk++){
     //          if(kk != j){
     //               LL += squaremod(occupieds[kk]);
     //          }
     //     }
     //     LL = apply(op,LL);
     //     LL += Vnuc;
     //     complex_function_3d Lx = one*Dx(LL);
     //     complex_function_3d Ly = one*Dy(LL);
     //     complex_function_3d Lz = one*Dz(LL);
     //     Fcwf temp(world);

     //     temp[0] = Lz*occupieds[j][0] + (Lx - myi*Ly)*occupieds[j][1];
     //     temp[1] =  (Lx + myi*Ly)*occupieds[j][0] - Lz*occupieds[j][1];
     //     temp[2] = Lz*occupieds[j][2] + (Lx - myi*Ly)*occupieds[j][3];
     //     temp[2].scale(-1.0);
     //     temp[3] = Lz*occupieds[j][3] - (Lx + myi*Ly)*occupieds[j][2];

     //     Eeff += std::real(inner(occupieds[0],temp));
     //}
     //if(world.rank()==0) print("Eeff = ", Eeff);


}

void DF::solve(World& world){

     // Start timer
     start_timer(world);

     //Begin calculation
     if(world.rank() == 0){
        print("\n\n   Dirac Fock Calculation");
        print("   ------------------------");
     }
    
     if(not DFparams.no_compute){
          if(DFparams.job == 0){
               solve_occupied(world);
          }
          else{
               if(world.rank()==0) print("Specify a better job parameter.");
          }
     }
     else{
          if(world.rank()==0) print("Requested no computation.");
     }
    
     // Report calculation time
     // Precision is set to 10 coming in, drop it to 2
     std::cout.precision(2);
     std::cout << std::fixed;
     Tensor<double> times = end_timer(world);
     if(world.rank() == 0) print("\n   Calculation time:", times[0],"\n");

     //Make density lineplots
     if(DFparams.lineplot){
          start_timer(world);
          if(world.rank()==0) print("***Making lineplots***");
          //make_density_lineplots(world, "density_lineplots", 100000, 0.005);
          //make_component_lineplots(world, "large_component_lineplots", "small_component_lineplots", 100000, 5);
          make_component_logplots(world, "large_component_lineplots", "small_component_lineplots", 1000, -6, 1);
          times = end_timer(world);
          if(world.rank()==0) print("     ", times[0]);
     }
     
}

//Print the number of coefficients being used by each orbital.
//This has been useful in the past for debugging
void DF::print_sizes(World& world, bool individual=false){
     if(world.rank()==0) print("\nPrinting orbital sizes:\n");
     int n = Init_params.num_occupied;
     double a,b1,b2,b3,b4;
     for(unsigned int j=0; j < n; j++){
          b1 = occupieds[j][0].size();
          b2 = occupieds[j][1].size();
          b3 = occupieds[j][2].size();
          b4 = occupieds[j][3].size();
          a = b1+b2+b3+b4;
          if(world.rank()==0){
               if(individual){
                    print("Orbital; ",j+1," size: ",b1,"+",b2,"+",b3,"+",b4,"=",a);
               }
               else{
                    print("Orbital; ",j+1," size: ",a);
               }
          }
     }
}

//kthxbye



