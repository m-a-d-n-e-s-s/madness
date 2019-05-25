/*
 *
 *   Written by: bsundahl and jscanderson
 *   Date: A long time ago...
 *
 * TODO:
 *  - Use energies calculated in diagonalization and calculate various totals at the end
 *  - Use total of orbital energies as convergence criterion, rather than recalculating
 *  - 
 *
 */ 

#include "DF.h"
//#include "Plot_VTK.h"
#include "fcwf.h"
#include "../chem/potentialmanager.h"

using namespace madness;

// Needed for rebalancing
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

// Needed for timers
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

//reports back with time spent in calculation (current time - the first time in ttt,sss)
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
     //if(world.rank()==0) print("Initial expnts are: \n",expnt, "\n");
     //if(world.rank()==0) print("Initial coeffs are: \n",coeff, "\n");
    int rank = coeff.dim(0);
    double Cdelta = 0.0;
    double max_kept;
    int max_j;
    for(int j = 0; j < rank; j++){
         if(expnt[j] > max_expnt){
               Cdelta += coeff[j]*std::pow(constants::pi/expnt[j],1.5);
         }
         else{
               max_kept = expnt[j];
               max_j = j;
               break;
         } 
    }
    //if(world.rank()==0) print("limit exponent is: ", max_expnt, " and max kept is: ", max_kept);
     coeff = coeff(Slice(max_j,-1));
     expnt = expnt(Slice(max_j,-1));
     //if(world.rank()==0) print("new expnts are: \n", expnt, "\n");
     //if(world.rank()==0) print("new coeffs are: \n", coeff, "\n");

     // Then calculate what the new coefficient needs to be out front
     coeff[0] = coeff[0] + Cdelta * std::pow(expnt[0]/constants::pi,1.5); 

     //reset rank because we use it below
     rank = coeff.dim(0);
     //----------------------------------------------------------------------------------------------

     
     
     if(world.rank()==0) print("rank: ", rank, "\nexponents:\n", expnt, "\ncoefficients:\n", coeff);


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


//Stolen from SCF.cc to aid in orthonormalization
Tensor<std::complex<double>> Q2(const Tensor<std::complex<double>>& s) {
    Tensor<std::complex<double>> Q = -0.5*s;
    for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.5;
    return Q;
}



//Functor to make the (Gaussian) nuclear potential
class GaussianNucleusFunctor : public FunctionFunctorInterface<double,3> {
     private:
          std::vector<int> m_Zlist;
          std::vector<coord_3d> m_Rlist;
          std::vector<int> m_Alist;
          std::vector<double> m_xi;
     public:
          // Constructor 
          GaussianNucleusFunctor(Molecule& molecule, double bohr_rad){

               //get atom coordinates
               m_Rlist = molecule.get_all_coords_vec();
               
               //get atomic numbers
               for(unsigned int i = 0; i < m_Rlist.size(); i++){
                    m_Zlist.push_back(molecule.get_atom_number(i));
               }
               
               //find atomic mass numbers for each atom
               int Alist[116] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64,69,74,75,80,79,84,85,88,89,90,93,98,98,102,103,106,107,114,115,120,121,130,127,132,133,138,139,140,141,144,145,152,153,158,159,162,162,168,169,174,175,180,181,184,187,192,193,195,197,202,205,208,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,265,266,264,272,277,284,289,288,292};
               for(unsigned int i = 0; i < m_Zlist.size(); i++){
                    m_Alist.push_back(Alist[m_Zlist[i]-1]);
               }

               //calculate factors necessary for the potential
               for(unsigned int i = 0; i < m_Zlist.size(); i++){
                    m_xi.push_back(sqrt(3.0/2.0)/(0.836*pow(m_Alist[i],1.0/3.0)+0.57)*bohr_rad);
                    //std::cout << m_xi << std::endl;
               }
          }
          
          //overload () operator 
          double operator() (const coord_3d&r) const {
               double result = 0.0;
               int n = m_Zlist.size();
               for(int i = 0; i < n; i++){
                    double x = r[0] - m_Rlist[i][0];
                    double y = r[1] - m_Rlist[i][1];
                    double z = r[2] - m_Rlist[i][2];
                    double r = sqrt(x*x+y*y+z*z);
                    result += -m_Zlist[i]*erf(m_xi[i]*r)/r;
               }
               return result;
          }

          //some getters
          std::vector<int> get_Zlist(){
               return m_Zlist;
          }
          std::vector<coord_3d> get_Rlist(){
               return m_Rlist;
          }
          std::vector<int> get_Alist(){
               return m_Alist;
          }
          std::vector<double> get_xi(){
               return m_xi;
          }
          int get_num_atoms(){
               return m_Rlist.size();
          }
};


//Functor to make the fermi nuclear charge distribution (not normalized) for a given center
class FermiNucDistFunctor : public FunctionFunctorInterface<double,3> {
     private:
          int m_A;
          double m_T;
          double m_C;
          std::vector<coord_3d> m_R;
     public:
          // Constructor 
          FermiNucDistFunctor(int& Z, coord_3d R, double bohr_rad){
               //m_T = 0.000043463700858425666; //2.3 fm in bohr
               m_T = 2.3/bohr_rad;
               m_R.push_back(R);
               
               //find atomic mass numbers for each atom
               int Alist[116] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64,69,74,75,80,79,84,85,88,89,90,93,98,98,102,103,106,107,114,115,120,121,130,127,132,133,138,139,140,141,144,145,152,153,158,159,162,162,168,169,174,175,180,181,184,187,192,193,195,197,202,205,208,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,265,266,264,272,277,284,289,288,292};
               m_A = Alist[Z-1];

               double PI = constants::pi;
               if(m_A < 5){
                    m_C = 0.000022291*pow(m_A, 1.0/3.0) - 0.0000090676;
               }
               else{
                    m_C = sqrt(5.0/3.0*pow((0.836*pow(m_A,1.0/3.0)+0.570)/bohr_rad,2) - 7.0/3.0*pow(PI*m_T/4.0/log(3.0),2));
               }
          }
          
          //overload () operator 
          double operator() (const coord_3d&r) const {
               double x = r[0] - m_R[0][0];
               double y = r[1] - m_R[0][1];
               double z = r[2] - m_R[0][2];
               double rr = sqrt(x*x+y*y+z*z);
               double result = 1.0/(1.0+exp(4.0*log(3.0)*(rr-m_C)/m_T));
               return result;
          }

          //Because the distribution is only nonzero in a small window around the center, need to create a special point
          std::vector<coord_3d> special_points() const {
               return m_R;
          }

          madness::Level special_level() {
               return 18;
          }

          void print_details(World& world){
               int Alist[116] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,63,64,69,74,75,80,79,84,85,88,89,90,93,98,98,102,103,106,107,114,115,120,121,130,127,132,133,138,139,140,141,144,145,152,153,158,159,162,162,168,169,174,175,180,181,184,187,192,193,195,197,202,205,208,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,265,266,264,272,277,284,289,288,292};
               double T = 2.3/52917.72490083583;
               double PI = constants::pi;
               
               if(world.rank()==0){
                    for(int i = 0; i < 116; i++){
                         double RMS = (0.836*pow(Alist[i],1.0/3.0)+0.570)/52917.72490083583;
                         double C;
                         if(Alist[i] < 5){
                              C = 0.000022291*pow(Alist[i], 1.0/3.0) - 0.0000090676;
                         }
                         else{
                              C = sqrt(5.0/3.0*pow(RMS,2)-7.0/3.0*pow(PI*T/4.0/log(3.0),2));
                         }
                         double xi = 3.0/2.0/pow(RMS,2);
                         printf("Z: %3i,  A: %3i,  RMS: %.10e,  C: %.10e,  xi: %.10e\n", i+1, Alist[i], RMS, C, xi);
                    }
               }
          }
};

//generic f(r)=||r|| function for calculation of the radial expectation value
double myr(const coord_3d& r){
     return std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

//Creates the fermi nuclear potential. Also calculates the nuclear repulsion energy
void DF::make_fermi_potential(World& world, real_convolution_3d& op, real_function_3d& potential, double& nuclear_repulsion_energy){
     if(world.rank()==0) print("\n***Making a Fermi Potential***");
     
     //Get list of atom coordinates
     std::vector<coord_3d> Rlist = Init_params.molecule.get_all_coords_vec();
     std::vector<int> Zlist(Rlist.size());
     unsigned int num_atoms = Rlist.size();
     real_function_3d temp;
     double tempnorm;

     for(unsigned int i = 0; i < num_atoms; i++){
          Zlist[i] = Init_params.molecule.get_atom_number(i);
          FermiNucDistFunctor rho(Zlist[i], Rlist[i],DFparams.bohr_rad);
          temp = real_factory_3d(world).functor(rho).truncate_mode(0);
          tempnorm = temp.trace();
          temp.scale(-Zlist[i]/tempnorm);
          if(i == 0){
               potential = temp;
               //rho.print_details(world);
          }
          else{
               potential += temp;
          }
     }

     potential = apply(op,potential);


     nuclear_repulsion_energy = 0.0;
     double rr;
     for(unsigned int m = 0; m < num_atoms; m++){
          for(unsigned int n = m+1; n < num_atoms; n++){
               coord_3d dist = Rlist[m] - Rlist[n];
               rr = std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
               nuclear_repulsion_energy += Zlist[m]*Zlist[n]/rr;
          }
     }

}

// Collective constructor
DF::DF(World & world, const char* filename) : DF(world, (world.rank() == 0 ? std::make_shared<std::ifstream>(filename) : nullptr))
{}

// Constructor that actually does stuff
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

     // Read in archive
     if(DFparams.nwchem){
          Init_params.readnw(world, DFparams.archive);
     }
     else{
          Init_params.read(world, DFparams.archive, DFparams.restart);
     }
     if(world.rank() == 0){
          Init_params.print_params();
          print_molecule(world);
     }   

     // Set some function defaults   
     FunctionDefaults<3>::set_thresh(DFparams.thresh); //Always use user-specified thresh
     FunctionDefaults<3>::set_truncate_mode(1);   

     //If user requests different k, then project functions
     if(DFparams.k != Init_params.order){
          FunctionDefaults<3>::set_k(DFparams.k);
          for(unsigned int i = 0; i < Init_params.num_occupied; i++){
               for(unsigned int j = 0; j < 4; j++){
                    Init_params.orbitals[i][j] = project(Init_params.orbitals[i][j], FunctionDefaults<3>::get_k(), DFparams.thresh, false);
               }
               world.gop.fence();
               Init_params.orbitals[i].truncate();
          }
          if(DFparams.job == 1){
               for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
                    for(unsigned int j = 0; j < 4; j++){
                         Init_params.virtuals[i][j] = project(Init_params.virtuals[i][j], FunctionDefaults<3>::get_k(), DFparams.thresh, false);
                    }
                    world.gop.fence();
                    Init_params.virtuals[i].truncate();
               }
          }
     }

     //Set local orbitals and energies to those from the archive
     energies = Init_params.energies;
     v_energies = Init_params.v_energies;
     occupieds = Init_params.orbitals;
     virtuals = Init_params.virtuals;
     total_energy = Init_params.Init_total_energy;

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
     Tpsi[2] = (myc*myc)*(psiz[0] + psix[1] - myi*psiy[1] - 2*myi*psi[2]);
     Tpsi[3] = (myc*myc)*(psix[0] + myi*psiy[0] - psiz[1] - 2*myi*psi[3]);

     return Tpsi * (-myi);
}

//function to calculate the kinetic + rest energy expectation value using Dirac Hamiltonian c*\alpha*p+\Beta*m*c*c
double DF::rele(World& world, Fcwf& psi){
     Fcwf Tpsi = apply_T(world, psi);
     
     std::complex<double> energy  = inner(psi, Tpsi);
     //if(world.rank()==0) print("   in rele: ", energy.real());

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
     for(unsigned int i = 0; i < n; i++){

          std::vector<complex_function_3d> temp(n-i);
          for(unsigned int j = 0; j < n-i; j++){
               temp[j] = complex_factory_3d(world);    
          }
          compress(world, temp);

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

          gaxpy(world, 1.0, temp, 1.0, occupieds[i][0]*conj(world,temp0));
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][1]*conj(world,temp1));
          gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][2]*conj(world,temp2));
          gaxpy(world, 1.0, temp, 1.0/(myc*myc), occupieds[i][3]*conj(world,temp3));

          truncate(world, temp);

          if(world.rank()==0) print(i, "Starting apply phase in K");

          temp = apply(world, op, temp);
          
          truncate(world,temp);

          if(world.rank()==0) print(i, "Exiting apply phase in K");

          Kpsis[i][0] += sum(world, mul(world, temp, temp0));
          Kpsis[i][1] += sum(world, mul(world, temp, temp1));
          Kpsis[i][2] += sum(world, mul(world, temp, temp2));
          Kpsis[i][3] += sum(world, mul(world, temp, temp3));
          
          if(world.rank()==0) print(i, "Exiting sum block in K");

          temp = conj(world, temp);

          temp0 = occupieds[i][0]*temp;
          temp1 = occupieds[i][1]*temp;
          temp2 = occupieds[i][2]*temp;
          temp3 = occupieds[i][3]*temp;

          if(world.rank()==0) print(i, "Entering final loop in K");

          for(unsigned int j = i+1; j < n; j++){
               Kpsis[j][0] += temp0[j-i];
               Kpsis[j][1] += temp1[j-i];
               Kpsis[j][2] += temp2[j-i];
               Kpsis[j][3] += temp3[j-i];
          }
          
          if(world.rank()==0) print(i, "Exiting final loop in K");

          Kpsis[i].truncate();
     }

     //Report time
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
}

//This function is old an only used in the virtual solver (that doesn't work)
Fcwf DF::apply_K(World& world, real_convolution_3d& op, Fcwf& phi){
     throw;
     complex_function_3d temp(world);
     Fcwf result(world);

     unsigned int n = Init_params.num_occupied;
     for(unsigned int i = 0; i < n; i++){
          temp = inner_func(world, occupieds[i], phi);
          temp.truncate();
          temp = apply(op,temp);
          if(i == 0){
               result = occupieds[i]*temp;
          }
          else{
               result += occupieds[i]*temp;
          }
     }
     result.truncate();
     return result;
}

//Diagonalize psis in the Fock space. psis is modified in place. requires Kpsis to be precomputed
//Kpsis are transformed in place.
void DF::diagonalize(World& world, real_function_3d& myV, real_convolution_3d& op, std::vector<Fcwf>& Kpsis){

     if(world.rank()==0) print("\n***Diagonalizing***");
     start_timer(world);

     unsigned int n = Init_params.num_occupied;
     Tensor<std::complex<double>> fock(n, n);
     Tensor<std::complex<double>> overlap(n, n);
     Tensor<std::complex<double>> U(n,n);
     Tensor<double> evals(n);
     std::vector<Fcwf> temp_orbitals;

     if(world.rank()==0) print("     Forming Matrices");
     start_timer(world);
     
     ////Form the Fock Matrix

     //calculate coulomb part
     if(world.rank() == 0) print("          Adding (V+J)psi");
     real_function_3d rho = real_factory_3d(world);
     for(unsigned int j = 0; j < n; j++){
          rho += squaremod(occupieds[j]);
     }
     
     //TODO: Here try moving the apply out to operate on the sum of the nuclear and electronic charge distributions
     real_function_3d potential = myV + apply(op,rho);
     potential.truncate();

     //add in coulomb parts to neworbitals
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals.push_back(occupieds[j]*potential); //add in coulomb term
     }

     if(world.rank() == 0) print("          Subtracting K*psi");
     //Move Kpsis to new orbitals, as they are part of the fock operator
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] -= Kpsis[j]; //yes this needs to be subtraction. exchange function doesn't include the negative.
     }

     //add in T_psi
     if(world.rank()==0) print("          Adding T*psi");
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] += apply_T(world, occupieds[j]);  //add in "kinetic" term
     }

     if(world.rank()==0) print("          Integrating to form Fock Matrix");
     start_timer(world);

     //Now compute the fock matrix
     fock = matrix_inner(world, occupieds, temp_orbitals);

     //symmetrize
     fock = (1.0/2.0)*(fock + conj_transpose(fock));
     
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("               ", times[0]); 


     ////and the overlap matrix
     if(world.rank()==0) print("          Integrating to form Overlap Matrix");
     start_timer(world);
     overlap = matrix_inner(world,occupieds,occupieds);
     
     //symmetrize
     overlap = (1.0/2.0)*(overlap + conj_transpose(overlap));

     times = end_timer(world);
     if(world.rank()==0) print("               ", times[0]);
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     //debugging: print fock and overlap matrices
     //if(world.rank()==0){
     //     print("fock:\n", fock);
     //     print("\noverlap:\n", overlap);
     //}
     
     if(world.rank()==0) print("     Eigensolver");
     start_timer(world);

     //Diagonalize
     sygv(fock, overlap, 1, U, evals);
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     //debugging: print matrix of eigenvectors
     //if(world.rank()==0) print("U:\n", U);

     //Before applying the transformation, fix arbitrary rotations introduced by the eigensolver. 
     if(world.rank()==0) print("     Removing Rotations");
     start_timer(world);

     double thresh_degenerate = DFparams.thresh*10.0;
     double csquared = 137.0359895*137.0359895; //electron rest energy

     //swap columns for a diagonally dominant matrix
     bool switched = true;
     while (switched) {
          switched = false;
          for (unsigned int kk = 0; kk < Init_params.num_occupied; kk++) {
               for (unsigned int j = kk + 1; j < Init_params.num_occupied; j++) {
                    double sold = std::real(U(kk,kk)*std::conj(U(kk,kk))) + std::real(U(j,j)*std::conj(U(j,j)));
                    double snew = std::real(U(kk,j)*std::conj(U(kk,j))) + std::real(U(j,kk)*std::conj(U(j,kk)));
                    //if (snew > sold and not ((evals[j] - evals[kk]) > thresh_degenerate * std::max(std::fabs(evals[kk])-csquared,1.0)) ) {
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
     for (unsigned int kk = 0; kk < Init_params.num_occupied; ++kk)
          U(_, kk).scale(std::conj(U(kk,kk))/std::abs(U(kk,kk)));

     
     //Find clusters of degenerate eigenvalues and rotate eigenvectors to maximize overlap with previous ones
     unsigned int ilo = 0; // first element of cluster
     if(world.rank()==0) print("          Degeneracy threshold: ",thresh_degenerate);
     
     while (ilo < Init_params.num_occupied - 1) {
         unsigned int ihi = ilo;
         while (fabs(evals[ilo] - evals[ihi + 1])
                < thresh_degenerate * std::fabs(evals[ilo])){// pow(10,floor(log10(std::fabs(evals[ilo]))))) { 
                //< thresh_degenerate * std::max(std::fabs(evals[ilo]-csquared),1.0)){// pow(10,floor(log10(std::fabs(evals[ilo]))))) { 
             ++ihi;
             if (ihi == Init_params.num_occupied - 1)
                 break;
         }
         unsigned int nclus = ihi - ilo + 1;
         if (nclus > 1) {
              if(world.rank()==0){
                    print("          found cluster from ", ilo + 1, " to " , ihi + 1);
                    for(unsigned int kk = ilo; kk <= ihi; kk++){
                         print("               ",evals[kk]);
                    }
              }


              //Use the polar decomposition to undo rotations:
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

     //Debugging: Print transformation matrix after rotation removal
     //if(world.rank()==0) print("U:\n", U);
     
     
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);

     if(world.rank()==0) print("     Applying Transformation");
     start_timer(world);

     ////Apply the transformation to the Exchange
     transform(world, Kpsis, U);


     ////Apply the transformation to the orbitals 
     transform(world, occupieds, U);

     //truncate
     for(int kk = 0; kk < Init_params.num_occupied; kk++){
           Kpsis[kk].truncate();
           occupieds[kk].truncate();
     }


     //truncate

     //debugging
     //for(unsigned int j=0; j < Init_params.num_occupied; j++){
     //     double tempdouble = rele(world, occupieds[j]);
     //     if(world.rank()==0) print("   after diag, rele ",j," = ",tempdouble);
     //}
     

     //Set energies = evals, and fix the energy printing stage.
     energies = evals;

     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

}

//returns a vector of orthonormal Fcwfs constructed from the input vector of Fcwfs
std::vector<Fcwf> orthogonalize(World& world, std::vector<Fcwf> orbitals){
     throw;
     int n = orbitals.size();
     std::vector<Fcwf> result;
     for(int i = 0; i < n; i++){
          result.push_back(copy(orbitals[i]));
     }
     std::complex<double> r(0.0,0.0);

     for(int i = 0; i < n; i++){
          result[i].normalize();
          for(int j = i+1; j < n; j++){
               r = inner(result[i],result[j]);
               result[j] -= result[i]*r;
          }
     }
     return result;

}

//faster orthogonalize that modifies the input functions in place
//TODO: The function below mimics one from SCF.cc. In the future we will probably want a different variant (which also should exist in SCF.cc or somwhere like that) that treats core and valence orbitals differently: i.e. doesn't mix valence orbitals into the core orbitals
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
          transform(world, occupieds, Q);
     } while (maxq>0.01);

     //normalize afterward
     for(unsigned int i = 0; i < n; i++){
          occupieds[i].normalize();
     }

}

//This is an old function for virtual orbitals that probably doesn't work anymore
//Different diagonalize function that doesn't try to calculate a coulomb part internally
Tensor<double> DF::diagonalize_virtuals(World& world, real_function_3d& JandV,real_convolution_3d& op, std::vector<Fcwf>& Kpsis){

     if(world.rank()==0) print("\n***Diagonalizing***");

     unsigned int n = Init_params.num_virtuals;
     Tensor<std::complex<double>> fock(n, n);
     Tensor<std::complex<double>> overlap(n, n);
     Tensor<std::complex<double>> U(n,n);
     Tensor<double> evals(n);
     std::vector<Fcwf> temp_orbitals;

     if(world.rank()==0) print("     Forming Matrices");
     
     ////Form the Fock Matrix
     //first apply Fock operator to all orbitals
     
     if(world.rank() == 0) print("          Moving K*psi");
     //Move Kpsis to new orbitals, as they are part of the fock operator
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals.push_back(Kpsis[j]);
          temp_orbitals[j].scale(-1.0);
     }

     //add in coulomb parts to neworbitals
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] += virtuals[j]*JandV; //add in coulomb term
     }

     //add in T_psi
     if(world.rank()==0) print("          Adding T*psi");
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] += apply_T(world, virtuals[j]);  //add in "kinetic" term
     }

     //We now have a vector of F*psi for each psi (in temp_orbitals).
     //Calculate each element of the fock matrix
     if(world.rank()==0) print("          Integrating to form Fock Matrix");
     for(unsigned int j = 0; j < n; j++){
          for(unsigned int k = 0; k < n; k++){
               fock(j,k) = inner(virtuals[j],temp_orbitals[k]); 
          }
     }

     ////Form the overlap matrix
     if(world.rank()==0) print("          Integrating to form Overlap Matrix");
     for(unsigned int j = 0; j < n; j++){
          for(unsigned int k = 0; k < n; k++){
               overlap(j,k) = inner(virtuals[j],virtuals[k]);
          }
     }
     
     //debugging: print fock and overlap matrices
     //if(world.rank()==0){
     //     print("real(fock):\n", real(fock));
     //     print("imag(fock):\n", imag(fock));
     //     print("\nreal(overlap):\n", real(overlap));
     //     print("\nimag(overlap):\n", imag(overlap));
     //}
     

     if(world.rank()==0) print("     Eigensolver");

     //Diagonalize
     sygv(fock, overlap, 1, U, evals);

     //debugging: print matrix of eigenvectors and eigenvalues
     //if(world.rank()==0) print("U:\n", U);
     //if(world.rank()==0) print("\nevals:\n",evals);

     //Before applying the transformation, fix arbitrary rotations introduced by the eigensolver. 
     if(world.rank()==0) print("     Removing Rotations");

     double thresh_degenerate = DFparams.thresh*100.0;
     double csquared = 137.0359895*137.0359895; //electron rest energy
     //swap columns for a diagonally dominant matrix
     bool switched = true;
     while (switched) {
          switched = false;
          for (unsigned int kk = 0; kk < n; kk++) {
               for (unsigned int j = kk + 1; j < n; j++) {
                    //double sold = std::norm(U(kk, kk)) + std::norm(U(j, j));
                    //double snew = std::norm(U(kk, j)) + std::norm(U(j, kk));
                    double sold = std::real(U(kk,kk)*std::conj(U(kk,kk))) + std::real(U(j,j)*std::conj(U(j,j)));
                    double snew = std::real(U(kk,j)*std::conj(U(kk,j))) + std::real(U(j,kk)*std::conj(U(j,kk)));
                    if (snew > sold and not ((evals[j] - evals[kk]) > thresh_degenerate * std::max(std::fabs(evals[kk])-csquared,1.0)) ) {
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
     for (unsigned int kk = 0; kk < n; ++kk)
          U(_, kk).scale(std::conj(U(kk,kk))/std::abs(U(kk,kk)));

     
     //Find clusters of degenerate eigenvalues and rotate eigenvectors to maximize overlap with previous ones
     unsigned int ilo = 0; // first element of cluster
     if(world.rank()==0) print("          Degeneracy threshold: ",thresh_degenerate);
     
     while (ilo < n - 1) {
         unsigned int ihi = ilo;
         while (fabs(evals[ilo] - evals[ihi + 1])
                < thresh_degenerate * std::max(std::fabs(evals[ilo]-csquared),1.0)){// pow(10,floor(log10(std::fabs(evals[ilo]))))) { 
             ++ihi;
             if (ihi == n - 1)
                 break;
         }
         unsigned int nclus = ihi - ilo + 1;
         if (nclus > 1) {
              if(world.rank()==0){
                    print("          found cluster from ", ilo + 1, " to " , ihi + 1);
                    for(unsigned int kk = ilo; kk <= ihi; kk++){
                         print("               ",evals[kk] - csquared);
                    }
              }


              //Use the polar decomposition to undo rotations:
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


     //if(world.rank()==0) print("U:\n", U);
     
     
     fock = inner(conj_transpose(U), inner(fock, U));
     //if(world.rank()==0)print("\n U^* F U:\n",fock,"\n");

     if(world.rank()==0) print("     Applying Transformation");

     ////Apply the transformation to the Exchange
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] = Fcwf(world);
          for(unsigned int k = 0; k < n; k++){
               temp_orbitals[j] += Kpsis[k]*U(k,j);
          }
     }
     for(unsigned int m = 0; m < n; m++){
           Kpsis[m] = temp_orbitals[m];
     }

     ////Apply the transformation to the orbitals
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals[j] = Fcwf(world);
          for(unsigned int k = 0; k < n; k++){
               temp_orbitals[j] += virtuals[k]*U(k,j); 
          }
     }
     for(unsigned int j = 0; j < n; j++){
          virtuals[j] = temp_orbitals[j];
     }

     return evals;

}

//Apply's Green's function to Vpsi (a Fcwf). Overwrites Vpsi with new Fcwf
void apply_BSH(World& world, Fcwf& Vpsi, double& eps, double& small, double& thresh){


     //necessary constants
     double myc = 137.0359895; //speed of light
     double c2 = myc*myc;
     std::complex<double> myi(0,1); //imaginary number
     std::complex<double> ic = myi*myc;
    
     //calculate exponent for equivalent BSH operator
     double mu = std::sqrt((myc*myc*myc*myc-eps*eps)/myc/myc);

     //if(world.rank() == 0) print("Hi, this is apply_BSH! mu is: ", mu);

     //create gradient BSH operators
     world.gop.fence();

     //create BSH operator
     if(world.rank()==0) print("mu: ", mu);
     std::shared_ptr<real_convolution_3d> op = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world, mu,small,thresh)); 
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator(world, mu, small, thresh); 
     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator_Joel(world, mu, 1e-8, thresh); 
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> allops(16);
     for(unsigned int i = 0; i < 4; i++){
          allops[i] = op;
          allops[4+i] = op3[0];
          allops[8+i] = op3[1];
          allops[12+i] = op3[2];
     }

     //create intermediate functions necessary to compute new components
     //ttt = wall_time();
     std::vector<complex_function_3d> temp(16);
     for(unsigned int i = 0; i < 4; i++){
          temp[i] = Vpsi[i];
          temp[4+i] = Vpsi[i];
          temp[8+i] = Vpsi[i];
          temp[12+i] = Vpsi[i];
     }

     temp = apply(world, allops, temp);

     //build transformation tensor
     Tensor<std::complex<double>> U(16,4);
     U(0,0) = c2+eps; U(14,0) = -ic; U(7,0) = -ic; U(11,0) = -myc;
     U(1,1) = c2+eps; U(6,1) = -ic; U(10,1) = myc; U(15,1) = ic;
     U(12,2) = -ic; U(5,2) = -ic; U(9,2) = -myc; U(2,2) = eps-c2;
     U(4,3) = -ic; U(8,3) = myc; U(13,3) = ic; U(3,3) = eps-c2;
     U *= (1.0/c2);

     //Apply transformation tensor and store back in Vpsi
     temp = transform(world, temp, U);
     Vpsi[0] = temp[0];
     Vpsi[1] = temp[1];
     Vpsi[2] = temp[2];
     Vpsi[3] = temp[3];

}

//Apply's Green's function to Vpsi (a Fcwf). Overwrites Vpsi with new Fcwf
//
//Instead of using the derivative of the Green's function, first applies the
//nonrelativistic Green's function, then applied (H_D + eps) to the result
void apply_BSH_new(World& world, Fcwf& Vpsi, double& eps, double& small, double& thresh){


     //necessary constants
     double myc = 137.0359895; //speed of light
     double c2 = myc*myc;
     std::complex<double> myi(0,1); //imaginary number
     std::complex<double> ic = myi*myc;
    
     //calculate exponent for equivalent BSH operator
     double mu = std::sqrt(-(2*eps*c2+eps*eps)/c2);

     //if(world.rank()==0) print("    mu = ", mu);

     world.gop.fence();

     //create BSH operator
     real_convolution_3d op = BSHOperator3D(world, mu,small,thresh); // Finer length scale and accuracy control

     //Apply BSH operator to Vpsi
     Vpsi = apply(world, op, Vpsi);

     //mu = (apply_T(world,Vpsi)*(1.0/c2)).norm2();
     //if(world.rank()==0) print("    after BSH and T= ", mu);

     //Apply (1/c^2)(H_D + eps) to Vpsi. Using apply_T for convenience, but this requires adding 2c^2Vpsi
     Vpsi = apply_T(world, Vpsi)*(1.0/c2) + Vpsi * ((eps+2*c2)/c2);

     //mu = Vpsi.norm2();
     //if(world.rank()==0) print("    after full apply = ", mu);

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


void DF::saveDF(World& world){
     Tensor<double> times = get_times(world);
     start_timer(world);
     if(world.rank()==0) print("\n***Saving at time: ",times[0]," ***");
     try{
          archive::ParallelOutputArchive output(world, DFparams.savefile.c_str(), 1);
          output & total_energy & Init_params.spinrestricted & Init_params.num_occupied & energies & Init_params.L & Init_params.order & Init_params.molecule;
          for(unsigned int i = 0; i < Init_params.num_occupied; i++){
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

//Creates the nuclear potential from the molecule object
void DF::make_gaussian_potential(World& world, real_function_3d& potential){
     if(world.rank()==0) print("\n***Making a Gaussian Potential***");
     GaussianNucleusFunctor Vfunctor(Init_params.molecule, DFparams.bohr_rad);
     potential = real_factory_3d(world).functor(Vfunctor).truncate_mode(0).truncate_on_project();
}

//Creates the nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
void DF::make_gaussian_potential(World& world, real_function_3d& potential, double& nuclear_repulsion_energy){
     if(world.rank()==0) print("\n***Making a Gaussian Potential***");
     GaussianNucleusFunctor Vfunctor(Init_params.molecule,DFparams.bohr_rad);
     potential = real_factory_3d(world).functor(Vfunctor).truncate_mode(0).truncate_on_project();
     std::vector<coord_3d> Rlist = Vfunctor.get_Rlist();
     std::vector<int> Zlist = Vfunctor.get_Zlist();
     nuclear_repulsion_energy = 0.0;
     double rr;
     int num_atoms = Rlist.size();
     for(int m = 0; m < num_atoms; m++){
          for(int n = m+1; n < num_atoms; n++){
               coord_3d dist = Rlist[m] - Rlist[n];
               rr = std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
               nuclear_repulsion_energy += Zlist[m]*Zlist[n]/rr;
          }
     }

}


void DF::DF_load_balance(World& world, real_function_3d& Vnuc){
     if(world.rank()==0) print("\n***Load Balancing***");
     start_timer(world);
     LoadBalanceDeux<3> lb(world);
     lb.add_tree(Vnuc, lbcost<double,3>(12.0,96.0),true);
     //Commenting out below block to test memory issues
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          for(int kk = 0; kk < 4; kk++){
               lb.add_tree(occupieds[j][kk], lbcost<std::complex<double>,3>(24.0,192.0),true);
          }
     }
     FunctionDefaults<3>::redistribute(world, lb.load_balance(2), true);
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
     
}

void DF::make_component_lineplots(World& world, const char* filename1, const char* filename2, int npt, double endpnt){
     std::vector<real_function_3d> large_densities;
     std::vector<real_function_3d> small_densities;
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities.push_back(squaremod_large(occupieds[i]));
          small_densities.push_back(squaremod_small(occupieds[i]));
     }

     //double h = Init_params.L*(1.0/(npt-1));
     double h = endpnt*(1.0/(npt-1));
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities[i].reconstruct();
          small_densities[i].reconstruct();
     }
     if(world.rank()==0){
          FILE* file1 = fopen(filename1,"w");
          FILE* file2 = fopen(filename2,"w");
          if(!file1) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          if(!file2) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               for(int j=0; j<npt; ++j){
                    coord_3d r({j*h,0.0,0.0});
                    fprintf(file1,"%.14e ", j*h);
                    fprintf(file2,"%.14e ", j*h);
                    plot_line_print_value(file1, large_densities[i].eval(r));
                    plot_line_print_value(file2, small_densities[i].eval(r));
                    fprintf(file1,"\n");
                    fprintf(file2,"\n");
               }
               fprintf(file1,"\n");
               fprintf(file2,"\n");
          }
          fclose(file1);
          fclose(file2);
     }
     world.gop.fence();

}

void DF::make_component_logplots(World& world, const char* filename1, const char* filename2, int npt, int startpnt, int endpnt){
     std::vector<real_function_3d> large_densities;
     std::vector<real_function_3d> small_densities;
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities.push_back(squaremod_large(occupieds[i]));
          small_densities.push_back(squaremod_small(occupieds[i]));
     }

     //double h = Init_params.L*(1.0/(npt-1));
     double h = (endpnt-startpnt)*(1.0/(npt-1));
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          large_densities[i].reconstruct();
          small_densities[i].reconstruct();
     }
     if(world.rank()==0){
          FILE* file1 = fopen(filename1,"w");
          FILE* file2 = fopen(filename2,"w");
          if(!file1) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          if(!file2) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               for(int j=0; j<npt; ++j){
                    coord_3d r({std::pow(10,startpnt + j*h),0.0,0.0});
                    fprintf(file1,"%.14e ", r[0]);
                    fprintf(file2,"%.14e ", r[0]);
                    plot_line_print_value(file1, large_densities[i].eval(r));
                    plot_line_print_value(file2, small_densities[i].eval(r));
                    fprintf(file1,"\n");
                    fprintf(file2,"\n");
               }
               fprintf(file1,"\n");
               fprintf(file2,"\n");
          }
          fclose(file1);
          fclose(file2);
     }
     world.gop.fence();

}

void DF::make_density_lineplots(World& world, const char* filename, int npt, double endpnt){
     std::vector<real_function_3d> densities;
     real_function_3d zero(world);
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          densities.push_back(squaremod(occupieds[i]));
     }

     double h = endpnt*(1.0/(npt-1));
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          densities[i].reconstruct();
     }
     if(world.rank()==0){
          FILE* file = fopen(filename,"w");
          if(!file) MADNESS_EXCEPTION("DF density lineplots: failed to open the plot file",0);
          for(unsigned int i=0; i<Init_params.num_occupied; i++){
               for(int j=0; j<npt; ++j){
                    coord_3d r({j*h,0.0,0.0});
                    fprintf(file,"%.14e ", j*h);
                    plot_line_print_value(file, densities[i].eval(r));
                    fprintf(file,"\n");
               }
               fprintf(file,"\n");
          }
          fclose(file);
     }
     world.gop.fence();

}

//One iteration
bool DF::iterate(World& world, real_function_3d& V, real_convolution_3d& op, real_function_3d& JandV, std::vector<Fcwf>& Kpsis, XNonlinearSolver<std::vector<Fcwf>, std::complex<double>, Fcwf_vector_allocator>& kainsolver, double& tolerance, int& iteration_number, double& nuclear_repulsion_energy){

     Tensor<double> times = get_times(world);

     if(world.rank()==0) print("\n\n\nIteration: ", iteration_number, " at ",times[0]);
     if(world.rank()==0) print("--------------");
     start_timer(world);

     std::vector<Fcwf> Residuals;
     Fcwf temp_function(world);
     double residualnorm;
     real_function_3d rho = real_factory_3d(world);

     bool iterate_again = false; //Assume iterations will stop

     //Diagonalize
     diagonalize(world, V, op, Kpsis);

     //Diagonalization forces us to recompute J
     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
          rho += squaremod(occupieds[kk]);
     }
     JandV = V + apply(op,rho); 
     JandV.truncate();



     //Apply BSH to each psi
     if(world.rank()==0) print("\n***Applying BSH operator***");
     start_timer(world);
     double maxresidual = -1.0;
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){

          //construct the function to which we will apply the BSH
          //occupieds[j].truncate();
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

          //Print residual norm to keep track
          if(world.rank()==0) printf("     Orbital: %3i,  Resid: %.10e\n",j+1, residualnorm);

          //Compare so that at the end we know the max residual
          maxresidual = std::max(maxresidual, residualnorm);

          //If the norm is big enough, we'll need to iterate again.
          if(residualnorm > tolerance) iterate_again = true;

          //Store residual function if we're using KAIN. Not necessary if we're not using kain
          if(iteration_number != 1 and DFparams.kain){
               Residuals.push_back(occupieds[j] - temp_function);
          }
          else{
               occupieds[j] = temp_function; //if not using KAIN, then just use the new orbital
          }
     }
     if(world.rank()==0) printf("                max Resid: %.10e\n",maxresidual);
     if(world.rank()==0) printf("                tolerance: %.10e\n",tolerance);
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //debugging
     //for(unsigned int j=0; j < Init_params.num_occupied; j++){
     //     double tempdouble = rele(world, occupieds[j]);
     //     if(world.rank()==0) print("   after BSH, rele ",j," = ",tempdouble);
     //}
     

     //Apply the kain solver, if called for
     if(iteration_number != 1 and DFparams.kain){
          if(world.rank()==0) print("\n***Applying KAIN Solver***");
          start_timer(world);

          //occupieds = kainsolver.update(occupieds, Residuals);

          //Replace above line with kain + step restriction
          Residuals = kainsolver.update(occupieds, Residuals); //Using Residuals for updated orbitals to save storage
          for(unsigned int i=0; i < Init_params.num_occupied; i++){

               //see how big of a step KAIN took for each orbital
               residualnorm = (occupieds[i]-Residuals[i]).norm2();

               //Restrict the step taken by KAIN if it's took big
               if(residualnorm > DFparams.maxrotn){
                    double s = DFparams.maxrotn / residualnorm;
                    if(world.rank()==0) print("     restricting step for orbital: ", i+1);
                    occupieds[i] = Residuals[i]*s + occupieds[i]*(1.0-s);
               }
               else{
                    occupieds[i] = Residuals[i];
               }
          }

          times = end_timer(world);
          if(world.rank()==0) print("     ", times[0]);
     }

     //truncate
     for(unsigned int i = 0; i < Init_params.num_occupied; i++) occupieds[i].truncate();

     //orthogonalize
     if(world.rank()==0) print("\n***Orthonormalizing***");
     start_timer(world);

     orthogonalize_inplace(world);

     //truncate here and normalize again
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
           occupieds[i].truncate();
           occupieds[i].normalize();
     }

     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //calculate new exchange. Has timer built in.
     if(world.rank()==0) print("\n***Recalculating Exchange***");
     exchange(world, op, Kpsis);

     //Calculate new J+V term
     if(world.rank()==0) print("\n***Recalculating Coulomb***");
     start_timer(world);
     rho = real_factory_3d(world);
     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
          rho += squaremod(occupieds[kk]);
     }
     JandV = V + apply(op,rho); 
     JandV.truncate();
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Calculate and print total energy
     //Simultaneously update eps
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

     real_function_3d Jop = apply(op,rho);

     //Compute kinetic energy contribution
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          //energies[j] = rele(world,occupieds[j]);
          //kinetic_energy += (energies[j]);
          kinetic_energy += rele(world, occupieds[j]);
     }
          
     //Compute electron-nuclear attraction energy contribution, taking advantage of vmra's inner
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
     nuclear_attraction_energy = nuclear_attraction_tensor.sum();
     
     //Compute electron-electron repulsion energy contribution, again using vmra
     coulomb_tensor = real(inner(world,occupieds1,mul(world,Jop,occupieds1)));
     coulomb_tensor += real(inner(world,occupieds2,mul(world,Jop,occupieds2)));
     coulomb_tensor += real(inner(world,occupieds3,mul(world,Jop,occupieds3)))*(1.0/(myc*myc));
     coulomb_tensor += real(inner(world,occupieds4,mul(world,Jop,occupieds4)))*(1.0/(myc*myc));
     coulomb_energy = 0.5*coulomb_tensor.sum();
     
     //Calculate Exchange energy contribution
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
     exchange_energy = 0.5*exchange_tensor.sum();

     //Loop through energies and calculate their new values using the individual contributions
     //(inside the "tensor" variables
     //for(unsigned int i = 0; i < Init_params.num_occupied; i++){
     //     energies[i] += nuclear_attraction_tensor[i]+coulomb_tensor[i] - exchange_tensor[i];
     //}

     //compute total energy
     total_energy = kinetic_energy + coulomb_energy - exchange_energy + nuclear_attraction_energy + nuclear_repulsion_energy;

     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Need r function to print r expectation values
     real_function_3d rfunc = real_factory_3d(world).f(myr);

     //Loop through occupied orbitals. Calculate <r>, number of coefficients, max depth.
     //Print these along with updated energies
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){

          //calculate <r>
          double r_expec = std::real(inner(occupieds[j], occupieds[j]*rfunc));

          //find number of coefficients
          int numcoeffs = occupieds[j][0].size() + occupieds[j][1].size() + occupieds[j][2].size() + occupieds[j][3].size();

          //find maximum depth
          int maxdepth = std::max(occupieds[j][0].max_depth(), occupieds[j][1].max_depth());
          maxdepth = std::max(int(occupieds[j][2].max_depth()), maxdepth);
          maxdepth = std::max(int(occupieds[j][3].max_depth()), maxdepth);

          //Print everything
          if(world.rank()==0){
               printf("                Orbital: %3i, Energy: %.10e, <r>: %8e, No. coeffs: %7i, Max depth: %3i\n",j+1, energies[j], r_expec, numcoeffs, maxdepth);

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
     
     //check total energy for convergence.
     //if(std::fabs(total_energy-old_total_energy) > DFparams.thresh*pow(10,floor(log10(std::fabs(total_energy))))){
     //     iterate_again = true; 
     //}

     //truncate
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          occupieds[j].truncate();
     }

     times = end_timer(world);
     if(world.rank()==0) print("     Iteration time:", times[0]);
 
     return iterate_again;
}

// Solves the for occupied orbitals
void DF::solve_occupied(World & world)
{

     //State what we're doing here
     if(world.rank()==0) print("\nSolving for ", Init_params.num_occupied, " occupied orbitals\n-----------------------------------\n");

     //Will need a coulomb operator
     real_convolution_3d op = CoulombOperator(world,DFparams.small,DFparams.thresh);

     //allocator is useful to have 
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

     //Initialize vector of Fcwfs to hold exchange applied to Psis
     std::vector<Fcwf> Kpsis = allocator();

     //Calculate initial exchange
     if(world.rank()==0) print("\n***Calculating Initial Exchange***");
     exchange(world, op, Kpsis);

     //Calculate initial J+V term
     if(world.rank()==0) print("\n***Calculating Initial Coulomb***");
     start_timer(world);
     real_function_3d rho = real_factory_3d(world);
     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
          rho += squaremod(occupieds[kk]);
     }
     real_function_3d JandV = Vnuc + apply(op,rho); 
     JandV.truncate();
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Set tolerance for residuals
     //double tol = pow(10,floor(0.5*log10(DFparams.thresh)));
     double tol = 50.0*DFparams.thresh; 
     

     //Now time to start iterating
     bool keep_going = true;
     int iteration_number = 1;
     //while(iteration_number < DFparams.max_iter){
     while((keep_going and iteration_number <= DFparams.max_iter) or iteration_number <= DFparams.min_iter){
          keep_going = iterate(world, Vnuc, op, JandV, Kpsis, kainsolver, tol, iteration_number, nuclear_repulsion_energy);
          
          //Load balance and save between iterations
          if(keep_going and iteration_number <= DFparams.lb_iter) DF_load_balance(world, Vnuc);
          if(DFparams.do_save) saveDF(world);

          //Increment iteration counter
          iteration_number++;
     }



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
          else if(DFparams.job == 1){
               solve_virtuals1(world);
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

void DF::solve_virtuals1(World& world){

     //First, move final occupied orbital to be first virtual orbital
     Fcwf final_occupied = copy(occupieds.back());
     occupieds.pop_back(); 
     Init_params.num_occupied -= 1;
     virtuals.emplace(virtuals.begin(),final_occupied);
     Init_params.num_virtuals += 1;

     //Do the same for energies
     double temp_energy = energies(energies.dim(0)-1);
     energies = energies(Slice(0,energies.dim(0)-2));
     Tensor<double> temp_v_energies(Init_params.num_virtuals);
     temp_v_energies[0] = temp_energy;
     for(unsigned int i = 1; i < Init_params.num_virtuals; i++){
          temp_v_energies(i) = v_energies(i-1);
     }
     v_energies = temp_v_energies;
     
     //First, perform Dirac Fock on the n-1 occupied orbitals
     solve_occupied(world);

     //Next, calculate the last occupied orbital simultaneously with the virtuals
     if(world.rank()==0) print("***Calculating Final Occupied and Virtual Orbitals***");

     //For now, recalculate JandV. In the future, maybe get this out of solve_occupied?
     //Will need a coulomb operator
     real_convolution_3d op = CoulombOperator(world,DFparams.small,DFparams.thresh);
     //Form nuclear potential
     real_function_3d Vnuc;
     double nuclear_repulsion_energy;
     if(DFparams.nucleus == 1){
          make_fermi_potential(world, op, Vnuc, nuclear_repulsion_energy);
     }
     else{
          make_gaussian_potential(world, Vnuc, nuclear_repulsion_energy);
     }
     //Calculate initial J+V term
     real_function_3d rho = real_factory_3d(world);
     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
          rho += squaremod(occupieds[kk]);
     }
     real_function_3d JandV = Vnuc + apply(op,rho); 

     //Set tolerance for residuals
     double tol = pow(10,floor(0.5*log10(DFparams.thresh)));

     //Orthogonalize the functions against the occupieds
     std::complex<double> tempcomplex;
     for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
          for(unsigned int j = 0; j < Init_params.num_occupied; j++){
               tempcomplex = inner(virtuals[i], occupieds[j]);
               virtuals[i] -= occupieds[j]*tempcomplex;
          }
          virtuals[i].normalize();
     }

     //orthonormalize the functions internally
     orthogonalize(world,virtuals);

     //Initialize exchange vector
     Fcwf_vector_allocator allocator(world,Init_params.num_virtuals);
     std::vector<Fcwf> Kpsis = allocator();
     for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
          Kpsis[i] = apply_K(world, op, virtuals[i]);
     }

     bool keep_going = true;
     int iteration_number = 1;
     double csquared = 137.0359895*137.0359895; //electron rest energy
     std::vector<Fcwf> Residuals = allocator();
     XNonlinearSolver<std::vector<Fcwf>, std::complex<double>, Fcwf_vector_allocator> kainsolver(allocator);
     kainsolver.set_maxsub(DFparams.maxsub);
     Fcwf tempfcwf(world);
     while(keep_going and iteration_number < DFparams.max_iter){
          
          if(world.rank()==0) print("Iteration:", iteration_number);

          keep_going = false;
          
          //diagonalize
          v_energies = diagonalize_virtuals(world, JandV, op, Kpsis);

          //Apply BSH to each function individually
          for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
               
               //Calculate function to apply BSH to
               tempfcwf = virtuals[i] * JandV;
               tempfcwf.scale(-1.0);
               tempfcwf += Kpsis[i];

               //if(world.rank()==0) print("Going into BSH, energy is:",v_energies[i] - csquared);
               apply_BSH(world,tempfcwf,v_energies[i],DFparams.small,DFparams.thresh);
               Residuals[i] = virtuals[i] - tempfcwf;
               double residualnorm = Residuals[i].norm2();
               if(world.rank()==0) print("Virtual ", i, "     Residual: ", residualnorm);
               if(residualnorm > tol) keep_going = true;
          }    

          //Apply KAIN, if called for
          if(DFparams.kain){
               virtuals = kainsolver.update(virtuals, Residuals);
          }
          else{
               virtuals = virtuals - Residuals;
          }

          //Orthogonalize the functions against the occupieds
          for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
               for(unsigned int j = 0; j < Init_params.num_occupied; j++){
                    tempcomplex = inner(virtuals[i], occupieds[j]);
                    virtuals[i] -= occupieds[j]*tempcomplex;
               }
               virtuals[i].normalize();
          }

          //orthonormalize the functions internally
          orthogonalize(world,virtuals);

          //Recalculate exchange
          for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
               Kpsis[i] = apply_K(world, op, virtuals[i]);
          }

          //Recalculate and print energies for the output
          if(world.rank()==0){
               print("\n***Printing Current Energies***");
          }

          for(unsigned int j = 0; j < Init_params.num_virtuals; j++){
               v_energies[j] = rele(world,virtuals[j]);
               v_energies[j] += real(inner(virtuals[j],virtuals[j]*JandV));
               v_energies[j] -= real(inner(virtuals[j],Kpsis[j]));
          }
          
          for(unsigned int j = 0; j < Init_params.num_virtuals; j++){
               if(world.rank()==0){
                    printf("                Virtual: %3i, Energy: %.10e\n",j+1, v_energies[j]-csquared);
               }
          }

     }


     //Fix the location of the last occupied orbital
     occupieds.push_back(copy(virtuals[0]));
     for(unsigned int i = 0; i < Init_params.num_virtuals-1; i++){
          virtuals[i] = virtuals[i+1];
     }
     virtuals.pop_back();
     Tensor<double> tempenergies(Init_params.num_occupied + 1);
     tempenergies(Init_params.num_occupied) = v_energies[0];
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          tempenergies(i) = energies(i);
     }
     energies = tempenergies;
     v_energies = v_energies(Slice(1,Init_params.num_virtuals - 1));
     Init_params.num_occupied += 1;
     Init_params.num_virtuals -= 1;

     //Recalculate and print total energy for the output
     if(world.rank()==0){
          print("\n***Printing Current Energies***");
     }
     double kinetic_energy = 0.0;
     double coulomb_energy = 0.0;
     double exchange_energy = 0.0;
     double nuclear_attraction_energy = 0.0;
     double myc = 137.0359895; //speed of light

     rho = real_factory_3d(world);
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          rho += squaremod(occupieds[j]);
     }
     real_function_3d Jop = apply(op,rho);

     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          energies[j] = rele(world,occupieds[j]);
          kinetic_energy += (energies[j] - myc*myc);
          double nuclear_attraction_energy_correction = real(inner(occupieds[j],occupieds[j]*Vnuc));
          nuclear_attraction_energy += nuclear_attraction_energy_correction;
          double coulomb_energy_correction = real(inner(occupieds[j],occupieds[j]*Jop));
          energies[j] += coulomb_energy_correction + nuclear_attraction_energy_correction;
          coulomb_energy += 0.5*coulomb_energy_correction;
          double exchange_energy_correction = real(inner(occupieds[j],apply_K(world,op,occupieds[j])));
          energies[j] -= exchange_energy_correction;
          exchange_energy += 0.5*exchange_energy_correction;
     }
     total_energy = kinetic_energy + coulomb_energy - exchange_energy + nuclear_attraction_energy + nuclear_repulsion_energy;
     
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          if(world.rank()==0){
               printf("                Orbital: %3i, Energy: %.10e\n",j+1, energies[j]-myc*myc);
          }
     }
     if(world.rank()==0){
          print("\n              Kinetic Energy: ",kinetic_energy);
          print("             Coulomb  Energy: ",coulomb_energy);
          print("             Exchange Energy: ",exchange_energy);
          print("   Nuclear Attraction Energy: ",nuclear_attraction_energy);
          print("    Nuclear Repulsion Energy: ",nuclear_repulsion_energy);
          print("                Total Energy: ",total_energy);
     }

     //Below is where we can now solve for virtual orbitals



}

//kthxbye



