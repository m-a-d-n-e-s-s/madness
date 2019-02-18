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
          GaussianNucleusFunctor(Molecule& molecule){

               //get atom coordinates
               m_Rlist = molecule.get_all_coords_vec();
               
               //get atomic numbers
               for(unsigned int i = 0; i < m_Rlist.size(); i++){
                    m_Zlist.push_back(molecule.get_atom_number(i));
               }
               
               //find atomic mass numbers for each atom
               int Alist[116] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59,63,65,70,72,75,79,80,84,85,87,89,91,93,96,98,101,103,106,108,112,115,119,122,127,127,131,133,137,139,140,141,144,145,150,152,157,159,162,165,167,169,173,175,178,181,184,186,190,192,195,197,200,204,207,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,264,266,264,272,277,284,289,288,292};
               for(unsigned int i = 0; i < m_Zlist.size(); i++){
                    m_Alist.push_back(Alist[m_Zlist[i]-1]);
               }

               //calculate factors necessary for the potential
               for(unsigned int i = 0; i < m_Zlist.size(); i++){
                    m_xi.push_back(sqrt(3.0/2.0)/(0.836*pow(m_Alist[i],1.0/3.0)+0.57)*52917.7211);
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
          FermiNucDistFunctor(int& Z, coord_3d R){
               m_T = 0.000043463700858425666; //2.3 fm in bohr
               m_R.push_back(R);
               
               //find atomic mass numbers for each atom
               int Alist[116] = {1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59,63,65,70,72,75,79,80,84,85,87,89,91,93,96,98,101,103,106,108,112,115,119,122,127,127,131,133,137,139,140,141,144,145,150,152,157,159,162,165,167,169,173,175,178,181,184,186,190,192,195,197,200,204,207,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,264,266,264,272,277,284,289,288,292};
               m_A = Alist[Z-1];

               double PI = constants::pi;
               if(m_A < 5){
                    m_C = 0.000022291*pow(m_A, 1.0/3.0) - 0.0000090676;
               }
               else{
                    m_C = sqrt(5.0/3.0*pow((0.836*pow(m_A,1.0/3.0)+0.570)/52917.7211,2) - 7.0/3.0*pow(PI*m_T/4.0/log(3.0),2));
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
          FermiNucDistFunctor rho(Zlist[i], Rlist[i]);
          temp = real_factory_3d(world).functor(rho).truncate_mode(0);
          tempnorm = temp.trace();
          temp.scale(-Zlist[i]/tempnorm);
          if(i == 0){
               potential = temp;
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
     //FunctionDefaults<3>::set_cubic_cell(-Init_params.L, Init_params.L);
     //FunctionDefaults<3>::set_k(Init_params.order);
     FunctionDefaults<3>::set_thresh(DFparams.thresh); //Always use user-specified thresh
     FunctionDefaults<3>::set_truncate_mode(1);   

     //If user requests different k, then project functions
     if(DFparams.k != Init_params.order){
          FunctionDefaults<3>::set_k(DFparams.k);
          for(unsigned int i = 0; i < Init_params.num_occupied; i++){
               for(unsigned int j = 0; j < 4; j++){
                    //if(world.rank()==0) print(i,j);
                    Init_params.orbitals[i][j] = project(Init_params.orbitals[i][j], FunctionDefaults<3>::get_k(), DFparams.thresh, false);
               }
               world.gop.fence();
               Init_params.orbitals[i].truncate();
          }
          if(DFparams.job == 1){
               for(unsigned int i = 0; i < Init_params.num_virtuals; i++){
                    for(unsigned int j = 0; j < 4; j++){
                         //if(world.rank()==0) print(i,j);
                         Init_params.virtuals[i][j] = project(Init_params.virtuals[i][j], FunctionDefaults<3>::get_k(), DFparams.thresh, false);
                    }
                    world.gop.fence();
                    Init_params.virtuals[i].truncate();
               }
          }
     }

     energies = Init_params.energies;
     v_energies = Init_params.v_energies;
     occupieds = Init_params.orbitals;
     virtuals = Init_params.virtuals;
     total_energy = Init_params.Init_total_energy;

     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("Preparation complete: ", times[0]);

}

//function to calculate the kinetic + rest energy expectation value using Dirac Hamiltonian c*\alpha*p+\Beta*m*c*c
double DF::rele(World& world, const Fcwf& psi){
     std::complex<double> energy = 0.0;
     double myc = 137.0359895; //speed of light in atomic units
     std::complex<double> myi(0,1);
     complex_derivative_3d Dx(world,0);
     complex_derivative_3d Dy(world,1);
     complex_derivative_3d Dz(world,2);
     Fcwf Tpsi(world);
     complex_function_3d temp(world);
     
     Tpsi[0] = Dz(psi[2]) + Dx(psi[3]) - myi*Dy(psi[3]) + myc*myi*psi[0];
     Tpsi[1] = Dx(psi[2]) + myi*Dy(psi[2]) - Dz(psi[3]) + myc*myi*psi[1];
     Tpsi[2] = Dz(psi[0]) + Dx(psi[1]) - myi*Dy(psi[1]) - myc*myi*psi[2];
     Tpsi[3] = Dx(psi[0]) + myi*Dy(psi[0]) - Dz(psi[1]) - myc*myi*psi[3];

     //The giant block below this exists because I was debugging
     //some hung queues. I still don't know what's causing them
     //but they occur in this function somewhere, somehow
     
     //if(world.rank()==0) print("a");
     //temp = Dz(psi[2]);
     //world.gop.fence();
     //if(world.rank()==0) print("b");
     //temp += Dx(psi[3]);
     //world.gop.fence();
     //if(world.rank()==0) print("c");
     //temp -= myi*Dy(psi[3]);
     //world.gop.fence();
     //if(world.rank()==0) print("d");
     //temp += myc*myi*psi[0];
     //world.gop.fence();
     //if(world.rank()==0) print("e");
     //Tpsi[0] = copy(temp);
     //world.gop.fence();

     //if(world.rank()==0) print("f");
     //temp = Dx(psi[2]);
     //world.gop.fence();
     //if(world.rank()==0) print("g");
     //temp += myi*Dy(psi[2]);
     //world.gop.fence();
     //if(world.rank()==0) print("h");
     //temp -= Dz(psi[3]);
     //world.gop.fence();
     //if(world.rank()==0) print("i");
     //temp += myc*myi*psi[1];
     //world.gop.fence();
     //if(world.rank()==0) print("j");
     //Tpsi[1] = copy(temp);
     //world.gop.fence();
     //
     //if(world.rank()==0) print("k");
     //temp = Dz(psi[0]);
     //world.gop.fence();
     //if(world.rank()==0) print("l");
     //temp += Dx(psi[1]);
     //world.gop.fence();
     //if(world.rank()==0) print("m");
     //temp -= myi*Dy(psi[1]);
     //world.gop.fence();
     //if(world.rank()==0) print("n");
     //temp -= myc*myi*psi[2];
     //world.gop.fence();
     //if(world.rank()==0) print("o");
     //Tpsi[2] = copy(temp);
     //world.gop.fence();
     //
     //if(world.rank()==0) print("p");
     //temp = Dx(psi[0]);
     //world.gop.fence();
     //if(world.rank()==0) print("q");
     //temp += myi*Dy(psi[0]);
     //world.gop.fence();
     //if(world.rank()==0) print("r");
     //temp -= Dz(psi[1]);
     //world.gop.fence();
     //if(world.rank()==0) print("s");
     //temp -= myc*myi*psi[3];
     //world.gop.fence();
     //if(world.rank()==0) print("t");
     //Tpsi[3] = copy(temp);
     //world.gop.fence();
     //if(world.rank()==0) print("u");

     energy = inner(psi, Tpsi);
     energy *= -1.0*myi*myc;

     return energy.real();
}

//returns a new Fcwf that is the result of applying the Dirac free-particle hamiltonian on psi
Fcwf apply_T(World& world, const Fcwf& psi){
     double myc = 137.0359895; //speed of light in atomic units
     std::complex<double> myi(0,1);
     complex_derivative_3d Dx(world,0);
     complex_derivative_3d Dy(world,1);
     complex_derivative_3d Dz(world,2);
     Fcwf Tpsi(world);

     Tpsi[0] = Dz(psi[2]) + Dx(psi[3]) - myi*Dy(psi[3]) + myc*myi*psi[0];
     Tpsi[1] = Dx(psi[2]) + myi*Dy(psi[2]) - Dz(psi[3]) + myc*myi*psi[1];
     Tpsi[2] = Dz(psi[0]) + Dx(psi[1]) - myi*Dy(psi[1]) - myc*myi*psi[2];
     Tpsi[3] = Dx(psi[0]) + myi*Dy(psi[0]) - Dz(psi[1]) - myc*myi*psi[3];

     return Tpsi * (-myi*myc);
}

////Calculates K*psi for each psi in the orbitals vector and stores them in result
//void DF::exchange(World& world, real_convolution_3d& op, std::vector<Fcwf>& Kpsis){
//     //start timer
//     start_timer(world);
//
//     complex_function_3d temp(world);
//
//     unsigned int n = Init_params.num_occupied;
//     for(unsigned int i = 0; i < n; i++){
//          for(unsigned int j = 0; j < n ; j++){
//               temp = inner_func(world,occupieds[j],occupieds[i]);
//               temp.truncate();
//
//               temp = apply(op,temp);
//               if(j == 0){
//                    Kpsis[i] = occupieds[j]*temp;
//               }
//               else{
//                    Kpsis[i] += occupieds[j]*temp;
//               }
//          }
//
//          //here too
//          Kpsis[i].truncate();
//     }
//
//     //Report time
//     Tensor<double> times = end_timer(world);
//     if(world.rank()==0) print("     ", times[0]);
//}
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

     //form occupied orbitals into 4 vectors of components
     unsigned int n = Init_params.num_occupied;
     //std::vector<complex_function_3d> occupieds0(n);
     //std::vector<complex_function_3d> occupieds1(n);
     //std::vector<complex_function_3d> occupieds2(n);
     //std::vector<complex_function_3d> occupieds3(n);
     //for(unsigned int i = 0; i < n; i++){
     //     occupieds0[i] = occupieds[i][0];
     //     occupieds1[i] = occupieds[i][1];
     //     occupieds2[i] = occupieds[i][2];
     //     occupieds3[i] = occupieds[i][3];
     //}

     //Calculate and accumulate exchange contributions
     for(unsigned int i = 0; i < n; i++){

          std::vector<complex_function_3d> temp(n-i);
          for(unsigned int j = 0; j < n-i; j++){
               temp[j] = complex_factory_3d(world);    
          }

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

          gaxpy(world, 1.0, temp, 1.0, occupieds[i][0]*conj(world, temp0));
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][1]*conj(world, temp1));
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][2]*conj(world, temp2));
          gaxpy(world, 1.0, temp, 1.0, occupieds[i][3]*conj(world, temp3));

          temp = apply(world, op, temp);

          Kpsis[i][0] += sum(world, mul(world, temp, temp0));
          Kpsis[i][1] += sum(world, mul(world, temp, temp1));
          Kpsis[i][2] += sum(world, mul(world, temp, temp2));
          Kpsis[i][3] += sum(world, mul(world, temp, temp3));
          
          temp = conj(world, temp);

          temp0 = occupieds[i][0]*temp;
          temp1 = occupieds[i][1]*temp;
          temp2 = occupieds[i][2]*temp;
          temp3 = occupieds[i][3]*temp;
          for(unsigned int j = i+1; j < n; j++){
               Kpsis[j][0] += temp0[j-i];
               Kpsis[j][1] += temp1[j-i];
               Kpsis[j][2] += temp2[j-i];
               Kpsis[j][3] += temp3[j-i];
          }
          





          //for(unsigned int j = i; j < n ; j++){ //parallelize this loop using 4-vectors approach

          //     /*TODO (vector inner_func)
          //      * loop over 4 component indices
          //      *   gather the n-i+1 (here, index j) functions for the component index
          //      *   vector multiply by occupieds[i] component
          //      *   accumulate result
          //      * truncate
          //      */


          //     //load balance like in SCF.cc KE routine


          //     temp = inner_func(world,occupieds[j],occupieds[i]);
          //     temp.truncate();

          //     temp = apply(op,temp);
          //     if(i == 0 && j == 0){
          //          Kpsis[i] = occupieds[j]*temp;
          //     }
          //     else if(i == 0){
          //          Kpsis[i] += occupieds[j]*temp;
          //          temp = temp.conj();
          //          Kpsis[j] = occupieds[i]*temp;
          //     }
          //     else if(i == j){
          //          Kpsis[i] += occupieds[j]*temp;
          //     }
          //     else{
          //          Kpsis[i] += occupieds[j]*temp;
          //          temp = temp.conj();
          //          Kpsis[j] += occupieds[i]*temp;
          //     }
          //}

          Kpsis[i].truncate();
     }

     //Report time
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);
}

Fcwf DF::apply_K(World& world, real_convolution_3d& op, Fcwf& phi){
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
     //first apply Fock operator to all orbitals
     //
     
     //calculate coulomb part
     if(world.rank() == 0) print("          Adding (V+J)psi");
     real_function_3d rho = real_factory_3d(world);
     for(unsigned int j = 0; j < n; j++){
          rho += squaremod(occupieds[j]);
     }
     
     //----------------------------------------------------------------------------------------------------------
     //DEBUGGING: Print each matrix that goes into the Fock matrix

     //real_function_3d fred = apply(op,rho);
     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals.push_back(occupieds[j]*fred); //add in coulomb term
     //}
     //fock = matrix_inner(world, occupieds, temp_orbitals);
     //if(world.rank()==0) print("J:\n",fock);

     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals[j] = occupieds[j]*myV;
     //}
     //fock = matrix_inner(world, occupieds, temp_orbitals);
     //if(world.rank()==0) print("V:\n",fock);

     //
     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals[j] = Kpsis[j]*(-1.0);
     //}
     //fock = matrix_inner(world, occupieds, temp_orbitals);
     //if(world.rank()==0) print("-K:\n",fock);

     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals[j] = apply_T(world, occupieds[j]);
     //}
     //fock = matrix_inner(world, occupieds, temp_orbitals);
     //if(world.rank()==0) print("T:\n",fock);
     //END DEBUGGING
     //----------------------------------------------------------------------------------------------------------

     //TODO: Here try moving the apply out to operate on the sum of the nuclear and electronic charge distributions
     real_function_3d potential = myV + apply(op,rho);
     potential.truncate();

     //add in coulomb parts to neworbitals
     for(unsigned int j = 0; j < n; j++){
          temp_orbitals.push_back(occupieds[j]*potential); //add in coulomb term
          //temp_orbitals[j] = occupieds[j]*potential;
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

     //We now have a vector of F*psi for each psi (in neworbitals).
     //Calculate each element of the fock matrix
     //
     //New method: Make a matrix for each component and add the matrices. This facilitates the
     //use of the new inner product
     //
     //1) Create four std::vectors each from occupieds and temp_orbitals, one for each component
     //2) Call the new matrix inner for each component. Result in each case gives a matrix
     //3) Sum all the matrices
     //std::vector<complex_function_3d> occupieds_1;
     //std::vector<complex_function_3d> occupieds_2;
     //std::vector<complex_function_3d> occupieds_3;
     //std::vector<complex_function_3d> occupieds_4;
     //std::vector<complex_function_3d> temp_orbitals_1;
     //std::vector<complex_function_3d> temp_orbitals_2;
     //std::vector<complex_function_3d> temp_orbitals_3;
     //std::vector<complex_function_3d> temp_orbitals_4;
     //for(unsigned int i = 0; i < n; i++){
     //     occupieds_1.push_back(occupieds[i][0]);
     //     occupieds_2.push_back(occupieds[i][1]);
     //     occupieds_3.push_back(occupieds[i][2]);
     //     occupieds_4.push_back(occupieds[i][3]);
     //     temp_orbitals_1.push_back(temp_orbitals[i][0]);
     //     temp_orbitals_2.push_back(temp_orbitals[i][1]);
     //     temp_orbitals_3.push_back(temp_orbitals[i][2]);
     //     temp_orbitals_4.push_back(temp_orbitals[i][3]);
     //}
     //Tensor<std::complex<double>> component1 = matrix_inner(world,occupieds_1,temp_orbitals_1);
     //Tensor<std::complex<double>> component2 = matrix_inner(world,occupieds_2,temp_orbitals_2);
     //Tensor<std::complex<double>> component3 = matrix_inner(world,occupieds_3,temp_orbitals_3);
     //Tensor<std::complex<double>> component4 = matrix_inner(world,occupieds_4,temp_orbitals_4);

     //fock = component1+component2+component3+component4;
     fock = matrix_inner(world, occupieds, temp_orbitals);
     
     //Old Method here:
     //for(unsigned int j = 0; j < n; j++){
     //     for(unsigned int k = 0; k < n; k++){
     //          fock(j,k) = inner(occupieds[j],temp_orbitals[k]); 
     //     }
     //}

     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("               ", times[0]); 


     ////Form the overlap matrix
     if(world.rank()==0) print("          Integrating to form Overlap Matrix");
     start_timer(world);
     //for(unsigned int j = 0; j < n; j++){
     //     for(unsigned int k = 0; k < n; k++){
     //          overlap(j,k) = inner(occupieds[j],occupieds[k]);
     //     }
     //}
     //component1 = matrix_inner(world,occupieds_1,occupieds_1);
     //component2 = matrix_inner(world,occupieds_2,occupieds_2);
     //component3 = matrix_inner(world,occupieds_3,occupieds_3);
     //component4 = matrix_inner(world,occupieds_4,occupieds_4);
     //overlap = component1 + component2 + component3 + component4;
     overlap = matrix_inner(world,occupieds,occupieds);
     
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
     //if(world.rank()==0) {
     //     print("evals:");
     //     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
     //          print(evals[j] - 137.0359895*137.0359895);    
     //     }
     //}

     //Before applying the transformation, fix arbitrary rotations introduced by the eigensolver. 
     if(world.rank()==0) print("     Removing Rotations");
     start_timer(world);

     double thresh_degenerate = DFparams.thresh*100.0;
     double csquared = 137.0359895*137.0359895; //electron rest energy
     //swap columns for a diagonally dominant matrix
     bool switched = true;
     while (switched) {
          switched = false;
          for (unsigned int kk = 0; kk < Init_params.num_occupied; kk++) {
               for (unsigned int j = kk + 1; j < Init_params.num_occupied; j++) {
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
     for (unsigned int kk = 0; kk < Init_params.num_occupied; ++kk)
          U(_, kk).scale(std::conj(U(kk,kk))/std::abs(U(kk,kk)));

     
     //Find clusters of degenerate eigenvalues and rotate eigenvectors to maximize overlap with previous ones
     unsigned int ilo = 0; // first element of cluster
     if(world.rank()==0) print("          Degeneracy threshold: ",thresh_degenerate);
     
     while (ilo < Init_params.num_occupied - 1) {
         unsigned int ihi = ilo;
         while (fabs(evals[ilo] - evals[ihi + 1])
                < thresh_degenerate * std::max(std::fabs(evals[ilo]-csquared),1.0)){// pow(10,floor(log10(std::fabs(evals[ilo]))))) { 
             ++ihi;
             if (ihi == Init_params.num_occupied - 1)
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

     //Debugging: Print transformation matrix after rotation removal
     //if(world.rank()==0) print("U:\n", U);
     
     
     fock = inner(conj_transpose(U), inner(fock, U));
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);
     //if(world.rank()==0)print("\n U^* F U:\n",fock,"\n");

     if(world.rank()==0) print("     Applying Transformation");
     start_timer(world);

     ////Apply the transformation to the Exchange
     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals[j] = Fcwf(world);
     //     for(unsigned int k = 0; k < n; k++){
     //          temp_orbitals[j] += Kpsis[k]*U(k,j);
     //     }
     //}
     //for(unsigned int m = 0; m < n; m++){
     //      Kpsis[m] = temp_orbitals[m];
     //}
     
     //Trying a different method of applying the transformation, using vmra's transform
     //for(unsigned int i = 0; i < n; i++){
     //     temp_orbitals_1[i] = Kpsis[i][0];
     //     temp_orbitals_2[i] = Kpsis[i][1];
     //     temp_orbitals_3[i] = Kpsis[i][2];
     //     temp_orbitals_4[i] = Kpsis[i][3];
     //}
     //temp_orbitals_1 = transform(world, temp_orbitals_1, U);
     //temp_orbitals_2 = transform(world, temp_orbitals_2, U);
     //temp_orbitals_3 = transform(world, temp_orbitals_3, U);
     //temp_orbitals_4 = transform(world, temp_orbitals_4, U);
     //for(unsigned int i = 0; i < n; i++){
     //     Kpsis[i][0] = temp_orbitals_1[i];
     //     Kpsis[i][1] = temp_orbitals_2[i];
     //     Kpsis[i][2] = temp_orbitals_3[i];
     //     Kpsis[i][3] = temp_orbitals_4[i];
     //}

     transform(world, Kpsis, U);

     ////Apply the transformation to the orbitals
     //for(unsigned int j = 0; j < n; j++){
     //     temp_orbitals[j] = Fcwf(world);
     //     for(unsigned int k = 0; k < n; k++){
     //          temp_orbitals[j] += occupieds[k]*U(k,j); 
     //     }
     //}
     //for(unsigned int j = 0; j < n; j++){
     //     occupieds[j] = temp_orbitals[j];
     //}
     
     //Trying a different method of applying the transformation, using vmra's transform
     //for(unsigned int i = 0; i < n; i++){
     //     temp_orbitals_1[i] = occupieds[i][0];
     //     temp_orbitals_2[i] = occupieds[i][1];
     //     temp_orbitals_3[i] = occupieds[i][2];
     //     temp_orbitals_4[i] = occupieds[i][3];
     //}
     //temp_orbitals_1 = transform(world, temp_orbitals_1, U);
     //temp_orbitals_2 = transform(world, temp_orbitals_2, U);
     //temp_orbitals_3 = transform(world, temp_orbitals_3, U);
     //temp_orbitals_4 = transform(world, temp_orbitals_4, U);
     //for(unsigned int i = 0; i < n; i++){
     //     occupieds[i][0] = temp_orbitals_1[i];
     //     occupieds[i][1] = temp_orbitals_2[i];
     //     occupieds[i][2] = temp_orbitals_3[i];
     //     occupieds[i][3] = temp_orbitals_4[i];
     //}
     
     transform(world, occupieds, U);
     
     
     times = end_timer(world);
     if(world.rank()==0) print("          ", times[0]);
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Debugging? Let's try setting energies to evals like we're supposed to.
     //energies = evals;

}

//returns a vector of orthonormal Fcwfs constructed from the input vector of Fcwfs
std::vector<Fcwf> orthogonalize(World& world, std::vector<Fcwf> orbitals){
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
//Make this a member function of DF
//TODO: The function below mimics one from SCF.cc. In the future we will probably want a different variant (which also should exist in SCF.cc or somwhere like that) that treats core and valence orbitals differently: i.e. doesn't mix valence orbitals into the core orbitals
void DF::orthogonalize_inplace(World& world){

     unsigned int n = occupieds.size();
     double maxq;

     //Debugging: original overlap matrix for comparison
     //Tensor<std::complex<double>> overlap = matrix_inner(world, occupieds, occupieds);
     //if(world.rank()==0) print("\noriginal overlap:\n", overlap);

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

     //Debugging: print new overlap matrix
     //overlap = matrix_inner(world, occupieds, occupieds);
     //if(world.rank()==0) print("\nnew overlap:\n", overlap);

}

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
//double iterate(World& world, complex_function_3d& V, Fcwf& psi, double& eps){
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
     //TODO: make my own GradBSH that first computes BSH with a ridiculous lo, and then accumulates to CDelta, then need something intelligent for the derivative of that result
     world.gop.fence();
     //double ttt = wall_time();
     //create BSH operator
     if(world.rank()==0) print("mu: ", mu);
     std::shared_ptr<real_convolution_3d> op = std::shared_ptr<real_convolution_3d>(BSHOperatorPtr3D(world, mu,small,thresh)); // Finer length scale and accuracy control
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator(world, mu, small, thresh); 
     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator_Joel(world, mu, 1e-8, thresh); 
     std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> allops(16);
     for(unsigned int i = 0; i < 4; i++){
          allops[i] = op;
          allops[4+i] = op3[0];
          allops[8+i] = op3[1];
          allops[12+i] = op3[2];
     }

    
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("              create operators:", ttt);


     //create intermediate functions necessary to compute new components
     //ttt = wall_time();
     std::vector<complex_function_3d> temp(16);
     for(unsigned int i = 0; i < 4; i++){
          temp[i] = Vpsi[i];
          temp[4+i] = Vpsi[i];
          temp[8+i] = Vpsi[i];
          temp[12+i] = Vpsi[i];
     }

     //Fcwf oppsi = apply(world, op, Vpsi);             //BSH(Vpsi)
     //Fcwf oppsix = apply(world, *op3[0], Vpsi);       //GradBSH_x(Vpsi)
     //Fcwf oppsiy = apply(world, *op3[1], Vpsi);       //GradBSH_y(Vpsi)
     //Fcwf oppsiz = apply(world, *op3[2], Vpsi);       //GradBSH_z(Vpsi)
     temp = apply(world, allops, temp);
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("            create derivatives:", ttt);
     //manually compute new components of wavefunction. See Blackledge paper and any definition of the Dirac single-particle Hamiltonian
     //ttt = wall_time();



     //build transformation tensor
     Tensor<std::complex<double>> U(16,4);
     U(0,0) = c2+eps; U(14,0) = -ic; U(7,0) = -ic; U(11,0) = -myc;
     U(1,1) = c2+eps; U(6,1) = -ic; U(10,1) = myc; U(15,1) = ic;
     U(12,2) = -ic; U(5,2) = -ic; U(9,2) = -myc; U(2,2) = eps-c2;
     U(4,3) = -ic; U(8,3) = myc; U(13,3) = ic; U(3,3) = eps-c2;
     U *= (1.0/c2);

     ////build vector to be transformed
     //std::vector<complex_function_3d> temp(16);
     //temp[0] = oppsi[0]; temp[1] = oppsi[1]; temp[2] = oppsi[2]; temp[3] = oppsi[3];
     //temp[4] = oppsix[0]; temp[5] = oppsix[1]; temp[6] = oppsix[2]; temp[7] = oppsix[3];
     //temp[8] = oppsiy[0]; temp[9] = oppsiy[1]; temp[10] = oppsiy[2]; temp[11] = oppsiy[3];
     //temp[12] = oppsiz[0]; temp[13] = oppsiz[1]; temp[14] = oppsiz[2]; temp[15] = oppsiz[3];

     temp = transform(world, temp, U);
     Vpsi[0] = temp[0];
     Vpsi[1] = temp[1];
     Vpsi[2] = temp[2];
     Vpsi[3] = temp[3];

     //Vpsi[0] = myc*myc*oppsi[0]-myi*myc*oppsiz[2]-myi*myc*(oppsix[3]-myi*oppsiy[3]);
     //Vpsi[1] = myc*myc*oppsi[1]-myi*myc*(oppsix[2]+myi*oppsiy[2])+myi*myc*oppsiz[3];
     //Vpsi[2] = -myi*myc*oppsiz[0]-myi*myc*(oppsix[1]-myi*oppsiy[1])-myc*myc*oppsi[2];
     //Vpsi[3] = -myi*myc*(oppsix[0]+myi*oppsiy[0])+myi*myc*oppsiz[1]-myc*myc*oppsi[3];
     ////finish up application of the dirac green's function
     //Vpsi = (Vpsi+(oppsi*eps))*(1.0/(myc*myc));




     //Vpsi.normalize(); //Don't want to do this because it messes with the equation KAIN is trying to solve.
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("          apply transformation:", ttt);

}

//Apply's Green's function to Vpsi (a Fcwf). Overwrites Vpsi with new Fcwf
//double iterate(World& world, complex_function_3d& V, Fcwf& psi, double& eps){
void apply_BSH_new(World& world, Fcwf& Vpsi, double& eps, double& small, double& thresh){


     //necessary constants
     double myc = 137.0359895; //speed of light
     double c2 = myc*myc;
     std::complex<double> myi(0,1); //imaginary number
     std::complex<double> ic = myi*myc;
    
     //calculate exponent for equivalent BSH operator
     double mu = std::sqrt((myc*myc*myc*myc-eps*eps)/myc/myc);

     //if(world.rank() == 0) print("Hi, this is apply_BSH! mu is: ", mu);

     //create gradient BSH operators
     //TODO: make my own GradBSH that first computes BSH with a ridiculous lo, and then accumulates to CDelta, then need something intelligent for the derivative of that result
     world.gop.fence();
     //double ttt = wall_time();
     //create BSH operator
     //if(world.rank()==0) print("mu: ", mu);

     real_convolution_3d op = BSHOperator3D(world, mu,small,thresh); // Finer length scale and accuracy control

     Vpsi = apply(world, op, Vpsi);
     Vpsi = (apply_T(world, Vpsi) + Vpsi * eps) * (1.0/c2);

     //Vpsi = apply(world, op, apply_T(world, Vpsi) + Vpsi * eps) * (1.0/c2);

     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator(world, mu, 1e-8, thresh); 
     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> op3 = GradBSHOperator_Joel(world, mu, 1e-8, thresh); 
     //std::vector<std::shared_ptr<SeparatedConvolution<double,3>>> allops(16);
     //for(unsigned int i = 0; i < 4; i++){
     //     allops[i] = op;
     //     allops[4+i] = op3[0];
     //     allops[8+i] = op3[1];
     //     allops[12+i] = op3[2];
     //}

    
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("              create operators:", ttt);


     //create intermediate functions necessary to compute new components
     //ttt = wall_time();
     //std::vector<complex_function_3d> temp(16);
     //for(unsigned int i = 0; i < 4; i++){
     //     temp[i] = Vpsi[i];
     //     temp[4+i] = Vpsi[i];
     //     temp[8+i] = Vpsi[i];
     //     temp[12+i] = Vpsi[i];
     //}

     //Fcwf oppsi = apply(world, op, Vpsi);             //BSH(Vpsi)
     //Fcwf oppsix = apply(world, *op3[0], Vpsi);       //GradBSH_x(Vpsi)
     //Fcwf oppsiy = apply(world, *op3[1], Vpsi);       //GradBSH_y(Vpsi)
     //Fcwf oppsiz = apply(world, *op3[2], Vpsi);       //GradBSH_z(Vpsi)
     //temp = apply(world, allops, temp);
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("            create derivatives:", ttt);
     //manually compute new components of wavefunction. See Blackledge paper and any definition of the Dirac single-particle Hamiltonian
     //ttt = wall_time();



     //build transformation tensor
     //Tensor<std::complex<double>> U(16,4);
     //U(0,0) = c2+eps; U(14,0) = -ic; U(7,0) = -ic; U(11,0) = -myc;
     //U(1,1) = c2+eps; U(6,1) = -ic; U(10,1) = myc; U(15,1) = ic;
     //U(12,2) = -ic; U(5,2) = -ic; U(9,2) = -myc; U(2,2) = eps-c2;
     //U(4,3) = -ic; U(8,3) = myc; U(13,3) = ic; U(3,3) = eps-c2;
     //U *= (1.0/c2);

     ////build vector to be transformed
     //std::vector<complex_function_3d> temp(16);
     //temp[0] = oppsi[0]; temp[1] = oppsi[1]; temp[2] = oppsi[2]; temp[3] = oppsi[3];
     //temp[4] = oppsix[0]; temp[5] = oppsix[1]; temp[6] = oppsix[2]; temp[7] = oppsix[3];
     //temp[8] = oppsiy[0]; temp[9] = oppsiy[1]; temp[10] = oppsiy[2]; temp[11] = oppsiy[3];
     //temp[12] = oppsiz[0]; temp[13] = oppsiz[1]; temp[14] = oppsiz[2]; temp[15] = oppsiz[3];

     //temp = transform(world, temp, U);
     //Vpsi[0] = temp[0];
     //Vpsi[1] = temp[1];
     //Vpsi[2] = temp[2];
     //Vpsi[3] = temp[3];

     //Vpsi[0] = myc*myc*oppsi[0]-myi*myc*oppsiz[2]-myi*myc*(oppsix[3]-myi*oppsiy[3]);
     //Vpsi[1] = myc*myc*oppsi[1]-myi*myc*(oppsix[2]+myi*oppsiy[2])+myi*myc*oppsiz[3];
     //Vpsi[2] = -myi*myc*oppsiz[0]-myi*myc*(oppsix[1]-myi*oppsiy[1])-myc*myc*oppsi[2];
     //Vpsi[3] = -myi*myc*(oppsix[0]+myi*oppsiy[0])+myi*myc*oppsiz[1]-myc*myc*oppsi[3];
     ////finish up application of the dirac green's function
     //Vpsi = (Vpsi+(oppsi*eps))*(1.0/(myc*myc));




     //Vpsi.normalize(); //Don't want to do this because it messes with the equation KAIN is trying to solve.
     //world.gop.fence();
     //ttt = wall_time() - ttt;
     //if(world.rank()==0) print("          apply transformation:", ttt);

}

// Prints norms of the given vector of functions
void DF::print_norms(World & world,
                     std::vector<Fcwf> f)
{
   // Container
   Tensor<double> norms(f.size());

   // Calc the norms
   for(unsigned int i = 0; i < f.size(); i++){
       norms(i) = f[i].norm2();
   }

   // Print em in a smart way
   if(world.rank() == 0) print(norms);
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

// Generates initial guess functions from nonrelativistic orbitals. Note that
// The length of the returned vector is twice that of the input vector
std::vector<Fcwf> to_fcwfs(World & world, std::vector<real_function_3d> & orbitals){
     // Get size
     unsigned int n = orbitals.size();

     complex_function_3d complexreader;
     std::vector<Fcwf> result;
     Fcwf spinup(world);
     Fcwf spindown(world);

     //Loop over input orbitals
     for(unsigned int i = 0; i < n; i++){
          complexreader = function_real2complex(orbitals[i]);
          spinup = Fcwf(copy(complexreader), complex_factory_3d(world), complex_factory_3d(world), complex_factory_3d(world));
          spindown = Fcwf(complex_factory_3d(world), copy(complexreader), complex_factory_3d(world), complex_factory_3d(world));
          result.push_back(spinup);
          result.push_back(spindown);
          
     }

     return result;

}

//Duplicates the energies in a vector of energies and adds the electron rest energy to each
Tensor<double> make_initial_energies(Tensor<double> energies){
     int n = energies.size();
     Tensor<double> result(2*n);
     double csquared = 137.0359895*137.0359895; //electron rest energy
     double temp;
     for(int i = 0; i < n; i++){
          temp = energies(i) + csquared;
          result(2*i) = temp;
          result(2*i+1) = temp; //energies get repeated because we double the number of orbitals
     }
     return result;
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
     GaussianNucleusFunctor Vfunctor(Init_params.molecule);
     potential = real_factory_3d(world).functor(Vfunctor).truncate_mode(0).truncate_on_project();
}

//Creates the nuclear potential from the molecule object. Also calculates the nuclear repulsion energy
void DF::make_gaussian_potential(World& world, real_function_3d& potential, double& nuclear_repulsion_energy){
     if(world.rank()==0) print("\n***Making a Gaussian Potential***");
     GaussianNucleusFunctor Vfunctor(Init_params.molecule);
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
     //for(unsigned int j = 0; j < Init_params.num_occupied; j++){
     //     for(int kk = 0; kk < 4; kk++){
     //          lb.add_tree(occupieds[j][kk], lbcost<std::complex<double>,3>(24.0,192.0),true);
     //     }
     //}
     FunctionDefaults<3>::redistribute(world, lb.load_balance(2), false);
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

//Heavily editing this to debug
void DF::make_density_lineplots(World& world, const char* filename, int npt, double endpnt){
     std::vector<real_function_3d> densities;
     real_function_3d zero(world);
     for(unsigned int i=0; i < Init_params.num_occupied; i++){
          densities.push_back(squaremod(occupieds[i]));
          //densities.push_back(copy(zero));
     }

     //double h = Init_params.L*(1.0/(npt-1));
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

     //Diagonalize, but not on the first iteration, unless it's a restarted job
     //if(iteration_number != 1 or DFparams.restart){
     //     diagonalize(world, V, op, Kpsis);
     //     //update JandV
     //     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
     //          rho += squaremod(occupieds[kk]);
     //     }
     //     JandV = V + apply(op,rho); 
     //}

     //Actually let's see what happens if I always diagonalize. Answer: still messes up first iteration
     //Returning to this idea since it turns out I wasn't updating JandV...
     diagonalize(world, V, op, Kpsis);
     for(unsigned int kk = 0; kk < Init_params.num_occupied; kk++){
          rho += squaremod(occupieds[kk]);
     }
     JandV = V + apply(op,rho); 



     //Apply BSH to each psi
     if(world.rank()==0) print("\n***Applying BSH operator***");
     start_timer(world);
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){

          //Debugging
          //std::vector<double> temporbitalnorms(4);
          //if(world.rank()==0) print("\nbefore:");
          //for(unsigned int l = 0; l < 4; l++) temporbitalnorms[l] = occupieds[j][l].norm2();
          //if(world.rank()==0) {
          //     print("          ", temporbitalnorms[0]);
          //     print("          ", temporbitalnorms[1]);
          //     print("          ", temporbitalnorms[2]);
          //     print("          ", temporbitalnorms[3]);
          //}

          //construct the function to which we will apply the BSH
          occupieds[j].truncate();
          temp_function = occupieds[j]*JandV;
          temp_function.scale(-1.0);
          temp_function += Kpsis[j];
          temp_function.truncate();

          //temp_function now holds (K-V-J)psi, so apply the BSH
          //apply_BSH(world,  temp_function, energies[j], DFparams.small, DFparams.thresh);
          apply_BSH_new(world,  temp_function, energies[j], DFparams.small, DFparams.thresh);
          

          //if(world.rank()==0) print("after:");
          //for(unsigned int l = 0; l < 4; l++) temporbitalnorms[l] = temp_function[l].norm2();
          //if(world.rank()==0) {
          //     print("          ", temporbitalnorms[0]);
          //     print("          ", temporbitalnorms[1]);
          //     print("          ", temporbitalnorms[2]);
          //     print("          ", temporbitalnorms[3]);
          //     print("\n");
          //}

          //debugging: Look at size of function before and after truncating here to see what it's doing
          //std::size_t funcsize = temp_function[0].size();
          //funcsize += temp_function[1].size();
          //funcsize += temp_function[2].size();
          //funcsize += temp_function[3].size();
          //if(world.rank()==0) print("               size after BSH:", funcsize);

          temp_function.truncate(); //try truncating here

          //funcsize = temp_function[0].size();
          //funcsize += temp_function[1].size();
          //funcsize += temp_function[2].size();
          //funcsize += temp_function[3].size();
          //if(world.rank()==0) print("          size after truncate:", funcsize);
          
          //Now calculate the residual
          //temp_function = occupieds[j] - temp_function; 
          residualnorm = (occupieds[j] - temp_function).norm2();

          //Print residual norm to keep track
          if(world.rank()==0) printf("     Orbital: %3i,  Resid: %.10e\n",j+1, residualnorm);

          //If the norm is big enough, we'll need to iterate again.
          if(residualnorm > tolerance) iterate_again = true;

          //Store residual function if we're using KAIN. Not necessary if we're not using kain
          if(iteration_number != 1 and DFparams.kain){
               Residuals.push_back(occupieds[j] - temp_function);
          }
          else{
               occupieds[j] = temp_function;
          }
     }
     if(world.rank()==0) printf("                tolerance: %.10e\n",tolerance);
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Apply the kain solver, if called for
     if(iteration_number != 1 and DFparams.kain){
          if(world.rank()==0) print("\n***Applying KAIN Solver***");
          start_timer(world);

          //occupieds = kainsolver.update(occupieds, Residuals);

          //Replace above line with kain + step restriction
          Residuals = kainsolver.update(occupieds, Residuals); //Using Residuals for new orbitals to save storage
          for(unsigned int i=0; i < Init_params.num_occupied; i++){
               residualnorm = (occupieds[i]-Residuals[i]).norm2();
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

     //truncate here
     for(unsigned int i = 0; i < Init_params.num_occupied; i++) occupieds[i].truncate();

     //orthogonalize
     if(world.rank()==0) print("\n***Orthonormalizing***");
     start_timer(world);
     //occupieds = orthogonalize(world,occupieds);
     orthogonalize_inplace(world);
     //truncate here and normalize again
     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
           occupieds[i].truncate();
           occupieds[i].normalize();
     }
     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //calculate new exchange
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

     //if(world.rank()==0) print("     Recalculating Energies");
     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          
          //if(world.rank()==0) print("     Adding terms for orbital ", j);
          energies[j] = rele(world,occupieds[j]);
          
          kinetic_energy += (energies[j] - myc*myc);
          //if(world.rank()==0) print("          Kinetic done.");
     }
          
     //temp_function = occupieds[j]*V;
     //nuclear_attraction_energy_correction = real(inner(occupieds[j],temp_function));
     //nuclear_attraction_energy += nuclear_attraction_energy_correction;
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
     nuclear_attraction_tensor += real(inner(world,occupieds3,mul(world,V,occupieds3)));
     nuclear_attraction_tensor += real(inner(world,occupieds4,mul(world,V,occupieds4)));
     nuclear_attraction_energy = nuclear_attraction_tensor.sum();
     //if(world.rank()==0) print("          Nuclear Attraction done.");
     
     //temp_function = occupieds[j]*Jop;
     //coulomb_energy_correction = real(inner(occupieds[j],temp_function));
     //energies[j] += coulomb_energy_correction + nuclear_attraction_energy_correction;
     //coulomb_energy += 0.5*coulomb_energy_correction;
     coulomb_tensor = real(inner(world,occupieds1,mul(world,Jop,occupieds1)));
     coulomb_tensor += real(inner(world,occupieds2,mul(world,Jop,occupieds2)));
     coulomb_tensor += real(inner(world,occupieds3,mul(world,Jop,occupieds3)));
     coulomb_tensor += real(inner(world,occupieds4,mul(world,Jop,occupieds4)));
     coulomb_energy = 0.5*coulomb_tensor.sum();
     //if(world.rank()==0) print("          Coulomb Repulsion done.");
     
     //exchange_energy_correction = real(inner(occupieds[j],Kpsis[j]));
     //energies[j] -= exchange_energy_correction;
     //exchange_energy += 0.5*exchange_energy_correction;
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
     exchange_tensor += real(inner(world,occupieds3,Kpsis3));
     exchange_tensor += real(inner(world,occupieds4,Kpsis4));
     exchange_energy = 0.5*exchange_tensor.sum();
     //if(world.rank()==0) print("          Exchange done.");

     for(unsigned int i = 0; i < Init_params.num_occupied; i++){
          energies[i] += nuclear_attraction_tensor[i]+coulomb_tensor[i] - exchange_tensor[i];
     }
     total_energy = kinetic_energy + coulomb_energy - exchange_energy + nuclear_attraction_energy + nuclear_repulsion_energy;

     times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Need r function to print r expectation values
     real_function_3d rfunc = real_factory_3d(world).f(myr);

     for(unsigned int j = 0; j < Init_params.num_occupied; j++){
          double r_expec = std::real(inner(occupieds[j], occupieds[j]*rfunc));
          if(world.rank()==0){
               printf("                Orbital: %3i, Energy: %.10e, <r>: %8e\n",j+1, energies[j]-myc*myc, r_expec);
          }
     }
     if(world.rank()==0){
          print("\n              Kinetic Energy: ",kinetic_energy);
          print("             Coulomb  Energy: ",coulomb_energy);
          print("             Exchange Energy: ",exchange_energy);
          print("   Nuclear Attraction Energy: ",nuclear_attraction_energy);
          print("    Nuclear Repulsion Energy: ",nuclear_repulsion_energy);
          print("                Total Energy: ",total_energy);
          print("       Total Energy Residual: ", std::fabs(total_energy - old_total_energy));
          print("Energy Convergence Threshold: " , DFparams.thresh*pow(10,floor(log10(std::fabs(total_energy)))),"\n");
     }
     
     //check total energy for convergence. the global variable thresh determines how many significant figures we look for
     if(std::fabs(total_energy-old_total_energy) > DFparams.thresh*pow(10,floor(log10(std::fabs(total_energy))))){
          iterate_again = true; 
     }

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
     Tensor<double> times = end_timer(world);
     if(world.rank()==0) print("     ", times[0]);

     //Set tolerance for residuals
     double tol = pow(10,floor(0.5*log10(DFparams.thresh)));
     

     //Now time to start iterating
     bool keep_going = true;
     int iteration_number = 1;
     //while(iteration_number < DFparams.max_iter){
     while(keep_going and iteration_number < DFparams.max_iter){
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



