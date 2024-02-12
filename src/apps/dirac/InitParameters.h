
/// \file InitParameters
/// \brief Input parameters for a Dirac Fock calculation, read from a specified archive resulting from a nonrelativistic moldft calculation or restarted from a previous Dirac Fock calculation.


#ifndef MADNESS_APPS_DFGUESSPARAMS_H_INCLUDED
#define MADNESS_APPS_DFGUESSPARAMS_H_INCLUDED

#include "fcwf.h"
#include <madness/chem/NWChem.h>

Function<std::complex<double>,3> function_real2complex(const Function<double,3>& r);
double myxfunc(const madness::coord_3d& r);
double myyfunc(const madness::coord_3d& r);

namespace madness{


     struct InitParameters{
          // Ground state parameters that are read in from archive
          std::string inFile;                      ///< Name of input archive to read in
          double Init_total_energy;                ///< Total energy of the nonrelativistic ground state
          bool spinrestricted;                     ///< Indicates if input calc. was spin-restricted
          bool closed_shell;
          unsigned int num_occupied;               ///< Number of orbitals
          Tensor<double> energies;                 ///< Energies of input orbitals
          Tensor<double> occ;                      ///< Occupancy of input orbitals
          double L;                                ///< Box size of input - Dirac Fock calcluation is in same box
          int order;                               ///< Order of polynomial used in input
          Molecule molecule;                       ///< The molecule used in input calculation
          std::vector<Fcwf> orbitals;              ///< The occupied orbitals 

          // Default constructor
          InitParameters() {}

          // Initializes InitParameters using the contents of file \c filename
          void read(World& world, const std::string& filename, bool restart, bool Krestricted){ 
               // Save the filename
               inFile = filename;
  
               //First check to see if we're starting the job from a saved DF calculation rather than a moldft calculation
               if(restart){
                    if(world.rank()==0) print("\n Reading initial data from restarted DF calculation");
                    archive::ParallelInputArchive input(world, filename.c_str());
                    input & Init_total_energy;
                    input & spinrestricted;
                    input & closed_shell;
                    input & num_occupied;
                    input & energies;
                    input & L;
                    input & order;
                    input & molecule;

                    //Code breaks if spinrestricted (state of archive) and Krestricted (requested by user) don't match
                    //This functionality could probably be added at some point, but it's a bit of work
                    MADNESS_CHECK(spinrestricted == Krestricted);

                    // Set this so we can read in whats
                    // written in the archive 
                    FunctionDefaults<3>::set_k(order);
                    FunctionDefaults<3>::set_cubic_cell(-L, L);

                    //Now we just have to unpack the orbitals
                    for(unsigned int i=0; i < num_occupied; i++){
                         Fcwf reader(world);
                         for(int j=0; j < 4; j++){
                              input & reader[j];
                         }
                         orbitals.push_back(copy(reader));
                    }
               }
               else{ //If we're not reading in from DF, then we're reading in from moldft

                    //some dummy variables for reading/computing
                    std::vector<int> dummy2;
                    Tensor<double> temp_energies;

                    //read in what's in the archive. See SCF::save_mos for how these archives are stored
                    archive::ParallelInputArchive input(world, filename.c_str());
                    unsigned int version=0;
                    std::string xc;
                    std::string localize_method;

                    input & version;
                    input & Init_total_energy;              // double
                    input & spinrestricted;      // bool
                    input & L;                   // double            box size
                    input & order;               // int               wavelet order
                    input & molecule;            // Molecule   
                    input & xc;
                    input & localize_method;

                    input & num_occupied;        // int
                    input & temp_energies;       // Tensor<double>    orbital energies
                    input & occ;                 // Tensor<double>    orbital occupations
                    input & dummy2;              // std::vector<int>  sets of orbitals(?)

                    //For now assume spin-restricted means closed shell in moldft
                    closed_shell = spinrestricted;

                    // Check that order is positive and less than 30
                    if (order < 1 or order > 30){
                         if(world.rank() == 0) print("\n   ***PLEASE NOTE***\n   Invalid wavelet order read from archive, setting to 8.\n   This seems to happen when the default wavelet order is used in moldft."); 
                         order = 8;
                    }

                    // Set this so we can read in whats
                    // written in the archive 
                    FunctionDefaults<3>::set_k(order);
                    FunctionDefaults<3>::set_cubic_cell(-L, L);
                     
                    //Now read in the orbitals and construct Fcwfs from them. 
                    complex_derivative_3d Dx(world,0);
                    complex_derivative_3d Dy(world,1);
                    complex_derivative_3d Dz(world,2);
                    //double myc = 137.0359895; //speed of light in atomic units
                    std::complex<double> myi(0,1);
                    if(spinrestricted){ 
                         //If the calculation was spin-restricted in moldft, then we only have "spin-up" orbitals

                         //Initialize some functions for reading in the orbitals
                         real_function_3d reader;
                         complex_function_3d complexreader;
                         Fcwf spinup(world);
                         Fcwf spindown(world); //used if !Krestricted
                         real_function_3d xfunc = real_factory_3d(world).f(myxfunc);
                         real_function_3d yfunc = real_factory_3d(world).f(myyfunc);
                         
                         //Handle Kramers-restricted and unrestricted cases differently
                         if(Krestricted){
                              //Loop over the occupied orbitals and convert
                              for(unsigned int i = 0; i < num_occupied; i++){
                                   //read in orbital
                                   input & reader;

                                   //change to a complex function
                                   complexreader = function_real2complex(reader);

                                   //build up the corresponding fcwfs, using kinetic balance to define the small component
                                   spinup[0] = complexreader;
                                   spinup[1] = complex_factory_3d(world);
                                   spinup[2] = (-myi) * Dz(complexreader);
                                   spinup[2].scale(0.5);
                                   spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                                   spinup[3].scale(0.5);
                                   spinup.normalize();
                                   orbitals.push_back(spinup);
                              }

                              //Update energies
                              energies = temp_energies;
                         }
                         else{
                              //Loop over the occupied orbitals and convert
                              for(unsigned int i = 0; i < num_occupied; i++){
                                   //read in orbital
                                   input & reader;

                                   //change to a complex function
                                   complexreader = function_real2complex(reader);

                                   //build up the corresponding fcwfs, using kinetic balance to define the small component
                                   spinup[0] = complexreader;
                                   spinup[1] = complex_factory_3d(world);
                                   spinup[2] = (-myi) * Dz(complexreader);
                                   spinup[2].scale(0.5);
                                   spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                                   spinup[3].scale(0.5);
                                   spinup.normalize();
                                   spindown[0] = complex_factory_3d(world);
                                   spindown[1] = complexreader;
                                   spindown[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                                   spindown[2].scale(0.5);
                                   spindown[3] = (myi) * Dz(complexreader);
                                   spindown[3].scale(0.5);
                                   spindown.normalize();
                                   orbitals.push_back(spinup);
                                   orbitals.push_back(spindown);
                              }

                              //Double length of energies tensor and fill in as needed.
                              energies = Tensor<double>(2*num_occupied);
                              for(unsigned int i = 0; i < num_occupied; i++){
                                   energies[2*i] = temp_energies[i];
                                   energies[2*i+1] = temp_energies[i];
                              }
                              num_occupied *= 2;
                         }
                    }
                    else{

                         if(world.rank()==0) print("number of alpha read in from moldft is:" ,num_occupied);

                         // Read in alpha ground state orbitals
                         real_function_3d reader;
                         complex_function_3d complexreader;
                         Fcwf fcwfreader(world);
                         for(unsigned int i = 0; i < num_occupied; i++){
                              input & reader;
                              complexreader = function_real2complex(reader);
                              fcwfreader[0] = complexreader;
                              fcwfreader[1] = complex_factory_3d(world);
                              fcwfreader[2] = (-myi) * Dz(complexreader);
                              fcwfreader[2].scale(0.5);
                              fcwfreader[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                              fcwfreader[3].scale(0.5);
                              fcwfreader.normalize();
                              orbitals.push_back(fcwfreader);
                         }

                         if(!Krestricted){
                              // Read in beta quantities
                              unsigned int num_betas=0;
                              input & num_betas;

                              Tensor<double> beta_energies;
                              input & beta_energies;

                              Tensor<double> dummy3;
                              input & dummy3;

                              std::vector<int> dummy4;
                              input & dummy4;

                              //read in beta ground state orbitals
                              for(unsigned int i = 0; i < num_betas; i++){
                                   input & reader;
                                   complexreader = function_real2complex(reader);
                                   fcwfreader[0] = complex_factory_3d(world);
                                   fcwfreader[1] = complexreader;
                                   fcwfreader[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                                   fcwfreader[2].scale(0.5);
                                   fcwfreader[3] = (myi) * Dz(complexreader);
                                   fcwfreader[3].scale(0.5);
                                   fcwfreader.normalize();
                                   orbitals.push_back(fcwfreader);
                              }

                              //Handle energy tensor and number of occupied orbitals
                              energies = Tensor<double>(num_occupied + num_betas);
                              for(unsigned int i = 0; i < num_occupied; i++){
                                   energies(i) = temp_energies(i);
                              }
                              for(unsigned int i = 0; i < num_betas; i++){
                                   energies(num_occupied + i) = beta_energies(i);
                              }
                              num_occupied += num_betas;

                         }
                         else{
                              energies = temp_energies;
                         }



                    }
                    
                    //reorder orbitals and energies in ascending order, if necessary.
                    double tempdouble;
                    Fcwf fcwfreader(world);
                    for(unsigned int i = 0; i < num_occupied; i++){
                         for(unsigned int j = i+1; j < num_occupied; j++){
                              if(energies(j) < energies(i)){
                                   if(world.rank()==0) print("swapping orbitals", i, " and ", j);
                                   tempdouble = energies(j);
                                   energies(j) = energies(i);
                                   energies(i) = tempdouble;
                                   fcwfreader = orbitals[j];
                                   orbitals[j] = orbitals[i];
                                   orbitals[i] = fcwfreader;
                              }
                         }
                    }
               }
          }

          //This function no longer works
          //TODO: Update this function before using it
          void readnw(World& world, const std::string& filename, bool Krestricted){
               //Called to read in initial parameters from an nwchem output file
               
               //For now just use default values for L and order
               order = 6;
               L = 50.0;
               FunctionDefaults<3>::set_k(order);
               FunctionDefaults<3>::set_cubic_cell(-L, L);

               //Need to set this to something...
               Init_total_energy = 0.0;
               
               //Construct interface object from slymer namespace
               slymer::NWChem_Interface nwchem(filename,std::cout);

               //For parallel runs, silencing all but 1 slymer instance
               if(world.rank() != 0) {
                    std::ostream dev_null(nullptr);
                    nwchem.err = dev_null;
               }

               //Read in basis set
               nwchem.read(slymer::Properties::Basis);

               //Read in the molecular orbital coefficients, energies, and occupancies
               nwchem.read(slymer::Properties::Energies | slymer::Properties::MOs | slymer::Properties::Occupancies);

               //Need to construct a molecule object by ourselves
               molecule = Molecule();
               unsigned int anum;
               double x,y,z,q;
               for(unsigned int i=0; i < nwchem.atoms.size(); i++){
                    anum = symbol_to_atomic_number(nwchem.atoms[i].symbol);
                    q = anum*1.0;
                    x = nwchem.atoms[i].position[0];
                    y = nwchem.atoms[i].position[1];
                    z = nwchem.atoms[i].position[2];
                    molecule.add_atom(x,y,z,q,anum);
               }

               //Find out how many orbitals we're dealing with by looking at the occupancies
               unsigned int numalpha(0), numbeta(0);

               bool have_beta(false);
               for(unsigned int i = 0; i < nwchem.beta_occupancies.size(); i++){
                    if(nwchem.beta_occupancies[i] > 0.0) have_beta = true;
               }
               
               if(have_beta){
                    //we're reading from an unrestricted calculation
                    //and for now we will assume this is an open shell calculation
                    for(unsigned int i = 0; i < nwchem.occupancies.size(); i++){
                         if(nwchem.occupancies[i] == 1.0) numalpha+=1;
                    }
                    for(unsigned int i = 0; i < nwchem.beta_occupancies.size(); i++){
                         if(nwchem.beta_occupancies[i] == 1.0) numbeta+=1;
                    }

                    //Right now DF can only handle a single unpaired electron
                    MADNESS_CHECK(numalpha-1 == numbeta);
               }
               else{
                    for(unsigned int i = 0; i < nwchem.occupancies.size(); i++){
                         if(nwchem.occupancies[i] == 2.0) numalpha += 1;
                    }
                    numbeta = numalpha;
               }
               closed_shell = !have_beta;
               
               //correctly set the number of occupied orbitals for the DF calculation
               if(Krestricted){
                    num_occupied = numalpha;
               }
               else{
                    num_occupied = numalpha+numbeta;
               }


               //Let's print everything so we have a visual check on what we're working with (for now)
               if(world.rank()==0) print("\nalpha occupancies:\n",nwchem.occupancies);
               if(world.rank()==0) print("\nbeta occupancies:\n",nwchem.beta_occupancies);
               if(world.rank()==0) print("\nenergies:\n",nwchem.energies);
               if(world.rank()==0) print("\nbeta energies:\n",nwchem.beta_energies);
               if(world.rank()==0) print("num alpha",numalpha);
               if(world.rank()==0) print("num beta",numbeta);


               //Now that we know how many orbitals we have. initialize and fill energy tensor
               energies = Tensor<double>(num_occupied);
               if(Krestricted){
                    for(unsigned int i=0; i < numalpha; i++){
                         energies[i] = nwchem.energies[i];
                    }
               }
               else{
                    if(closed_shell){
                         for(unsigned int i=0; i < numalpha; i++){
                              energies[2*i] = nwchem.energies[i];
                              energies[2*i+1] = nwchem.energies[i];
                         }
                    }
                    else{
                         for(unsigned int i=0; i < numalpha-1; i++){
                              energies[2*i] = nwchem.energies[i];
                              energies[2*i+1] = nwchem.energies[i];
                         }
                         energies[2*(numalpha-1)] = nwchem.energies[numalpha-1];
                    }
               }

               //Cast the 'basis set' into a Gaussian basis and iterate over it
               vector_real_function_3d temp1;
               int ii = 0;
               for(auto basis : slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set)) {
                    //Get the center of gaussian as its special point
                    std::vector<coord_3d> centers;
                    coord_3d r;
                    r[0] = basis.get().center[0]; r[1] = basis.get().center[1]; r[2] = basis.get().center[2];
                    centers.push_back(r);

                    //Now make the function
                    temp1.push_back(FunctionFactory<double,3>(world).functor(std::shared_ptr<FunctionFunctorInterface<double,3>>(new slymer::Gaussian_Functor(basis.get(), centers))));
                    double norm2 = temp1[ii].norm2();
                    if(world.rank() == 0) print("function", ii, "has norm", norm2);
                    ii++;
               }

               //Normalize aos
               normalize(world, temp1);

               //Transform aos now to get alpha mos
               vector_real_function_3d temp = transform(world, temp1, nwchem.MOs , true);

               //Convert and store alpha occupied MOs.
               complex_function_3d complexreader(world);
               Fcwf spinup(world);
               Fcwf spindown(world);
               complex_derivative_3d Dx(world,0);
               complex_derivative_3d Dy(world,1);
               complex_derivative_3d Dz(world,2);
               std::complex<double> myi(0,1);
               for(unsigned int i = 0; i < numalpha-1; i++){
                    complexreader = function_real2complex(temp[i]);
                    spinup[0] = complexreader;
                    spinup[1] = complex_factory_3d(world);
                    spinup[2] = (-myi) * Dz(complexreader);
                    spinup[2].scale(0.5);
                    spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                    spinup[3].scale(0.5);
                    spinup.normalize();
                    orbitals.push_back(spinup);
                    if(!Krestricted){
                         spindown[0] = complex_factory_3d(world);
                         spindown[1] = complexreader;
                         spindown[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                         spindown[2].scale(0.5);
                         spindown[3] = (myi) * Dz(complexreader);
                         spindown[3].scale(0.5);
                         spindown.normalize();
                         orbitals.push_back(spindown);
                    }
               }
               complexreader = function_real2complex(temp[numalpha-1]);
               spinup[0] = complexreader;
               spinup[1] = complex_factory_3d(world);
               spinup[2] = (-myi) * Dz(complexreader);
               spinup[2].scale(0.5);
               spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
               spinup[3].scale(0.5);
               spinup.normalize();
               orbitals.push_back(spinup);
               if(closed_shell and not Krestricted){
                    spindown[0] = complex_factory_3d(world);
                    spindown[1] = complexreader;
                    spindown[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                    spindown[2].scale(0.5);
                    spindown[3] = (myi) * Dz(complexreader);
                    spindown[3].scale(0.5);
                    spindown.normalize();
                    orbitals.push_back(spindown);
               }

               //Assure that the numbers line up
               MADNESS_ASSERT(num_occupied == orbitals.size());

          }

          // Prints all information
          void print_params() const
          {
               madness::print("\n     Input Parameters");
               madness::print("   -------------------------");
               madness::print("         Input Archive:", inFile);
               madness::print("       Spin Restricted:", spinrestricted);
               madness::print("  No. of occ. orbitals:", num_occupied);
               madness::print("                     L:", L);
               madness::print("         Wavelet Order:", order);
               madness::print("  Initial Total Energy:", Init_total_energy);
               madness::print("      Orbital Energies:", energies);
          }
     };
}
#endif

//kthxbye
