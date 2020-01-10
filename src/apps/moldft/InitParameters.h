
/// \file InitParameters
/// \brief Input parameters for a Dirac Fock calculation, read from a specified archive resulting from a nonrelativistic moldft calculation or restarted from a previous Dirac Fock calculation.


#ifndef MADNESS_APPS_DFGUESSPARAMS_H_INCLUDED
#define MADNESS_APPS_DFGUESSPARAMS_H_INCLUDED

#include "fcwf.h"
#include "NWChem.h"

Function<std::complex<double>,3> function_real2complex(const Function<double,3>& r);
double myxfunc(const madness::coord_3d& r);
double myyfunc(const madness::coord_3d& r);

namespace madness{


     struct InitParameters{
          // Ground state parameters that are read in from archive
          std::string inFile;                      ///< Name of input archive to read in
          double Init_total_energy;                ///< Total energy of the nonrelativistic ground state
          bool spinrestricted;                     ///< Indicates if input calc. was spin-restricted
          unsigned int num_occupied;               ///< Number of orbitals
          unsigned int num_virtuals;
          Tensor<double> energies;                 ///< Energies of input orbitals
          Tensor<double> v_energies;
          Tensor<double> occ;                      ///< Occupancy of input orbitals
          double L;                                ///< Box size of input - Dirac Fock calcluation is in same box
          int order;                               ///< Order of polynomial used in input
          Molecule molecule;                       ///< The molecule used in input calculation
          std::vector<Fcwf> orbitals;              ///< The occupied orbitals 
          std::vector<Fcwf> virtuals;

          // Default constructor
          InitParameters() {}

          // Initializes InitParameters using the contents of file \c filename
          void read(World& world, const std::string& filename, bool restart){ 
               // Save the filename
               inFile = filename;
  
               if(restart){
                    if(world.rank()==0) print("\n Reading initial data from restarted DF calculation");
                    archive::ParallelInputArchive input(world, filename.c_str());
                    input & Init_total_energy;
                    input & spinrestricted;
                    input & num_occupied;
                    input & energies;
                    input & L;
                    input & order;
                    input & molecule;

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
               else{
                    std::vector<int> dummy2;
                    Tensor<double> temp_energies;

                    archive::ParallelInputArchive input(world, filename.c_str());
                    input & Init_total_energy;              // double
                    input & spinrestricted;      // bool
                    input & num_occupied;        // int
                    input & energies;            // Tensor<double>    orbital energies
                    input & occ;                 // Tensor<double>    orbital occupations
                    input & dummy2;              // std::vector<int>  sets of orbitals(?)
                    input & L;                   // double            box size
                    input & order;               // int               wavelet order
                    input & molecule;            // Molecule   

                    // Check that order is positive and less than 30
                    if (order < 1 or order > 30){
                         if(world.rank() == 0) print("\n   ***PLEASE NOTE***\n   Invalid wavelet order read from archive, setting to 8.\n   This seems to happen when the default wavelet order is used in moldft."); 
                         order = 8;
                    }

                    // Set this so we can read in whats
                    // written in the archive 
                    FunctionDefaults<3>::set_k(order);
                    FunctionDefaults<3>::set_cubic_cell(-L, L);
                     
                     
                    complex_derivative_3d Dx(world,0);
                    complex_derivative_3d Dy(world,1);
                    complex_derivative_3d Dz(world,2);
                    double myc = 137.0359895; //speed of light in atomic units
                    std::complex<double> myi(0,1);
                    if(spinrestricted){ 
                         // Read in nonrelativistic ground state orbitals, changing them to Fcwfs
                         real_function_3d reader;
                         complex_function_3d complexreader;
                         Fcwf spinup(world);
                         Fcwf spindown(world);
                         real_function_3d xfunc = real_factory_3d(world).f(myxfunc);
                         real_function_3d yfunc = real_factory_3d(world).f(myyfunc);
                         for(unsigned int i = 0; i < num_occupied; i++){
                              input & reader;
                              complexreader = function_real2complex(reader);
                              spinup[0] = complexreader;
                              spinup[1] = complex_factory_3d(world);
                              spinup[2] = (-myi) * Dz(complexreader);
                              spinup[2].scale(0.5);
                              spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                              spinup[3].scale(0.5);
                              spinup.normalize();
                              //spindown[0] = complex_factory_3d(world);
                              //spindown[1] = complexreader;
                              //spindown[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                              //spindown[2].scale(0.5);
                              //spindown[3] = (myi) * Dz(complexreader);
                              //spindown[3].scale(0.5);
                              //spindown.normalize();
                              //spinup = Fcwf(complexreader, complex_factory_3d(world), complex_factory_3d(world), complex_factory_3d(world));
                              //spindown = Fcwf(complex_factory_3d(world), copy(complexreader), complex_factory_3d(world), complex_factory_3d(world));
                              orbitals.push_back(spinup);
                              //orbitals.push_back(spindown);
                         }

                         ////Bring in p block, construct 4-component wavefunctions, and diagonalize under j_z
                         //std::vector<Fcwf> jzorbitals;
                         //std::vector<Fcwf> justorbitals;
                         //for(unsigned int i=0; i<3; i++){
                         //     input & reader;
                         //     complexreader = function_real2complex(reader);
                         //     spinup[0] = complexreader;
                         //     spinup[1] = complex_factory_3d(world);
                         //     spinup[2] = (-myi) * Dz(complexreader);
                         //     spinup[2].scale(0.5);
                         //     spinup[3] = (-myi) * (Dx(complexreader) + myi * Dy(complexreader));
                         //     spinup[3].scale(0.5);
                         //     spinup.normalize();
                         //     spindown[0] = complex_factory_3d(world);
                         //     spindown[1] = complexreader;
                         //     spindown[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                         //     spindown[2].scale(0.5);
                         //     spindown[3] = (myi) * Dz(complexreader);
                         //     spindown[3].scale(0.5);
                         //     spindown.normalize();
                         //     ////spinup = Fcwf(complexreader, complex_factory_3d(world), complex_factory_3d(world), complex_factory_3d(world));
                         //     ////spindown = Fcwf(complex_factory_3d(world), copy(complexreader), complex_factory_3d(world), complex_factory_3d(world));
                         //     //spinup = spinup + spindown;
                         //     //orbitals.push_back(spinup);
                         //     //orbitals.push_back(spindown);
                         //     justorbitals.push_back(spinup);
                         //     justorbitals.push_back(spindown);
                         //}


                         ////for(unsigned int i=0; i<6; i++){
                         ////     spinup[0] = -myi*(xfunc*Dy(justorbitals[i][0])-yfunc*Dx(justorbitals[i][0])) + 0.5*justorbitals[i][0];
                         ////     spinup[1] = -myi*(xfunc*Dy(justorbitals[i][1])-yfunc*Dx(justorbitals[i][1])) - 0.5*justorbitals[i][1];
                         ////     spinup[2] = -myi*(xfunc*Dy(justorbitals[i][2])-yfunc*Dx(justorbitals[i][2])) + 0.5*justorbitals[i][2];
                         ////     spinup[3] = -myi*(xfunc*Dy(justorbitals[i][3])-yfunc*Dx(justorbitals[i][3])) - 0.5*justorbitals[i][3];
                         ////     jzorbitals.push_back(spinup);
                         ////}

                         ////Tensor<std::complex<double>> jztensor = matrix_inner(world,justorbitals,jzorbitals);
                         ////Tensor<std::complex<double>> overlap = matrix_inner(world,justorbitals,justorbitals);

                         ////Tensor<std::complex<double>> U(6,6);
                         ////Tensor<double> evals(6);

                         ////if(world.rank()==0){
                         ////     print("\njz:\n", jztensor);
                         ////     print("\noverlap:\n", overlap);
                         ////}
                         ////
                         ////sygv(jztensor,overlap,1,U,evals);

                         ////if(world.rank()==0) print("U:\n",U);
                         ////if(world.rank()==0) print("evals:\n",evals);

                         ////transform(world,justorbitals,U);
                         ////for(unsigned int i = 0; i < 6; i++){
                         ////     if(evals[i] > 0){
                         ////          if(world.rank()==0) print("Add in state with eval: ", evals[i]);
                         ////          orbitals.push_back(justorbitals[i]);
                         ////     }
                         ////}
                         //orbitals.push_back(justorbitals[0]);
                         //orbitals.push_back(justorbitals[2]);
                         //orbitals.push_back(justorbitals[4]);
                         ////duplicate the energies
                         ////energies = Tensor<double>(2*num_occupied);
                         ////double csquared = 137.0359895*137.0359895;
                         ////double temp;
                         ////for(unsigned int i = 0; i < num_occupied; i++){
                         ////     temp = temp_energies(i);//+csquared;
                         ////     energies(2*i) = temp;
                         ////     energies(2*i+1) = temp;
                         ////}

                         ////correct the number of orbitals
                         //num_occupied *= 2;
                    }
                    else{

                         if(world.rank()==0) print("number of alpha read in is:" ,num_occupied);

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

                         //// Read in beta quantities
                         //unsigned int num_betas;
                         //input & num_betas;

                         //Tensor<double> beta_energies;
                         //input & beta_energies;
                         ////NEED TO ADD THESE INTO TOTAL ENERGIES MATRIX

                         //Tensor<double> dummy3;
                         //input & dummy3;

                         //std::vector<int> dummy4;
                         //input & dummy4;

                         //if(world.rank()==0) print("made it here: ", num_betas);
                         //
                         ////read in beta ground state orbitals
                         //for(unsigned int i = 0; i < num_betas; i++){
                         //     input & reader;
                         //     complexreader = function_real2complex(reader);
                         //     fcwfreader[0] = complex_factory_3d(world);
                         //     fcwfreader[1] = complexreader;
                         //     fcwfreader[2] = (-myi) * (Dx(complexreader) - myi * Dy(complexreader));
                         //     fcwfreader[2].scale(0.5);
                         //     fcwfreader[3] = (myi) * Dz(complexreader);
                         //     fcwfreader[3].scale(0.5);
                         //     fcwfreader.normalize();
                         //     orbitals.push_back(fcwfreader);
                         //}

                         ////fix up energies tensor and num_occupied
                         //energies = Tensor<double>(num_occupied + num_betas);
                         //double csquared = 137.0359895*137.0359895;
                         //for(unsigned int i = 0; i < num_occupied; i++){
                         //     energies(i) = temp_energies(i);// + csquared;
                         //}
                         //for(unsigned int i = 0; i < num_betas; i++){
                         //     energies(num_occupied + i) = beta_energies(i);// + csquared;
                         //}
                         //num_occupied += num_betas;

                    }

                    //For convenience, reorder orbitals and energies in ascending order.
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

          void readnw(World& world, const std::string& filename){
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
               num_virtuals = 0;

               bool have_beta(false);
               for(unsigned int i = 0; i < nwchem.beta_occupancies.size(); i++){
                    if(nwchem.beta_occupancies[i] > 0.0) have_beta = true;
               }
               
               if(have_beta){
                    //we're reading from an open-shell calculation
                    for(unsigned int i = 0; i < nwchem.occupancies.size(); i++){
                         (nwchem.occupancies[i] == 1.0) ? numalpha+=1 : num_virtuals+=1;
                    }
                    for(unsigned int i = 0; i < nwchem.beta_occupancies.size(); i++){
                         (nwchem.beta_occupancies[i] == 1.0) ? numbeta+=1 : num_virtuals+=1;
                    }
               }
               else{
                    for(unsigned int i = 0; i < nwchem.occupancies.size(); i++){
                         (nwchem.occupancies[i] == 2.0) ? numalpha += 1 : num_virtuals += 2;
                    }
                    numbeta = numalpha;
               }
               num_occupied = numalpha+numbeta;

               //Let's print everything so we have a visual check on what we're working with (for now)
               if(world.rank()==0) print("\nalpha occupancies:\n",nwchem.occupancies);
               if(world.rank()==0) print("\nbeta occupancies:\n",nwchem.beta_occupancies);
               if(world.rank()==0) print("\nenergies:\n",nwchem.energies);
               if(world.rank()==0) print("\nbeta energies:\n",nwchem.beta_energies);
               if(world.rank()==0) print("num alpha",numalpha);
               if(world.rank()==0) print("num beta",numbeta);
               if(world.rank()==0) print("num virtuals",num_virtuals);


               //Now that we know how many orbitals we have. initialize energy tensors
               energies = Tensor<double>(num_occupied);
               v_energies = Tensor<double>(num_virtuals);

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

               //Keep track of how many energies i've stored
               int energy_index = 0;
               int v_energy_index = 0;

               //nonrelativistic energies need to be adjusted for relativistic calculations
               double csquared = 137.0359895*137.0359895;

               //Convert and store alpha occupied MOs and virtuals. If closed shell, do betas too.
               complex_function_3d complexreader(world);
               Fcwf fcwfreader(world);
               for(unsigned int i = 0; i < temp.size(); i++){
                    complexreader = function_real2complex(temp[i]);
                    fcwfreader = Fcwf(complexreader, complex_factory_3d(world), complex_factory_3d(world), complex_factory_3d(world));
                    if(have_beta){
                         if(nwchem.occupancies[i] == 1.0){
                              orbitals.push_back(fcwfreader);
                              energies[energy_index] = nwchem.energies[i] + csquared;
                              energy_index++;
                         }
                         else{
                              virtuals.push_back(fcwfreader);
                              v_energies[v_energy_index] = nwchem.energies[i] + csquared;
                              v_energy_index++;
                         }

                    }
                    else{
                         if(nwchem.occupancies[i] == 2.0){
                              orbitals.push_back(fcwfreader);
                              energies[energy_index] = nwchem.energies[i] + csquared;
                              energy_index++;
                              fcwfreader = Fcwf(complex_factory_3d(world), complexreader, complex_factory_3d(world), complex_factory_3d(world));
                              orbitals.push_back(fcwfreader);
                              energies[energy_index] = nwchem.energies[i] + csquared;
                              energy_index++;

                         }
                         else{
                              virtuals.push_back(fcwfreader);
                              v_energies[v_energy_index] = nwchem.energies[i] + csquared;
                              v_energy_index++;
                              fcwfreader = Fcwf(complex_factory_3d(world), complexreader, complex_factory_3d(world), complex_factory_3d(world));
                              virtuals.push_back(fcwfreader);
                              v_energies[v_energy_index] = nwchem.energies[i] + csquared;
                              v_energy_index++;
                         }
                    }
               }


               //If we're doing an open shell calculation we need to read in the beta orbitals explicitly
               if(have_beta){
                    //Transform aos again to get beta mos
                    temp = transform(world, temp1, nwchem.beta_MOs, true);

                    //Convert and store beta occupied MOs and virtuals
                    for(unsigned int i = 0; i < temp.size(); i++){
                         complexreader = function_real2complex(temp[i]);
                         fcwfreader = Fcwf(complex_factory_3d(world), complexreader, complex_factory_3d(world), complex_factory_3d(world));
                         if(nwchem.beta_occupancies[i] == 1.0){
                              orbitals.push_back(fcwfreader);
                              energies[energy_index] = nwchem.beta_energies[i] + csquared;
                              energy_index++;
                         }
                         else{
                              virtuals.push_back(fcwfreader);
                              v_energies[v_energy_index] = nwchem.beta_energies[i] + csquared;
                              v_energy_index++;
                         }
                    }
               }
               
               //Assure that the numbers line up
               MADNESS_ASSERT(num_occupied == orbitals.size());
               MADNESS_ASSERT(num_virtuals == virtuals.size());

               //Need to sort
               double tempdouble;
               //Sort occupied
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
               //sort virtual
               for(unsigned int i = 0; i < num_virtuals; i++){
                    for(unsigned int j = i+1; j < num_virtuals; j++){
                         if(v_energies(j) < v_energies(i)){
                              if(world.rank()==0) print("swapping virtuals", i, " and ", j);
                              tempdouble = v_energies(j);
                              v_energies(j) = v_energies(i);
                              v_energies(i) = tempdouble;
                              fcwfreader = virtuals[j];
                              virtuals[j] = virtuals[i];
                              virtuals[i] = fcwfreader;
                         }
                    }
               }

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
