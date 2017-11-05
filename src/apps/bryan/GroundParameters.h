
/// \file GroundParameters
/// \brief Input parameters for a response calculation, read from a specified archive.


#ifndef MADNESS_APPS_GROUNDPARAMS_H_INCLUDED
#define MADNESS_APPS_GROUNDPARAMS_H_INCLUDED

namespace madness 
{

   struct GroundParameters
   {
      // Ground state parameters that are read in from archive
      std::string inFile;                      ///< Name of input archive to read in ground state
      bool spinrestricted;                     ///< Indicates if ground state calc. was open or closed shell
      unsigned int num_orbitals;               ///< Number of orbitals in ground state
      Tensor<double> energies;                 ///< Energy of ground state orbitals
      Tensor<double> occ;                      ///< Occupancy of ground state orbitals
      double L;                                ///< Box size of ground state - response calcluation is in same box
      int order;                               ///< Order of polynomial used in ground state
      Molecule molecule;                       ///< The molecule used in ground state calculation
      std::vector<real_function_3d> orbitals;  ///< The ground state orbitals 

      // Default constructor
      GroundParameters() {}

      // Initializes ResponseParameters using the contents of file \c filename
      void read(World& world, const std::string& filename)
      { 
         // Save the filename
         inFile = filename;

         double dummy1;
         std::vector<int> dummy2;

         archive::ParallelInputArchive input(world, filename.c_str());
         input & dummy1;              // double
         input & spinrestricted;      // bool
         input & num_orbitals;        // int
         input & energies;            // Tensor<double>    orbital energies
         input & occ;                 // Tensor<double>    orbital occupations
         input & dummy2;              // std::vector<int>  sets of orbitals(?)
         input & L;                   // double            box size
         input & order;               // int               wavelet order
         input & molecule;            // Molecule   

         // Check that order is positive and less than 30
         if (order < 1 or order > 30)
         {
            if(world.rank() == 0) print("\n   ***PLEASE NOTE***\n   Invalid wavelet order read from archive, setting to 8.\n   This seems to happen when the default wavelet order is used in moldft."); 
            order = 8;
         }

         // Set this so we can read in whats
         // written in the archive 
         //FunctionDefaults<3>::set_k(order);

         // Read in ground state orbitals
         for(unsigned int i = 0; i < num_orbitals; i++)
         {
              real_function_3d reader;
              input & reader;
              orbitals.push_back(reader); 
         }
      }

      // Prints all information
      void print_params() const
      {
         madness::print("\n   Ground State Parameters");
         madness::print("   -------------------------");
         madness::print("  Ground State Archive:", inFile);
         madness::print("       Spin Restricted:", spinrestricted);
         madness::print("    Number of orbitals:", num_orbitals);
         madness::print("                     L:", L);
         madness::print("         Wavelet Order:", order);
         madness::print("      Orbital Energies:", energies);
      }
   };
}
#endif
