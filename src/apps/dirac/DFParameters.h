
/// \file DFParameters
/// \brief Input parameters for a Dirac Fock calculation.

#ifndef MADNESS_APPS_DFPARAMS_H_INCLUDED
#define MADNESS_APPS_DFPARAMS_H_INCLUDED

#include <madness/chem/molecule.h>

namespace madness {

     struct DFParameters{
          // List of input parameters
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //               If you add something here, don't forget to add it to serializable!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          std::string archive;         ///< Name of input archive to read in ground state
          int job;                     ///< Indicates the type of job to do
                                       ///<   Values |   Job
                                       ///<   --------------------------------------------------------------
                                       ///<     0    |   Dirac Fock on occupied orbitals only (Default)
                                       ///<   --------------------------------------------------------------
                                       ///<     Add more values here
                                       ///<
          int max_iter;                ///< Maximum number of iterations
          double small;                ///< Minimum length scale to be resolved
          double thresh;               ///< Accuracy criterion when truncating
          double dconv;                ///< Accuracy criterion for charge density. Defaults to thresh.
          int k;                       ///< Number of legendre polynomials in scaling basis
          bool kain;                   ///< Turns on KAIN nonlinear solver 
          int maxsub;                  ///< Sets maximum subspace size for KAIN
          double maxrotn;              ///< maximum step allowed by kain
          bool restart;                ///< Indicates this is a restarted DF job
          int nucleus;                 ///< Indicates which nucleus model to use (1 for fermi, anything else for Gaussian)
          bool do_save;                ///< Whether or not to save after each iteration. Defaults to true. Turn off with 'no_save'
          std::string savefile;        ///< Gives the file to save the archive each iteration Default: DFrestartdata (in working directory)
          int lb_iter;                 ///< How many iterations to load balance (after the initial load balancing)
          bool nwchem;                 ///< Indicates archive given is actually an nwchem file for starting the job
          bool lineplot;               ///< Whether or not to make lineplots at the end of the job
          bool no_compute;             ///< If true, will skip all computation
          double bohr_rad;             ///< bohr radius in fm (default: 52917.7211)
          int min_iter;                ///< minimum number of iterations (default: 2)
          bool Krestricted;            ///< Calculation should be performed in Kramers-restricted manner (default: false)
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //               If you add something here, don't forget to add it to serializable!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          template<typename Archive>
          void serialize(Archive& ar){
               ar & archive & job & max_iter & small & thresh & k & kain & maxsub & maxrotn & restart & nucleus & do_save & savefile & lb_iter & nwchem & lineplot & no_compute & bohr_rad & min_iter & Krestricted;
          }

          // Default constructor
          DFParameters()
          : job(0)
          , max_iter(20)
          , small(1e-5)
          , thresh(1e-6)
          , dconv(1e-6)
          , k(8)
          , kain(false)
          , maxsub(10)
          , maxrotn(0.25)
          , restart(false)
          , nucleus(0)
          , do_save(true)
          , savefile("DFrestartdata")
          , lb_iter(20)
          , nwchem(false)
          , lineplot(false)
          , no_compute(false)
          , bohr_rad(52917.7210544)  // bohr radius in fm from CODATA 2022
          , min_iter(2)
          , Krestricted(false)
          {}

          // Initializes DFParameters using the contents of file \c filename
          void read_file(const std::string& filename){
               std::ifstream f(filename.c_str());
               read(f);
          }

          // Initializes DFParameters using the contents of stream \c f
          void read(std::istream& f){
               position_stream(f, "DiracFock");
               std::string s;

               while(f >> s){
                    if(s == "end"){
                         break;
                    }
                    else if (s == "archive"){
                         f >> archive;
                    }
                    else if (s == "job"){
                         f >> job;
                    }
                    else if (s == "max_iter"){
                         f >> max_iter;
                    }
                    else if (s == "small"){
                         f >> small;
                    }
                    else if (s == "thresh"){
                         f >> thresh;
                         dconv = thresh;
                    }
                    else if (s == "dconv"){
                         f >> dconv;
                    }
                    else if (s == "k"){
                         f >> k;
                    }
                    else if (s == "kain"){
                         kain = true;
                    }
                    else if (s == "maxsub"){
                         f >> maxsub;
                    }
                    else if (s == "maxrotn"){
                         f >> maxrotn;
                    }
                    else if (s == "restart"){
                         restart = true;
                    }
                    else if (s == "nucleus"){
                         f >> nucleus;
                    }
                    else if (s == "no_save"){
                         do_save = false;
                    }
                    else if (s == "savefile"){
                         f >> savefile;
                    }
                    else if (s == "lb_iter"){
                         f >> lb_iter;
                    }
                    else if (s == "nwchem"){
                         nwchem = true;
                    }
                    else if (s == "lineplot"){
                         lineplot = true;
                    }
                    else if (s == "no_compute"){
                         no_compute = true;
                    }
                    else if (s == "bohr_rad"){
                         f >> bohr_rad;
                    }
                    else if (s == "min_iter"){
                         f >> min_iter;
                    }
                    else if (s == "Krestricted"){
                         Krestricted = true;
                    }
                    else{
                       std::cout << "Dirac Fock: unrecognized input keyword " << s << std::endl;
                       MADNESS_EXCEPTION("input error", 0); 
                    }
               }
          } // end read()

          // Prints all information
          void print_params() const{
               madness::print("\n   Input Dirac Fock Parameters");
               madness::print("   -------------------------");
               madness::print("            Initial Guess File:", archive);
               madness::print("                           Job:", job);
               madness::print("          Refinement Threshold:", thresh);
               madness::print("                             k:", k);
               madness::print("Smallest Resolved Length Scale:", small);
               madness::print("                Max Iterations:", max_iter);
               madness::print("               Use KAIN Solver:", kain);
               if(kain) madness::print("     KAIN Solver Subspace Size:", maxsub);
               madness::print("                          Save:", do_save);
               madness::print("                     save file:", savefile);
               if(nucleus == 1){
                    madness::print("                       Nucleus: fermi");
               }
               else if (nucleus == 2) {
                    madness::print("                       Nucleus: point");
               }
               else{
                    madness::print("                       Nucleus: gaussian");
               }
               madness::print("           Kramers restriction:", Krestricted);
               madness::print("                  Do Lineplots:", lineplot);
          }
     };

}
#endif

//kthxbye
