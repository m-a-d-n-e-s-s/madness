
/// \file DFParameters
/// \brief Input parameters for a Dirac Fock calculation.

#ifndef MADNESS_APPS_DFPARAMS_H_INCLUDED
#define MADNESS_APPS_DFPARAMS_H_INCLUDED

#include <chem/molecule.h>

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
                                       ///<     1    |   "Method 1" calculation of virtuals of a
                                       ///<          |   single-valence state.
                                       ///<
          int print_level;             ///< Controls the amount and style of printing. Higher values print more
                                       ///<   Values |   What gets printed
                                       ///<   ----------------------------
                                       ///<     1    |   Print out each step in the calculation,
                                       ///<          |   along with timings
                                       ///<   ----------------------------
                                       ///<     2    |   Debug level. Prints EVERYTHING!!!

          bool plot;                   ///< Turn on plotting of final orbitals. Output format is .vts 
          bool plot_range;             ///< Controls which orbitals will be plotted 
          std::vector<int> plot_data;  ///< Orbitals to plot
          int max_iter;                ///< Maximum number of iterations
          //double econv;                ///< Convergence criterion for the orbital energies
          double small;                ///< Minimum length scale to be resolved
          double thresh;               ///< Accuracy criterion when truncating
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
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //               If you add something here, don't forget to add it to serializable!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          // NOT YET IMPLEMENTED
          std::vector<double> protocol_data;
          bool localized;              ///< Flag to use localized orbitals or not. MUST BE TRUE IF USING LOCALIZED GROUND STATE ORBITALS!!

          template<typename Archive>
          void serialize(Archive& ar){
               ar & archive & job & print_level & max_iter & small & thresh & k & kain & maxsub & maxrotn & restart & nucleus & do_save & savefile & lb_iter & nwchem & lineplot;
          }

          // Default constructor
          DFParameters()
          : job(0)
          , print_level(0)
          , plot(false)
          , max_iter(20)
          , small(1e-5)
          , k(8)
          , thresh(1e-6)
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
                    else if (s == "print_level"){
                         f >> print_level;
                    }
                    else if (s == "plot"){
                         plot = true;
                         std::string buf;
                         std::getline(f,buf);
                         plot_data = std::vector<int>();
                         std::string d;
                         std::stringstream t(buf);
                         t >> d;
                         if (d == "range"){
                              plot_range = true;
                              t >> d;
                              plot_data.push_back(std::stoi(d));
                              t >> d;
                              for(int z = plot_data[0]; z < std::stoi(d); z++) plot_data.push_back(z);
                         }
                         else{
                              // Add in the one we just read
                              plot_data.push_back(std::stoi(d));
                              while(t >> d)  plot_data.push_back(std::stoi(d));
                         }
                    }
                    else if (s == "max_iter"){
                         f >> max_iter;
                    }
                    else if (s == "small"){
                         f >> small;
                    }
                    else if (s == "thresh"){
                         f >> thresh;
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
               else{
                    madness::print("                       Nucleus: gaussian");
               }
               madness::print("                  Do Lineplots:", lineplot);
               madness::print("           Plot Final Orbitals:", plot);
               if(plot and plot_range) madness::print("                 Plot Start:", plot_data[0]);
               if(plot and plot_range) madness::print("                   Plot End:", plot_data[1]);
               if(plot and not plot_range) madness::print("        Orbitals to be Plotted:", plot_data);
               madness::print("                   Print Level:", print_level);
          }
     };

}
#endif

//kthxbye
