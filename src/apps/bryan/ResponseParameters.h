
/// \file ResponseParameters
/// \brief Input parameters for a response calculation.

#ifndef MADNESS_APPS_RESPONSEPARAMS_H_INCLUDED
#define MADNESS_APPS_RESPONSEPARAMS_H_INCLUDED

#include <chem/molecule.h>

namespace madness 
{

   struct ResponseParameters
   {
      // List of input parameters
      std::string archive;               ///< Name of input archive to read in ground state
      std::string nwchem;                ///< Root name of nwchem files for intelligent starting guess
      int states;                        ///< Number of excited states requested
      int print_level;                   ///< Controls the amount and style of printing. Higher values print more
                                         ///<   Values |   What gets printed
                                         ///<   ----------------------------
                                         ///<     1    |   Print out each step in the calculation,
                                         ///<          |   along with timings
                                         ///<   ----------------------------
                                         ///<     2    |   Debug level. Prints EVERYTHING!!!

      bool tda;                          ///< Turn on Tam-Danchof approximation (only calculate excitations)
      bool plot;                         ///< Turn on plotting of final orbitals. Output format is .vts 
      bool plot_range;                   ///< Controls which orbitals will be plotted 
      std::vector<int> plot_data;        ///< Orbitals to plot
      double plot_L;                     ///< Controls the plotting box size 
      int plot_pts;                      ///< Controls number of points in plots
      int max_iter;                      ///< Maximum number of iterations
      double dconv;                      ///< Convergence criterion for the orbital density 
      double small;                      ///< Minimum length scale to be resolved
      std::vector<double> protocol_data; ///< Different thresholds for truncation
      int larger_subspace;               ///< Number of iterations to diagonalize in a subspace consisting of old and new vectors
      int k;                             ///< Polynomial order to use in calculation
      bool random;                       ///< Use a random guess for initial response functions
      bool store_potential;              ///< Store the potential instead of computing each iteration
      bool e_window;                     ///< Use an energy window to excite from
      double range_low;                  ///< Energy range (lower end) for orbitals to excite from
      double range_high;                 ///< Energy range (upper end) for orbitals to excite from
      bool plot_initial;                 ///< Flag to plot the ground state orbitals read in from archive
      bool localized;                    ///< Flag to use localized orbitals or not. MUST BE TRUE IF USING LOCALIZED
      bool restart;                      ///< Flag to restart from file
      std::string resp_archive;          ///< Response restart archive

      // NOT YET IMPLEMENTED
      std::string xc_data;
      bool kain;  

      // Used to broadcast data to all mpi ranks
      template<typename Archive>
      void serialize(Archive& ar)
      {
         ar & archive
            & nwchem 
            & states 
            & print_level 
            & tda 
            & plot 
            & plot_range 
            & plot_data 
            & plot_L 
            & plot_pts 
            & max_iter 
            & dconv
            & small
            & protocol_data 
            & larger_subspace
            & k
            & random 
            & store_potential 
            & e_window 
            & range_low 
            & range_high 
            & plot_initial
            & localized
            & restart 
            & resp_archive;
      }

      // Default constructor
      ResponseParameters()
      : states(1)
      , nwchem("")
      , print_level(1)
      , tda(false)
      , plot(false)
      , plot_L(-1.0)
      , plot_pts(201)
      , max_iter(20)
      , dconv(1e-3)
      , small(1e-6)
      , protocol_data(madness::vector_factory(1e-4, 1e-6))
      , larger_subspace(0)
      , k(0)
      , random(false)
      , store_potential(false)
      , e_window(false)
      , range_low(0.0)
      , range_high(1.0)
      , plot_initial(false)
      , localized(false)
      , restart(false)
      {}

      // Initializes ResponseParameters using the contents of file \c filename
      void read_file(const std::string& filename)
      {
         std::ifstream f(filename.c_str());
         read(f);
      }

      // Initializes ResponseParameteres using the contents of stream \c f
      void read(std::istream& f)
      {
         position_stream(f, "response");
         std::string s;
         xc_data = "hf";
         protocol_data = madness::vector_factory(1e-4, 1e-6);

         while(f >> s)
         {
            if(s == "end")
            {
               break;
            }
            else if (s == "archive")
            {
               f >> archive;
            }
            else if (s == "nwchem")
            {
               f >> nwchem;
            }
            else if (s == "restart")
            {
               restart = true;
               f >> resp_archive;
            }
            else if (s == "states")
            {
               f >> states;
            }
            else if (s == "print_level")
            {
               f >> print_level;
            }
            else if (s == "tda")
            {
               tda = true;
            }
            else if (s == "larger_subspace")
            {
               f >> larger_subspace;
            }
            else if (s == "k")
            {
               f >> k;
            }
            else if (s == "plot")
            {
               plot = true;
               std::string buf;
               std::getline(f,buf);
               plot_data = std::vector<int>();
               std::string d;
               std::stringstream t(buf);
               t >> d;
               if (d == "range")
               {
                  plot_range = true;
                  t >> d;
                  plot_data.push_back(std::stoi(d));
                  t >> d;
                  for(int z = plot_data[0]+1; z < std::stoi(d); z++) plot_data.push_back(z);
               }
               else
               {
                  // Add in the one we just read
                  plot_data.push_back(std::stoi(d));
                  while(t >> d)  plot_data.push_back(std::stoi(d));
               }
            }
            else if (s == "plot_pts")
            {
               f >> plot_pts;
            }
            else if (s == "plot_L")
            {
               f >> plot_L;
            }
            else if (s == "max_iter")
            {
               f >> max_iter;
            }
            else if (s == "dconv")
            {
               f >> dconv;
            }
            else if (s == "small")
            {
               f >> small;
            }
            else if (s == "plot_initial")
            {
               plot_initial = true;
            }
            else if (s == "random")
            {
               random = true;
            }
            else if (s == "store_potential")
            {
               store_potential = true;
            }
            else if (s == "localized")
            {
               localized = true;
            }
            else if (s == "range")
            {
               e_window = false;
               f >> range_low;
               f >> range_high;
            }
            else if (s == "xc")
            {
               char buf[1024];
               f.getline(buf,sizeof(buf));
               xc_data = buf;
            }
            else if (s == "protocol")
            {
               std::string buf;
               std::getline(f,buf);
               protocol_data = std::vector<double>();
               double d;
               std::stringstream s(buf);
               while (s >> d) protocol_data.push_back(d);
            }
            else 
            {
               std::cout << "response: unrecognized input keyword " << s << std::endl;
               MADNESS_EXCEPTION("input error", 0); 
            }
         }
      } // end read()

      // Prints all information
      void print_params() const
      {
         madness::print("\n   Input Response Parameters");
         madness::print("   -------------------------");
         madness::print("                XC Functional:", xc_data);
         madness::print("            Ground State File:", archive);
         if(nwchem != "") madness::print("                  NWChem File:", nwchem);
         madness::print("             States Requested:", states);
         madness::print("            TDA Approximation:", tda);
         madness::print("           Localized Orbitals:", localized);
         madness::print("                Energy Window:", e_window, " (Not yet implemented)");
         if(e_window) madness::print("           Energy Range Start:", range_low);
         if(e_window) madness::print("             Energy Range End:", range_high);
         if(k>0) madness::print("                            k:", k);
         madness::print("     Use Random Initial Guess:", random);
         madness::print("              Store Potential:", store_potential);
         madness::print("               Max Iterations:", max_iter);
         madness::print("   Larger Subspace Iterations:", larger_subspace);
         madness::print("Density Convergence Threshold:", dconv);
         madness::print("                     Protocol:", protocol_data);
         if(plot_initial) madness::print("        Plot Initial Orbitals:", plot_initial);
         if(plot) madness::print("          Plot Final Orbitals:", plot);
         if(plot and plot_pts != 201) madness::print("          Plot Num. of Points:", plot_pts);
         if(plot and plot_L > 0.0) madness::print("                Plot Box Size:", plot_L);
         if(plot and plot_range) madness::print("                   Plot Start:", plot_data[0]);
         if(plot and plot_range) madness::print("                     Plot End:", plot_data.back());
         if(plot and not plot_range) madness::print("       Orbitals to be Plotted:", plot_data);
         madness::print("                  Print Level:", print_level);
      }
   };

} // namespace madness
#endif
