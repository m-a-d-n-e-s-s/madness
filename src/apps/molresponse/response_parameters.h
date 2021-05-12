// Copyright 2021 Adrian Hurtado

/// \file ResponseParameters
/// \brief Input parameters for a response calculation.

#ifndef SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_
#define SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_

#include <chem/QCCalculationParametersBase.h>
#include <chem/molecule.h>
#include <chem/xcfunctional.h>
#include <madness/mra/mra.h>
#include <madness/world/parallel_archive.h>

#include <functional>
#include <numeric>
#include <string>
#include <vector>

using namespace madness;

class ResponseParameters : public QCCalculationParametersBase {
 public:
  ResponseParameters() : QCCalculationParametersBase() {
    // initialize with: key, value, comment (optional), allowed values (optional)
    initialize<std::string>("archive", "restartdata", "file to read ground parameters from");
    initialize<std::string>("nwchem", "", "Root name of nwchem files for intelligent starting guess");
    initialize<size_t>("states", 1, "Number of excited states requested");
    initialize<int>("print_level", 3, "0: no output; 1: final energy; 2: iterations; 3: timings; 10: debug");
    initialize<bool>("tda", false, "turn on Tam-Danchof approximation (excitations energy");
    initialize<bool>("plot", false, "turn on plotting of final orbitals. Output format is .vts");
    initialize<bool>("plot_range", false, "controls which orbitals will be plotted");
    initialize<std::vector<int>>("plot_data", std::vector<int>{0}, "Orbitals to plot");
    initialize<std::vector<double>>(
        "plot_cell", std::vector<double>(), "lo hi in each dimension for plotting (default is all space)");
    initialize<double>("plot_L", -1.0, "Controls the plotting box size");
    initialize<size_t>("plot_pts", 201, "Controls number of points in plots");
    initialize<bool>("plot_all_orbitals", false, "Turn on 2D plotting of response orbitals ");

    initialize<size_t>("maxiter", 25, "maximum number of iterations");

    initialize<double>("dconv", 3.e-4, "recommended values: 1.e-4 < dconv < 1.e-8");
    initialize<bool>("dconv_set", false, "Convergence flage for the orbtial density");

    initialize<bool>("dconv_set", false, "Convergence flage for the orbtial density");

    initialize<bool>("guess_xyz", false, "TODO : check what this is for");

    initialize<double>("small", 1.e10, "smallest length scale we need to resolve");
    initialize<std::vector<double>>("protocol_data", {1.e-4, 1.e-6}, "calculation protocol");

    initialize<size_t>(
        "larger_subspace", 0, "Number of iterations to diagonalize in a subspace consisting of old and new vectors");
    initialize<int>("k", 7, "polynomial order");

    initialize<bool>("random", false, "Use random guess for initial response functions");
    initialize<bool>("store_potential", true, "Store the potential instead of computing each iteration");
    initialize<bool>("e_range", false, "Use an energy range to excite from");
    initialize<double>("e_range_lo", 0, "Energy range (lower end) for orbitals to excite from");
    initialize<double>("e_range_hi", 1, "Energy range (upper end) for orbitals to excite from");
    initialize<bool>("plot_initial", false, "Flag to plot the ground state orbitals read in from archivie");
    // Restart Parameters
    initialize<bool>("restart", false, "Flag to restart scf loop from file");
    initialize<std::string>("restart_file", "", "file to read ground parameters from");
    // kain
    initialize<bool>("kain", false, "Turn on Krylov Accelarated Inexact Newton Solver");
    initialize<double>("maxrotn", 1.0, "Max orbital rotation per iteration");
    initialize<size_t>("maxsub", 10, "size of iterative subspace ... set to 0 or 1 to disable");
    initialize<std::string>("xc", "hf", "XC input line");
    initialize<bool>("save", false, "if true save orbitals to disk");
    initialize<std::string>("save_file", "", "File name to save orbitals for restart");
    initialize<bool>("save_density", false, "Flag to save density at each iteration");
    initialize<std::string>("save_density_file", "", "File name to save density for restart");
    initialize<bool>("load_density", false, "Flag to load density for restart");
    initialize<std::string>("load_density_file", "", "File name to load density for restart");
    initialize<size_t>("guess_max_iter", 5, "maximum number of guess iterations");
    // properties
    initialize<bool>("property", false, "Flag to turn on frequency dependent property calc");
    initialize<std::string>("response_type", "excited_state", "dipole,nuclear,order2,order3");
    initialize<bool>("dipole", false, "Flag to turn on frequency dependent property calc");
    initialize<bool>("nuclear", false, "Flag to turn on frequency dependent property calc");
    initialize<bool>("order2", false, "Flag to turn on frequency dependent property calc");
    initialize<bool>("order3", false, "Flag to turn on frequency dependent property calc");
    initialize<vector<std::string>>("response_types", {"", "", ""}, "response types make up order2 or order3 ");
    initialize<double>("omega", 0.0, "Incident energy for dynamic response");
  }

  void read_and_set_derived_values(World &world, std::string inputfile, std::string tag) {
    // read the parameters from file and brodcast
    // tag
    read(world, inputfile, tag);
    set_derived_value("dconv", sqrt(get<double>("econv")) * 0.1);
  }

  // convenience getters
  double econv() const { return get<double>("econv"); }
  double dconv() const { return get<double>("dconv"); }
  bool localize() const { return get<bool>("localize"); }
  std::string local() const { return get<std::string>("local"); }
  std::pair<std::string, double> ncf() const { return get<std::pair<std::string, double>>("ncf"); }
  int maxiter() const { return get<int>("maxiter"); }
};

namespace madness {

struct ResponseParameters {
  // List of input parameters
  std::string archive;  ///< Name of input archive to read in ground state
  std::string nwchem;   ///< Root name of nwchem files for intelligent starting guess
  size_t states;        ///< Number of excited states requested
  size_t print_level;   ///< Controls the amount and style of printing. Higher
                        ///< values prsize_t more
                        ///<   Values |   What gets printed
                        ///<   ----------------------------
                        ///<     1    |   Print out timing of each step in the
                        ///<     calculation,
                        ///<          |   along with energy, energy res., and
                        ///<          func. res.
                        ///<   ----------------------------
                        ///<     2    |   Debug level. Prints EVERYTHING!!!

  bool tda;                           ///< Turn on Tam-Danchof approximation (only calculate excitations)
  bool plot;                          ///< Turn on plotting of final orbitals. Output format is .vts
  bool plot_range;                    ///< Controls which orbitals will be plotted
  std::vector<int> plot_data;         ///< Orbitals to plot
  double plot_L;                      ///< Controls the plotting box size
  size_t plot_pts;                    ///< Controls number of points in plots
  bool plot_all_orbitals;             // Turn on 2D plotting of response orbitals
  size_t max_iter;                    ///< Maximum number of iterations
  double dconv;                       ///< Convergence criterion for the orbital density
  bool dconv_set;                     ///< Convergence flag for the orbital density
  bool guess_xyz;                     ///< Convergence flag for the orbital density
  double small;                       ///< Minimum length scale to be resolved
  std::vector<double> protocol_data;  ///< Different thresholds for truncation
  size_t larger_subspace;             ///< Number of iterations to diagonalize in a
                                      ///< subspace consisting of old and new vectors
  size_t k;                           ///< Polynomial order to use in calculation
  bool random;                        ///< Use a random guess for initial response functions
  bool store_potential;               ///< Store the potential instead of computing each
                                      ///< iteration
  bool e_window;                      ///< Use an energy window to excite from
  double range_low;                   ///< Energy range (lower end) for orbitals to excite from
  double range_high;                  ///< Energy range (upper end) for orbitals to excite from
  bool plot_initial;                  ///< Flag to plot the ground state orbitals read in from
                                      ///< archive
  bool restart;                       ///< Flag to restart from file
  std::string restart_file;           ///< Flag to restart from file
  bool kain;                          ///< Flag to use KAIN solver
  double maxrotn;
  size_t maxsub;          ///< How many previous iterations KAIN will store
  std::string xc;         ///< Controls the HF or DFT switch, as well as which DFT
                          ///< functional is used
  bool save;              ///< Controls if orbitals will be saved each iteration
  std::string save_file;  ///< Flag to save to file

  bool save_density;
  std::string save_density_file;  ///< Flag to save to file

  bool load_density;
  std::string load_density_file;  ///< Flag to save to file

  size_t guess_max_iter;  ///< Maximum number of iterations for guess functions

  // Start of properties
  bool property;  ///< Flag that this is a properties calculation

  std::string response_type;  //

  bool dipole;   ///< Flag that this is a properties calculation
  bool nuclear;  ///< Flag that this is a properties calculation
  bool order2;   ///< Flag that this is a properties calculation
  bool order3;   ///< Flag that this is a properties calculation

  vector<std::string> response_types;  // string holding response type 1

  double omega;  ///< Incident energy for polarizability

  // TESTING
  bool old;
  bool old_two_electron;
  // END TESTING

  // Used to broadcast data to all mpi ranks
  template <typename Archive>
  void serialize(Archive &ar) {
    ar &archive &nwchem &states &print_level &tda &plot &plot_range &plot_data &plot_L &plot_pts &plot_all_orbitals
        &max_iter &dconv &dconv_set &guess_xyz &small &protocol_data &larger_subspace &k &random &store_potential
            &e_window &range_low &range_high &plot_initial &restart &restart_file &kain &maxrotn &maxsub &xc &save
                &save_file &guess_max_iter &property &response_type &dipole &nuclear &order2 &order3 &response_types
                    &omega &old &old_two_electron;
  }

  // Default constructor
  ResponseParameters()
      : archive("restartdata"),
        nwchem(""),
        states(1),
        print_level(1),
        tda(false),
        plot(false),
        plot_range(false),
        plot_data(std::vector<int>{0}),
        plot_L(-1.0),
        plot_pts(201),
        plot_all_orbitals(false),
        max_iter(40),
        dconv(0),
        dconv_set(false),
        guess_xyz(false),
        small(1e-6),
        protocol_data(madness::vector_factory(1e-6, 1e-8)),
        larger_subspace(0),
        k(0),
        random(false),
        store_potential(true),
        e_window(false),
        range_low(0.0),
        range_high(1.0),
        plot_initial(false),
        restart(false),
        restart_file(""),
        kain(false),
        maxrotn(1.0),
        maxsub(10),
        xc("hf"),
        save(false),
        save_file(""),
        save_density(false),
        save_density_file(""),
        load_density(false),
        load_density_file(""),
        guess_max_iter(5),
        property(false),
        response_type(),
        dipole(false),
        nuclear(false),
        order2(false),
        order3(false),
        response_types({"", "", ""}),
        omega(0.0),
        old(true),
        old_two_electron(false) {}

  // Initializes ResponseParameters using the contents of file \c filename
  void read_file(const std::string &filename) {
    std::ifstream f(filename.c_str());
    read(f);
  }

  // Initializes ResponseParameteres using the contents of stream \c f
  void read(std::istream &f) {
    position_stream(f, "response");
    std::string s;
    xc = "hf";
    protocol_data = madness::vector_factory(1e-6);

    while (f >> s) {
      if (s == "end") {
        break;
      } else if (s == "archive") {
        f >> archive;
      } else if (s == "nwchem") {
        f >> nwchem;
      } else if (s == "restart") {
        restart = true;
        f >> restart_file;
      } else if (s == "print_level") {
        f >> print_level;
      } else if (s == "tda") {
        tda = true;
      } else if (s == "larger_subspace") {
        f >> larger_subspace;
      } else if (s == "k") {
        f >> k;
      } else if (s == "plot") {
        plot = true;
        std::string buf;
        std::getline(f, buf);
        plot_data = std::vector<int>();
        std::string d;
        std::stringstream t(buf);
        t >> d;
        if (d == "range") {
          plot_range = true;
          t >> d;
          plot_data.push_back(std::stoi(d));
          t >> d;
          for (int z = plot_data[0] + 1; z <= std::stoi(d); z++) plot_data.push_back(z);
        } else {
          // Add in the one we just read
          plot_data.push_back(std::stoi(d));
          while (t >> d) plot_data.push_back(std::stoi(d));
        }
      } else if (s == "plot_pts") {
        f >> plot_pts;
      } else if (s == "plot_L") {
        f >> plot_L;
      } else if (s == "plot_all_orbitals") {
        plot_all_orbitals = true;
      } else if (s == "max_iter") {
        f >> max_iter;
      } else if (s == "dconv") {
        f >> dconv;
        dconv_set = true;
      } else if (s == "small") {
        f >> small;
      } else if (s == "plot_initial") {
        plot_initial = true;
      } else if (s == "random") {
        random = true;
      } else if (s == "guess_xyz") {
        guess_xyz = true;
      } else if (s == "store_potential_off") {
        store_potential = false;
      } else if (s == "range") {
        e_window = false;
        f >> range_low;
        f >> range_high;
      } else if (s == "xc") {
        f >> xc;

        if (not(xc == "hf" or xc == "lda" or xc == "pbe0" or xc == "b3lyp")) {
          MADNESS_EXCEPTION("Unsupported DFT functional requested.", 0);
        }
      } else if (s == "protocol") {
        std::string buf;
        std::getline(f, buf);
        protocol_data = std::vector<double>();
        double d;
        std::stringstream s(buf);
        while (s >> d) protocol_data.push_back(d);
      } else if (s == "kain") {
        kain = true;
      } else if (s == "maxrotn") {
        f >> maxrotn;
      } else if (s == "maxsub") {
        f >> maxsub;
      } else if (s == "save") {
        save = true;
        f >> save_file;
      } else if (s == "save_density") {
        save_density = true;
        f >> save_density_file;
      } else if (s == "load_density") {
        load_density = true;
        f >> load_density_file;
      } else if (s == "guess_iter") {
        f >> guess_max_iter;
      } else if (s == "property") {
        property = true;
        f >> response_type;
        if (response_type.compare("dipole") == 0) {
          f >> omega;
          dipole = true;
        } else if (response_type.compare("nuclear") == 0) {
          f >> omega;
          nuclear = true;
        } else if (response_type.compare("2ndOrder") == 0) {
          f >> response_types[0];
          f >> response_types[1];
          f >> omega;
          order2 = true;
        } else if (response_type.compare("3rdOrder") == 0) {
          f >> response_types[0];
          f >> response_types[1];
          f >> response_types[2];
          f >> omega;
          order3 = true;
        } else if (response_type.compare("excited_state") == 0) {
          f >> states;
        } else {
          MADNESS_EXCEPTION("Not a an avaible response type", 0);
        }
      } else if (s == "old_two_electron") {
        // Use potential manager, for debugging. Remove
        // after it works
        old_two_electron = true;
      } else {
        std::cout << "response: unrecognized input keyword " << s << std::endl;
        MADNESS_EXCEPTION("input error", 0);
      }
    }
  }  // end read()
  // sets the number of states based on the property type and
  // number of molecules in the atoms
  // for order calculations the number of states is the
  // multiplication of states that the calculation is derived from
  void SetNumberOfStates(Molecule &molecule) {
    vector<std::string> calculation_type;
    vector<bool> calc_flags;
    if (dipole) {
      states = 3;
    } else if (nuclear) {
      states = 3 * molecule.natom();
    } else if (order2) {
      vector<int> nstates;  // states 1
      for (size_t i = 0; i < 2; i++) {
        if (response_types[i] == "dipole") {
          nstates.push_back(3);
        } else if (response_types[i] == "nuclear") {
          nstates.push_back(3 * molecule.natom());
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      }
      states = std::accumulate(nstates.begin(), nstates.end(), 1, std::multiplies<>());
    } else if (order3) {
      vector<int> nstates;  // states 1
      for (size_t i = 0; i < 3; i++) {
        if (response_types[i] == "dipole") {
          nstates.push_back(3);
        } else if (response_types[i] == "nuclear") {
          nstates.push_back(3 * molecule.natom());
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      }
      states = std::accumulate(nstates.begin(), nstates.end(), 1, std::multiplies<>());
    }
  }

  // Prints all information
  void print_params() const {
    madness::print("\n   Input Response Parameters");
    madness::print("   -------------------------");
    madness::print("       Response XC Functional:", xc);
    madness::print("            Ground State File:", archive);
    if (nwchem != "") madness::print("                  NWChem File:", nwchem);
    if (restart) madness::print("                 Restart File:", restart_file);
    madness::print("             States Requested:", states);
    if (!property) madness::print("            TDA Approximation:", tda);
    if (e_window and !property) madness::print("                Energy Window:", e_window, " (Not yet implemented)");
    if (e_window and !property) madness::print("           Energy Range Start:", range_low);
    if (e_window and !property) madness::print("             Energy Range End:", range_high);
    if (k > 0) madness::print("                            k:", k);
    if (!property and random) madness::print("                Initial Guess: Random");
    if (!property and !random and nwchem == "") madness::print("                Initial Guess: Solid Harmonics * MOs");
    madness::print("              Store Potential:", store_potential);
    madness::print("         Max Guess Iterations:", guess_max_iter);
    madness::print("               Max Iterations:", max_iter);
    if (!property) madness::print("   Larger Subspace Iterations:", larger_subspace);
    madness::print("                     Use KAIN:", kain);
    if (kain) madness::print("           KAIN Subspace Size:", maxsub);
    madness::print("                Save orbitals:", save);
    if (dconv != 0.0) madness::print("Density Convergence Threshold:", dconv);
    if (!property) madness::print("                     Protocol:", protocol_data);
    if (plot_initial and !property) madness::print("        Plot Initial Orbitals:", plot_initial);
    if (plot and !property) madness::print("          Plot Final Orbitals:", plot);
    if (plot and plot_pts != 201 and !property) madness::print("          Plot Num. of Points:", plot_pts);
    if (plot and plot_L > 0.0 and !property) madness::print("                Plot Box Size:", plot_L);
    if (plot and plot_range and !property) madness::print("                   Plot Start:", plot_data[0]);
    if (plot and plot_range and !property) madness::print("                     Plot End:", plot_data.back());
    if (plot and not plot_range and !property) madness::print("       Orbitals to be Plotted:", plot_data);
    madness::print("                  Print Level:", print_level);

    // Start of property printing
    vector<std::string> calculation_type;
    calculation_type.push_back("Linear Dipole Perturbation");
    calculation_type.push_back("Linear Nuclear Displacement");
    calculation_type.push_back("2nd Order");
    calculation_type.push_back("3rd Order");
    vector<bool> calc_flags;
    calc_flags.push_back(dipole);
    calc_flags.push_back(nuclear);
    calc_flags.push_back(order2);
    calc_flags.push_back(order3);

    for (std::size_t i = 0; i < calc_flags.size(); i++) {
      if (calc_flags[i]) {
        madness::print("                     Property: ", calculation_type[i]);
        madness::print("           Incident Frequency: ", omega);
      }
    }
  }
};  // namespace madness
}  // namespace madness
#endif  // SRC_APPS_MOLRESPONSE_RESPONSE_PARAMETERS_H_
