#include <chem/SCF.h>
#include <madness/world/worldmem.h>
#include <molresponse/density.h>
#include <molresponse/ground_parameters.h>
#include <molresponse/response_parameters.h>
#include <stdlib.h>

#include "TDDFT.h" // All response functions/objects enter through this
//#include "molresponse/density.h"
#include "molresponse/global_functions.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static inline int file_exists(const char* inpname) {
  struct stat buffer;
  size_t rc = stat(inpname, &buffer);
  return (rc == 0);
}
#endif

template <typename T> void test_same(const T& t1, const T& t2) {
  if (t1 != t2) {
    print("t1, t2", t1, t2);
    using madness::operators::operator<<;
    std::cout << "++" << t1 << "++" << std::endl;
    std::cout << "++" << t2 << "++" << std::endl;

    throw std::runtime_error("failure in test");
    ;
  }
}

struct inputfile {
  std::string fname;
  inputfile(const std::string filename, const std::string lines) {
    fname = filename;
    std::ofstream myfile;
    myfile.open(fname);
    myfile << lines << std::endl;
    myfile.close();
  }

  ~inputfile() { remove(fname.c_str()); }
};
void run_density(World& world, density_vector& rho) {
  // Create the TDDFT object
  print("Computing Density");
  TDDFT calc(world, rho);
  if (calc.r_params.response_type().compare("excited_state") == 0) {
    print("Entering Excited State Response Runner");
    calc.solve_excited_states(world);
  } else {
    print("Entering Frequency Response Runner");
    calc.compute_freq_response(world);
  }
  //
  // densityTest.PlotResponseDensity(world);

  if (calc.r_params.response_type().compare("dipole") == 0) { //
    print("Computing Alpha");
    Tensor<double> alpha = rho.ComputeSecondOrderPropertyTensor(world);
    print("Second Order Analysis");
    rho.PrintSecondOrderAnalysis(world, alpha);
  }
}
void read_and_create_density(World& world, inputfile ifile, std::string tag) {
  GroundParameters g_params;
  ResponseParameters r_params;
  r_params.read_and_set_derived_values(world, "input1", tag);
  std::string ground_file = r_params.archive();
  g_params.read(world, ground_file);
  g_params.print_params();
  g_params.molecule().print();
  r_params.print();
  density_vector d1 = set_density_type(world, r_params, g_params);
  d1.PrintDensityInformation();
  run_density(world, d1);
}
bool test_create_dipole_save(World& world) {
  print("entering test_dipole");
  std::string inputlines = R"input(dipole_test
			archive restartdata
			first_order True
			dipole True
			save_density True
			save_density_file "restart"
			maxiter 5# asd
			end
			)input";
  inputfile ifile("input1", inputlines);
  read_and_create_density(world, ifile, "dipole_test");
  return true;
}
bool test_create_nuclear(World& world) {
  print("entering test_create_nuclear");
  std::string inputlines = R"input(nuclear_test
			archive restartdata
			nuclear True
			end
			)input";
  inputfile ifile("input1", inputlines);
  read_and_create_density(world, ifile, "nuclear_test");

  return true;
}
bool test_create_order2_dd(World& world) {
  print("entering test_create_order2_dd");
  std::string inputlines = R"input(order2_dd
			archive restartdata
			order2 True
			d2_types dd
			end
			)input";
  inputfile ifile("input1", inputlines);
  read_and_create_density(world, ifile, "order2_dd");
  return true;
}
bool test_create_order2_dn(World& world) {
  print("entering test_create_order2_dn");
  std::string inputlines = R"input(order2_dn
			archive restartdata
			order2 True
			d2_types dn
			end
			)input";
  inputfile ifile("input1", inputlines);
  read_and_create_density(world, ifile, "order2_dn");
  return true;
}

int main(int argc, char** argv) {
  // Initialize MADNESS mpi
  initialize(argc, argv);
  int success = 0;
  {
    //{  // limite lifetime of world so that finalize() can execute cleanly
    World world(SafeMPI::COMM_WORLD);
    molresponse::start_timer(world);
    startup(world, argc, argv, true);
    print_meminfo(world.rank(), "startup");
    std::cout.precision(6);
    //      FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3> >(world)));
    test_create_dipole_save(world);
    // test_create_nuclear(world);
    // test_create_order2_dd(world);
    // test_create_order2_dn(world);
    if (world.rank() == 0)
      printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    world.gop.fence();
  }
  finalize();

  return success;
}
