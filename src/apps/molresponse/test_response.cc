#include <chem/SCF.h>
#include <madness/world/worldmem.h>
#include <molresponse/density.h>
#include <molresponse/ground_parameters.h>
#include <molresponse/response_parameters.h>
#include <stdlib.h>

#include "TDDFT.h"  // All response functions/objects enter through this
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

template <typename T>
void test_same(const T& t1, const T& t2) {
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
void read_and_create_density(World& world, inputfile ifile, std::string tag) {
  GroundParameters g_params;
  ResponseParameters param1;
  param1.read_and_set_derived_values(world, "input1", tag);
  std::string ground_file = param1.archive();
  g_params.read(world, ground_file);
  g_params.print_params();
  param1.print();
  density_vector d1 = set_density_type(world, param1, g_params);
  d1.PrintDensityInformation();
}
bool test_create_dipole(World& world) {
  print("entering test_dipole");
  std::string inputlines = R"input(dipole_test
			archive restartdata
			dipole True
			#dconv 1.e-4
			maxiter 12# asd
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
    try {
      test_create_dipole(world);
      test_create_nuclear(world);
      // test_create_order2_dd(world);
      // test_create_order2_dn(world);
    } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
      success = 1;
    } catch (const madness::MadnessException& e) {
      print(e);
      error("caught a MADNESS exception");
      success = 1;
    } catch (const madness::TensorException& e) {
      print(e);
      error("caught a Tensor exception");
      success = 1;
    } catch (const char* s) {
      print(s);
      error("caught a string exception");
      success = 1;
    } catch (const std::string& s) {
      print(s);
      error("caught a string (class) exception");
      success = 1;
    } catch (std::exception& e) {
      print("\n\tan error occurred .. ");
      print(e.what());
      success = 1;
    } catch (...) {
      print("\n\tan unknown error occurred .. ");
      success = 1;
    }
    if (world.rank() == 0) printf("\nfinished at time %.1fs\n\n", wall_time());
    world.gop.fence();
    world.gop.fence();
  }
  finalize();

  return success;
}
