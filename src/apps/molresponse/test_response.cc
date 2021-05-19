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

bool test_derived(World& world) {
  print("entering test_derived");
  std::string inputlines = R"input(mp3
			archive restartdata
			dipole True
			#dconv 1.e-4
			maxiter 12# asd
			end)input";
  inputfile ifile("input1", inputlines);

  ResponseParameters param;
  param.read_and_set_derived_values(world, "input1", "mp3");
  param.print();

  // test_same(param.econv(), 1.e-4);
  // test_same(param.dconv(), sqrt(param.econv()) * 0.1);
  return true;
}
bool test_read_archive_set_states(World& world) {
  print("entering test_derived");
  std::string inputlines = R"input(dipole_test
			archive restartdata
			dipole True
			#dconv 1.e-4
			maxiter 12# asd
			end
			nuclear_test
			nuclear True
			end
			order2_dd
			print_level 10
			order2 True
			d2_types dd
			end
			order2_dn
			order2 True
			d2_types dn
			end
			order2_nn
			order2 True
			d2_types nn
			end
			)input";
  inputfile ifile("input1", inputlines);

  GroundParameters g_params;
  ResponseParameters param1;
  std::string ground_file = param1.archive();

  g_params.read(world, ground_file);
  g_params.print_params();
  param1.read_and_set_derived_values(world, "input1", "dipole_test");
  param1.print();
  density_vector d1 = set_density_type(world, param1, g_params);
  d1.PrintDensityInformation();

  ResponseParameters param2;
  param2.read_and_set_derived_values(world, "input1", "nuclear_test");
  param2.print();
  density_vector d2 = set_density_type(world, param2, g_params);
  d2.PrintDensityInformation();

  ResponseParameters param3;
  param3.read_and_set_derived_values(world, "input1", "order2_dd");
  param3.print();
  density_vector d3 = set_density_type(world, param3, g_params);
  d3.PrintDensityInformation();

  ResponseParameters param4;
  param4.read_and_set_derived_values(world, "input1", "order2_dn");
  param4.print();

  ResponseParameters param5;
  param5.read_and_set_derived_values(world, "input1", "order2_nn");
  param5.print();

  // test_same(param.econv(), 1.e-4);
  // test_same(param.dconv(), sqrt(param.econv()) * 0.1);
  return true;
}

bool test_create_density(World& world) {
  print("entering test_derived");
  std::string inputlines = R"input(dipole
			archive restartdata
			dipole True
			#dconv 1.e-4
			maxiter 12# asd
			end
			)input";

  inputfile ifile("input1", inputlines);

  ResponseParameters param1;
  param1.read_and_set_derived_values(world, "input1", "dipole_test");
  param1.print();

  // test_same(param.econv(), 1.e-4);
  // test_same(param.dconv(), sqrt(param.econv()) * 0.1);
  return true;
}

int main(int argc, char** argv) {
  // Initialize MADNESS mpi
  initialize(argc, argv);
  int success = 0;
  {
    //{  // limite lifetime of world so that finalize() can execute cleanly
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv, true);
    print_meminfo(world.rank(), "startup");

    ResponseParameters r_params;
    r_params.print();

    const std::string ground_file = "restartdata";
    GroundParameters g_params;
    g_params.read(world, ground_file);
    g_params.print_params();

    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3> >(world)));

    std::cout.precision(6);
    // This makes a default input file name of 'input'
    try {
      test_derived(world);
      test_read_archive_set_states(world);
    } catch (const SafeMPI::Exception& e) {
      print(e);
      error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
      print(e);
      error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
      print(e);
      error("caught a Tensor exception");
    } catch (const char* s) {
      print(s);
      error("caught a string exception");
    } catch (const std::string& s) {
      print(s);
      error("caught a string (class) exception");
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
