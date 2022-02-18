//
// Created by adrianhurtado on 2/11/22.
//
#define CATCH_CONFIG_RUNNER
#include <xc.h>

#include "ExcitedResponse.hpp"
#include "ResponseExceptions.hpp"
#include "TDDFT.h"
#include "apps/chem/SCF.h"
#include "apps/external_headers/catch.hpp"
#include "apps/external_headers/tensor_json.hpp"
#include "madness/world/worldmem.h"
#include "response_functions.h"
#include "runners.hpp"
#include "string"
#include "timer.h"
#include "write_test_input.h"
#include "x_space.h"

#if defined(HAVE_SYS_TYPES_H) && defined(HAVE_SYS_STAT_H) && defined(HAVE_UNISTD_H)

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static inline int file_exists(const char *input_name) {
    struct stat buffer {};
    size_t rc = stat(input_name, &buffer);
    return (rc == 0);
}

#endif

int main(int argc, char *argv[]) {
    World &world = madness::initialize(argc, argv);
    int result = 0;
    world.gop.fence();
    startup(world, argc, argv);
    { result = Catch::Session().run(argc, argv); }

    return result;

    // print_meminfo(world.rank(), "startup");
    // std::cout.precision(6);
    // print_stats(world);
}

// Run a single chosen molecule and check results
TEST_CASE("Run MOLDFT/RESPONSE") {

    xc_func_type func;
    double rho[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    double sigma[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
    double exc[5];
    // everything is 1 here so i =1
    // all version stuff
    int i, vmajor, vminor, vmicro, func_id = 1;
    /* Get the libxc version */
    xc_version(&vmajor, &vminor, &vmicro);
    printf("Libxc version: %d.%d.%d\n", vmajor, vminor, vmicro);
    /* Initialize the functional */
    if (xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0) {
        fprintf(stderr, "Functional '%d' not found\n", func_id);
    }
    /* Evaluate the energy density, depending on the family */
    switch (func.info->family) {
        case XC_FAMILY_LDA:
            xc_lda_exc(&func, 5, rho, exc);
            break;
        case XC_FAMILY_GGA:
        case XC_FAMILY_HYB_GGA:
            xc_gga_exc(&func, 5, rho, sigma, exc);
            break;
    }
    /* Print out density and energy density per particle */
    for (i = 0; i < 5; i += 1) { printf("%lf %lf\n", rho[i], exc[i]); }

    print("The functional", func.info->name, " is ");
    switch(func.info->kind){
        case(XC_EXCHANGE):
            print("an exchange functional");
            break;
        case(XC_CORRELATION):
            print("a correlation functional");
            break;
        case(XC_EXCHANGE_CORRELATION):
            print("a exchange-correlation functional");
            break;
        case(XC_KINETIC):
            printf("a kinetic energy functional");
            break;
        default:
            printf(" of unkown kind");
            break;
    }

    /* Deallocate memory */
    xc_func_end(&func);
}