//
// Created by adrianhurtado on 2/3/22.
//

#include "FrequencyResponse.hpp"

#include "property.h"

void FrequencyResponse::initialize(World &world) {

    X_space PQ;



    if (r_params.dipole()) {
        if (world.rank() == 0) print("creating dipole property operator");
        PQ.X = PropertyRHS(world, DipoleVector(world));
        PQ.Y = PQ.X.copy();
        print("P: ", PQ.X.norm2());
        print("Q: ", PQ.Y.norm2());
    } else if (r_params.nuclear()) {
        if (world.rank() == 0) print("creating nuclear property operator");
        PQ.X = PropertyRHS(world, NuclearVector(world, molecule));
        PQ.Y = PQ.X.copy();
    }
    print("Pre iteration Information");
    print("Number of Response States: ", r_params.num_states());
    print("Number of Ground States: ", r_params.num_orbitals());
    print("k = ", FunctionDefaults<3>::get_k());
    print("protocol threshold = ", FunctionDefaults<3>::get_k());

    print("Property rhs func k = ", PQ.X[0][0].k());
    print("Property func k thresh= ", PQ.X[0][0].thresh());

    print("Property rhs func Q k = ", PQ.Y[0][0].k());
    print("Property func Q k thresh = ", PQ.Y[0][0].thresh());

    print("Property rhs func P norms", PQ.X.norm2());
    print("Property rhs func Q norms", PQ.Y.norm2());
}
//
// Here i should print some information about the calculation we are
// about to do
