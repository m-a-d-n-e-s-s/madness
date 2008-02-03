#include <examples/molecularbasis.h>

using namespace madness;


int main() {
    AtomicBasisSet g("sto-3g");

    g.print_all();

    g.print(Molecule("input"));

    return 0;
}
