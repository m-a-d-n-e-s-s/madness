#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/world/world_object.h>

using namespace std;
using namespace madness;

class Array : public WorldObject<Array> {
    vector<double> v;
public:
    /// Make block distributed array with size elements
    Array(World& world, size_t size) 
        : WorldObject<Array>(world), v((size-1)/world.size()+1)
    {
        process_pending();
    };

    /// Return the process in which element i resides
    ProcessID owner(size_t i) const {return i/v.size();};

    Future<double> read(size_t i) const {
        if (owner(i) == get_world().rank())
            return Future<double>(v[i-get_world().rank()*v.size()]);
        else
            return send(owner(i), &Array::read, i);
    };

    void write(size_t i, double value) {
        if (owner(i) == get_world().rank())
            v[i-get_world().rank()*v.size()] = value;
        else
            send(owner(i), &Array::write, i, value);
    };
};

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    Array a(world, 10000), b(world, 10000);

    // Without regard to locality, initialize a and b
    for (int i=world.rank(); i<10000; i+=world.size()) {
        a.write(i, 10.0*i);
        b.write(i,  7.0*i);
    }
    world.gop.fence();

    // All processes verify 100 random values from each array
    for (int j=0; j<100; j++) {
        size_t i = world.rand()%10000;
        Future<double> vala = a.read(i);
        Future<double> valb = b.read(i);
        // Could do work here until results are available
        MADNESS_ASSERT(vala.get() == 10.0*i);
        MADNESS_ASSERT(valb.get() ==  7.0*i);
    }
    world.gop.fence();

    if (world.rank() == 0) print("OK!");
    finalize();
    return 0;
}
