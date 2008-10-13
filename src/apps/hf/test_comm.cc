#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>

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

    /// Read element i
    Future<double> read(size_t i) const {
        if (owner(i) == world.rank())
            return Future<double>(v[i-world.rank()*v.size()]);
        else
            return send(owner(i), &Array::read, i);
    };

    /// Write element i
    Void write(size_t i, double value) {
        if (owner(i) == world.rank())
            v[i-world.rank()*v.size()] = value;
        else
            send(owner(i), &Array::write, i, value);
        return None;
    };
};


int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    madness::World world(MPI::COMM_WORLD);

    Array a(world, 10000), b(world, 10000);

    // Without regard to locality, initialize a[i]=i*10, b[i]=i*7
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
    MPI::Finalize();
}

