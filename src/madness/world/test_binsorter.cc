#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/binsorter.h>

using namespace madness;

typedef std::pair<int,double> valueT;

// Too lazy to wrap this stuff into a class
int P;
int me;
double local_sorted_sum;

void inserter(const std::pair<int,double>& t) {
    // For testing verify destination ... normally don't need to pass owner in the value
    if ((t.first%P) != me) throw "bad index";
    // Checksum on values ... normally would insert into local container
    local_sorted_sum += t.second;
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(SafeMPI::COMM_WORLD);
    world.args(argc,argv);

    P = world.size();
    me = world.rank();
    local_sorted_sum = 0.0;
    const unsigned int N = 10000*P;

    double local_sum = 0.0;
    
    // Deliberately use a small binsize to stress flushing and messaging
    BinSorter<valueT,void(*)(const valueT&)> sorter(world,inserter,11);

    // Compute checksum locally and remotely to verify results
    for (unsigned int i=0; i<N; i++) {
        const ProcessID owner = (9973u*i)%P; // crude multiplicative ran# generator
        const double value = 1.0/(i+1);
        sorter.insert(owner, valueT(owner,value));
        local_sum += value;
    }

    sorter.finish();
    
    world.gop.sum(local_sum);

    world.gop.sum(local_sorted_sum);

    if (world.rank() == 0) print(local_sum, local_sorted_sum);

    bool OK = (std::abs(local_sum-local_sorted_sum) < 1e-14*local_sum);

    if (world.rank() == 0) print("OK?", OK);

    world.gop.fence();
    finalize();

    return 0;
}
    

