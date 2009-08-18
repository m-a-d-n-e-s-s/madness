#include <world/world.h>
#include <string>

using namespace madness;
using namespace std;

class Fred {
public:
    string a(const string& input) const {
        return input + string("a");
    }

    string b(const string& input) const {
        return input + string("b");
    }
};


int main(int argc, char** argv) {
    madness::initialize(argc,argv);
    madness::World world(MPI::COMM_WORLD);

    TaskAttributes attr;
    attr.set_stealable(true);

    Fred fred;
    
    Future<string> r = world.taskq.add(fred, &Fred::a, string("Hello"),attr);
    Future<string> s = world.taskq.add(fred, &Fred::b, r, attr);

    Future<string> fff;
    Future<string> ggg = world.taskq.add(fred, &Fred::a, fff, attr);

    world.taskq.steal(100);

    fff.set("done");

    print(s.get(), ggg.get());

    madness::finalize();
    return 0;
}
