#include <world/world.h>
using namespace madness;

class A {
public:
    void member() {
        PROFILE_MEMBER_FUNC(A);
    }
};

void b() {
    PROFILE_FUNC;

    A a;
    a.member();
}

void a() {
    PROFILE_FUNC;

    b();
}

int main(int argc, char** argv) {
    initialize(argc, argv);

    World world(MPI::COMM_WORLD);
    world.args(argc,argv);

    {
        PROFILE_BLOCK(main);
        a();
    }

    ThreadPool::end();
    print_stats(world);
    finalize();
    return 0;
}
