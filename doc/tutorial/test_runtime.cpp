#include <madness/world/worldgop.h>
#include <cassert>

int task (int i) { return i + 1; }

int main(int argc, char* argv[]) {
  using namespace madness;

  World &world = initialize(argc, argv);
  const auto me = world.rank();
  const auto nproc = world.size();

  if (me == 0) {
    for (auto t = 0; t != 1000; ++t) {
      const auto proc = t % nproc;
      auto f = world.taskq.add(proc, &task, t);
      assert(f.get() == t+1);
    }
  }
  world.gop.fence();

  finalize();
  return 0;
}

