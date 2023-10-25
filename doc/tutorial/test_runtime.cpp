#include <madness/world/worldgop.h>

int task (int i) { return i + 1; }

int main(int argc, char* argv[]) {
  using namespace madness;

  World &world = initialize(argc, argv);
  const auto me = world.rank();
  const auto nproc = world.size();

  Future<int> f;
  if (me == 0) {
    f = world.taskq.add(&task, 0);
    for (auto t = 1; t != 1000; ++t) {
      const auto proc = t % nproc;
      f = world.taskq.add(
          proc, &task, f);
    }
  }
  world.gop.fence();
  if (f.get() != 1000)
    abort();

  finalize();
  return 0;
}

