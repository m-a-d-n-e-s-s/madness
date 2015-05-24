#include <madness/world/MADworld.h>

// This program is used to do a simple test of the task queue.

const int NGEN=100;
const int NTASK=100000;
madness::AtomicInt total_count;

// Global variables for thread local storage
pthread_key_t thread_key;
madness::AtomicInt thread_index;
unsigned long* thread_counters = nullptr;


// Access thread local storage
unsigned long* get_tls() {
    unsigned long* thread_counter =
            reinterpret_cast<unsigned long*>(pthread_getspecific(thread_key));
    if(! thread_counter) {
        thread_counter = thread_counters + (thread_index++);
        pthread_setspecific(thread_key, thread_counter);
    }

    return thread_counter;
}

// Initialize thread local storage
void init_tls(const unsigned long nthreads) {
    thread_index = 0;
    thread_counters = new unsigned long[nthreads];
    for(unsigned int i = 0; i < nthreads; ++i)
        thread_counters[i] = 0;
    pthread_key_create(&thread_key, nullptr);

    get_tls();
}

// Cleanup thread local storage
void cleanup_tls() {
    delete [] thread_counters;
    pthread_key_delete(thread_key);
}

class Task : public madness::TaskInterface {
    const int gen;
    bool done;
public:
    Task(int gen) : gen(gen), done(false) {}

    virtual void run(madness::World& world) {
        MADNESS_ASSERT(!done);

        unsigned long* thread_counter = get_tls();

        ++(*thread_counter);
        total_count++;
        if (gen > 1) world.taskq.add(new Task(gen-1));

        done = true;
    }

    static bool finished() {return total_count==(NGEN*NTASK);}
};

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    madness::World world(SafeMPI::COMM_WORLD);

    total_count = 0;
    init_tls(madness::ThreadPool::size() + 1);

    // Get start time.
    double start = madness::wall_time();

    // Add tasks to the task queue
    for (int i=0; i<NTASK; ++i)
        world.taskq.add(new Task(NGEN));

    // Wait for tasks to finish
    world.await(& Task::finished);

    // Get finish time.
    double finish = madness::wall_time();


    std::cout << "Total tasks = " << total_count
            << "\nTotal runtime = " << finish - start
            << " (s)\nTasks per thread:\n";
    for (unsigned long i = 0; i < (madness::ThreadPool::size() + 1); ++i)
        std::cout << i << " " << thread_counters[i] << "\n";

    cleanup_tls();
    madness::finalize();

    return 0;
}
