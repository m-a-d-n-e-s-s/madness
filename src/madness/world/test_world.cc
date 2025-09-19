/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#include <vector>
#include <numeric>
#include <algorithm>

#include <madness/world/MADworld.h>
#include <madness/world/world_object.h>
#include <madness/world/worlddc.h>

#if MADNESS_CATCH_SIGNALS
# include <csignal>
#endif

#ifdef HAVE_PARSEC
# include <parsec.h>
# ifdef PARSEC_HAVE_CUDA
#  include <cuda_runtime.h>
# endif
#endif

using namespace madness;
using namespace std;

void test0(World& world) {
    PROFILE_FUNC;
    const size_t n=555;
    char buf[n+2];
    buf[0] = buf[n+1] = 127;

    archive::BufferOutputArchive arout(buf+1,n-2);

    arout & 1L & 99.9;
    arout.store("hello",6);
    arout & 7L & 77.7;
    MADNESS_CHECK(buf[0]==127 && buf[n+1]==127);

    long i = 0l;
    double a = 0.0;

    archive::BufferInputArchive arin(buf+1,n-2);
    char s[8];
    s[0] = s[7] = 66;
    arin & i & a;
    MADNESS_CHECK(i==1 && a==99.9);
    arin.load(s+1,6);
    MADNESS_CHECK(s[0]==66 && s[1]=='h' && s[5]=='o' && s[7]==66);
    arin & i & a;
    MADNESS_CHECK(i==7 && a==77.7);

    if (world.rank() == 0) print("test0 (serialization to/from buf) seems to be working");
}

class B {
    long b;
public:
    B(long b=0) : b(b) {};
    void set(long value) {
        b=value;
    };
    long get() const {
        return b;
    };
    ~B() {
        print("B destructor");
    };
    template <typename Archive> void serialize(Archive &ar) {
        ar&b;
    }
};


#include <complex>
typedef std::complex<double> double_complex;

class TestTask : public TaskInterface {
public:
    void run(World& world) {
        print("Hi, I am running!");
    }
};

class TTT {
private:
    int state;
public:
    TTT() : state(0) {};

    static void fred() {
        print("Oops-a-daisy!");
    };

    static int mary() {
        return 99;
    };

    static int carl() {
        return 88;
    };

    static int dave(World* world) {
        return world->rank();
    };

    static int bert(int input) {
        return input+1;
    };

    static double_complex sara(double a, const double_complex& b) {
        return a*b;
    };

    static string kate(World* world, const string& msg, double d) {
        ostringstream s;
        s << "Process " << world->rank() << " says '" << msg << "' and "
        << d << " right back at you!";
        return s.str();
    };

    double jody(double a, double b, double c) {
        return a+b+c + state;
    };

    double hugh(vector< Future<int> >& a) {
        double sum = 0.0;
        for (int i=0; i<(int)a.size(); ++i) sum += a[i].get();
        return sum;
    };
};

double dumb(int a1, int a2, int a3, int a4, int a5, int a6, int a7) {
    return a1+a2+a3+a4+a5+a6+a7;
}



void test5(World& world) {
    PROFILE_FUNC;
    int nproc = world.size();
    ProcessID me = world.rank();
    ProcessID right = (me+1)%nproc;
    TaskInterface* task = new TestTask();
    task->inc();
    task->inc();
    world.taskq.add(task);
    print("added the task ... about to dec");
    task->dec();
    task->dec();
    print("fencing");
    world.gop.fence();
    print("done with fence");

    world.taskq.add(TTT::fred);
    print("MAKING MARY");
    Future<int> mary = world.taskq.add(TTT::mary);
    print("MADE MARY");
    Future<int> carl = world.taskq.add(right,TTT::carl);
    Future<int> dave = world.taskq.add(right,TTT::dave, &world);
    Future<int> bert_input;
    Future<int> bert = world.taskq.add(TTT::bert,bert_input);
    Future<double> sara1;
    Future<double_complex> sara2;
    Future<double_complex> sara = world.taskq.add(TTT::sara,sara1,sara2);
    Future<string> kate2;
    Future<double> kate3;
    Future<string> kate = world.taskq.add(TTT::kate,&world,kate2,kate3);
    Future<string> katy2;
    Future<double> katy3;
    Future<string> katy = world.taskq.add(TTT::kate,&world,&katy2,&katy3);
    Future<string> cute = world.taskq.add(right,TTT::kate,&world,string("Boo!"),-42.0);
    TTT ttt;
    Future<double> jody = world.taskq.add(ttt,&TTT::jody,1.0,2.0,3.0);
    Future<double> duh = world.taskq.add(me,dumb,0,1,2,3,4,5,6);
    print("done with making futs");

    bert_input.set(7);
    sara1.set(3.0);
    sara2.set(double_complex(2.1,1.2));
    kate2.set(string("Who's your daddy?"));
    kate3.set(3.14);
    katy2.set(string("Your momma"));
    katy3.set(6.28);

    vector< Future<int> > futv = future_vector_factory<int>(7);
    Future<double> hugh = world.taskq.add(ttt,&TTT::hugh,futv);
    for (int i=0; i<7; ++i) {
        print("assigning",i,futv[i]);
        futv[i].set(i);
    }

    print("about to fence again");
    world.gop.fence();
    print("finished fence again");

    MADNESS_CHECK(mary.probe());
    MADNESS_CHECK(carl.probe());
    MADNESS_CHECK(dave.probe());
    MADNESS_CHECK(bert.probe());
    MADNESS_CHECK(sara.probe());
    MADNESS_CHECK(kate.probe());
    MADNESS_CHECK(katy.probe());
    MADNESS_CHECK(cute.probe());
    MADNESS_CHECK(jody.probe());
    MADNESS_CHECK(hugh.probe());
    MADNESS_CHECK(duh.probe());

    MADNESS_CHECK(mary.get() == 99);
    MADNESS_CHECK(carl.get() == 88);
    MADNESS_CHECK(dave.get() == right);
    MADNESS_CHECK(bert.get() == 8);
    MADNESS_CHECK(hugh.get() == 21.0);
    MADNESS_CHECK(duh.get() == 21.0);
    print("Sara says",sara.get().real(),sara.get().imag());
    print("Kate says",kate.get());
    print("Katy says",katy.get());
    print("Cute says",cute.get());
    print("Jody says",jody.get());

    if (me == 0) print("test5 (tasks and futures) OK");
}

class TestBarrier : public TaskInterface {
    int count; // does not need to be volatile since barrier includes necessary memory fences
public:

    TestBarrier(const madness::TaskAttributes& attr)
        : TaskInterface(attr)
        , count(0)
    {
        print("Testing barrier with nthread", attr.get_nthread());
    }

#if defined(__INTEL_COMPILER) || defined(__PGI)
  using madness::TaskInterface::run;
#endif

    void run(World& world, const TaskThreadEnv& env) {
        // Using the barrier each thread takes turns to update
        // the shared counter.

        env.barrier();

        int nthread = env.nthread();
        int id = env.id();
        for (int i=0; i<100; ++i) {
            for (int p=0; p<nthread; ++p) {
                env.barrier();
                if (p == id) count += (p+1);
            }
        }
        env.barrier();
        if (id == 0)
            print("     result from sum", count, "expected", 100*nthread*(nthread+1)/2);
    }
};

class TimeBarrier : public TaskInterface {
public:

    TimeBarrier(const madness::TaskAttributes& attr)
        : TaskInterface(attr)
    {
        print("Timing barrier with nthread", attr.get_nthread());
    }

#if defined(__INTEL_COMPILER) || defined(__PGI)
  using madness::TaskInterface::run;
#endif

    void run(World& world, const TaskThreadEnv& env) {
        // Barrier a zillion times

		for (int i=0; i<1000000; ++i) {
	        env.barrier();
		}
    }
};


// test multithreaded tasks
void test_multi(World& world) {
    // Test the correctness and performance of the barrier
    for (unsigned int i=1; i<=ThreadPool::size()+1; ++i) {
        world.taskq.add(new TestBarrier(TaskAttributes::multi_threaded(i)));
        double start = cpu_time();
        world.taskq.add(new TimeBarrier(TaskAttributes::multi_threaded(i)));
        double used = cpu_time()-start;
        print("barrier took", used*10.0,"micro seconds per call");
        world.gop.fence();
    }
}


class Foo : public WorldObject<Foo> {
    int a;
    std::vector<double> dbuf_short_;
    std::vector<double> dbuf_long_;
public:
    Foo(World& world, int a)
            : WorldObject<Foo>(world)
            , a(a) {
      process_pending();
      dbuf_short_.reserve((world.nproc() > 1 ? RMI::max_msg_len() : 1024)/sizeof(double)-5);
      dbuf_long_.reserve((world.nproc() > 1 ? RMI::max_msg_len() : 1024)/sizeof(double)+5);

      // make sure values are integer so equality correctness test is robust with compiler optimization
      std::generate_n(std::back_inserter(dbuf_short_), dbuf_short_.capacity(), [&world]() { return world.rand(); } );
      std::generate_n(std::back_inserter(dbuf_long_), dbuf_long_.capacity(), [&world]() { return world.rand(); } );
    }

    virtual ~Foo() { }

    int get0() {
        return a;
    }
    int get1(int a1) {
        return a+a1;
    }
    int get2(int a1, char a2) {
        return a+a1+a2;
    }
    int get3(int a1, char a2, short a3) {
        return a+a1+a2+a3;
    }
    int get4(int a1, char a2, short a3, long a4) {
        return a+a1+a2+a3+a4;
    }
    int get5(int a1, char a2, short a3, long a4, short a5) {
        return a+a1+a2+a3+a4+a5;
    }
    double getbuf0(const std::vector<double>& buf) {
        return (double)a + std::accumulate(buf.begin(), buf.end(), 0.0);
    }

    int get0c() const {
        return a;
    }
    int get1c(int a1) const {
        return a+a1;
    }
    int get2c(int a1, char a2) const {
        return a+a1+a2;
    }
    int get3c(int a1, char a2, short a3) const {
        return a+a1+a2+a3;
    }
    int get4c(int a1, char a2, short a3, long a4) const {
        return a+a1+a2+a3+a4;
    }
    int get5c(int a1, char a2, short a3, long a4, short a5) const {
        return a+a1+a2+a3+a4+a5;
    }
    double getbuf0c(const std::vector<double>& buf) const {
        return (double)a + std::accumulate(buf.begin(), buf.end(), 0.0);
    }

    Future<int> get0f() {
        return Future<int>(a);
    }

    const std::vector<double>& dbuf_short() const { return dbuf_short_; }
    const std::vector<double>& dbuf_long() const { return dbuf_long_; }
    const std::vector<double>& dbuf() const { return dbuf_long(); }

    // ping-pong via AMs
    void ping_am(int from, int speed) {
      madness::print("got an AM ping from proc ", from, " speed=", speed);
      if (speed < 10)
        this->send(from, &Foo::pong_am, this->get_world().rank(), speed + 1);
    }
    void pong_am(int from, int speed) {
      madness::print("got an AM pong from proc ", from, " speed=", speed);
      if (speed < 10)
        this->send(from, &Foo::ping_am, this->get_world().rank(), speed + 1);
    }

    // ping-pong via tasks
    void ping(int from, int speed) {
      madness::print("got a ping from proc ", from, " speed=", speed);
      if (speed < 10)
        this->task(from, &Foo::pong, this->get_world().rank(), speed + 1);
    }
    void pong(int from, int speed) {
      madness::print("got a pong from proc ", from, " speed=", speed);
      if (speed < 10)
        this->task(from, &Foo::ping, this->get_world().rank(), speed + 1);
    }

};

void test6(World& world) {
    PROFILE_FUNC;
    uniqueidT id;
    {
      ProcessID me = world.rank();
      ProcessID nproc = world.nproc();
      world.srand(73); // Everyone needs the same seed
      Foo a(world, me * 100);
      id = a.id();
      MADNESS_CHECK(world.ptr_from_id<Foo>(id) && world.ptr_from_id<Foo>(id).value() == &a);
      MADNESS_CHECK(world.id_from_ptr(&a) && world.id_from_ptr(&a).value() == id);
      const auto dbuf_sum =
          std::accumulate(a.dbuf().begin(), a.dbuf().end(), 0.0);

      if (me == 0) {
        print(a.id());
        for (ProcessID p = 0; p < nproc; ++p) {
          MADNESS_CHECK(a.send(p, &Foo::get0).get() == p * 100);
          MADNESS_CHECK(a.task(p, &Foo::get0).get() == p * 100);

          MADNESS_CHECK(a.send(p, &Foo::get0f).get() == p * 100);
          MADNESS_CHECK(a.task(p, &Foo::get0f).get() == p * 100);

          MADNESS_CHECK(a.send(p, &Foo::get1, 1).get() == p * 100 + 1);
          MADNESS_CHECK(a.task(p, &Foo::get1, Future<int>(1)).get() ==
                        p * 100 + 1);

          MADNESS_CHECK(a.send(p, &Foo::get2, 1, 2).get() == p * 100 + 3);
          MADNESS_CHECK(a.task(p, &Foo::get2, 1, 2).get() == p * 100 + 3);

          MADNESS_CHECK(a.send(p, &Foo::get3, 1, 2, 3).get() == p * 100 + 6);
          MADNESS_CHECK(a.task(p, &Foo::get3, 1, 2, 3).get() == p * 100 + 6);

          MADNESS_CHECK(a.send(p, &Foo::get4, 1, 2, 3, 4).get() ==
                        p * 100 + 10);
          MADNESS_CHECK(a.task(p, &Foo::get4, 1, 2, 3, 4).get() ==
                        p * 100 + 10);

          MADNESS_CHECK(a.send(p, &Foo::get5, 1, 2, 3, 4, 5).get() ==
                        p * 100 + 15);
          MADNESS_CHECK(a.task(p, &Foo::get5, 1, 2, 3, 4, 5).get() ==
                        p * 100 + 15);

          MADNESS_CHECK(a.task(p, &Foo::getbuf0, a.dbuf()).get() ==
                        p * 100 + dbuf_sum);

          MADNESS_CHECK(a.send(p, &Foo::get0c).get() == p * 100);
          MADNESS_CHECK(a.task(p, &Foo::get0c).get() == p * 100);

          MADNESS_CHECK(a.send(p, &Foo::get1c, 1).get() == p * 100 + 1);
          MADNESS_CHECK(a.task(p, &Foo::get1c, 1).get() == p * 100 + 1);

          MADNESS_CHECK(a.send(p, &Foo::get2c, 1, 2).get() == p * 100 + 3);
          MADNESS_CHECK(a.task(p, &Foo::get2c, 1, 2).get() == p * 100 + 3);

          MADNESS_CHECK(a.send(p, &Foo::get3c, 1, 2, 3).get() == p * 100 + 6);
          MADNESS_CHECK(a.task(p, &Foo::get3c, 1, 2, 3).get() == p * 100 + 6);

          MADNESS_CHECK(a.send(p, &Foo::get4c, 1, 2, 3, 4).get() ==
                        p * 100 + 10);
          MADNESS_CHECK(a.task(p, &Foo::get4c, 1, 2, 3, 4).get() ==
                        p * 100 + 10);

          MADNESS_CHECK(a.send(p, &Foo::get5c, 1, 2, 3, 4, 5).get() ==
                        p * 100 + 15);
          MADNESS_CHECK(a.task(p, &Foo::get5c, 1, 2, 3, 4, 5).get() ==
                        p * 100 + 15);

          MADNESS_CHECK(a.task(p, &Foo::getbuf0c, a.dbuf()).get() ==
                        p * 100 + dbuf_sum);
        }
      } // me == 0

      for (ProcessID p = 0; p != nproc; ++p) {
        a.send(p, &Foo::ping_am, me, 1);
        a.task(p, &Foo::ping, me, 1);
      }

#ifdef MADNESS_WORLDOBJECT_FUTURE_TRACE
      for (ProcessID p = 0; p != nproc; ++p) {
        auto f = a.task(p, &Foo::get0);
        a.trace(f);
      }
#endif

      world.gop.fence();

#ifdef MADNESS_WORLDOBJECT_FUTURE_TRACE
      MADNESS_CHECK(a.trace_status_nfuture_registered() == (a.trace_futures() ? nproc : 0));
      MADNESS_CHECK(decltype(a)::trace_status_nfuture_assigned(a.id()) ==
                    (decltype(a)::trace_futures(a.id()) ? nproc : 0));
#endif

      // stress the large message protocol ... off by default
      if (0) {
        const auto dbuf_sum_long =
            std::accumulate(a.dbuf_long().begin(), a.dbuf_long().end(), 0.0);
        const auto dbuf_sum_short =
            std::accumulate(a.dbuf_short().begin(), a.dbuf_short().end(), 0.0);
#if 0 // uncomment to STRESS the large msg protocol
      const size_t nmsg = 128;
#else
        const size_t nmsg = 1;
#endif
        std::vector<Future<double>> results;
        std::vector<double> results_ref;
        for (size_t m = 0; m != nmsg; ++m) {
          for (ProcessID p = 0; p < nproc; ++p) {
            results.push_back(a.task(p, &Foo::getbuf0c, a.dbuf_long()));
            results_ref.push_back(p * 100 + dbuf_sum_long);
            results.push_back(a.task(p, &Foo::getbuf0c, a.dbuf_short()));
            results_ref.push_back(p * 100 + dbuf_sum_short);
          }
        }
        world.gop.fence();
        for (size_t r = 0; r != results.size(); r += 2) {
          MADNESS_CHECK(results[r].get() == results_ref[r]);
        }
      }
    }

    // test that the object is gone
    auto ptr_opt = world.ptr_from_id<Foo>(id);
#ifndef NDEBUG
    MADNESS_CHECK(ptr_opt && *ptr_opt == nullptr);
#else
    MADNESS_CHECK(!ptr_opt);
#endif

    print("test 6 (world object active message and tasks) seems to be working");
}


class TestFutureForwarding : public WorldObject<TestFutureForwarding> {
public:
    TestFutureForwarding(World& world)
            : WorldObject<TestFutureForwarding>(world) {
        this->process_pending();
    }

    virtual ~TestFutureForwarding() { }

    Future<int> test(int state) {
        if (state < 99) {
            return send(get_world().random_proc(), &TestFutureForwarding::test, state+1);
        }
        else {
            return Future<int>(state+1);
        }
    }
};

void test6a(World& world) {
    PROFILE_FUNC;

    if (world.size() < 2) return;

    TestFutureForwarding t(world);
    if (world.rank() == 0) {
        Future<int> fred = t.test(0);
        world.gop.fence();
        MADNESS_CHECK(fred.get() == 100);
    }
    else {
        world.gop.fence();
    }
    if (world.rank() == 0) {
        print("If got here test6a is OK!");
    }
}


void test7(World& world) {
    PROFILE_FUNC;
    int nproc = world.size();
    ProcessID me = world.rank();
    WorldContainer<int,double> c(world);

    typedef WorldContainer<int,double>::iterator iterator;
    typedef WorldContainer<int,double>::const_iterator const_iterator;
    typedef WorldContainer<int,double>::futureT futureT;

    // Everyone inserts distinct values 0..1000 into the container,
    // fences, and then tries to read all values back

    // Note that insertion with key or accessor should be safe
    for (int i=me; i<1000; i+=nproc) c.replace(i,(double) i);
    world.gop.fence();

    for (int i=999; i>=0; --i) {
        futureT fut = c.find(i);
        iterator it = fut.get();
        MADNESS_CHECK(it != c.end());
        double j = it->second;
        MADNESS_CHECK(j == i);
    }
    world.gop.fence();

    // Check that unset keys return end correctly
    for (int i=10001; i<10020; ++i) {
        MADNESS_CHECK(c.find(i).get() == c.end());
    }

    // Check that other iterators compare correctly
    MADNESS_CHECK(c.find(10).get() == c.find(10).get());
    MADNESS_CHECK(c.find(11).get() != c.find(12).get());
    MADNESS_CHECK(c.end() == c.end());
    MADNESS_CHECK(c.find(12).get() != c.end());

    // Loop thru local stuff
    for (iterator it=c.begin(); it != c.end(); ++it) {
        MADNESS_CHECK(it->first == it->second);
    };

    // Check shallow copy and const iterator
    const WorldContainer<int,double> d(c);

    // Loop thru local stuff with a const iterator
    for (const_iterator it=d.begin(); it != d.end(); ++it) {
        MADNESS_CHECK(it->first == it->second);
    };

    world.gop.fence();
    if (me == 0) print("test7 (world container basics) OK");
}

void test8(World& world) {
    PROFILE_FUNC;
    vector<unsigned char> v;
    archive::VectorOutputArchive arout(v);
    arout & &world;

    World* p = nullptr;
    archive::VectorInputArchive arin(v);
    arin & p;
    MADNESS_CHECK(p==&world);
    if (world.rank() == 0) print("test8 (serializing world pointer) OK");
}

void null_func() { }

int val_func() {
    return 1;
}

int val1d_func(int input) {
    return input+1;
}

void test9(World& world) {
    PROFILE_FUNC;
    const int ntask = 100000;

    double used = -cpu_time();
    for (int i=0; i<ntask; ++i) world.taskq.add(null_func);
    used += cpu_time();
    print("Time to add",ntask,"null, local tasks",used,"time/task",used/ntask);

    used = -cpu_time();
    world.taskq.fence();
    used += cpu_time();
    print("Time to run",ntask,"null, local tasks",used,"time/task",used/ntask);

    vector< Future<int> > v = future_vector_factory<int>(ntask);
    used = -cpu_time();
    for (int i=0; i<ntask; ++i) v[i] = world.taskq.add(val_func);
    used += cpu_time();
    print("Time to add",ntask,"value, local tasks",used,"time/task",used/ntask);

    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA0");
    std::cout.flush();
    world.taskq.fence();
    print("AAAAAAAAAAAAAAAA1");
    std::cout.flush();
    used += cpu_time();
    print("Time to run",ntask,"value, local tasks",used,"time/task",used/ntask);
    v.clear();
    print("AAAAAAAAAAAAAAAA");
    std::cout.flush();
    Future<int> input;
    Future<int> result = input;
    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA2");
    std::cout.flush();
    for (int i=0; i<ntask; ++i) {
        result = world.taskq.add(val1d_func,result);
    }
    used += cpu_time();
    print("AAAAAAAAAAAAAAAA3");
    std::cout.flush();
    print("Time to make",ntask,"chain of tasks",used,"time/task",used/ntask);
    input.set(0);
    used = -cpu_time();
    print("AAAAAAAAAAAAAAAA4");
    std::cout.flush();
    world.taskq.fence();
    print("AAAAAAAAAAAAAAAA5");
    std::cout.flush();
    used += cpu_time();
    print("Time to  run",ntask,"chain of tasks",used,"time/task",used/ntask);
    MADNESS_CHECK(result.get() == ntask);
    if (world.rank() == 0) print("test9 (time task creation and processing) OK");
}


class Mary {
private:
    mutable uint64_t val;
public:
    Mary() : val(0) {}

    void inc() const {
        val++;
    }
    void add(int i) {
        val += i;
    }
    void fred(int i, double j) {
        val += i*(int)j;
    };

    string cary0() {
        return string("Cary0 sends greetings");
    };

    string cary(int i) {
        ostringstream s;
        val += i;
        s << "Cary sends greetings: " << i << " " << val << endl;
        return s.str();
    };

    string alan(int i, int j, int proc) {
        ostringstream s;
        val += i*j;
        s << "Alan sends greetings: " << i << " " << j << " " << val << endl;
        return s.str();
    };

    double galahad(const string& str, int j, double z) {
        istringstream s(str);
        int i;
        s >> i;
        //val += i*j*z;
        print("Galahad",str,i,j,z,val);
        return val;
    };


    uint64_t get() const {
        return val;
    };

    bool get_me_twice(World* world, const WorldContainer<int,Mary>& d) {
        return true;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & val;
    }
};

void pounder(const WorldContainer<int,Mary>& m, int ind) {
    for (int i=0; i<1000/m.get_world().size()+10; ++i) 
        m.send(ind, &Mary::inc);
    print("pounder task finished sending");
}

void test10(World& world) {
    PROFILE_FUNC;
    // test forwarding methods to an item
    ProcessID me = world.rank();
    int nproc = world.size();
    WorldContainer<int,Mary> m(world);
    typedef WorldContainer<int,Mary>::iterator iterator;
    //world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::inc);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_CHECK(int(it->second.get()) == nproc);
    }
    world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::add,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_CHECK(long(it->second.get()) == nproc*(nproc+1)/2);
    }
    world.gop.fence();

    for (int i=0; i<nproc; ++i)
        m.send(i,&Mary::fred,2,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_CHECK(long(it->second.get()) == nproc*(3*nproc-1)/2);
    }
    world.gop.fence();
    print("finished forwarding test, starting pounder");

    // Test that item methods are executed atomically by having
    // everyone pound on one item

    const int ind = 9999999;
    if (world.rank() == 0) m.replace(std::pair<int,Mary>(ind,Mary()));
    world.gop.fence();
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    world.taskq.add(pounder, m, ind);
    print("pounding tasks submitted");
    world.gop.fence();
    print("finished pounding");
    if (world.rank() == 0)
      MADNESS_CHECK(long(m.find(ind).get()->second.get()) == nproc * (1000/world.size()+10) * 7);

    world.gop.fence();

    Future<double>  galahad = m.task(ProcessID(0),&Mary::galahad,string("1"),me,3.14);
    world.gop.fence();
    print("result of galahad",galahad.get());

    print("main making vector of results");
    //vector< Future<string> > results(nproc,Future<string>::default_initializer());
    vector< Future<string> > results = future_vector_factory<string>(nproc);
    vector< Future<bool> > b = future_vector_factory<bool>(nproc);
    print("main finished making vector of results");
    for (int i=0; i<nproc; ++i) {
        print("main making task",i);
        results[i] = m.task(i,&Mary::alan,3,4,world.rank());
        b[i] = m.send(i,&Mary::get_me_twice,&world,m);
        print("main finished making task",i);
    }
    print("about to fence");
    world.gop.fence();

    for (int i=0; i<nproc; ++i) {
        MADNESS_CHECK(results[i].probe());
        MADNESS_CHECK(b[i].probe());
        print("results",i,results[i].get(),b[i].get());
    };

    world.gop.fence();

    if (me == 0) print("test10 (messaging to world container items) OK");
}


struct Key {
    typedef unsigned long ulong;
    ulong n, i, j, k;
    hashT hashval;

    Key() {};  // Empty default constructor for speed - but is therefore junk

    Key(ulong n, ulong i, ulong j, ulong k)
            : n(n), i(i), j(j), k(k), hashval(0)
    {
        madness::hash_combine(hashval, n);
        madness::hash_combine(hashval, i);
        madness::hash_combine(hashval, j);
        madness::hash_combine(hashval, k);
    }

    hashT hash() const {
        return hashval;
    }

    template <typename opT>
    void foreach_child(const opT& op) const {
        ulong n2 = n+1;
        ulong i2 = i<<1;
        ulong j2 = j<<1;
        ulong k2 = k<<1;
        for (int p=0; p<2; ++p)
            for (int q=0; q<2; ++q)
                for (int r=0; r<2; ++r)
                    op(Key(n2,i2+p,j2+q,k2+r));
    }

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & n & i & j & k & hashval;
    }

    bool operator==(const Key& b) const {
        // It's extremely probable that different keys will have a different hash
        return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
    }
};


ostream& operator<<(ostream& s, const Key& key) {
    s << "Key(" << key.n << "," << key.i << "," << key.j << "," << key.k << "," << key.hash() << ")";
    return s;
}

struct Node {
    typedef WorldContainer<Key,Node> dcT;
    Key key;
    double value;
    bool isleaf;
    Node() : value(0.0), isleaf(true) {};
    Node(double value) : value(value), isleaf(true) {};
    Node(const Node& node) : key(node.key), value(node.value), isleaf(node.isleaf) {};

    struct do_random_insert {
        dcT& d;
        double value;
        do_random_insert(dcT& d, double value)
                : d(d), value(value) {}
        void operator()(const Key& key) const {
            d.task(key,&Node::random_insert,d, key, value);
        }
    };

    void random_insert(const dcT& constd, const Key& keyin, double valin) {
        dcT& d = const_cast<dcT&>(constd);
        MADNESS_CHECK(valin<=1.0);
        //print("inserting",keyin,valin);
        key = keyin;
        value = valin;
        isleaf = true;
        if (value>0.25 && key.n<9 && d.size()<4000) {
            isleaf = false;
            World& world = d.get_world();
            double ran = world.drand();
            MADNESS_CHECK(ran>=0 && ran<1.0);
            key.foreach_child(do_random_insert(d,value*ran));
        }
    }

    template <class Archive>
    void serialize(Archive& ar) {
        ar & key & value & isleaf;
    }

    bool is_leaf() const {
        return isleaf;
    }

    double get() const {
        return value;
    }

    void set(double v) {
        value = v;
    }
};

ostream& operator<<(ostream& s, const Node& node) {
    s << "Node(" << node.get() << "," << node.is_leaf() << ")" << endl;
    return s;
}


void walker1(WorldContainer<Key,Node>& d, const Key& key);

struct Walker1 {
    WorldContainer<Key,Node>& d;
    Walker1(WorldContainer<Key,Node>& d) : d(d) {};
    void operator()(const Key& key) const {
        walker1(d,key);
    };
};


void walker1(WorldContainer<Key,Node>& d, const Key& key) {
    static double counter = 0;
    WorldContainer<Key,Node>::iterator it = d.find(key).get();
    if (it != d.end()) {
        Node node = it->second;
        node.set(++counter);
        d.erase(key);
        d.replace(key,node);
        it = d.find(key).get();
        MADNESS_CHECK(it != d.end());
        MADNESS_CHECK(it->second.get() == counter);
        if (!node.is_leaf()) {
            key.foreach_child(Walker1(d));
        }
    }
}

void walker2(WorldContainer<Key,Node>& d, const Key& key) {
    static double counter = 1;
    WorldContainer<Key,Node>::iterator it = d.find(key).get();
    if (it != d.end()) {
        Node node = it->second;
        node.set(++counter);
        d.replace(key,node);
        it = d.find(key).get();
        MADNESS_CHECK(it != d.end());
        if (it->second.get() != counter) {
            print("failing",it->second.get(),counter,key,d.owner(key));
        }
        MADNESS_CHECK(it->second.get() == counter);
        if (!node.is_leaf()) {
            key.foreach_child(Walker1(d));
        }
    }
}

void test11(World& world) {
    PROFILE_FUNC;
    // Test the various flavours of erase
    ProcessID me = world.rank();
    WorldContainer<Key,Node> d(world);

    // First build an oct-tree with random depth
    world.srand(); // Each process will have a different random number
    print("first ran#",world.drand());
    world.gop.fence();
    if (me == 0) {
        Key root = Key(0,0,0,0);
        d.task(root,&Node::random_insert,d,root,1.0);
    }
    world.gop.fence();

    print("size before erasing",d.size());
    //d.clear();
    d.erase(d.begin(),d.end());
    print("size after erasing",d.size());
    world.srand(); // Each process will have a different random number
    print("first ran#",world.drand());
    world.gop.fence();
    // rebuild the tree in the same container
    if (me == 0) {
        Key root = Key(0,0,0,0);
        d.task(root,&Node::random_insert,d,root,1.0);
    }
    world.gop.fence();
    print("size after rebuilding",d.size());

    // Test get, erase, and re-insert of nodes with new value by node 0
    if (me == 0) {
        Key root = Key(0,0,0,0);
        walker1(d,root);
    }
    world.gop.fence();
    print("walker1 done");


    // Test get and re-insert of nodes with new value by node 0
    if (me == 0) {
        Key root = Key(0,0,0,0);
        walker2(d,root);
    }
    world.gop.fence();
    print("walker2 done");

    print("size before clearing",d.size());
    d.clear();
    print("size after clearing",d.size());
    if (me == 0) print("test11 (erasing and inserting in world containers) OK");
}


void test12(World& world) {
    PROFILE_FUNC;
    if (world.size() != 1) return;
    // Test file IO
    ProcessID me = world.rank();
    WorldContainer<int,double> d(world);

    // Everyone puts 100 distinct entries in the container
    for (int i=0; i<100; ++i) d.replace(me*100 + i, me*100+i);

    world.gop.fence();

    archive::BinaryFstreamOutputArchive out("testme.ar");
    out & d;
    out.close();

    world.gop.fence();

    archive::BinaryFstreamInputArchive in("testme.ar");
    WorldContainer<int,double> c(world);
    in & c;

    world.gop.fence();

    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        MADNESS_CHECK(c.probe(key));
        MADNESS_CHECK(c.find(key).get()->second == key);
    }

    world.gop.fence();

    if (world.rank() == 0) print("test12 (container archive I/O) OK");
}

void test13(World& world) {
    PROFILE_FUNC;
    // Basic functionality with 1 (default) writer
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> fout(world, "fred");
    fout & 1.0 & "hello";
    fout.close();

    double v=0.0;
    char s[6];
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> fin(world, "fred");
    fin & v & s;
    fin.close();
    fin.remove();

    print("This is what I read", v, s);
    world.gop.fence();


    // Store and load an archive with multiple writers

    int nio = (world.size()-1)/2 + 1;
    if (nio > 10) nio = 10;

    print("nio",nio);

    ProcessID me = world.rank();
    WorldContainer<int,double> d(world);
    // Everyone puts 100 distinct entries in the container
    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        d.replace(key, double(key));
    }

    world.gop.fence();

    fout.open(world,"fred",nio);
    fout & d;
    fout.close();

    WorldContainer<int,double> c(world);
    fin.open(world,"fred");
    fin & c;

    for (int i=0; i<100; ++i) {
        int key = me*100+i;
        MADNESS_CHECK(c.find(key).get()->second == key);
    }

    fin.close();
    archive::ParallelOutputArchive<>::remove(world, "fred");

    print("Test13 OK");
    world.gop.fence();
}

void test14(World& world) {

  if (world.size() > 1) {
    static size_t call_counter = 0;
    ++call_counter;

    const auto n = 1 + std::numeric_limits<int>::max()/sizeof(int);
    //const auto n = 1000000;

    auto iarray = std::make_unique<int[]>(n);
    iarray[0] = -1;
    iarray[n-1] = -1;

    world.gop.set_max_reducebcast_msg_size(std::numeric_limits<int>::max()/(std::min(10ul,call_counter)));
    world.gop.broadcast(iarray.get(), n, 0);

    if (world.rank() == 1) {
      MADNESS_CHECK(iarray[0] == -1);
      MADNESS_CHECK(iarray[n-1] == -1);
    }

    print("Test14 OK");
  }
  world.gop.fence();
}

void test15(World& world) {

  if (world.size() > 1) {
    //const auto n = 1 + std::numeric_limits<int>::max()/sizeof(int);
    const auto n = 1000000;
    auto iarray = std::make_unique<int[]>(n);

    if (world.rank() == 1)
      std::iota(iarray.get(), iarray.get()+n, 0);
    else
      std::fill(iarray.get(), iarray.get()+n, 0);

    world.gop.max(iarray.get(), n);

    if (world.rank() == 1) {
      MADNESS_CHECK(iarray[0] == 0);
      MADNESS_CHECK(iarray[n-1] == n-1);
    }

    print("Test15 OK");
  }
  world.gop.fence();
}

inline bool is_odd(int i) {
    return i & 0x1;
}

inline bool is_even(int i) {
    return !is_odd(i);
}

void work_odd(World& world) {
    test5(world);
    test6(world);
    test6a(world);
    test7(world);
    test8(world);
    test9(world);
    test10(world);
    //test11(world);
    // test12(world); cannot run due to filename collision
    // test13(world);
    test14(world);
    test15(world);
    world.gop.fence();
}

void work_even(World& world) {
    test5(world);
    test6(world);
    test6a(world);
    test7(world);
    test8(world);
    test9(world);
    test10(world);
    //test11(world);
    // test12(world); cannot run due to filename collision
    // test13(world);
    test14(world);
    test15(world);
    world.gop.fence();
}

void test_multi_world(World& world) {
    if (world.size() < 2) return;

    // Make two more worlds: odd processes, even processes

    // a) make list of ranks of processes in the subgroups
    //
    // Only process belonging to the subgroups participate
    // in the next steps (this does not mean that some processes
    // can be left out here! Intracomm::Create() only works
    // if all processes are engaged (blame MPI_Comm_create))
    //
    // b) make MPI group and hence new MPI sub-communcator
    //
    // c) make new worlds and do work

    std::cout << "\n\nREPEATING TESTS IN MULTI-WORLD\n\n" << std::endl;

    std::cout << "== multiple worlds created with Intracomm::Create()==" << std::endl;
    std::vector<int> odd, even;
    for (int i=0; i<world.size(); ++i) {
        if (is_odd(i))
            odd.push_back(i);
        else
            even.push_back(i);
    }

    const int color = world.rank() % 2;

    if (color) {   // Odd processes
      SafeMPI::Group g_odd = world.mpi.comm().Get_group().Incl(odd.size(), &odd[0]);
      SafeMPI::Intracomm comm_odd = world.mpi.comm().Create(g_odd);
      {
        World world_odd(comm_odd);
        work_odd(world_odd);
      }
    }
    else {                      // Even processes
      SafeMPI::Group g_even = world.mpi.comm().Get_group().Incl(even.size(),&even[0]);
      SafeMPI::Intracomm comm_even = world.mpi.comm().Create(g_even);
      {
        World world_even(comm_even);
        work_even(world_even);
      }
    }
    world.gop.fence();

    // try to do the same but use MPI_Comm_split
    {
      std::cout << "== multiple worlds created with Intracomm::Split()==" << std::endl;
      SafeMPI::Intracomm comm = world.mpi.comm().Split(color, world.rank());
      std::cout << "split_comm.size() = " << comm.Get_size() << std::endl;
      World subworld(comm);
      if (color == 1)
        work_odd(subworld);
      else
        work_even(subworld);
    }
    world.gop.fence();

    // now split world by host and run tasks on each host's world
    {
      std::cout << "== multiple worlds created with Intracomm::Split_type()==" << std::endl;
      SafeMPI::Intracomm comm = world.mpi.comm().Split_type(SafeMPI::Intracomm::SHARED_SPLIT_TYPE, world.rank());
      std::cout << "split_comm.size() = " << comm.Get_size() << std::endl;
      World subworld(comm);
      work_even(subworld);
    }
    world.gop.fence();

}

#ifdef HAVE_PARSEC
# ifdef PARSEC_HAVE_CUDA

extern void __cuda_hello_world(); // in hello_world.cu
class GPUHelloWorldTask : public TaskInterface {
public:
    void run(World& world) {
      __cuda_hello_world();
    }
};

void test_cuda0(World& world) {
  world.taskq.add(new GPUHelloWorldTask());
}
# endif
#endif

#if  MADNESS_CATCH_SIGNALS
void mad_signal_handler( int signum ) {
  // announce the signal
  std::cerr << "MADNESS caught signal " << signum << " will wait for you to attach a debugger" << std::endl;

  bool DebugWait = true;
  while (DebugWait) {
  }
}
#endif

int main(int argc, char** argv) {

#if  MADNESS_CATCH_SIGNALS
    signal(SIGSEGV, mad_signal_handler);
#endif
    World& world = initialize(argc,argv);

    // redirectio(world);
    print("The processor frequency is",cpu_frequency());
    print("There are",world.size(),"processes and I am process",world.rank(),"with",ThreadPool::size(),"threads");

    world.args(argc,argv);

    world.gop.fence();
    test_multi_world(world);

    try {
        PROFILE_BLOCK(main_program);

        test0(world);
        test5(world);
        test6(world);
        test6a(world);
        test7(world);
        test8(world);
        test9(world);
        test10(world);
        //test11(world);
        test12(world);
        test13(world);
        test14(world);
        test15(world);

        for (int i=0; i<10; ++i) {
          print("REPETITION",i);
          test_multi_world(world);
        }

#ifdef PARSEC_HAVE_CUDA
        test_cuda0(world);
#endif
    }
    catch (const SafeMPI::Exception& e) {
        print(e);
        error("caught an MPI exception");
    }
    catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    }
    catch (const char* s) {
        print(s);
        error("caught a string exception");
    }
    catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");

    //print_stats(world);
    finalize();
    return 0;
}
