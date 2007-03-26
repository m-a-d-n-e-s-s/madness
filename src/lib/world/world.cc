#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>



//////////////////////////////
/*    Test code below here  */
//////////////////////////////

using namespace madness;
using namespace std;

class B {
    long b;
public:
    B(long b=0) : b(b) {};
    void set(long value) {b=value;};
    long get() const {return b;};
    ~B(){print("B destructor");};
    template <typename Archive> void serialize(Archive &ar) {ar&b;}
};

void handler(World& world, ProcessID from, const AmArg& arg) {
    world.mpi.Send(arg.buf[0]+1, from, 33);
}

void hello(World& world, ProcessID from, const AmArg& arg) {
    print(world.mpi.rank(),"got hello from",from);
    cout.flush();
}

void test1(World& world) {
    ProcessID me = world.mpi.rank();
    long nproc = world.mpi.nproc();

    //world.mpi.set_debug(true);

    if (nproc == 1) throw "Gimme someone to talk to!";

    print(me,"ENTERING FENCE 0");
    world.gop.fence();
    print(me,"EXITING FENCE 0");

    for (int i=0; i<20; i++) {
      long reply = -1;
      ProcessID p = world.mpi.random_proc_not_me();
      world.am.send_recv(p,handler,AmArg(me+1000),&reply,sizeof(reply),p,33);
      if (reply != me+1001) {
	print("bad reply",reply,me+1001);
        throw "Ooops ...";
      }
    }
    print(me,"ENTERING FENCE 1");
    world.gop.fence();
    print(me,"EXITING FENCE 1");

    if (nproc < 16) {
      world.am.broadcast(hello,AmArg());  // Everyone says hello to everyone else
      world.gop.fence();
    }

    if (me == 0) print("AM test1 seems to be working");
}

void test2_handler(World& world, ProcessID src, void *buf, size_t len) {
    short *s = (short *) buf;
    short lens = (len-16)/sizeof(short);
    if ((lens*sizeof(short)) != (len-16))
        error("test2: msg not a multiple of sizeof(short)");

    for (short i=0; i<lens; i++)
        if (s[i+8] != i) throw "test2: bad contents";

    if (lens < 10000) {
        lens++;
	s = new short[lens+8];
        for (short i=0; i<lens; i++) s[i+8] = i;
        ProcessID dest = world.mpi.random_proc();
        world.am.send_long(dest, test2_handler, s, lens*sizeof(short)+16);
        world.am.fence();
        delete [] s;
    }
}


void test2(World& world) {
    ProcessID me = world.mpi.rank();
    short buf[] = {0,0,0,0,0,0,0,0,0}; // First 8 shorts for system header, 9th for data

    // Each process starts sending to a random process an array of
    // shorts of length 1 initialized to 0.  The handler checks the
    // contents, adds another entry on the end, and forwards to
    // another random process unless it already has more than 1000
    // elements.  On average, each process will send/recv 1000
    // messages.

    //world.set_debug(true);

    test2_handler(world, me, buf, sizeof(buf));
    world.gop.fence();
    if (me == 0) print("AM test2 seems to be working");
}

void test3_handler(World& world, ProcessID src, void *buf, size_t len) {
    short *s = (short *) buf;
    short lens = (len-16)/sizeof(short);
    if ((lens*sizeof(short)) != (len-16))
        error("test2: msg not a multiple of sizeof(short)");

    for (short i=0; i<lens; i++)
        if (s[i+8] != i) throw "test2: bad contents";

    if (lens < 10000) {
        lens++;
	s = new short[lens+8];
        for (short i=0; i<lens; i++) s[i+8] = i;
        ProcessID dest = world.mpi.random_proc();
        world.am.send_long_managed(dest, test3_handler, s, lens*sizeof(short)+16);
    }
}


void test3(World& world) {
    ProcessID me = world.mpi.rank();
    short buf[] = {0,0,0,0,0,0,0,0,0}; // First 8 shorts for system header, 9th for data

    // Same as test2 but using managed long buffers

    test3_handler(world, me, buf, sizeof(buf));
    world.gop.fence();
    if (me == 0) print("AM test3 seems to be working");
}

void testdumbser() {
    const size_t n=555;
    char buf[n+2];
    buf[0] = buf[n+1] = 127;

    BufferOutputArchive arout(buf+1,n-2);

    arout & 1L & 99.9;
    arout.store("hello",6);
    arout & 7L & 77.7;
    MADNESS_ASSERT(buf[0]==127 && buf[n+1]==127);

    long i;
    double a;

    BufferInputArchive arin(buf+1,n-2);
    char s[8];
    s[0] = s[7] = 66;
    arin & i & a; 
    MADNESS_ASSERT(i==1 && a==99.9);
    arin.load(s+1,6);
    MADNESS_ASSERT(s[0]==66 && s[1]=='h' && s[5]=='o' && s[7]==66);
    arin & i & a; 
    MADNESS_ASSERT(i==7 && a==77.7);
}

void test4(World& world) {
    int nproc = world.mpi.nproc();
    ProcessID me = world.mpi.rank();
    ProcessID left = (me+nproc-1) % nproc;
    ProcessID right = (me+1) % nproc;

    // First check local futures
    Future<int> a;
    a.set(1);
    MADNESS_ASSERT(a.get() == 1);

    Future<int> b;
    RemoteReference< FutureImpl<int> > rb=b.remote_ref(world), rc;
    world.mpi.Sendrecv(&rb,sizeof(rb),MPI::BYTE,left,1,
                       &rc,sizeof(rc),MPI::BYTE,right,1);
    Future<int> c(rc);
    c.set(me);
    world.gop.fence();
    MADNESS_ASSERT(b.get() == left);
    print("about to enter final barrier");
    world.gop.barrier();
    print("leaving final barrier");
    
    if (me == 0) print("test4 OK");
}

class Foo : public WorldObject<Foo> {
    int a;
public:
    Foo(World& world, int a) 
        : WorldObject<Foo>(world)
        , a(a) 
    {
        process_pending();
    };

    int get0() {return a;};
    int get1(int a1) {return a+a1;};
    int get2(int a1, char a2) {return a+a1+a2;};
    int get3(int a1, char a2, short a3) {return a+a1+a2+a3;};
    int get4(int a1, char a2, short a3, long a4) {return a+a1+a2+a3+a4;};
    int get5(int a1, char a2, short a3, long a4, short a5) {return a+a1+a2+a3+a4+a5;};

    int get0c() const {return a;};
    int get1c(int a1) const {return a+a1;};
    int get2c(int a1, char a2) const {return a+a1+a2;};
    int get3c(int a1, char a2, short a3) const {return a+a1+a2+a3;};
    int get4c(int a1, char a2, short a3, long a4) const {return a+a1+a2+a3+a4;};
    int get5c(int a1, char a2, short a3, long a4, short a5) const {return a+a1+a2+a3+a4+a5;};


};

void test5(World& world) {
    ProcessID me = world.rank();
    ProcessID nproc = world.nproc();
    Foo a(world, me*100);

    if (me == 0) {
        print(a.id());
        for (ProcessID p=0; p<nproc; p++) {
            MADNESS_ASSERT(a.send(p,&Foo::get0).get() == p*100);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get0).get() == p*100);

            MADNESS_ASSERT(a.send(p,&Foo::get1,1).get() == p*100+1);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get1,1).get() == p*100+1);

            MADNESS_ASSERT(a.send(p,&Foo::get2,1,2).get() == p*100+3);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get2,1,2).get() == p*100+3);

            MADNESS_ASSERT(a.send(p,&Foo::get3,1,2,3).get() == p*100+6);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get3,1,2,3).get() == p*100+6);            

            MADNESS_ASSERT(a.send(p,&Foo::get4,1,2,3,4).get() == p*100+10);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get4,1,2,3,4).get() == p*100+10);            

            MADNESS_ASSERT(a.send(p,&Foo::get5,1,2,3,4,5).get() == p*100+15);
            if (p) MADNESS_ASSERT(a.send(p,&Foo::get5,1,2,3,4,5).get() == p*100+15);            
        }
    }
    world.gop.fence();
    print("test 5 seems to be working");
}



void test6(World& world) {
    int nproc = world.mpi.nproc();
    ProcessID me = world.mpi.rank();
    World::poll_all();
    WorldContainer<int,double> c(world);
    //world.gop.fence();  // Currently needed until we can handle pending

    typedef WorldContainer<int,double>::iterator iterator;
    typedef WorldContainer<int,double>::const_iterator const_iterator;
    typedef WorldContainer<int,double>::future future;

    // Everyone inserts distinct values 0..1000 into the container,
    // fences, and then tries to read all values back
    for (int i=me; i<100; i+=nproc) c.insert(i,(double) i);
    print("finished inserting");
    world.gop.fence();

    for (int i=0; i<100; i++) {
        future fut = c.find(i);
        iterator it = fut.get();
	if (it == c.end()) {
	  print("failing",i);
	}
	MADNESS_ASSERT(it != c.end());
        double j = it->second;
        MADNESS_ASSERT(j == i);
    }
    world.gop.fence();
    print("finished finding");
    
    // Check that unset keys return end correctly
    for (int i=10001; i<10020; i++) {
        MADNESS_ASSERT(c.find(i).get() == c.end());
    }

    // Check that other iterators compare correctly
    MADNESS_ASSERT(c.find(10).get() == c.find(10).get());
    MADNESS_ASSERT(c.find(11).get() != c.find(12).get());
    MADNESS_ASSERT(c.end() == c.end());
    MADNESS_ASSERT(c.find(12).get() != c.end());

    // Loop thru local stuff
    for (iterator it=c.begin(); it != c.end(); ++it) {
        MADNESS_ASSERT(it->first == it->second);
    };


    // Check shallow copy and const iterator
    const WorldContainer<int,double> d(c);

    // Loop thru local stuff with a const iterator
    for (const_iterator it=d.begin(); it != d.end(); ++it) {
        MADNESS_ASSERT(it->first == it->second);
    };

    world.gop.fence();
    if (me == 0) print("test6 OK");
}

void test7(World& world) {
    vector<unsigned char> v;
    VectorOutputArchive arout(v);
    arout & &world;

    World* p;
    VectorInputArchive arin(v);
    arin & p;
    MADNESS_ASSERT(p==&world);
    if (world.rank() == 0) print("test7 OK");
}



#include <complex>
typedef std::complex<double> double_complex;

class TestTask : public TaskInterface {
 public:
    void run(World& world) {print("Hi, I am running!");}
};

class TTT {
private:
    int state;
public:
    TTT() : state(0) {};

    static void fred(){print("Oops-a-daisy!");};

    static int mary() {return 99;};

    static int carl() {return 88;};

    static int dave(World* world) {return world->mpi.rank();};

    static int bert(int input) {return input+1;};

    static double_complex sara(double a, const double_complex& b) {
        return a*b;
    };

    static string kate(World* world, const string& msg, double d) {
        ostringstream s;
        s << "Process " << world->mpi.rank() << " says '" << msg << "' and " 
          << d << " right back at you!";
        return s.str();
    };

    double jody(double a, double b, double c) {
        return a+b+c + state;
    };

    double hugh(vector< Future<int> >& a) {
        double sum = 0.0;
        for (int i=0; i<(int)a.size(); i++) sum += a[i].get();
        return sum;
    };
};

double dumb(int a1, int a2, int a3, int a4, int a5, int a6, int a7) {
    return a1+a2+a3+a4+a5+a6+a7;
}



void test8(World& world) {
    int nproc = world.mpi.nproc();
    ProcessID me = world.mpi.rank();
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
    Future<int> mary = world.taskq.add(TTT::mary);
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
    Future<string> cute = world.taskq.add(right,TTT::kate,&world,"Boo!",-42.0);
    TTT ttt;
    Future<double> jody = world.taskq.add(ttt,&TTT::jody,1.0,2.0,3.0);
    Future<double> duh = world.taskq.add(me,dumb,0,1,2,3,4,5,6);
    print("done with making futs");

    bert_input = 7;
    sara1 = 3.0;
    sara2 = double_complex(2.1,1.2);
    kate2 = string("Who's your daddy?");
    kate3 = 3.14;

    vector< Future<int> > futv = future_vector_factory<int>(7);
    Future<double> hugh = world.taskq.add(ttt,&TTT::hugh,futv);
    for (int i=0; i<7; i++) {
        print("assigning",i,futv[i]);
        futv[i] = i;
    }
    
    print("about to fence again");
    world.gop.fence();
    print("finished fence again");

    MADNESS_ASSERT(mary.probe());
    MADNESS_ASSERT(carl.probe());
    MADNESS_ASSERT(dave.probe());
    MADNESS_ASSERT(bert.probe());
    MADNESS_ASSERT(sara.probe());
    MADNESS_ASSERT(kate.probe());
    MADNESS_ASSERT(cute.probe());
    MADNESS_ASSERT(jody.probe());
    MADNESS_ASSERT(hugh.probe());
    MADNESS_ASSERT(duh.probe());

    MADNESS_ASSERT(mary.get() == 99);
    MADNESS_ASSERT(carl.get() == 88);
    MADNESS_ASSERT(dave.get() == right);
    MADNESS_ASSERT(bert.get() == 8);
    MADNESS_ASSERT(hugh.get() == 21.0);
    MADNESS_ASSERT(duh.get() == 21.0);
    print("Sara says",sara.get().real(),sara.get().imag());
    print("Kate says",kate.get());
    print("Cute says",cute.get());
    print("Jody says",jody.get());
    
    if (me == 0) print("test9 OK");
}

void null_func(){};

int val_func() {return 1;};

int val1d_func(int input) {
    return input+1;
}

void test9(World& world) {
    const int ntask = 100000;

    double used = -cpu_time();
    for (int i=0; i<ntask; i++) world.taskq.add(null_func);
    used += cpu_time();
    print("Time to add",ntask,"null, local tasks",used,"time/task",used/ntask);
    
    used = -cpu_time();
    world.taskq.fence();
    used += cpu_time();
    print("Time to run",ntask,"null, local tasks",used,"time/task",used/ntask);

    vector< Future<int> > v = future_vector_factory<int>(ntask);
    used = -cpu_time();
    for (int i=0; i<ntask; i++) v[i] = world.taskq.add(val_func);
    used += cpu_time();
    print("Time to add",ntask,"value, local tasks",used,"time/task",used/ntask);
    
    used = -cpu_time();
    world.taskq.fence();
    used += cpu_time();
    print("Time to run",ntask,"value, local tasks",used,"time/task",used/ntask);
    v.clear();

    Future<int> input;
    Future<int> result = input;
    used = -cpu_time();
    for (int i=0; i<ntask; i++) {
        result = world.taskq.add(val1d_func,result);
    }
    used += cpu_time();
    print("Time to make",ntask,"chain of tasks",used,"time/task",used/ntask);
    input = 0;
    used = -cpu_time();
    world.taskq.fence();
    used += cpu_time();
    print("Time to  run",ntask,"chain of tasks",used,"time/task",used/ntask);
    MADNESS_ASSERT(result.get() == ntask);
}


class Mary {
private:
    int val;
public:
    Mary() : val(0) {
        print("MAKING Mary",(void*)this);
    };
    void inc() {
        print("INC Mary",(void*)this,val);
        val++;
        print("Mary was just incremented",val);
    };
    void add(int i) {
        print("ADD Mary",(void*)this,val,i);
        val += i;
        print("Mary was just added",i,val);
    };
    void fred(int i, double j) {
        print("FRED Mary",(void*)this,val,i,j);
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

    string alan(int i, int j) {
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
        

    int get() const {return val;};

    bool get_me_twice(World* world, const WorldContainer<int,Mary>& d) {
        return true;
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
        ar & val;
    }
};

void test10(World& world) {
    // test forwarding methods to an item
    ProcessID me = world.mpi.rank();
    int nproc = world.mpi.nproc();
    World::poll_all();
    WorldContainer<int,Mary> m(world);
    typedef WorldContainer<int,Mary>::iterator iterator;
 
    for (int i=0; i<nproc; i++) 
        m.send(i,&Mary::inc);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(it->second.get() == nproc);
    }
    world.gop.barrier();

    for (int i=0; i<nproc; i++) 
        m.send(i,&Mary::add,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(it->second.get() == nproc*(nproc+1)/2);
    }
    world.gop.fence();

    for (int i=0; i<nproc; i++) 
        m.send(i,&Mary::fred,2,me);
    world.gop.fence();

    for (iterator it=m.begin(); it!=m.end(); ++it) {
        print("mary",it->first,it->second.get());
        MADNESS_ASSERT(it->second.get() == nproc*(3*nproc-1)/2);
    }

    Future<double>  galahad = m.task(ProcessID(0),&Mary::galahad,string("1"),me,3.14);
    world.gop.fence();
    print("result of galahad",galahad.get());


    print("main making vector of results");
    //vector< Future<string> > results(nproc,Future<string>::default_initializer());
    vector< Future<string> > results = future_vector_factory<string>(nproc);
    vector< Future<bool> > b = future_vector_factory<bool>(nproc);
    print("main finished making vector of results");
    for (int i=0; i<nproc; i++) {
        print("main making task",i);
        results[i] = m.task(i,&Mary::alan,3,4);
        b[i] = m.send(i,&Mary::get_me_twice,&world,m);
        print("main finished making task",i);
    }
    print("about to fence");
    world.gop.fence();

    for (int i=0; i<nproc; i++) {
        MADNESS_ASSERT(results[i].probe());
        MADNESS_ASSERT(b[i].probe());
        print("results",i,results[i].get(),b[i].get());
    };

    if (me == 0) print("test8 OK");
}


struct Key {
    typedef unsigned long ulong;
    ulong n, i, j, k;
    hashT hashval;

    Key() {};  // Empty default constructor for speed - but is therefore junk

    Key(ulong n, ulong i, ulong j, ulong k)
        : n(n), i(i), j(j), k(k), hashval(madness::hash(&this->n,4,0)) {};

    hashT hash() const {
        return hashval;
    };

    template <typename opT> 
    void foreach_child(const opT& op) {
        ulong n2 = n+1;
        ulong i2 = i<<1;
        ulong j2 = j<<1;
        ulong k2 = k<<1;
        for (int p=0; p<2; p++)
            for (int q=0; q<2; q++)
                for (int r=0; r<2; r++)
                    op(Key(n2,i2+p,j2+q,k2+r));
    };

    template <typename Archive>
    void serialize(const Archive& ar) {
	ar & n & i & j & k & hashval;
    }

    bool operator==(const Key& b) const {
        // It's extremely probable that different keys will have a different hash
        return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
    };
};

ostream& operator<<(ostream& s, const Key& key) {
    s << "Key(" << key.n << "," << key.i << "," << key.j << "," << key.k << "," << key.hash() << ")";
    return s;
}

// struct Node {
//     typedef WorldContainer<Key,Node> dcT;
//     Key key;
//     double value;
//     bool isleaf;
//     Node() : value(0.0), isleaf(true) {};
//     Node(double value) : value(value), isleaf(true) {};

//     struct do_random_insert {
//         dcT& d;
//         double value;
//         do_random_insert(dcT& d, double value) 
//             : d(d), value(value) {};
//         void operator()(const Key& key) const {
//             d.task(key,&Node::random_insert,d, key, value);
//         };
//     };
 
//     void random_insert(const dcT& constd, const Key& keyin, double valin) {
//         dcT& d = const_cast<dcT&>(constd);
//         //print("inserting",keyin,valin);
//         key = keyin;
//         value = valin;
//         isleaf = true;
//         if (value > 0.25) {
//             isleaf = false;
//             World& world = d.world();
//             double ran = world.mpi.drand();
//             key.foreach_child(do_random_insert(d,value*ran)); 
//         }
//     };

//     template <class Archive>
//     void serialize(Archive& ar) {
//         ar & key & value & isleaf;
//     }
// };


// void test11(World& world) {
//     // Test the various flavours of erase
//     ProcessID me = world.rank();
//     WorldContainer<Key,Node> d(world);

//     // First build an oct-tree with random depth
//     if (me == 0) {
//         Key root = Key(0,0,0,0);
//         d.task(root,&Node::random_insert,d,root,1.0);
//     }
//     world.gop.fence();
//     print("size before erasing",d.size());
//     //d.clear();
//     d.erase(d.begin(),d.end());
//     print("size after erasing",d.size());
//     world.gop.fence();
// }



int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    
    World world(MPI::COMM_WORLD);
    redirectio(world);

    print("The processor frequency is",cpu_frequency());
    print("there are",world.mpi.nproc(),"processes and I am process",world.mpi.rank());

    for (int arg=1; arg<argc; arg++) {
        if (std::strcmp(argv[arg],"-dx")==0) xterm_debug("world", 0);
        if (std::strcmp(argv[arg],"-dam")==0) world.am.set_debug(true);
        if (std::strcmp(argv[arg],"-dmpi")==0) world.mpi.set_debug(true);
        if (std::strcmp(argv[arg],"-dref")==0) world.mpi.set_debug(true);
    }

    world.gop.fence();

    try {
//         test1(world);
//         test2(world);
//         test3(world);
        test4(world);
        test5(world);
        test6(world);
        test7(world);
        test8(world);
        test9(world);
        test10(world);
        //        test11(world);
    } catch (MPI::Exception e) {
        error("caught an MPI exception");
    } catch (madness::MadnessException e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const char* s) {
        print(s);
        error("caught a string exception");
    } catch (...) {
        error("caught unhandled exception");
    }

    print("entering final fence");
    world.gop.fence();
    print("done with final fence");
    MPI::Finalize();
    return 0;
}
