#include <iostream>
#include <tasks/tasks.h>
#include <tasks/sav.h>
#include <misc/vector_factory.h>
#include <misc/print.h>

using namespace std;
using namespace madness;



void doit(int a) {
    cout << comm_default->rank() << " Doing " << a << endl;
}

void setit(SAV<int>& a) {
    cout << comm_default->rank() << " Setting 99" << endl;
    a.set(99);
}

void chain(int a, SAV<int>& b) {
    cout << comm_default->rank() << " Chaining " << a << endl;
    b.set(a-1);
}

int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    Communicator comm;
    madness::comm_default = &comm;
    TaskQueue q;

    typedef Task1in< void (*)(int), SAV<int> > doitT;
    typedef Task1out< void (*)(SAV<int>&), SAV<int> > setitT;
    typedef Task1in1out< void (*)(int,SAV<int>&), SAV<int>, SAV<int> > chainT;

    cout << "Hello from " << comm.rank() << endl;


    // A bunch of tasks on node 0 all depending on the same
    // assignment at node 1;

    if (comm.rank() == 0) {
        SAV<int> a(1,1,true);
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.add(new doitT(&doit, a));
        q.wait();
    } else if(comm.rank() == 1) {
        SAV<int> a(0,1,false);
        a.set(33);
    }

    // A chain of tasks interleaved between 0 and 1
    if (comm.rank() == 0) {
        SAV<int> a(1,1,true),b(1,2,false);
        SAV<int> c(1,3,true),d(1,4,false);
        SAV<int> e(1,5,true),f(1,6,false);
        q.add(new chainT(&chain, a, b));
        q.add(new chainT(&chain, c, d));
        q.add(new chainT(&chain, e, f));
    } else {
        SAV<int> a(0,1,false),b(0,2,true);
        SAV<int> c(0,3,false),d(0,4,true);
        SAV<int> e(0,5,false),f(0,6,true);
        q.add(new chainT(&chain, b, c));
        q.add(new chainT(&chain, d, e));
        q.add(new doitT(&doit, f));
        a.set(10);
    }
    q.wait();

    // Send a vector from 0 to 1 ... send a Tensor back
    if (comm.rank() == 0) {
        SAV< vector<int> > v(1,44,false);
        v.set(vector_factory(3,2,1));
        SAV< Tensor<double> > t(1,45,true);
        while (!t.probe());
        std::cout << t.get();
    } else {
        SAV< vector<int> > v(0,44,true);
        while (!v.probe());
        const vector<int>& w = v.get();
        print(w.size(),w[0],w[1],w[2]);
        SAV< Tensor<double> > t(0,45,false);
        t.set(Tensor<double>(2,3).fillrandom());
    }
    print(comm.rank(),"HERE");
    // Same but using fixed size vectors and tensors.
    vector<long> dims = vector_factory(2L,4L);
    if (comm.rank() == 0) {
        SAV< vector<int> > v(1,44,false,6);
        v.set(vector_factory(5,4,3,2,1,0));
        SAV< Tensor<double> > t(1,45,true,dims);
        while (!t.probe());
        std::cout << t.get();
    } else {
        SAV< vector<int> > v(0,44,true,6);
        while (!v.probe());
        const vector<int>& w = v.get();
        print(w.size(),w[0],w[1],w[2]);
        SAV< Tensor<double> > t(0,45,false,dims);
        t.set(Tensor<double>(dims).fillindex());
    }

    comm.close();
    MPI::Finalize();
    return 0;
}


