#include <iostream>
#include <tasks/tasks.h>
#include <misc/print.h>
#include <misc/misc.h>
#include <cstdlib>

#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

using namespace std;
using namespace madness;

void ping_pong_handler(Communicator& comm, ProcessID src, const AMArg& arg);    

class TaskPingPong : public TaskInterface {
    int a;
    ProcessID dest;
public:
    TaskPingPong(int a, ProcessID dest) : a(a+1), dest(dest) {};
    bool probe() const {return true;};
    void run() {
        print(comm_default->rank(),"executing",a);
        if (a<=10) madness::comm_default->am_send(dest,ping_pong_handler,AMArg(a));
    };
    ~TaskPingPong() {};
};   

void ping_pong_handler(Communicator& comm, ProcessID src, const AMArg& arg) {
    taskq.add_local(new TaskPingPong(arg.arg0,src));
} 

void pass_the_parcel(Communicator& comm, ProcessID src, VectorInputArchive& ar) {
    long i;
    ar & i;
    if (i<100) {
        ProcessID dest = rand()%comm.nproc();
        print(src,"just gave",comm.rank(),"this",i,"now sending to",dest);
        std::vector<unsigned char> v(128);
        VectorOutputArchive ar(v);
        i++;
        ar & i;
        taskq.add_generic(dest, pass_the_parcel,v);
    }
    else {
        print(src,"just gave",comm.rank(),"this",i,"finished");  
    }
}


void feedme(Communicator& comm, ProcessID src, VectorInputArchive& ar) {
    long i;
    double d;
    ar & i & d;
    print(src,"just gave",madness::comm_default->rank(),"this",i,d);
}

long have_fun_sum;
void have_fun(Communicator& comm, ProcessID src, const AMArg& arg) {
    // This will generate (p-1)^2 tasks on each process
    long i = arg.arg0;
    have_fun_sum++;
    if (have_fun_sum < 10) print(src,"asked me if I was having fun yet!");
    if (i < 2) taskq.add_am_broadcast(have_fun, AMArg(i+1));
}



int main(int argc, char* argv[]) {
    MPI::Init(argc, argv);
    Communicator comm;
    madness::comm_default = &comm;
    madness::redirectio(comm);
    print("ping_pong handler",comm_default->am_register(ping_pong_handler));
    print("task generic handler",comm_default->am_register(task_generic_handler));
    print("task generic broadcast handler",comm_default->am_register(task_generic_broadcast_handler));
    print("have_fun handler",comm_default->am_register(have_fun));
    taskq.register_generic_op(feedme);
    taskq.register_generic_op(pass_the_parcel);
    
    std::srand(comm.rank()*1023);
    std::vector<unsigned char> v;
    VectorOutputArchive ar(v);
    
    //madness::xterm_debug(*comm_default,0,0);
    //comm.set_debug(true);
    
    try {
        if (comm_default->rank() == 0) taskq.add_local(new TaskPingPong(0,1));
        taskq.global_fence();
        print(comm.rank(),"EXITING SERVER 1");
    
        ar.close();
        ar.open();
        ar & 0L;
        if (comm.rank() == 0) taskq.add_generic(0,pass_the_parcel,v);
        taskq.global_fence();
        print(comm.rank(),"EXITING SERVER 2");
        
        ar.close();
        ar.open();
        ar & 99L & 3.14159;
        taskq.add_generic(0,feedme,v);
        taskq.add_generic(1,feedme,v);
        taskq.global_fence();
        print(comm.rank(),"EXITING SERVER 3");
        ar.close();
        ar.open();
        ar & 42L & 2.71828;
    
        taskq.add_generic_broadcast(feedme,v);
        taskq.global_fence();
        print(comm.rank(),"EXITING SERVER 4");
        
        have_fun_sum = 0;
        taskq.global_fence();
        if (comm.rank()==0) taskq.add_am_broadcast(have_fun,AMArg(0));
        taskq.global_fence();
        long np1 = (comm.nproc()-1);
        if (have_fun_sum != np1*np1) print("Sorry ... fun was not had",have_fun_sum);
        print(comm.rank(),"EXITING SERVER 5");
    } catch (char const* msg) {
        std::cerr << "Exception (string): " << msg << std::endl;
        comm.Abort();
    } catch (std::exception& e) {
        std::cerr << "Exception (std): " << e.what() << std::endl;
        comm.Abort();
    } catch (MadnessException& e) {
        std::cerr << e << std::endl;
        comm.Abort();
    } catch (TensorException& e) {
        std::cerr << e << std::endl;
        comm.Abort();
    } catch (MPI::Exception& e) {
        std::cerr << "Exception (mpi): code=" << e.Get_error_code()
        << ", class=" << e.Get_error_class()
        << ", string=" << e.Get_error_string() << std::endl;
        comm.Abort();
    } catch (...) {
        std::cerr << "Exception (general)" << std::endl;
        comm.Abort();
    }
    
    
    
    comm.close();
    
    MPI::Finalize();
    return 0;
}


