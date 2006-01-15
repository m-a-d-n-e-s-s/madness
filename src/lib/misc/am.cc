#include <iostream>
using std::cout;
using std::endl;

#include <misc/communicator.h>
using madness::Communicator;
using madness::AMArg;

static unsigned long next = 1;

int myrand(void) {
  next = next * 1103515245 + 12345;
  return((unsigned)(next/65536) % 32768);
}

void mysrand(unsigned seed) {
  next = seed;
}

ProcessID random_proc(const Communicator& comm) {
  while (1) {
    ProcessID p = myrand()%comm.size();
    if (p != comm.rank()) return p;
  }
}

void handler(Communicator& comm, ProcessID from, const AMArg& arg) {
  comm.send(AMArg(arg.function,arg.arg0+1),from);
}

int main(int argc, char** argv) {
  MADMPIInit(argc, argv);
  Communicator comm;
  ProcessID me = comm.rank();
  long nproc = comm.nproc();

  if (nproc == 1) throw "Gimme someone to talk to!";

  mysrand(me);

  int handle = comm.am_register(handler);

  for (int i=0; i<2; i++) {
      long reply = 0;
      ProcessID p = random_proc(comm);
      cout << me << " sending " << p << endl;
      comm.am_send(AMArg(handle,me+1000),p);
      cout.flush();
      comm.recv(reply,p);
      cout << me << " received " << endl;
      cout.flush();
      if (reply != me+1001) throw "Ooops ...";
      comm.am_poll();
  }

  cout << me << "entering barrier" << endl;
  cout.flush();
  comm.am_barrier();

  MADMPIFinalize();
  return 0;
}

  
