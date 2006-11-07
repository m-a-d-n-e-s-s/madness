#include <vector>
#include <iostream>
using std::cout;
using std::endl;
#include <cstring>
using std::strcat;
#include <mra/mra.h>
#include <misc/communicator.h>
#include <serialize/textfsar.h>
using madness::archive::TextFstreamInputArchive;
using madness::archive::TextFstreamOutputArchive;
#include <serialize/binfsar.h>
using madness::archive::BinaryFstreamInputArchive;
using madness::archive::BinaryFstreamOutputArchive;
#include <serialize/vecar.h>
using madness::archive::VectorInputArchive;
using madness::archive::VectorOutputArchive;

namespace madness {

  template<typename T>
  class DiskDir : public Function<T>
  {
  private:
    char directoryname[256];
    char filename[256];
    Communicator comm;
    long partLevel;
  public:

    DiskDir(){};

    DiskDir(char directory[256], Communicator& comm) {
      strcpy(directoryname,directory);
      strcpy(filename, directory);
      partLevel = 0;
    };

    DiskDir(char directory[256], long Level, Communicator& comm) {
      strcpy(directoryname,directory);
      strcpy(filename, directory);
      partLevel = Level;
    };

    void operator<<(Function<T>& f) {
      f.save(filename, partLevel, comm);
    };

    void operator>>(Function<T>& f) {
      f.load(filename, comm);
    };

    void operator<<(vector< Function<T> >& f) {
      TextFstreamOutputArchive oar(comm.rank() ? 0 : filename);
      int nsize = f.size();
      oar & nsize;
      for(int i = 0; i < nsize; i++){
	if (comm.rank() == 0) {
	  f[i].saveMain(filename, oar, partLevel, comm);
	}
	else {
	  f[i].saveLoadWorker4DD(comm, true);
	}
	comm.Barrier();
      }
      oar.close();
    };

    void operator>>(vector< Function<T> >& f) {
      TextFstreamInputArchive iar(comm.rank() ? 0 : filename);
      long partLevel;
      int nsize;
      if (comm.rank() == 0) {
	iar & nsize;
      }
      comm.Barrier();
      comm.Bcast(nsize, 0);
      f.resize(nsize);
      for(int i = 0; i < nsize; i++) {
	if (comm.rank() == 0) {
	  iar & partLevel;
	  iar & FunctionDefaults::k;
	  iar & FunctionDefaults::thresh;
	  iar & FunctionDefaults::initial_level;
	  iar & FunctionDefaults::max_refine_level;
	  iar & FunctionDefaults::truncate_method;
	  iar & FunctionDefaults::autorefine;
	}
	comm.Bcast(FunctionDefaults::k, 0);
	comm.Bcast(FunctionDefaults::thresh, 0);
	comm.Bcast(FunctionDefaults::initial_level, 0);
	comm.Bcast(FunctionDefaults::max_refine_level, 0);
	comm.Bcast(FunctionDefaults::truncate_method, 0);
	comm.Bcast(FunctionDefaults::autorefine, 0);
	if (comm.rank() == 0) {
	  bool active_flag;
	  iar & active_flag;
	  if(active_flag) f[i].loadManager4DD(filename, iar, comm, active_flag, partLevel);
	  for ( int j = 1; j < comm.size(); j++) {
	    archive::MPIOutputArchive arout(comm, j);
	    arout & -1;
	  }
	}
	else {
	  f[i].saveLoadWorker4DD(comm, false);
	}
      }
      comm.Barrier();
      if(comm.rank() == 0) iar.close();
    };

    template<typename T2>
    void operator<<(T2& arbitrary) {
      TextFstreamOutputArchive oar(comm.rank() ? 0 : filename);
      if(comm.rank() == 0) oar & arbitrary;
    };

    template<typename T2>
    void operator>>(T2& arbitrary) {
      TextFstreamInputArchive iar(comm.rank() ? 0 : filename);
      if(comm.rank() == 0) iar & arbitrary;
    };


    void operator<<(char* arbitrary) {
      TextFstreamOutputArchive oar(comm.rank() ? 0 : filename);
      if(comm.rank() == 0) oar & arbitrary;
    };

    void operator>>(char* arbitrary) {
      TextFstreamInputArchive iar(comm.rank() ? 0 : filename);
      if(comm.rank() == 0) iar & arbitrary;
    };

    /*
    void operator<<(char[] arbitrary) {
      TextFstreamOutputArchive oar(comm.rank() ? 0 : filename);
      oar & arbitrary;
    };

    void operator>>(char[] arbitrary) {
      TextFstreamInputArchive iar(comm.rank() ? 0 : filename);
      iar & arbitrary;
    };
    */

    void operator << (vector<T>& arbitrary) {
      if(comm.rank() == 0) {
	TextFstreamOutputArchive oar(filename);
        oar & arbitrary.size();
        int nsize = arbitrary.size();
        for(int i = 0; i < nsize;i++) {
          oar & arbitrary[i];
        }
      }
      comm.Barrier();
    };

    void operator >> (vector<T>& arbitrary) {
      if(comm.rank() == 0) {
	TextFstreamInputArchive iar(filename);
	unsigned int n;
        iar & n;
	int nsize = n;
	arbitrary.resize(n);
        for(int i = 0; i < nsize;i++) {
          iar & arbitrary[i];
        }
      }
      comm.Barrier();
    };

    void operator[](char file[256]) {
      char newfilename[256];
      strcpy(newfilename, directoryname);
      strcat(newfilename, file);
      strcpy(filename, newfilename);
      cout << newfilename << endl;
    }

    void operator()(char file[256], long index) {
      char newfilename[256], tmp[16];
      strcpy(newfilename, directoryname);
      strcat(newfilename, file);
      sprintf(tmp, "%lu", index);
      strcat(newfilename, tmp);
      strcpy(filename, newfilename);
      cout << newfilename << endl;
    }

    void operator()(char file[256]) {
      char newfilename[256];
      strcpy(newfilename, directoryname);
      strcat(newfilename, file);
      strcpy(filename, newfilename);
      cout << newfilename << endl;
    }

    ~DiskDir(){};
  };

}
