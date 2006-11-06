#include <iostream>
using std::cout;
using std::endl;

#include <cmath>
using std::abs;

#include <algorithm>
using std::max;
using std::sort;

#include <complex>
#include <cstring>
using std::strcpy;

#include "mra.h"
#include "twoscale.h"
#include "legendre.h"

#include <misc/madexcept.h>
#include <misc/shared_ptr.h>

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

  template <typename T>
  template <class Archive>
  void Function<T>::save_local(Archive& ar)
  {
    ar & FunctionDefaults::k;
    ar & FunctionDefaults::thresh;
    ar & FunctionDefaults::initial_level;
    ar & FunctionDefaults::max_refine_level;
    ar & FunctionDefaults::truncate_method;
    ar & FunctionDefaults::autorefine;
    //_save_local(ar, tree());
    _save_local(ar, tree());
  };
  
  template <typename T>
  template <class Archive>
  void Function<T>::_save_local(Archive& ar, const OctTreeT* tree)
  {
    //    void Function<T>::_save_local(Archive& ar, const OctTreeT* tree) {
    ar & isactive(tree);
    ar & tree->n() & tree->x() & tree->y() & tree->z();
    if(isactive(tree)) {
      const TensorT *t = coeff(tree);
      ar & (t != 0);
      if(t) ar & *t;
    }
    FORIJK(
	   ar & (tree->child(i,j,k)!=0);
	   if (tree->child(i,j,k)) {
	     _save_local(ar, tree->child(i,j,k));
	   }
	   );
  }

  template <typename T>
  template <class Archive>
  void Function<T>::load_local(const Archive& ar)
  {
    ar & FunctionDefaults::k;
    ar & FunctionDefaults::thresh;
    ar & FunctionDefaults::initial_level;
    ar & FunctionDefaults::max_refine_level;
    ar & FunctionDefaults::truncate_method;
    ar & FunctionDefaults::autorefine;
    bool active_flag;
    ar & active_flag;
    //_load_local(ar, tree());
    _load_local(ar, tree());
  }

  template <typename T>
  template <class Archive>
  void Function<T>::_load_local(const Archive& ar, OctTreeT* tree)
  {
    //    void Function<T>::_load_local(const Archive& ar, OctTreeT* tree) {
    set_active(tree);
    Level n_local;
    Translation x_local, y_local, z_local;
    ar & n_local & x_local & y_local & z_local;
    bool have_coeffs;
    ar & have_coeffs;
    if(have_coeffs) {
      TensorT t;
      ar & t;
      set_coeff(tree, t);
    }
    FORIJK(
	   bool have_child;
	   ar & have_child;
	   if (have_child) {
	     OctTreeTPtr child = tree->child(i, j, k);
	     if (!child) child = tree->insert_local_child(i,j,k);
	     bool active_flag;
	     ar & active_flag;
	     if (active_flag) _load_local(ar, child);
	   }
	   );
  }

  template <typename T>
  void Function<T>::save(const char* f, const long partLevel, Communicator& comm)
  {
    if (comm.rank() == 0) {
      TextFstreamOutputArchive oar(f);
      saveMain(f, oar, partLevel, comm);
      oar.close();
    }
    else {
      saveLoadWorker(tree(), comm, true);
    }
  }

  template <typename T>
  template <class Archive>
  void Function<T>::saveMain(const char* f, Archive& ar,
			     const long partLevel, Communicator& comm)
  {
    ar & partLevel;
    ar & FunctionDefaults::k;
    ar & FunctionDefaults::thresh;
    ar & FunctionDefaults::initial_level;
    ar & FunctionDefaults::max_refine_level;
    ar & FunctionDefaults::truncate_method;
    ar & FunctionDefaults::autorefine;
    //saveManager(f, ar, tree(), partLevel, comm);
    saveManager(f, ar, tree(), partLevel, comm);
    if( comm.size() != 1) {
      for ( int i = 1; i < comm.size(); i++) {
	archive::MPIOutputArchive arout(comm, i);
	arout & -1;
      }
    }
  };

  template <typename T>
  template <class Archive>
  void Function<T>::saveManager(const char* f, Archive& ar, const OctTreeT* tree, const long partLevel, Communicator& commFunc)
  {
//    void Function<T>::saveManager(const char* f, Archive& ar, const OctTreeT* tree, const long partLevel, Communicator& commFunc) {
    if (isremote(tree)) {
      shadowManager_save(f, ar, tree, tree->n(), tree->x(), tree->y(), 
			 tree->z(), tree->rank(), partLevel, commFunc);
    }
    else {
      ar & isactive(tree);
      ar & tree->n() & tree->x() & tree->y() & tree->z();
      if(isactive(tree)) {
	const TensorT *t = coeff(tree);
	ar & (t != 0);
	if(t) ar & *t;
      }
      FORIJK(
	     OctTreeTPtr child = tree->child(i,j,k);
	     if(child) {
	       if(partLevel > 0 && child->n() == partLevel && child->islocal()) {
		 char ftest[256];
		 produceNewFilename(f, partLevel, child, ftest);
		 TextFstreamOutputArchive oar(ftest);
		 if(child->islocal()) {
		   oar & (child!=0);
		 }
		 saveManager(f, oar, child, partLevel, commFunc);
	       }
	       else {
		 if(child->islocal()) {
		   ar & (child!=0);
		 }
		 saveManager(f, ar, child, partLevel, commFunc);
	       }
	     }
	     else {
	       if(partLevel > 0 && tree->n() == (partLevel-1)) {
		 char ftest[256];
		 produceNewFilename(f, partLevel, child, ftest);
		 TextFstreamOutputArchive oar(ftest);
		 oar & (child!=0);
	       }
	       else{
		 ar & (child!=0);
	       }
	     }
	     );
    }
  }

  template <typename T>
  template <class Archive>
  void Function<T>::shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree, Level n, Translation x, Translation y, Translation z, ProcessID remoteRank, const long partLevel, Communicator& commFunc) 
//    void Function<T>::shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree, Level n, Translation x, Translation y, Translation z, ProcessID remoteRank, const long partLevel, Communicator& commFunc) 
  {
    int nRemoteBranch;
    madness::archive::MPIOutputArchive arout(commFunc, remoteRank);
    // Sending start branch data.
    arout & 1 & n & x & y & z;
    //      MPI_Barrier(MPI_COMM_WORLD);
    arout.flush();
    madness::archive::MPIInputArchive arin(commFunc, remoteRank);
    std::vector<localTreeMember> subtreeList; 
    std::vector<int> fileList;
    // Getting local subtree list.
    arin & nRemoteBranch;
    //      MPI_Barrier(MPI_COMM_WORLD);
    subtreeList.resize(nRemoteBranch);
    for (int i = 0; i < nRemoteBranch; i++) {
      arin & subtreeList[i];
    }
    //      MPI_Barrier(MPI_COMM_WORLD);
    int iter = 0;
    if(subtreeList[iter].n == partLevel)
      {
        char ftest[256];
        produceNewFilename2(f, partLevel, subtreeList[iter], ftest);
        TextFstreamOutputArchive oar1(ftest);
        dataSaveInShaManSave(f, oar1, tree, subtreeList[iter], remoteRank, partLevel, commFunc);
        iter += 1;
        _shadowManager_save(f, oar1, tree, subtreeList, remoteRank, partLevel, commFunc, nRemoteBranch, iter);
      }
    else{
      dataSaveInShaManSave(f, ar, tree, subtreeList[iter], remoteRank, partLevel, commFunc);
      iter += 1;
      _shadowManager_save(f, ar, tree, subtreeList, remoteRank, partLevel, commFunc, nRemoteBranch, iter);
    }
    //      MPI_Barrier(MPI_COMM_WORLD);
  }

  template <typename T>
  template <class Archive>
  void Function<T>::_shadowManager_save(const char* f, Archive& ar, const OctTreeT* tree, std::vector<localTreeMember>& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc, const int nRemoteBranch, int& iter) 
//    void Function<T>::dataSaveInShaManSave(const char* f, Archive& ar, const OctTreeT* tree, localTreeMember& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc) 
  {
    FORIJK(
	   if(partLevel > 0 && subtreeList[iter].n == partLevel && !subtreeList[iter].remote) {
	     char ftest[256];
	     produceNewFilename2(f, partLevel, subtreeList[iter], ftest);
	     TextFstreamOutputArchive oar2(ftest);
	     if(iter < nRemoteBranch) {
	       dataSaveInShaManSave(f, oar2, tree, subtreeList[iter], remoteRank, partLevel, commFunc);
	       iter += 1;
	       if(subtreeList[iter-1].have_child)
		 {
		   _shadowManager_save(f, oar2, tree, subtreeList, remoteRank, partLevel, commFunc, nRemoteBranch, iter);
		 }
	     }
	   }
	   else{
	     if(iter < nRemoteBranch) {
	       dataSaveInShaManSave(f, ar, tree, subtreeList[iter], remoteRank, partLevel, commFunc);
	       iter += 1;
	       if(subtreeList[iter-1].have_child)
		 {
		   _shadowManager_save(f, ar, tree, subtreeList, remoteRank, partLevel, commFunc, nRemoteBranch, iter);
		 }
	     }
	   }
	   );
  }

  template <typename T>
  template <class Archive>
  void Function<T>::dataSaveInShaManSave(const char* f, Archive& ar, const OctTreeT* tree, localTreeMember& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc) 
  {
    madness::archive::MPIInputArchive arin(commFunc, remoteRank);
    if(!subtreeList.remote) {
      ar & subtreeList.have_child;
    }
    if(subtreeList.have_child) {
      if(subtreeList.remote) {
	if(subtreeList.rank==0) {
          // Going to other local tree on master.
	  OctTreeT* treemaster = tree->find(subtreeList.n, subtreeList.x, subtreeList.y, subtreeList.z);
	  saveManager(f, ar, treemaster, partLevel, commFunc);
	}
	else {
          // Going to other local tree on client.
	  shadowManager_save(f, ar, tree, subtreeList.n, subtreeList.x, subtreeList.y, subtreeList.z, subtreeList.rank, partLevel, commFunc);
	}
      }
      else {
	// Storing active flag.
	ar & subtreeList.active;
	ar & subtreeList.n & subtreeList.x & subtreeList.y & subtreeList.z;
	if(subtreeList.active) {
	  // Inquiring Coefficients.
	  bool coeffsPointer;
	  commFunc.Recv(&coeffsPointer, 1, subtreeList.rank, SAVE_TAG);
	  ar & (coeffsPointer != 0);
	  // Storing Coefficients.
	  if(coeffsPointer) {
	    TensorT Coefficients(2*k, 2*k, 2*k);
	    arin & Coefficients;
	    ar & Coefficients;
	  }
	}
      }
    }
  }

  template <typename T>
  template <class Archive>
  void Function<T>::shadowManager_load(const char* f, const Archive& ar, OctTreeT* tree, Level n, Translation x, Translation y, Translation z, ProcessID remoteRank, Communicator& commFunc, const long partLevel, bool active_flag, bool have_child) 
  {
    //int nRemoteBranch;
    madness::archive::MPIOutputArchive arout(commFunc, remoteRank);
    // Send start branch data.
    arout & n & x & y & z;
    madness::archive::MPIInputArchive arin(commFunc, remoteRank);
    Level n_local;
    Level n_local2;
    Level n_local3;
    Translation x_local, y_local, z_local;
    Translation x_local2, y_local2, z_local2;
    Translation x_local3, y_local3, z_local3;
    ProcessID local_remoteRank;
    //bool remoteFlag;
    arout & have_child & active_flag;
    arout.flush();
    ar & n_local & x_local & y_local & z_local;
    bool inquireCoeffs;
    if(active_flag) { 
      ar & inquireCoeffs;
      arout & inquireCoeffs;
      arout.flush();
      if(inquireCoeffs) {
	Tensor<T> Coefficients(2*k, 2*k, 2*k);
	ar & Coefficients;
	arout & Coefficients;
	arout.flush();
      }
    }
    FORIJK( 
	   n_local2 = n_local + 1;
	   x_local2 = 2*x_local + i;
	   y_local2 = 2*y_local + j;
	   z_local2 = 2*z_local + k;
	   if(partLevel > 0 && n_local2 == partLevel) {
	     char ftest[256];
	     produceNewFilename3(f, partLevel, n_local2, x_local2, y_local2, z_local2, ftest);
	     TextFstreamInputArchive iar(ftest);
	     iar & have_child;
	     arout & have_child;
	     arout.flush();
	     if(have_child) {
	       iar & active_flag;
	       arin & local_remoteRank & n_local3 & x_local3 & y_local3 & z_local3;
	       if(local_remoteRank == 0) {
		 OctTreeT* treemaster = tree->find(n_local3, x_local3, y_local3, z_local3);
		 loadManager(f, iar, treemaster, commFunc, active_flag, 0, have_child);
	       }
	       else {
		 if(remoteRank != local_remoteRank) {
		   madness::archive::MPIOutputArchive arout2(commFunc, local_remoteRank);
		   arout2 & 1;
		   arout2.flush();
		 }
		 shadowManager_load(f, iar, tree, n_local3, x_local3, y_local3, z_local3, local_remoteRank, commFunc, partLevel, active_flag, have_child);
	       }
	       /*
        //OctTreeTPtrchild = tree -> child(i, j, k);
        //bool have_child;
        ar & have_child;
        arout & have_child;
        arout.flush();
        if(have_child) {
          ar & active_flag;
          arin & local_remoteRank & n_local & x_local & y_local & z_local;
          if(local_remoteRank == 0) {
            OctTreeT* treemaster = tree->find(n_local, x_local, y_local, z_local);
            loadManager(f, ar, treemaster, commFunc, active_flag, 0, have_child);
	       */
	     }
	   }
	   else {
	     ar & have_child;
	     arout & have_child;
	     arout.flush();
	     if(have_child) {
	       ar & active_flag;
	       arin & local_remoteRank & n_local3 & x_local3 & y_local3 & z_local3;
	       if(local_remoteRank == 0) {
		 OctTreeT* treemaster = tree->find(n_local3, x_local3, y_local3, z_local3);
		 loadManager(f, ar, treemaster, commFunc, active_flag, 0, have_child);
	       }
	       else {
		 if(remoteRank != local_remoteRank) {
		   madness::archive::MPIOutputArchive arout2(commFunc, local_remoteRank);
		   arout2 & 1;
		   arout2.flush();
		 }
		 shadowManager_load(f, ar, tree, n_local3, x_local3, y_local3, z_local3, local_remoteRank, commFunc, partLevel, active_flag, have_child);
	       }
	     }
	   }
	   );
  }

  template <typename T>
  void Function<T>::saveLoadWorker(OctTreeT* tree, Communicator& commFunc, bool save) 
  {
    int msg;
    madness::archive::MPIInputArchive arrecv(commFunc, 0);
    madness::archive::MPIOutputArchive arsend(commFunc, 0);
    while (1)
      {
        arrecv & msg;
        if (msg == -1) {
          break;
        }
        else if (msg == 1) {
          if(save) {
            Level n_c;
            Translation x_c, y_c, z_c;
            arrecv & n_c & x_c & y_c & z_c;
            //      MPI_Barrier(MPI_COMM_WORLD);
            OctTreeT* treeclient = tree->find(n_c, x_c, y_c, z_c);
            std::vector<localTreeMember> subtreeList; 
            localTreeList(subtreeList, treeclient);
            int nRemoteBranch = subtreeList.size();
            arsend & nRemoteBranch;
            arsend.flush();
	    //      MPI_Barrier(MPI_COMM_WORLD);
            for (int i = 0; i < nRemoteBranch; i++) {
              arsend & subtreeList[i];
            }
            arsend.flush();
	    //      MPI_Barrier(MPI_COMM_WORLD);
            sendRecvDataWorker_save(treeclient, commFunc);
	    //      MPI_Barrier(MPI_COMM_WORLD);
          }
          else {
            sendRecvDataWorker_load(tree, commFunc);
          }
        }
      }
  }
  
  template <typename T>
  void Function<T>::localTreeList(std::vector<localTreeMember> &subtreeList, const OctTreeT* tree)
  {
    localTreeMember branchList;
    if(tree) {
      branchList.have_child = true;
    }
    else {
      branchList.have_child = false;
    }
    if(branchList.have_child) {
      branchList.x      = tree->x();
      branchList.y      = tree->y();
      branchList.z      = tree->z();
      branchList.n      = tree->n();
      if(tree->rank() == -1) {
	branchList.rank   = comm()->rank();
      }
      else {
	branchList.rank   = tree->rank();
      }
      branchList.remote = tree->isremote();
      branchList.active = isactive(tree);
    }
    subtreeList.push_back(branchList);
    FORIJK( 
	   OctTreeTPtr child = tree->child(i,j,k);
	   if(child) {
	     localTreeList(subtreeList, child.get());
	   }
	   else {
	     if(tree->islocal()) {
	       localTreeMember branchList;
	       branchList.have_child = child;
	       branchList.x      = i;
	       branchList.y      = j;
	       branchList.z      = k;
	       branchList.rank   = 0;
	       branchList.remote   = false;
	       branchList.active   = false;
	       subtreeList.push_back(branchList);
	     }
	   }
	   );
  }

  template <typename T>
  void Function<T>::sendRecvDataWorker_save(OctTreeT* tree, Communicator& commFunc)
  {
    //madness::archive::MPIInputArchive arrecv(commFunc, 0);
    madness::archive::MPIOutputArchive arsend(commFunc, 0);
    if(isactive(tree)) {
      if ( tree->isremote() && !(tree->parent()) ) {
      }
      else {
	TensorT *t = coeff(tree);
	bool coeffsPointer = false;
	if(t) coeffsPointer = true;
	commFunc.Send(&coeffsPointer, 1, 0, SAVE_TAG);
	if(t) {
	  arsend & *t;
	  arsend.flush();
	}
      }
    }
    FOREACH_LOCAL_CHILD(OctTreeTPtr, tree, 
			//sendRecvDataWorker_save(child.get(), commFunc);
			sendRecvDataWorker_save(child, commFunc);
			);
  }

  template <typename T>
  void Function<T>::sendRecvDataWorker_load(const OctTreeT* treeclient, Communicator& commFunc)
  {
    madness::archive::MPIInputArchive arrecv(commFunc, 0);
    madness::archive::MPIOutputArchive arsend(commFunc, 0);
    Level n_c;
    Translation x_c, y_c, z_c;
    arrecv & n_c & x_c & y_c & z_c;
    OctTreeT* tree = treeclient->find(n_c, x_c, y_c, z_c);
    bool have_child, active_flag;
    arrecv & have_child & active_flag;
    if(active_flag) {
      set_active(tree);
      bool inquireCoeffs;
      arrecv & inquireCoeffs;
      if(inquireCoeffs) {
	Tensor<T> Coefficients(2*k, 2*k, 2*k);
	arrecv & Coefficients;
	set_coeff(tree, Coefficients);
      }
    }
    FORIJK( 
	   OctTreeTPtr child = tree->child(i, j, k);
	   arrecv & have_child;
	   if(!child) {
	     if(have_child) {
	       child = tree->insert_local_child(i,j,k);
	     }
	   }
	   if(child) {
	     if(child->rank() == -1) {
	       arsend & comm()->rank() & child->n() & child->x() & child->y() & child->z();
	     }
	     else {
	       arsend & child->rank() & child->n() & child->x() & child->y() & child->z();
	     }
	     arsend.flush();
	     if(islocal(child)) {
	       sendRecvDataWorker_load(child.get(), commFunc);
	     }
	   }
	   );
  }

  template <typename T>
  void Function<T>::load(const char* f, Communicator& comm)
  {
    TextFstreamInputArchive iar(comm.rank() ? 0 : f);
    long partLevel;
    if (comm.rank() == 0) {
      iar & partLevel;
      iar & FunctionDefaults::k;
      iar & FunctionDefaults::thresh;
      iar & FunctionDefaults::initial_level;
      iar & FunctionDefaults::max_refine_level;
      iar & FunctionDefaults::truncate_method;
      iar & FunctionDefaults::autorefine;
    }
    //comm.Bcast(partLevel, 0);
    comm.Bcast(FunctionDefaults::k, 0);
    comm.Bcast(FunctionDefaults::thresh, 0);
    comm.Bcast(FunctionDefaults::initial_level, 0);
    comm.Bcast(FunctionDefaults::max_refine_level, 0);
    comm.Bcast(FunctionDefaults::truncate_method, 0);
    comm.Bcast(FunctionDefaults::autorefine, 0);
    if (comm.rank() == 0) {
      bool active_flag;
      iar & active_flag;
      if(active_flag) loadManager(f, iar, tree(), comm, active_flag, partLevel, true);
      for ( int i = 1; i < comm.size(); i++) {
	archive::MPIOutputArchive arout(comm, i);
	arout & -1;
      }
      iar.close();
    }
    else {
      //saveLoadWorker(tree(), comm, false);
      saveLoadWorker(tree(), comm, false);
    }
  }
  
  template <typename T>
  template <class Archive>
  void Function<T>::loadManager(const char* f, const Archive& ar, OctTreeT* tree, Communicator& commFunc, bool active_flag, const long partLevel, bool have_child)
  {
    //    void Function<T>::loadManager(const char* f, const Archive& ar, OctTreeT* tree, Communicator& commFunc, bool active_flag, const long partLevel, bool have_child) {
    if (isremote(tree)) {
      madness::archive::MPIOutputArchive arout(commFunc, tree->rank());
      arout & 1;
      arout.flush();
      if(active_flag){
	set_active(tree);
      }
      else {
	set_inactive(tree);
      }
      shadowManager_load(f, ar, tree, tree->n(), tree->x(), tree->y(), tree->z(), tree->rank(), commFunc, partLevel, active_flag, have_child);
    }
    else {
      if(active_flag){
	set_active(tree);
      }
      else {
	set_inactive(tree);
      }
      Level n_local;
      Translation x_local, y_local, z_local;
      ar & n_local & x_local & y_local & z_local;
      if(active_flag) {
	bool inquireCoeffs;
	ar & inquireCoeffs;
	if(inquireCoeffs) {
	  TensorT t;
	  ar & t;
	  set_coeff(tree, t);
	}
      }
      FORIJK( 
	     //bool have_child;
	     OctTreeTPtr child = tree->child(i, j, k);
	     if(partLevel > 0 && tree->n() == (partLevel-1)) {
	       char ftest[256];
	       produceNewFilename(f, partLevel, child, ftest);
	       TextFstreamInputArchive iar(ftest);
	       iar & have_child;
	       if(have_child) {
		 //OctTreeTPtr child = tree->child(i, j, k);
		 if(!child) child = tree->insert_local_child(i,j,k);
		 //bool active_flag;
		 iar & active_flag;
		 loadManager(f, iar, child, commFunc, active_flag, partLevel, have_child);
	       }
	     }
	     else {
	       ar & have_child;
	       if(have_child) {
		 //OctTreeTPtr child = tree->child(i, j, k);
		 if(!child) child = tree->insert_local_child(i,j,k);
		 //bool active_flag;
		 ar & active_flag;
		 loadManager(f, ar, child, commFunc, active_flag, partLevel, have_child);
	       }
	     }
	     );
    }
  }
  
  template <typename T>
  void Function<T>::produceNewFilename(const char* f, const long partLevel, const OctTreeTPtr& tree, char newfilename[256])
  {
    char tmp[16];
    strcpy(newfilename, f);
    strcat(newfilename, "_");
    sprintf(tmp, "%d", tree->n());
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", tree->x());
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", tree->y());
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", tree->z());
    strcat(newfilename, tmp);
  }
  
  template <typename T>
  void Function<T>::produceNewFilename2(const char* f, const long partLevel, localTreeMember& subtreeList, char newfilename[256])
  {
    char tmp[16];
    strcpy(newfilename, f);
    strcat(newfilename, "_");
    sprintf(tmp, "%d", subtreeList.n);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", subtreeList.x);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", subtreeList.y);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", subtreeList.z);
    strcat(newfilename, tmp);
  }
  
  template <typename T>
  void Function<T>::produceNewFilename3(const char* f, const long partLevel, Level n, Translation x, Translation y, Translation z, char newfilename[256])
  {
    char tmp[16];
    strcpy(newfilename, f);
    strcat(newfilename, "_");
    sprintf(tmp, "%d", n);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", x);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", y);
    strcat(newfilename, tmp);
    strcat(newfilename, "_");
    sprintf(tmp, "%lu", z);
    strcat(newfilename, tmp);
  }
  template void Function<double>::save(const char* f, const long partLevel, Communicator& comm);

  template void Function<double>::load(const char* f, Communicator& comm);

  template void Function<double>::save_local<TextFstreamOutputArchive>( TextFstreamOutputArchive& ar);

  template void Function<double>::_save_local<TextFstreamOutputArchive>( TextFstreamOutputArchive& ar, const OctTreeT* tree); 

  template void Function<double>::load_local<TextFstreamInputArchive>( const TextFstreamInputArchive& ar); 

  template void Function<double>::_load_local<TextFstreamInputArchive>( const TextFstreamInputArchive& ar, OctTreeT* tree); 

  template void Function<double>::saveMain<TextFstreamOutputArchive>( const char*f, TextFstreamOutputArchive& ar, const long partLevel, Communicator& comm); 

  template void Function<double>::saveManager<TextFstreamOutputArchive>( const char* f, TextFstreamOutputArchive& ar, const OctTreeT* tree, const long partLevel, Communicator& commFunc);

  template void Function<double>::shadowManager_save<TextFstreamOutputArchive>( const char* f, TextFstreamOutputArchive& ar, const OctTreeT* tree, Level n, Translation x, Translation y, Translation z, ProcessID remoteRank, const long partLevel, Communicator& commFunc); 

  template void Function<double>::shadowManager_load<TextFstreamInputArchive>( const char* f, const TextFstreamInputArchive& ar, OctTreeT* tree, Level n, Translation x, Translation y, Translation z, ProcessID remoteRank, Communicator& commFunc, const long partLevel, bool active_flag, bool have_child); 

  template void Function<double>::loadManager<TextFstreamInputArchive>( const char* f, const TextFstreamInputArchive& ar, OctTreeT* tree, Communicator& commFunc, bool active_flag, const long partLevel, bool have_child); 

  template void Function<double>::dataSaveInShaManSave<TextFstreamOutputArchive>( const char* f, TextFstreamOutputArchive& ar, const OctTreeT* tree, localTreeMember& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc); 

  template void Function<double>::_shadowManager_save<TextFstreamOutputArchive>( const char* f, TextFstreamOutputArchive& ar, const OctTreeT* tree, std::vector<localTreeMember>& subtreeList, ProcessID remoteRank, const long partLevel, Communicator& commFunc, const int nRemoteBranch, int& iter); 

  //    template void Function<double>::produceNewFilename(const char* f, const long partLevel, const OctTreeTPtr& tree, char newfilename);     
}
