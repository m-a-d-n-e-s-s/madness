#ifndef SENDRECV_H
#define SENDRECV_H

/// \file sendrecv.h
/// \brief Defines send/recv for OctTree

#include <iostream>
#include <algorithm>
#include <cmath>
#include <list>

#include <mad_types.h>
#include <misc/print.h>
#include <misc/communicator.h>
#include <misc/shared_ptr.h>
#include <misc/misc.h>
#include <serialize/mpiar.h>
#include <octtree/octtree.h>
using namespace madness;


namespace madness {
    
    template <class T>
    void sendSubtree(OctTree<T> *tree, ProcessID me, ProcessID dest)
    {
	Communicator comm;
	archive::MPIOutputArchive ardest(comm, dest);

//	bool debug = false;
	bool debug = true;
	
	if (debug)
	    std::cout << "sendSubtree: about to send tree to " << dest << std::endl;
	if (debug)
	{
	    std::cout << "sendSubtree: subtree looks like this:" << std::endl;
	    tree->depthFirstTraverseAll();
	}
	ardest & *tree;
	
	if (debug)
	    std::cout << "sendSubtree: all done" << std::endl;
    }


    template <class T>
    void recvSubtree(OctTree<T> *t, ProcessID me, ProcessID source)
    {
	Communicator comm;
	madness::redirectio(comm);
	comm.print();
	archive::MPIInputArchive arsource(comm, source);

 	bool debug = true;
// 	bool debug = false;

	if (debug)
	{
	    std::cout << "recvSubtree: waiting to receive subtree from " << source << std::endl;
	}
	arsource & *t;
	if (debug)
	{
	    int s = t->tallyNodes();
	    std::cout << "recvSubtree: received subtree of size " 
		<< s << std::endl;
	}

	if (debug)
	    std::cout << "recvSubtree: all done!" << std::endl;
    }


    template <class T>
    void sendMsg(T msg, ProcessID me, ProcessID dest)
    {
	Communicator comm;
	madness::redirectio(comm);
	archive::MPIOutputArchive ardest(comm, dest);
	
	ardest & msg;
    }

    template <class T>
    void recvMsg(T *msg, ProcessID me, ProcessID source)
    {
	Communicator comm;
	madness::redirectio(comm);
	archive::MPIInputArchive arsource(comm, source);

	arsource & *msg;
    }

    template <class T>
    struct less
    {
	bool operator() (const T p1, const T p2)
	{
	    if (!p1)
		return true;
	    if (!p2)
		return false;
	    return *p1 < *p2;
	}
    };

}

#endif
