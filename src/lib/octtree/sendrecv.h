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

#define EXIT_CODE -1
#define COMPUTE_LOCAL_COST 1


namespace madness {
    
    class RootList
    {
	public:
	    Translation x, y, z;
	    Level n;
	    ProcessID current_owner, future_owner;
	    RootList(): x(0), y(0), z(0), n(-1), current_owner(0), future_owner(0) {};
	    RootList(Translation x0, Translation y0, Translation z0, Level n0, ProcessID current,
		ProcessID future): x(x0), y(y0), z(z0), n(n0), current_owner(current),
		future_owner(future)
	    {};
	    template <class T>
	    RootList(OctTree<T> *t, ProcessID current, ProcessID future): x(t->x()), y(t->y()), 
		z(t->z()), n(t->n()), current_owner(current), future_owner(future)
	    {};

	    template <class Archive>
            inline void serialize(const Archive& ar) {ar & x & y & z & n & current_owner & future_owner;}

	    friend bool operator < (const RootList& t1, const RootList& t2)
	    {
		if (t1.n > t2.n)
		    return true;
		else if (t1.n < t2.n)
		    return false;
		else
		{
		    if (t1.x < t2.x)
			return true;
		    else if (t1.x > t2.x)
			return false;
		    else
		    {
			if (t1.y < t2.y)
			    return true;
			else if (t1.y > t2.y)
			    return false;
			else
			{
			    if (t1.z < t2.z)
				return true;
			    else
				return false;
			}
		    }
		}
	    }
    };

    template <class T>
    void parallelPartition(OctTree<T> *root)
    {
	Communicator comm;
	ProcessID me, np;
	std::vector<RootList> *globalList = new std::vector<RootList>(), *localList = new std::vector<RootList>();
	
	bool debug = true;

	me = comm.rank();
	np = comm.nproc();

	Cost cost, idealPartition, accumulate, sofar, subtotal;
	
    }

    template <class T>
    Cost computeGlobalCost(RootList root, OctTree<T> *tree, std::vector<OctTree<T>* > const treeList)
    {
	Communicator comm;
	int me = comm.rank(), np = comm.nproc();
	Cost total = 0;
	std::cout << "computeGlobalCost: beginning of function" << std::endl;
	if (me == 0)
	{
	    total = globalCostManager(root, tree, treeList);
	    // then, shut down infinite loop of worker
	    std::cout << "computeGlobalCost: sending exit code" << std::endl;
	    for (int i = 1; i < np; i++)
	    {
		archive::MPIOutputArchive arout(comm, i);
		arout & EXIT_CODE;
	    }
	}
	else
	{
	    total = globalCostWorker(root, tree, treeList);
	}
	return total;
    }

    template <class T>
    Cost globalCostManager(RootList root, OctTree<T> *tree, std::vector<OctTree<T>* > const treeList)
    {
	Communicator comm;
	int me = comm.rank(); 
	int rank, nRemoteKids, i;
	RootList remoteRoot;
	Cost total = 0, subtotal = 0;

	bool debug = true;
	
	if (debug)
	    std::cout << "globalCostManager: at beginning of function" << std::endl;
	
	if (root.current_owner == me)
	{
	    if (debug)
	    {
		std::cout << "globalCostManager: I own this root" << std::endl;
	    }
	    int npieces = treeList.size();
            for (i = 0; i < npieces; i++)
            {
                tree = treeList[i]->find(root.n, root.x, root.y, root.z);
                if ((tree) && (tree->n() == root.n) && (tree->x() == root.x) && 
			(tree->y() == root.y) && (tree->z() == root.z))
                {
                    break;
                }
            }
	    
	    if (tree->isParent())
	    {
		FOREACH_CHILD(OctTree<T>, tree,
		    if (child->isremote())
			rank = child->rank();
		    else
			rank = me;
		    if (debug)
		    {
			std::cout << "globalCostManager: sending child " <<
				"(" << child->x() << "," << child->y() << "," << child->z() << ")"
				<< "to globalCostManager function" << std::endl;
		    }
		    total += globalCostManager(RootList(child->x(), child->y(), child->z(),
				child->n(), rank, rank), child, treeList);
		);
	    }
	    subtotal = tree->getCost();
	    if (debug)
	    {
		std::cout << "globalCostManager: subtotal from local tree = "
			<< subtotal << std::endl;
	    }
//	    tree->setLocalSubtreeCost(subtotal);
	    total += subtotal;
	    subtotal = 0;
	}
	else
	{
	    if (debug)
	    {
		std::cout << "globalCostManager: I don't own this root; processor " <<
			root.current_owner << " does" << std::endl;
	    }
	    archive::MPIOutputArchive arout(comm, root.current_owner);
	    arout & COMPUTE_LOCAL_COST & root;
	    archive::MPIInputArchive arin(comm, root.current_owner);
	    arin & subtotal & nRemoteKids;
	    if (debug)
	    {
		std::cout << "globalCostManager: subtotal received from remote tree = "
			<< subtotal << std::endl;
	    }
	    total += subtotal;
	    for (i = 0; i < nRemoteKids; i++)
	    {
		arin & remoteRoot;
		total += globalCostManager(remoteRoot, tree, treeList);
	    }
	}
	return total;
    }
    
    template <class T>
    Cost globalCostWorker(RootList root, OctTree<T> *tree, std::vector<OctTree<T>* > const treeList)
    {
	Communicator comm;
	int me = comm.rank(), npieces = treeList.size();
	archive::MPIInputArchive arrecv(comm, 0);
	archive::MPIOutputArchive arsend(comm, 0);

	std::cout << "globalCostWorker: at beginning of function" << std::endl;

	while (1) // infinite loop
	{
	    int msg, i;
            OctTree<T> *t = new OctTree<T>();
            arrecv & msg;
            if (msg == EXIT_CODE)
            {
                std::cout << "Processor " << me << ": received exit code " << msg << std::endl;
                break;
            }
            else if (msg == COMPUTE_LOCAL_COST)
            {
		Cost cost = 0;
		std::vector<RootList> rootList = new std::vector<RootList>();
                std::cout << "Processor " << me << ": received compute local cost code " << msg <<
                        std::endl;
                arrecv & root;
                for (i = 0; i < npieces; i++)
                {
                    t = treeList[i]->find(root.n, root.x, root.y, root.z);
                    if ((t) && (t->n() == root.n) && (t->x() == root.x) && (t->y() == root.y) &&
                        (t->z() == root.z))
                    {
                        break;
                    }
                }
		std::cout << "Processor " << me << ": found tree: " <<
		    "(" << t->x() << "," << t->y() << "," << t->z() << "," << t->n() << ")"
		    << std::endl;
                cost += t->computeLocalCost(&rootList); 
		int rlsize = rootList.size();
                arsend & cost & rlsize; 
		for (i = 0; i < rlsize; i++)
		{
                    arsend & rootList[i];
		}
            }
            else
            {
                std::cout << "Processor " << me << ": received unknown code " << msg << std::endl;
            }
	}
	return 0;
    }

    template <class T>
    void exchangeTrees(std::vector<RootList> *globalList, std::vector<OctTree<T>* > *treeList)
    {
	std::vector<RootList> *localList = new std::vector<RootList>();
	Communicator comm;
	madness::redirectio(comm);
	int me = comm.rank(), nproc = comm.nproc();
	int glength = 0, tlength = 0, i, j;
	RootList root;

//	bool debug = true;
	bool debug = false;

	if (globalList)
	    glength = globalList->size();
	if (treeList)
	    tlength = treeList->size();
	
	if (debug)
	{
	    std::cout << "exchangeTrees: about to bcast glength = " << glength << std::endl;
	}
	comm.Bcast(&glength, 1, 0);
	if (debug)
	{
	    std::cout << "exchangeTrees: after bcast glength = " << glength << std::endl;
	    if (me == 0)
	    {
		for (i = 0; i < glength; i++)
		{
		    std::cout << "    n = " << (*globalList)[i].n << ", (" << (*globalList)[i].x <<
				"," << (*globalList)[i].y << "," << (*globalList)[i].z << ")" << std::endl;
		}
	    }
	}
	for (i = 0; i < glength; i++)
	{
	    if (me == 0)
	    {
		root = (*globalList)[i];
		if (debug)
		{
		    std::cout << "exchangeTrees: root " << i << ": n = " << root.n << ", " <<
			"(" << root.x << "," << root.y << "," << root.z << ")" << "current = "
			<< root.current_owner << ", future = " << root.future_owner << std::endl;
		    std::cout << "exchangeTrees: now root will be broadcast: " << std::endl;
		}
/*
		for (j = 1; j < nproc; j++)
		{
		    archive::MPIOutputArchive arout(comm, j);
		    arout & root;
		}
*/
	    }
	    else
	    {
//		archive::MPIInputArchive arin(comm, 0);
//		arin & root;
	    }
	    comm.Bcast(&root, 1, 0);
	    if (debug)
	    {
		std::cout << "exchangeTrees: root " << i << ": n = " << root.n << ", " <<
			"(" << root.x << "," << root.y << "," << root.z << ")" << "current = "
			<< root.current_owner << ", future = " << root.future_owner << std::endl;
	    	std::cout << "exchangeTrees: after bcast root number " << i << std::endl;
	    }
	    if (me != 0)
	    {
		globalList->push_back(root);
	    }
	    if ((root.current_owner == me) || (root.future_owner == me))
	    {
		if (debug)
		{
		    std::cout << "exchangeTrees: I have something to do with tree " <<
			"(" << root.x << "," << root.y << "," << root.z << "," << root.n << ")"
			<< std::endl;
		}
//		localList.push_back(root);
		localList->push_back(root);
	    }
	}
//	int llength = localList.size();
	int llength = localList->size();
	if (debug)
	{
	   std::cout << "exchangeTrees: my localList is of length " << llength << std::endl;
	}
	for (i = 0; i < llength; i++)
	{
	    if (debug)
		std::cout << std::endl;
//	    if (localList[i].current_owner == localList[i].future_owner)
	    if ((*localList)[i].current_owner == (*localList)[i].future_owner)
	    {
		// do nothing, unless localList[i] is not currently a local root
		if (debug)
		{
		    std::cout << "exchangeTrees: current_owner = future_owner" << std::endl;
		}
		OctTree<T> *t = new OctTree<T>();
		for (j = 0; j < tlength; j++)
		{
		    if (debug)
		    {
			std::cout << "exchangeTrees: looking for tree in treeList[" << j << "]"
				<< std::endl;
			(*treeList)[j]->depthFirstTraverseAll();
			std::cout << "exchangeTrees: will I find it?" << std::endl;
		    }
//		    t = (*treeList)[j]->find(localList[i].n, localList[i].x, localList[i].y, localList[i].z);
		    t = (*treeList)[j]->find((*localList)[i].n, (*localList)[i].x, (*localList)[i].y, 
				(*localList)[i].z);
//		    if ((t) && (t->n() == localList[i].n) && (t->x() == localList[i].x) &&
//				(t->y() == localList[i].y) && (t->z() == localList[i].z))
		    if ((t) && (t->n() == (*localList)[i].n) && (t->x() == (*localList)[i].x) &&
				(t->y() == (*localList)[i].y) && (t->z() == (*localList)[i].z))
			break;
		    else
			t = 0;
		}
		if (t != (*treeList)[j])
		{
		    OctTree<T> *p = new OctTree<T>();
		    p = t->parent();
//		    Translation x = localList[i].x, y = localList[i].y, z = localList[i].z;
		    Translation x = (*localList)[i].x, y = (*localList)[i].y, z = (*localList)[i].z;
		    x -= 2*(x/2); y -= 2*(y/2); z -= 2*(z/2);
		    if (debug)
		    {
			std::cout << "exchangeTrees: about to insert remote child " <<
			    "(" << x << "," << y << "," << z << ")" << std::endl;
		    }
//		    p->insert_remote_child(x, y, z, localList[i].future_owner);
		    p->insert_remote_child(x, y, z, (*localList)[i].future_owner);
		    treeList->push_back(t);
		    tlength++;
		}
		if (debug)
		{
		    std::cout << "exchangeTrees: end of current = future: tlength = " << tlength 
			<< std::endl;
		}
	    }
//	    else if (localList[i].current_owner == me)
	    else if ((*localList)[i].current_owner == me)
	    {
		// send it to its new owner
		// we assume that the element is indeed owned by this processor
		if (debug)
		{ std::cout << "exchangeTrees: I am current_owner" << std::endl;
		}
//		archive::MPIOutputArchive arsend(comm, localList[i].future_owner);
		OctTree<T> *t = new OctTree<T>();
		for (j = 0; j < tlength; j++)
		{
		    if (debug)
		    {
			std::cout << "exchangeTrees: looking for tree in treeList[" << j << "]"
				<< std::endl;
			(*treeList)[j]->depthFirstTraverseAll();
			std::cout << "exchangeTrees: will I find it?" << std::endl;
		    }
		    t = (*treeList)[j]->find((*localList)[i].n, (*localList)[i].x, (*localList)[i].y, 
						(*localList)[i].z);
		    if (debug)
		    {
			if (t)
			{
			    std::cout << "exchangeTrees: t = (" << t->x() << "," << t->y() << "," <<
				t->z() << ")" << std::endl;
			}
			else
			{
			    std::cout << "exchangeTrees: t not yet found" << std::endl;
			}
		    }
//		    if ((t) && (t->n() == localList[i].n) && (t->x() == localList[i].x) &&
//				(t->y() == localList[i].y) && (t->z() == localList[i].z))
		    if ((t) && (t->n() == (*localList)[i].n) && (t->x() == (*localList)[i].x) &&
				(t->y() == (*localList)[i].y) && (t->z() == (*localList)[i].z))
		    {
			if (debug)
			{
			    std::cout << "exchangeTrees: found tree in treeList[" << j << "]:"
				" (" << t->x() << "," << t->y() << "," << t->z() << ")"<< std::endl;
			}
			break;
		    }
		    else
		    {
			t = 0;
			if (debug)
			{
			    std::cout << "exchangeTrees: did not find tree in treeList[" << j << "]"
				<< std::endl;
			}
		    }
		}
		if (debug)
		{
		    std::cout << "exchangeTrees: I am about to send the tree to " <<
//			localList[i].future_owner << std::endl;
			(*localList)[i].future_owner << std::endl;
		    std::cout << "exchangeTrees: t = : " << (t) << std::endl;
		}
//		double x = 3.14159;
//		arsend & x;
//		arsend & *t;
		sendSubtree(t, me, (*localList)[i].future_owner);
//		sendSubtree(t, me, localList[i].future_owner);
		if (debug)
		{
		    std::cout << "exchangeTrees: sent tree t" << std::endl;
		}
//		OctTree<T> *p = t->parent();
		OctTree<T> *p = new OctTree<T>();
		p = t->parent();
		if ((p) && (p->islocal()))
		{
//		    Translation x = localList[i].x, y = localList[i].y, z = localList[i].z;
		    Translation x = (*localList)[i].x, y = (*localList)[i].y, z = (*localList)[i].z;
		    x -= 2*(x/2); y -= 2*(y/2); z -= 2*(z/2);
//		    p->insert_remote_child(x, y, z, localList[i].future_owner);
		    p->insert_remote_child(x, y, z, (*localList)[i].future_owner);
		}
		else
		{
		    treeList->erase((treeList->begin()+j));
		    tlength--;
		}
		if (debug)
		{
		    std::cout << "exchangeTrees: at end of I am current_owner, tlength = " <<
			tlength << std::endl;
		}
	    }
	    else
	    {
		// receive from current owner
		if (debug)
		{
		    std::cout << "exchangeTrees: I am future_owner" << std::endl;
		}
//		archive::MPIInputArchive arrecv(comm, localList[i].current_owner);
		OctTree<T> *t = new OctTree<T>();
//		double x;
//		arrecv & x;
//		arrecv & *t;
//		recvSubtree(t, me, localList[i].current_owner);
		recvSubtree(t, me, (*localList)[i].current_owner);
		if (debug)
		{
		    std::cout << "exchangeTrees: at least I received something!" << std::endl;
		}
//		treeList->push_back(t);
	    }
	
	}
    }


/*
    template <class T>
    Cost computeLocalCost(Translation x, Translation y, Translation z, Level n, std::vector<RootList> *roots)
    {
	// First, find the node we're talking about

	// Then, compute its local cost, returning that and a vector of roots
    }
*/

    template <class T>
    void sendSubtree(OctTree<T> *tree, ProcessID me, ProcessID dest)
    {
	Communicator comm;
	archive::MPIOutputArchive ardest(comm, dest);

	bool debug = false;
//	bool debug = true;
	
	if (debug)
	    std::cout << "sendSubtree: about to send tree to " << dest << std::endl;
	if (debug)
	{
	    std::cout << "sendSubtree: subtree looks like this:" << std::endl;
	    tree->depthFirstTraverseAll();
	    std::cout << "sendSubtree: done with depthFirstTraverseAll" << std::endl;
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

// 	bool debug = true;
 	bool debug = false;

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
