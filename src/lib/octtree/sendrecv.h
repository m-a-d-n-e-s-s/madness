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
    void sendSubtree(std::vector<OctTree<T>* > *treeList, ProcessID me, ProcessID dest)
    {
	int listsize = treeList->size();
	Communicator comm;
	archive::MPIOutputArchive ardest(comm, dest);

	bool debug = false;
	
	if (debug)
	    std::cout << "sendSubtree: listsize = " << listsize << std::endl;

	ardest & listsize;

	for (int i = 0; i < listsize; i++)
	{
	    if (debug)
		std::cout << "sendSubtree: about to send treeList[" << i << "]" << std::endl;
	    OctTree<T> *tree = (*treeList)[i];
	    OctTree<T> tree2;
	    /* Convert tree into form that can be passed to another process */
	    tree2 = copyTree(tree);
	    if (debug)
	    {
		std::cout << "sendSubtree: subtree looks like this:" << std::endl;
		tree2.depthFirstTraverseAll();
	    }
	    ardest & tree2;
	    if (debug)
		std::cout << "sendSubtree: sent treeList[" << i << "]" << std::endl;
	}
	if (debug)
	    std::cout << "sendSubtree: all done" << std::endl;
    }


    template <class T>
    void recvSubtree(std::vector<OctTree<T>* > *treeList, ProcessID me, ProcessID source)
    {
	Communicator comm;
	madness::redirectio(comm);
	comm.print();
	int vlen;
	archive::MPIInputArchive arsource(comm, source);

// 	bool debug = true;
 	bool debug = false;

	arsource & vlen;
	if (debug)
	   std::cout << "recvSubtree: vlen = " << vlen << std::endl;


	for (int i = 0; i < vlen; i++)
	{
	    OctTree<T> *t = new OctTree<T>();
	    if (debug)
		std::cout << "recvSubtree: about to receive subtree number " << i << std::endl;
	    arsource & *t;
	    if (debug)
	    {
		int s = t->tallyNodes();
		std::cout << "recvSubtree: received subtree number " << i << " of size " 
			<< s << std::endl;
	    }
	    treeList->push_back(t);
	    if (debug)
		std::cout << "recvSubtree: pushed subtree number " << i << std::endl;
	}

	if (debug)
	    std::cout << "recvSubtree: all done!" << std::endl;
    }

    template <class T>
    OctTree<T> copyTree(OctTree<T> *tree)
    {
//	bool debug = true;
	bool debug = false;
	OctTree<T> *treecopy = new OctTree<T>();
	OctTree<T> *current_tree;

	treecopy = new OctTree<T> (*tree);
	current_tree = treecopy;
	if (debug)
	{
	    std::cout << "copyTree: depthFirstTraverseAll() at very beginning" << std::endl;
	    current_tree->depthFirstTraverseAll();
	    std::cout << "copyTree: end of depthFirstTraverseAll()" << std::endl;
	}

	if (debug)
	{
	    std::cout << "copyTree: tree (" << tree->x() << "," << tree->y() << "," << tree->z() << ")"
		<< std::endl;
	    std::cout << "copyTree: treecopy (" << treecopy->x() << "," << 
		treecopy->y() << "," << treecopy->z() << ")" << std::endl;
	}
	current_tree->setRemote(false);
	current_tree->setRank(-1);
	if (debug)
	{
	    std::cout << "copyTree: rank of current_tree has been set" << std::endl;
	}
	FORIJK(
	    if (tree->child(i,j,k))
	    {
	        if (debug)
		    std::cout << std::endl;
		OctTree<T> *child = tree->child(i,j,k);
	    	if (debug)
	    	{
		    std::cout << "copyTree: made child(" << i << "," << j << "," << k << ")" << std::endl;
	    	}
		/* if the child is remote */
		if (child->rank() != tree->rank())
		{
		    if (debug)
		    {
			std::cout << "copyTree: child (" << child->x() << "," << child->y() << ","
				<< child->z() << ") is remote" << std::endl;
		    }
		    current_tree->insert_remote_child(i, j, k, child->rank());
		    if (debug)
		    {
			std::cout << "copyTree: inserted remote child" << std::endl;
		    }
		}
		/* if the child is local */
		else
		{
		    if (debug)
		    {
			std::cout << "copyTree: child (" << child->x() << "," << child->y() << ","
				<< child->z() << ") is local" << std::endl;
		    }
		    OctTree<T> current_child = copyTree(child);
		    if (debug)
		    {
			std::cout << "copyTree: after returning from copyTree local current child,"
				<< " current_child looks like:" << std::endl;
			current_child.depthFirstTraverseAll();
			std::cout << "copyTree: end of traversal of current_child" << std::endl;
		    }
		    current_tree->insert_local_child(i, j, k, current_child);
		    if (debug)
		    {
			std::cout << "copyTree: back from recursive call to copyTree" << std::endl;
		    }
		}
	    }
	);
	if ((tree->parent()) && (tree->parent()->rank() != tree->rank()))
	{
	    if (debug)
	    {
		std::cout << "copyTree: remote parent (" << tree->parent()->x() << "," 
			<< tree->parent()->y() << "," << tree->parent()->z() << ")" 
			<< " rank " << tree->parent()->rank() << std::endl;
	    }
/********************************************************************************/
/* This part is still a work in progress.  I need to figure out how to get the 	*/
/* remote parent node included in the tree that will be passed 			*/
/********************************************************************************/
/*
	    OctTree<T> *parent = treecopy->insert_remote_parent(tree->parent()->rank());
	    if (debug)
	    {
		std::cout << "copyTree: parent has been made: parent is (" << parent->x() <<
			"," << parent->y() << "," << parent->z() << "), layer " << 
			parent->n() << " Remote " << parent->isremote() << " Processor " << 
			parent->rank() << std::endl;
	    }
	    treecopy = parent;
*/
	}
	if (debug)
	{
	    std::cout << "copyTree: depthFirstTraverseAll() at very end" << std::endl;
	    treecopy->depthFirstTraverseAll();
	    std::cout << "copyTree: end of depthFirstTraverseAll()" << std::endl;
	}
	return *treecopy;
    }

}

#endif
