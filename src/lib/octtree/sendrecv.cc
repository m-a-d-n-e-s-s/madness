//#ifndef SENDRECV_H
//#define SENDRECV_H

/// \file sendrecv.cc
/// \brief Implements send/recv for OctTree

#include <octtree/sendrecv.h>
using namespace madness;

namespace madness {

    /********************************************************************************/
    /* parallelPartition: processor 0 is manager process that determines cost of 	*/
    /* tree globally, ideal partition size, and repartitions tree.  Uses producer	*/
    /* model for workers, where workers wait for message from manager about what to */
    /* do and respond accordingly.  In the end, this returns the list of trees to	*/
    /* be moved, globalList, which should then be fed into exchangeTrees in order 	*/
    /* to complete the process.							*/
    /********************************************************************************/

    template <class T>
    void parallelPartition(OctTree<T> *root, std::vector<RootList> *globalList,
                           std::vector<OctTree<T>* > *treeList) {
        Communicator comm;
        ProcessID me;

//	bool debug = true;
        bool debug = false;

        me = comm.rank();

// if I am processor 0, then I am the manager
        if (me == 0) {
            if (debug) {
                std::cout << "parallelPartition: before partitionManager" << std::endl;
            }
            partitionManager(root, globalList, treeList);
            if (debug) {
                std::cout << "parallelPartition: back from partitionManager" << std::endl;
                int size = globalList->size();
                for (int i = 0; i < size; i++) {
                    std::cout << "Subtree: " << std::endl;
                    std::cout << "Layer " << (*globalList)[i].n << ": " <<
                    "(" << (*globalList)[i].x << "," << (*globalList)[i].y << "," <<
                    (*globalList)[i].z << ")" << std::endl;

                }
                std::cout << "parallelPartition: all done" << std::endl << std::endl;
            }
        }
// otherwise, I am a worker
        else {
            if (debug) {
                std::cout << "parallelPartition: before partitionWorker" << std::endl;
            }
            partitionWorker(treeList);
            if (debug) {
                std::cout << "parallelPartition: back from partitionWorker" << std::endl;
            }
        }

    }

    template <class T>
    void partitionManager(OctTree<T> *root, std::vector<RootList> *globalList,
                          std::vector<OctTree<T>* > *treeList) {
        Communicator comm;
        ProcessID np = comm.nproc();
        std::vector<RootList> *localList = new std::vector<RootList>();
        RootList *rl = new RootList(root, 0, 0);
        Cost cost, idealPartition, spaceleft;
        cost = computeGlobalCost(*rl, root, *treeList);
        idealPartition = (Cost) floor((1.0*cost)/np);

//	bool debug = true;
        bool debug = false;

// we work down from processor np-1 to processor 0, to ensure that global root of tree
// will be on processor 0
        for (int p = np-1; p >= 0; p--) {
            if (debug) {
                std::cout << std::endl << std::endl;
                std::cout << "      BEGINNING OF PARTITION NUMBER " << p << std::endl;
            }
            if (p != 0) {
//  the idea here is that updating the partition size will smooth out any imbalances.
//  for example, imagine 4 processors and 7 cost.  ideal_0 = 1, resulting in loads of
//  1, 1, 1, 4.  If instead we update the ideal size when it becomes larger, we would have
//  a load of 1, 2, 2, 2 (much better balance).
                Cost iP = (Cost) floor((1.0*cost)/(p+1));
                if (iP > idealPartition) {
                    idealPartition = iP;
                }
                if (debug) {
                    std::cout << "partitionManager: ideal partition size = " << idealPartition
                    << std::endl;
                }
            }
            spaceleft = idealPartition;
            if (debug) {
                std::cout << "partitionManager: spaceleft = " << spaceleft << ", " <<
                "root->getLocalSubtreeCost = " << root->getLocalSubtreeCost() << std::endl;
            }
            if ((root->getLocalSubtreeCost() < spaceleft) || (p == 0)) {
                if (debug) {
                    std::cout << "partitionManager: enough space for subtree " <<
                    "(" << root->x() << "," << root->y() << "," << root->z() << ");"
                    << " put it on list" << std::endl;
                }
                int rank;
                if (root->islocal())
                    rank = 0;
                else
                    rank = root->rank();
                globalList->push_back(RootList(root, rank, p));
                root->setSendto(p);
                if (debug) {
                    std::cout << "partitionManager: subtree placed on list" << std::endl;
                }
                std::vector<RootList> remoteList;
                remoteList = findRemoteChildren(root, p);
                int rsize = remoteList.size();
                if (debug) {
                    std::cout << "partitionManager: found " << rsize << " remote subtrees for list"
                    << std::endl;
                }
                for (int i = 0; i < rsize; i++) {
                    int gsize = globalList->size();
                    bool already_there = false;
                    for (int j = 0; j < gsize; j++) {
                        if (remoteList[i].equals((*globalList)[j])) {
                            already_there = true;
                            break;
                        }
                    }
                    if (!(already_there)) {
                        globalList->push_back(remoteList[i]);
                        if (debug) {
                            std::cout << "partitionManager: put remote subtree n = " << remoteList[i].n
                            << ", (" << remoteList[i].x << "," << remoteList[i].y << "," <<
                            remoteList[i].z << "), " << remoteList[i].current_owner << " -> "
                            << remoteList[i].future_owner << " on globalList" << std::endl;
                        }
                    }
                }
            } else {
                if (debug) {
                    std::cout << "partitionManager: not enough space for subtree; calling "
                    << "globalPartitionSubtree" << std::endl;
                }
                *localList = globalPartitionSubtree(root, &spaceleft, p);
                if (debug) {
                    std::cout << "partitionManager: returned from globalPartitionSubtree" << std::endl;
                }
                int size = localList->size();
                for (int i = 0; i < size; i++) {
                    globalList->push_back((*localList)[i]);
                }
                if (debug) {
                    std::cout << "partitionManager: added " << size << " new entries to globalList"
                    << std::endl;
                }
            }
// recompute the remaining cost at every iteration
            cost = computeGlobalCost(*rl, root, *treeList);
        }
        for (int p = 1; p < np; p++) {
            archive::MPIOutputArchive arsend(comm, p);
            if (debug) {
                std::cout << "partitionManager: about to send exit code " << EXIT_CODE <<
                " to processor " << p << std::endl;
            }
            arsend & EXIT_CODE;
            if (debug) {
                std::cout << "partitionManager: sent exit code " << EXIT_CODE <<
                " to processor " << p << std::endl;
            }
        }
    }


    template <class T>
    std::vector<RootList> globalPartitionSubtree(OctTree<T> *root, Cost *spaceleft, int p) {
        Communicator comm;
        ProcessID me = comm.rank();
        std::vector<RootList> tmpList, localList;
        RootList tmproot;
        Cost costofroot;
        int success;
        bool keepgoing = true;
//	bool debug = true;
        bool debug = false;

        if (debug) {
            std::cout << "globalPartitionSubtree: at beginning of function, spaceleft = " <<
            *spaceleft << std::endl;
        }
        FOREACH_CHILD(OctTree<T>, root,
                      if ((*spaceleft == 0) || (child->getSendto() != -1)) {
                      keepgoing = false;
                      if (_debug) {
                              std::cout << "globalPartitionSubtree: don't bother with child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                      }
                      else {
                          keepgoing = true;
                      }
                      if ((child->islocal()) && (keepgoing)) {
                      if (_debug) {
                              std::cout << "globalPartitionSubtree: keep going with child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                          Cost cost = child->getLocalSubtreeCost();
                          if (_debug) {
                              std::cout << "globalPartitionSubtree: child (" << child->x() << "," <<
                              child->y() << "," << child->z() << "), with cost " << cost <<
                              " is local" << std::endl;
                          }
                          if (cost <= *spaceleft) {
                              *spaceleft -= cost;
                              tmpList.push_back(RootList(child, me, p));
                              child->setSendto(p);
                              if (_debug) {
                                  std::cout << "globalPartitonSubtree: added child to list" << std::endl;
                              }
                              std::vector <RootList> remoteList;
                              remoteList = findRemoteChildren(child, p);
                              int rsize = remoteList.size();
                              if (_debug) {
                                  std::cout << "globalPartitionSubtree: found " << rsize << "  new remote children"
                                  << std::endl;
                              }
                              for (int i = 0; i < rsize; i++) {
                                  tmpList.push_back(remoteList[i]);
                              }
                          } else {
                              if (_debug) {
                                  std::cout << "globalPartitonSubtree: child cost too big for partition"
                                  << std::endl;
                              }
                              if ((child->isParent()) && (keepgoing))
                                  tmpList = globalPartitionSubtree(child, spaceleft, p);
                          }
                          int size = tmpList.size();
                          for (int i = 0; i < size; i++) {
                              localList.push_back(tmpList[i]);
                          }
                          tmpList.clear();
                      } else if (keepgoing) {
                      if (_debug) {
                              std::cout << "globalPartitionSubtree: keep going with remote child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                          archive::MPIOutputArchive arsend(comm, child->rank());
                          archive::MPIInputArchive arrecv(comm, child->rank());
                          if (_debug) {
                              std::cout << "globalPartitionSubtree: about to send GLOB_PART_SUB_SIGNAL "
                              << GLOB_PART_SUB_SIGNAL << " to processor " << child->rank() <<
                              " concerning partition " << p << " and space left " << *spaceleft
                              << std::endl;
                          }
                          arsend & GLOB_PART_SUB_SIGNAL & RootList(child,child->rank(), p) & *spaceleft;
                          arrecv & success;
                          if (_debug) {
                              std::cout << "globalPartitionSubtree:  about to add " << success <<
                              " new roots to list" << std::endl;
                          }
                          if (success > 0) {
                              for (int i = 0; i < success; i++) {
                                  arrecv & tmproot & costofroot;
                                  *spaceleft -= costofroot;
                                  tmproot.future_owner = p;
                                  localList.push_back(tmproot);
                                  if (_debug) {
                                      std::cout << "globalPartitionSubtree: added root " <<
                                      "(" << tmproot.x << "," << tmproot.y << "," << tmproot.z << ")"
                                      << " with cost " << costofroot << " to the localList" << std::endl;
                                  }
                              }
                          }
                          int rsize;
                          RootList remtmp;
                          std::vector<RootList> remoteList;
                          arrecv & rsize;

                          for (int i = 0; i < rsize; i++) {
                              arrecv & remtmp;
                              remtmp.future_owner = p;
                              remoteList.push_back(remtmp);
                          }
                          for (int i = 0; i < rsize; i++) {
                              int lsize = localList.size();
                              bool already_there = false;
                              for (j = 0; j < lsize; j++) {
                                  // if this node is already on the list
                                  if (remoteList[i].equals(localList[j])) {
                                      already_there = true;
                                      break;
                                  }
                              }
                              if (already_there) {
                                  if (_debug) {
                                      std::cout << "globalPartitionSubtree: put remote subtree n = "
                                      << remoteList[i].n << ", (" << remoteList[i].x << "," <<
                                      remoteList[i].y << "," << remoteList[i].z << "), " <<
                                      remoteList[i].current_owner << " -> " <<
                                      remoteList[i].future_owner << " on globalList" << std::endl;
                                  }
                              } else {
                                  if (_debug) {
                                      std::cout << "globalPartitionSubtree: put remote subtree n = "
                                      << remoteList[i].n << ", (" << remoteList[i].x << "," <<
                                      remoteList[i].y << "," << remoteList[i].z << "), " <<
                                      remoteList[i].current_owner << " -> " <<
                                      remoteList[i].future_owner << " on globalList" << std::endl;
                                  }
                                  localList.push_back(remoteList[i]);
                              }
                          }
                      }
                     );
        if (debug) {
            std::cout << "globalPartitionSubtree: at end of function" << std::endl;
        }
        return localList;
    }


    template <class T>
    void partitionWorker(std::vector<OctTree<T>* > *treeList) {
        Communicator comm;
        int signal;
        int spaceleft;
        RootList rootToFind;
        std::vector<RootList> localList;
        std::vector<Cost> costList;
        bool keepGoing = true;
//	bool debug = true;
        bool debug = false;

        if (debug) {
            std::cout << "partitionWorker: before infinite loop" << std::endl;
        }
        archive::MPIInputArchive arrecv(comm, 0);
        archive::MPIOutputArchive arsend(comm, 0);

        while (keepGoing) // infinite loop
        {
            if (debug) {
                std::cout << "partitionWorker: at beginning of loop" << std::endl;
            }
            arrecv & signal;
            if (debug) {
                std::cout << "partitionWorker: after receiving signal " << signal << std::endl;
            }

// do various things based on which signal received
            switch (signal) {
            case GLOBAL_COST:
            case COMPUTE_LOCAL_COST:
                if (debug) {
                    std::cout << "partitionWorker: received signal GLOBAL_COST" << std::endl;
                }
                globalCostWorker(rootToFind, *treeList); // put the right stuff inside parens
                break;
            case GLOB_PART_SUB_SIGNAL:
                if (debug) {
                    std::cout << "partitionWorker: received signal GLOB_PART_SUB_SIGNAL" << std::endl;
                }
                arrecv & rootToFind & spaceleft;
                // find the root
                int lsize = treeList->size();
                OctTree<T> *tree = new OctTree<T>();
                for (int i = 0; i < lsize; i++) {
                    tree = (*treeList)[i]->findDown(rootToFind.n, rootToFind.x, rootToFind.y, rootToFind.z);
                    if ((tree->x() == rootToFind.x) && (tree->y() == rootToFind.y) &&
                            (tree->z() == rootToFind.z) && (tree->n() == rootToFind.n)) {
                        break;
                    }
                }
                if (tree) {
                    if (tree->getLocalSubtreeCost() <= spaceleft) {
                        if (debug) {
                            std::cout << "partitionWorker: found tree (" << tree->x() << "," <<
                            tree->y() << "," << tree->z() << ")" << std::endl;
                        }
                        arsend & 1 & rootToFind & tree->getLocalSubtreeCost();
                        tree->setSendto(rootToFind.future_owner);
                        std::vector<RootList> remoteList;
                        remoteList = findRemoteChildren(tree, rootToFind.future_owner);
                        int rsize = remoteList.size();
                        arsend & rsize;
                        if (debug) {
                            std::cout << "partitionWorker: found " << rsize << " remote children"
                            << std::endl;
                        }
                        for (int i = 0; i < rsize; i++) {
                            arsend & remoteList[i];
                        }
                    } else {
                        if (debug) {
                            std::cout << "partitionWorker: root (" << rootToFind.x << "," <<
                            rootToFind.y << "," << rootToFind.z << ") too expensive;"
                            << " try its children" << std::endl;
                        }
                        std::vector<RootList> remoteList;
                        for (int i = 0; i < 2; i++) {
                            if (spaceleft == 0) break;
                            for (int j = 0; j < 2; j++) {
                                if (spaceleft == 0) break;
                                for (int k = 0; k < 2; k++) {
                                    if (spaceleft == 0) break;
                                    OctTree<T> *child = tree->child(i,j,k);
                                    if (child) {
                                        partitionWorkerSub(child, &localList, &spaceleft, &costList,
                                                           &remoteList, rootToFind.future_owner);
                                    }
                                }
                            }
                        }

                        int size = localList.size();
                        if (debug) {
                            std::cout << "partitionWorker: have list of size " << size << " to send"
                            << std::endl;
                        }
                        arsend & size;
                        for (int i = 0; i < size; i++) {
                            arsend & localList[i] & costList[i];
                        }
                        int rsize = remoteList.size();
                        arsend & rsize;
                        if (debug) {
                            std::cout << "partitionWorker: Sub found " << rsize << " remote children"
                            << std::endl;
                        }
                        for (int i = 0; i < rsize; i++) {
                            arsend & remoteList[i];
                        }
                    }
                }
                break;
            case EXIT_CODE:
                // terminate
                keepGoing = false;
                if (debug) {
                    std::cout << "partitionWorker: received EXIT_CODE " << signal << std::endl;
                }
                break;
            default:
                std::cout << "partitionWorker: unknown signal " << signal << std::endl;
                break;
            }
        }
    }


    template <class T>
    void partitionWorkerSub(OctTree<T> *root, std::vector<RootList> *localList, Cost *spaceleft,
                            std::vector<Cost> *costList, std::vector<RootList> *remoteList, int p) {
        Communicator comm;
        int me = comm.rank();
//	bool debug = true;
        bool debug = false;

        if (debug) {
            std::cout << "partitionWorkerSub: at beginning of function, spaceleft = " <<
            *spaceleft << std::endl;
        }

        if (*spaceleft == 0) {
            return;
        }
        Cost rootcost = root->getLocalSubtreeCost();
        if (rootcost <= *spaceleft) {
            if (debug) {
                std::cout << "partitionWorkerSub: this subtree fits into partition" << std::endl;
            }
            localList->push_back(RootList(root, me, me));

            *remoteList = findRemoteChildren(root, p);
            root->setSendto(p);
            costList->push_back(rootcost);
            *spaceleft -= rootcost;
        } else {
            if (debug) {
                std::cout << "partitionWorkerSub: this subtree doesn't fit into partition "
                << "(cost = " << rootcost << "); try the children" << std::endl;
            }
            FOREACH_CHILD(OctTree<T>, root,
                          if (_debug) {
                          std::cout << "partitionWorkerSub: at beginning of FOREACH_CHILD"
                          << std::endl;
                      }
                      partitionWorkerSub(child, localList, spaceleft, costList, remoteList, p);
                      if (*spaceleft == 0) {
                          if (_debug) {
                                  std::cout << "partitionWorkerSub: no space left; we're outta here"
                                  << std::endl;
                              }
                              return;
                          }
                         );
        }
        if (debug) {
            std::cout << "partitionWorkerSub: at end of function" << std::endl;
        }
    }

    /********************************************************************************/
    /* computeGlobalCost: processor 0 is manager process that determines cost of 	*/
    /* tree globally.  Uses producer model for workers, where workers wait for 	*/
    /* message from manager about what to do and respond accordingly.  In the end,  */
    /* this returns the cost of the tree represented by RootList root on manager.	*/
    /********************************************************************************/

    template <class T>
    Cost computeGlobalCost(RootList root, OctTree<T> *tree, std::vector<OctTree<T>* > const treeList) {
        Communicator comm;
        int me = comm.rank(), np = comm.nproc();
        Cost total = 0;
//	bool debug = true;
        bool debug = false;
        if (debug)
            std::cout << "computeGlobalCost: beginning of function" << std::endl;
        if (me == 0) {
            for (int i = 1; i < np; i++) {
                archive::MPIOutputArchive arout(comm, i);
                if (debug)
                    std::cout << "computeGlobalCost: sending GLOBAL_COST code to " << i << std::endl;
                arout & GLOBAL_COST;
                if (debug)
                    std::cout << "computeGlobalCost: sent exit code to " << i << std::endl;
            }
            if (debug)
                std::cout << "computeGlobalCost: about to call globalCostManager" << std::endl;
            total = globalCostManager(root, treeList);
            // then, shut down infinite loop of worker
            if (debug) {
                std::cout << "computeGlobalCost: total = " << total << std::endl;
                std::cout << "computeGlobalCost: sending exit code" << std::endl;
            }
            for (int i = 1; i < np; i++) {
                archive::MPIOutputArchive arout(comm, i);
                if (debug)
                    std::cout << "computeGlobalCost: sending exit code to " << i << std::endl;
                arout & EXIT_CODE;
                if (debug)
                    std::cout << "computeGlobalCost: sent exit code to " << i << std::endl;
            }
        } else {
            archive::MPIInputArchive arin(comm, 0);
            arin & GLOBAL_COST;
            total = globalCostWorker(root, treeList);
        }
        return total;
    }

    template <class T>
    Cost globalCostManager(RootList root, std::vector<OctTree<T>* > const treeList) {
        Communicator comm;
        int me = comm.rank();
        int rank, nRemoteKids, i;
        RootList remoteRoot;
        Cost total = 0, subtotal = 0;
        OctTree<T> *tree = new OctTree<T>();

//	bool debug = true;
        bool debug = false;

        if (debug)
            std::cout << "globalCostManager: at beginning of function" << std::endl;

        if (root.current_owner == me) {
            if (debug) {
                std::cout << "globalCostManager: I own this root" << std::endl;
            }
            int npieces = treeList.size();
            if (debug) {
                std::cout << "globalCostManager: treeList.size() = " << npieces << std::endl;
            }
            for (i = 0; i < npieces; i++) {
                if (debug) {
                    std::cout << "globalCostManager: looking for tree in treeList" << std::endl;
                    std::cout << "globalCostManager: treeList[" << i << "], n = " <<  treeList[i]->n() <<
                    ", (" << treeList[i]->x() << "," << treeList[i]->y() << "," << treeList[i]->z()
                    << ")" << std::endl;
                }
                if ((treeList[i]->n() == root.n) && (treeList[i]->x() == root.x) &&
                        (treeList[i]->y() == root.y) && (treeList[i]->z() == root.z)) {
                    tree = treeList[i];
                } else {
                    tree = treeList[i]->findDown(root.n, root.x, root.y, root.z);
                }
                if (debug) {
                    std::cout << "globalCostManager: done looking for tree in treeList, tree = "
                    << tree << std::endl;
                }
                if ((tree) && (tree->n() == root.n) && (tree->x() == root.x) &&
                        (tree->y() == root.y) && (tree->z() == root.z)) {
                    if (debug) {
                        std::cout << "globalCostManager: found the right tree!" << std::endl;
                    }
                    break;
                }
            }

            if (tree->getSendto() != -1) {
                if (debug) {
                    std::cout << "globalCostManager: this tree is already accounted for " <<
                    "(" << tree->x() << "," << tree->y() << "," << tree->z() << ")" << std::endl;
                }
                return 0;
            } else if (tree->isParent()) {
                FOREACH_CHILD(OctTree<T>, tree,
                              if (child->isremote())
                              rank = child->rank();
                              else
                                  rank = me;
                                  if (debug) {
                                      std::cout << "globalCostManager: sending child " <<
                                      "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                                          << "to globalCostManager function" << std::endl;
                                      }
                              total += globalCostManager(RootList(child->x(), child->y(), child->z(),
                                                                  child->n(), rank, rank), treeList);
                             );
            }
            subtotal = tree->getCost();
            if (debug) {
                std::cout << "globalCostManager: subtotal from local tree = "
                << subtotal << std::endl;
            }
//	    tree->setLocalSubtreeCost(subtotal);
            total += subtotal;
            subtotal = 0;
        } else {
            if (debug) {
                std::cout << "globalCostManager: I don't own this root; processor " <<
                root.current_owner << " does" << std::endl;
            }
            archive::MPIOutputArchive arout(comm, root.current_owner);
            if (debug) {
                std::cout << "globalCostManager: about to send COMPUTE_LOCAL_COST " <<
                COMPUTE_LOCAL_COST << " and root " << root.x << " " << root.y <<
                " " << root.z << " " << root.n << " to processor " <<
                root.current_owner << std::endl;
            }
//	    arout & GLOBAL_COST & COMPUTE_LOCAL_COST & root;
            arout & COMPUTE_LOCAL_COST & root;
            if (debug) {
                std::cout << "globalCostManager: sent root to processor " << root.current_owner <<
                std::endl;
            }
            archive::MPIInputArchive arin(comm, root.current_owner);
            arin & subtotal & nRemoteKids;
            if (debug) {
                std::cout << "globalCostManager: received subtotal " << subtotal <<
                " and nRemoteKids " << nRemoteKids << std::endl;
            }
            total += subtotal;
            for (i = 0; i < nRemoteKids; i++) {
                arin & remoteRoot;
                total += globalCostManager(remoteRoot, treeList);
            }
        }
        tree->setLocalSubtreeCost(total);
        return total;
    }

    template <class T>
    Cost globalCostWorker(RootList root, std::vector<OctTree<T>* > const treeList) {
        Communicator comm;
        int me = comm.rank(), npieces = treeList.size();
        archive::MPIInputArchive arrecv(comm, 0);
        archive::MPIOutputArchive arsend(comm, 0);
//	bool debug = true;
        bool debug = false;

        if (debug)
            std::cout << "globalCostWorker: at beginning of function" << std::endl;

        while (1) // infinite loop
        {
            int msg, i;
            OctTree<T> *t = new OctTree<T>();
            arrecv & msg;
            if (msg == EXIT_CODE) {
                if (debug)
                    std::cout << "Processor " << me << ": received exit code " << msg << std::endl;
                break;
            } else if (msg == COMPUTE_LOCAL_COST) {
                Cost cost = 0;
//		std::vector<RootList> *rootList = new std::vector<RootList>();
                std::vector<RootList> rootList;
                if (debug) {
                    std::cout << "Processor " << me << ": received compute local cost code " << msg <<
                    std::endl;
                }
                arrecv & root;
                if (debug)
                    std::cout << "Processor " << me << ": received root" << std::endl;
                for (i = 0; i < npieces; i++) {
                    t = treeList[i]->findDown(root.n, root.x, root.y, root.z);
                    if ((t) && (t->n() == root.n) && (t->x() == root.x) && (t->y() == root.y) &&
                            (t->z() == root.z)) {
                        break;
                    }
                }
                if (debug) {
                    std::cout << "Processor " << me << ": found tree: " <<
                    "(" << t->x() << "," << t->y() << "," << t->z() << "," << t->n() << ")"
                    << std::endl;
                }
                cost += t->computeLocalCost(&rootList);
//		Cost tcost = t->computeLocalCost(&rootList);
//		cost += tcost;
                t->setLocalSubtreeCost(cost);
                int rlsize = rootList.size();
                if (debug)
                    std::cout << "about to send cost " << cost << " and rlsize " << rlsize << std::endl;
                arsend & cost & rlsize;
                for (i = 0; i < rlsize; i++) {
                    arsend & rootList[i];
                }
            } else {
                std::cout << "Processor " << me << ": received unknown code " << msg << std::endl;
            }
        }
        return 0;
    }


    template <class T>
    std::vector<RootList> findRemoteChildren(OctTree<T> *t, int p) {
        std::vector<RootList> remoteList;
        std::vector<RootList> subList;

        FOREACH_CHILD(OctTree<T>, t,
                      if (child->isremote()) {
                      remoteList.push_back(RootList(child, child->rank(), p));
                      }
                      else if (child->isParent()) {
                      subList = findRemoteChildren(child, p);
                          int size = subList.size();
                          for (int i = 0; i < size; i++) {
                              remoteList.push_back(subList[i]);
                          }
                          subList.clear();
                      }
                     );
        return remoteList;
    }


    template <class T>
    void exchangeTrees(std::vector<RootList> *globalList, std::vector<SharedPtr<OctTree<T> > > *treeList,
                       bool glue) {
        Communicator comm;
        int me = comm.rank();
        int glength = 0, tlength = 0, i, j;
        RootList root;

//	bool debug = true;
        bool debug = false;


        /* globalList: the list of trees that need to be exchanged */
        /* treeList: the list of subtrees that this processor owns */

        if (globalList)
	{
            glength = globalList->size();
	    sort(globalList->begin(), globalList->end());
	}
        if (treeList)
            tlength = treeList->size();

        /* First, we send out the list of trees that need to be exchanged */

        if (debug) {
            std::cout << "exchangeTrees: about to bcast glength = " << glength << std::endl;
        }
        comm.Bcast(&glength, 1, 0);
        if (debug) {
            std::cout << "exchangeTrees: after bcast glength = " << glength << std::endl;
            if (me == 0) {
                for (i = 0; i < glength; i++) {
                    std::cout << "    n = " << (*globalList)[i].n << ", (" << (*globalList)[i].x <<
                    "," << (*globalList)[i].y << "," << (*globalList)[i].z << ")" << std::endl;
                }
            }
        }
        for (i = 0; i < glength; i++) {
            if (me == 0) {
                root = (*globalList)[i];
                if (debug) {
                    std::cout << "exchangeTrees: root " << i << ": n = " << root.n << ", " <<
                    "(" << root.x << "," << root.y << "," << root.z << ")" << "current = "
                    << root.current_owner << ", future = " << root.future_owner << std::endl;
                    std::cout << "exchangeTrees: now root will be broadcast: " << std::endl;
                }
            } else {}
            comm.Bcast(&root, 1, 0);
            if (debug) {
//                std::cout << "exchangeTrees: root " << i << ": n = " << root.n << ", " <<
//                "(" << root.x << "," << root.y << "," << root.z << ")" << "current = "
//                << root.current_owner << ", future = " << root.future_owner << std::endl;
                std::cout << "exchangeTrees: after bcast root number " << i << std::endl;
            }
            if (me != 0) {
                globalList->push_back(root);
            }

            /* if this root is somehow connected with me, print this out for debug */
            if (debug) {
                if ((root.current_owner == me) || (root.future_owner == me)) {
                    std::cout << "exchangeTrees: I have something to do with tree " <<
                    "(" << root.x << "," << root.y << "," << root.z << "," << root.n << ")"
                    << std::endl;
                }
            }
        }
        MPI::COMM_WORLD.Barrier();
        /* Now, do something with each tree on the globalList */
        for (i = 0; i < glength; i++) {
            if (debug)
                std::cout << std::endl;
            /* if I own this root and will continue to own this root */
            if (((*globalList)[i].current_owner == me) &&
                    ((*globalList)[i].current_owner == (*globalList)[i].future_owner)) {
                // do nothing, unless localList[i] is not currently a local root
                if (debug) {
                    std::cout << "exchangeTrees: current_owner = future_owner" << std::endl;
                }
                OctTree<T> *t = new OctTree<T>();
                for (j = 0; j < tlength; j++) {
                    if (debug) {
                        std::cout << "exchangeTrees: looking for tree in treeList[" << j << "]"
                        << std::endl;
//                        (*treeList)[j]->depthFirstTraverseAll();
                        std::cout << "exchangeTrees: will I find it?" << std::endl;
                    }
                    t = (*treeList)[j]->findDown((*globalList)[i].n, (*globalList)[i].x, (*globalList)[i].y,
                                             (*globalList)[i].z);
                    if ((t) && (t->n() == (*globalList)[i].n) && (t->x() == (*globalList)[i].x) &&
                            (t->y() == (*globalList)[i].y) && (t->z() == (*globalList)[i].z)) {
                        if (debug) {
                            std::cout << "exchangeTrees: found tree n = " << t->n() << ", (" <<
                            t->x() << "," << t->y() << "," << t->z() << ")" << std::endl;
                        }
                        break;
                    } else
                        t = 0;
                }
                /* if the tree is not already a root */
                if (t != (*treeList)[j].get()) {
                    OctTree<T> *p = new OctTree<T>();
                    p = t->parent();
                    ProcessID future_owner = -1;
                    for (int k = i+1; k < glength; k++) {
                        if ((*globalList)[i].isDescendant((*globalList)[k])) {
                            future_owner = (*globalList)[k].future_owner;
                            break;
                        }
                    }
                    OctTree<T> *pprime = new OctTree<T>(p->n(), p->x(), p->y(), p->z(),true, 0,
                                                        future_owner, 0);

                    Translation x = (*globalList)[i].x, y = (*globalList)[i].y, z = (*globalList)[i].z;
                    x -= 2*(x/2);
                    y -= 2*(y/2);
                    z -= 2*(z/2);
                    SharedPtr<OctTree<T> > tptr = p->child(x,y,z);
                    treeList->push_back(tptr);
                    if (debug) {
                        std::cout << "exchangeTrees: pushed back tree n = " << t->n() << " " <<
                        (treeList->back())->n() << ",(" << t->x() << "," << t->y() << "," <<
                        t->z() << "), (" << (treeList->back())->x() << "," <<
                        (treeList->back())->y() << "," << (treeList->back())->z() << ")"
                        << std::endl;
                    }
                    tlength++;
                    if (debug) {
                        std::cout << "exchangeTrees: about to insert remote child " <<
                        "(" << x << "," << y << "," << z << ")" << " into parent " <<
                        "(" << p->x() << "," << p->y() << "," << p->z() << ")" << std::endl;
                    }
                    pprime->setChild(x,y,z, tptr);
                    if (debug) {
                        std::cout << "exchangeTrees: just set pprime's child" << std::endl;
                    }
		// Make sure this is right 7-03-06
                    t->setParent(pprime);
                    if (debug) {
                        std::cout << "exchangeTrees: just set t's parent" << std::endl;
                    }
/*
		    if (p->islocal())
		    {
*/
                        SharedPtr<OctTree<T> > q = SharedPtr<OctTree<T> >();
                        q = p->insert_remote_child(x, y, z, (*globalList)[i].future_owner);
                        if (debug) {
                            std::cout << "exchangeTrees: inserted remote child " <<
                            "(" << q->x() << "," << q->y() << "," << q->z() << ")" << " into parent " <<
                            "(" << p->x() << "," << p->y() << "," << p->z() << ")" << std::endl;
                            std::cout << "exchangeTrees: after inserting remote child, is t still OK?" <<
                            std::endl << "        " << "n = " << t->n() << " " <<
                            (treeList->back())->n() << ", (" << t->x() << "," << t->y() << "," <<
                            t->z() << "), (" << (treeList->back())->x() << "," <<
                            (treeList->back())->y() << "," << (treeList->back())->z() << ")"
                            << std::endl;
                        }
/*
		    }
		    else
		    {
			std::cout << "I bet this causes a segfault" << std::endl;
			p->setChild(x,y,z, SharedPtr<OctTree<T> >(0));
		    }
*/
                }
                if (debug) {
                    std::cout << "exchangeTrees: end of current = future: tlength = " << tlength
                    << std::endl;
                }
            }
            /* if I own this root but need to send it elsewhere */
            else if ((*globalList)[i].current_owner == me) {
                // send it to its new owner
                // we assume that the element is indeed owned by this processor
                if (debug) {
                    std::cout << "exchangeTrees: I am current_owner" << std::endl;
                }
                OctTree<T> *t = new OctTree<T>();
                for (j = 0; j < tlength; j++) {
                    if (debug) {
                        std::cout << "exchangeTrees: looking for tree in treeList[" << j << "]"
                        << std::endl;
//                        (*treeList)[j]->depthFirstTraverseAll();
                        std::cout << "exchangeTrees: will I find it?" << std::endl;
                    }
                    t = (*treeList)[j]->findDown((*globalList)[i].n, (*globalList)[i].x, (*globalList)[i].y,
                                             (*globalList)[i].z);
                    if (debug) {
                        if (t) {
                            std::cout << "exchangeTrees: t = (" << t->x() << "," << t->y() << "," <<
                            t->z() << ")" << std::endl;
                        } else {
                            std::cout << "exchangeTrees: t not yet found" << std::endl;
                        }
                    }
                    if ((t) && (t->n() == (*globalList)[i].n) && (t->x() == (*globalList)[i].x) &&
                            (t->y() == (*globalList)[i].y) && (t->z() == (*globalList)[i].z)) {
                        if (debug) {
                            std::cout << "exchangeTrees: found tree in treeList[" << j << "]:"
                            " (" << t->x() << "," << t->y() << "," << t->z() << ")"<< std::endl;
                        }
                        break;
                    } else {
                        t = 0;
                        if (debug) {
                            std::cout << "exchangeTrees: did not find tree in treeList[" << j << "]"
                            << std::endl;
                        }
                    }
                }
                if (debug) {
                    std::cout << "exchangeTrees: I am about to send the tree to " <<
                    (*globalList)[i].future_owner << std::endl;
                    std::cout << "exchangeTrees: t = : " << (t) << std::endl;
                }
                sendSubtree(t, me, (*globalList)[i].future_owner);
                if (debug) {
                    std::cout << "exchangeTrees: sent tree t" << std::endl;
                }
                OctTree<T> *p = new OctTree<T>();
                p = t->parent();
                if (debug) {
                    if (p) {
                        std::cout << "exchangeTrees: the parent exists and is (" << p->x() <<
                        "," << p->y() << "," << p->z() << ")" << std::endl;
                    } else {
                        std::cout << "exchangeTrees: the parent does not exist " << std::endl;
                    }
                }
                /* if the parent of the tree just sent is local, then insert a remote child in its place */
                if ((p) && (p->islocal())) {
                    if (debug) {
                        std::cout << "exchangeTrees: the parent of tree just sent is local" << std::endl;
                    }
                    Translation x = (*globalList)[i].x, y = (*globalList)[i].y, z = (*globalList)[i].z;
                    x -= 2*(x/2);
                    y -= 2*(y/2);
                    z -= 2*(z/2);
                    p->insert_remote_child(x, y, z, (*globalList)[i].future_owner);
                }
                /* otherwise, remove the tree from my treeList */
                else if (p) {
                    Translation x = (*globalList)[i].x, y = (*globalList)[i].y, z = (*globalList)[i].z;
                    x -= 2*(x/2);
                    y -= 2*(y/2);
                    z -= 2*(z/2);
                    p->insert_remote_child(x, y, z, (*globalList)[i].future_owner);
/*
		    std::cout << "I bet this causes a segfault too" << std::endl;
		    if (p)
                        p->setChild(x, y, z, SharedPtr<OctTree<T> >(0));
		    std::cout << "no, made it past" << std::endl;
*/
		    bool deleteit = true;
		    FOREACH_CHILD(SharedPtr<OctTree<T> >, p,
			if (child->islocal())
			    deleteit = false;
		    );
		    if (deleteit)
		    {
                    	if (debug) {
                            std::cout << "exchangeTrees: the parent of tree just sent is remote, "
                            << "so I'll delete this tree" << std::endl;
			}
                        treeList->erase((treeList->begin()+j));
                        tlength--;
                        if (debug) {
                            std::cout << "exchangeTrees: also need to delete the tree itself" << std::endl;
                        }
//                        delete t;
                    	if (debug) {
                            std::cout << "exchangeTrees: deleted the tree itself" << std::endl;
                    	}
                    }
		    else
		    {
			if (debug)
			{
			    std::cout << "exchangeTrees: this tree has other local children, " <<
				"don't delete it!" << std::endl;
			}
		    }
		}
            } else if ((*globalList)[i].future_owner == me) {
                // receive from current owner
                if (debug) {
                    std::cout << "exchangeTrees: I am future_owner" << std::endl;
                }
                OctTree<T> *t = new OctTree<T>();
                recvSubtree(t, me, (*globalList)[i].current_owner);
                if (debug) {
                    std::cout << "exchangeTrees: at least I received something!" << std::endl;
                }
                SharedPtr<OctTree<T> > tptr = SharedPtr<OctTree<T> >(t);
                int parent_loc = -1;
                for (int j = i+1; j < glength; j++) {
                    if ((*globalList)[i].isDescendant((*globalList)[j])) {
                        parent_loc = (*globalList)[j].future_owner;
                        break;
                    }
                }
                Translation x = t->x(), y = t->y(), z = t->z();
                OctTree<T> *p = new OctTree<T>(t->n()-1, x/2, y/2, z/2, true, 0, parent_loc, 0);
                x -= 2*(x/2);
                y -= 2*(y/2);
                z -= 2*(z/2);
                p->setChild(x,y,z, tptr);
                treeList->push_back(tptr);
                tlength++;
            }

            if (debug) {
                int tlength = treeList->size();
                std::cout << "exchangeTrees: at end of loop, list of length " <<
                tlength << ":" << std::endl;
                for (int i = 0; i < tlength; i++) {
                    std::cout << "tree " << i << " of " << tlength << ":" << std::endl;
//                    (*treeList)[i]->depthFirstTraverseAll();
                }
            }
        }
        if (debug) {
            int tlength = treeList->size();
            std::cout << "exchangeTrees: before sorting and gluing together, list of length " <<
            tlength << ":" << std::endl;
            for (int i = 0; i < tlength; i++) {
                std::cout << "tree " << i << " of " << tlength << ":" << std::endl;
//                (*treeList)[i]->depthFirstTraverseAll();
            }
        }

//	sort((*treeList).begin(), (*treeList).end(), less<OctTree<T>* > ());
	tlength = treeList->size();
//	sort((*treeList).begin(), (*treeList).end());
	sort((*treeList).begin(), (*treeList).end(), lessPtr<OctTree<T> > ());
	if (debug)
	{
	    std::cout << "exchangeTrees: before gluing, number of trees = " << tlength << std::endl;
	}
        if (glue)
            glueTrees(treeList);
	tlength = treeList->size();
	if (debug)
	{
	    std::cout << "exchangeTrees: about to return with " << tlength << " trees" << std::endl;
	}
    }

    template <class T>
    void sendSubtree(OctTree<T> *tree, ProcessID me, ProcessID dest) {
        Communicator comm;
        archive::MPIOutputArchive ardest(comm, dest);

//	bool debug = true;
        bool debug = false;

        if (debug) {
            std::cout << "sendSubtree: about to send tree to " << dest << std::endl;
            std::cout << "sendSubtree: subtree looks like this:" << std::endl;
            tree->depthFirstTraverseAll();
            std::cout << "sendSubtree: done with depthFirstTraverseAll" << std::endl;
        }
        ardest & *tree;

        if (debug)
            std::cout << "sendSubtree: all done" << std::endl;
    }

    template <class T>
    void sendSubtree(SharedPtr<OctTree<T> > tree, ProcessID me, ProcessID dest) {
        Communicator comm;
        archive::MPIOutputArchive ardest(comm, dest);

//	bool debug = true;
        bool debug = false;

        if (debug) {
	    std::cout << "sendSubtree: SharedPtr version" << std::endl;
            std::cout << "sendSubtree: about to send tree to " << dest << std::endl;
            std::cout << "sendSubtree: subtree looks like this:" << std::endl;
            tree->depthFirstTraverseAll();
            std::cout << "sendSubtree: done with depthFirstTraverseAll" << std::endl;
        }
        ardest & *(tree.get());

        if (debug)
            std::cout << "sendSubtree: all done" << std::endl;
    }


    template <class T>
    void recvSubtree(OctTree<T> *t, ProcessID me, ProcessID source) {
        Communicator comm;
        archive::MPIInputArchive arsource(comm, source);

// 	bool debug = true;
        bool debug = false;

        if (debug) {
            std::cout << "recvSubtree: waiting to receive subtree from " << source << std::endl;
        }
        arsource & *t;
//        arsource & t;
        if (debug) {
            int s = t->tallyNodes();
            std::cout << "recvSubtree: received subtree of size "
            << s << std::endl;
        }

        if (debug)
            std::cout << "recvSubtree: all done!" << std::endl;
    }

    template <class T>
    void recvSubtree(SharedPtr<OctTree<T> > p, ProcessID me, ProcessID source) {
        Communicator comm;
        archive::MPIInputArchive arsource(comm, source);
	OctTree<T> *t = new OctTree<T>();

// 	bool debug = true;
        bool debug = false;

        if (debug) {
	    std::cout << "recvSubtree: SharedPtr version" << std::endl;
            std::cout << "recvSubtree: waiting to receive subtree from " << source << std::endl;
        }
        arsource & *t;
        if (debug) {
            int s = t->tallyNodes();
            std::cout << "recvSubtree: received subtree of size "
            << s << std::endl;
        }

	p = SharedPtr<OctTree<T> > (t);

        if (debug)
            std::cout << "recvSubtree: all done!" << std::endl;
    }


    template <class T>
    void sendMsg(T msg, ProcessID me, ProcessID dest) {
        Communicator comm;
//        madness::redirectio(comm);
        archive::MPIOutputArchive ardest(comm, dest);

        ardest & msg;
    }

    template <class T>
    void recvMsg(T *msg, ProcessID me, ProcessID source) {
        Communicator comm;
//        madness::redirectio(comm);
        archive::MPIInputArchive arsource(comm, source);

        arsource & *msg;
    }

    /********************************************************************************/
    /* glueTrees: takes list of trees, finds any tree in list that is a child in	*/
    /* another of the trees, and sticks it in there.  E.g. if we had tree (2,6,6)   */
    /* and tree (1,3,3), we would glue the first tree into the second tree because  */
    /* (2,6,6) is the child of (1,3,3).						*/
    /********************************************************************************/

    template <class T>
    void glueTrees(std::vector<SharedPtr<OctTree<T> > > *treeList) {
        int N = treeList->size();

//        bool debug = true;
	bool debug = false;

	if (debug)
	{
	    std::cout << "glueTrees: beginning of function, N = " << N << std::endl;
	}

//	sort((*treeList).begin(), (*treeList).end());
	sort((*treeList).begin(), (*treeList).end(), lessPtr<OctTree<T> > ());

	if (debug)
	{
	    std::cout << "glueTrees: sorted treeList" << std::endl;
	}

	for (int i = N-1; i >= 0; i--)
	{
	    bool removeit = true;
	    SharedPtr<OctTree<T> > t = SharedPtr<OctTree<T> >();
	    t = (*treeList)[i];
	    if (t->islocal())
	    {
		removeit = false;
	    }
	    else
	    {
	        FOREACH_CHILD(SharedPtr<OctTree<T> >, t,
		    if (child->islocal())
		    {
		        removeit = false;
		    }
	        );
	    }
	    if (removeit)
	    {
//		if (debug)
		    std::cout << "glueTrees: removing tree " << i << " that is purely remote" << std::endl;
                treeList->erase((treeList->begin()+i));
		N = treeList->size();
	    }
	}

	if (debug)
	{
	    std::cout << "glueTrees: list of " << N << " local subtrees" << std::endl;
	    for (int q = 0; q < N; q++)
	    {	
//		std::cout << "    n = " << (*treeList)[q]->n() << ", (" << (*treeList)[q]->x() <<
//			"," << (*treeList)[q]->y() << "," << (*treeList)[q]->z() << ")" << std::endl;
		std::cout << "Subtree " << q << " of " << N << ":" << std::endl;
		(*treeList)[q]->depthFirstTraverseAll();
	    }
	}

	for (int i = 0; i < N-1; i++)
	{
	    if (debug)
	    {
		std::cout << "glueTrees: at top of i loop, i = " << i << ", and N = " << N <<
			std::endl;
	    	std::cout << "glueTrees: list of local subtrees" << std::endl;
	    	for (int q = 0; q < N; q++)
	    	{	
		    std::cout << "    n = " << (*treeList)[q]->n() << ", (" << (*treeList)[q]->x() <<
			"," << (*treeList)[q]->y() << "," << (*treeList)[q]->z() << ")" << std::endl;
	    	}
	    }
	    int j = i, tmpi = -1, tmpj = -1;
	    do
	    {
		j++;
		if (debug)
		{
		    std::cout << "glueTrees: at top of j loop (after j++), j = " << j << " and N = "
			<< N << std::endl;
/*
	    	    for (int q = 0; q < N; q++)
	    	    {	
		    	std::cout << "    n = " << (*treeList)[q]->n() << ", (" << (*treeList)[q]->x() 
				<< "," << (*treeList)[q]->y() << "," << (*treeList)[q]->z() << ")" << 
				std::endl;
	    	    }
*/
		}
		if ((*treeList)[i]->equals((*treeList)[j].get()))
		{
		    if (debug)
		    {
			std::cout << "glueTrees: treeList[" << i << "] == treeList[" << j << "]"
				<< std::endl;
			std::cout << "    [n = " << (*treeList)[i]->n() << ", (" << (*treeList)[i]->x()
				<< "," << (*treeList)[i]->y() << "," << (*treeList)[i]->z() << ") ]"
				<< std::endl;
		    }
		    // merge these trees together, deleting one of them
		    if ((*treeList)[i]->islocal()) // just on the off chance this is The Root
		    {
			tmpi = i; 
			tmpj = j;
			if (debug)
			{
			    std::cout << "glueTrees: this node, " << i << ", is The Root" << std::endl;
			}
		    }
		    else
		    {
			tmpi = j;
			tmpj = i;
		    }
		    FORIJK(
			SharedPtr<OctTree<T> > child = SharedPtr<OctTree<T> >();
			child = (*treeList)[tmpj]->child(i,j,k);
			if (child.get())
			{
			    if (debug)
			    {
				std::cout << "glueTrees: child(" << i << "," << j << "," << k << ") "
					<< "exists" << std::endl;
			    }
			    // if treeList[tmpi] doesn't already have a local node here
			    if (!((*treeList)[tmpi]->child(i,j,k).get()) || 
				((*treeList)[tmpi]->child(i,j,k)->isremote()))
			    {
				if (debug)
				{
				    std::cout << "glueTrees: child(" << i << "," << j << "," << k << ")"
					<< " either does not exist or is not local, so insert new node "
					<< "here" << std::endl;
				}
				if (child->islocal())
				    (*treeList)[tmpi]->insert_local_child(child);
			    }
			}
		    );
		    // remove treeList[tmpj]
		    if (tmpj > -1)
		    {
			if (debug)
			{
			    std::cout << "glueTrees: erase treeList[" << tmpj << "]" << std::endl;
			}
                	treeList->erase((treeList->begin()+tmpj));
		    }
		    N--;
		    j--;
		}
		else
		{
		    if (debug)
		    {
			std::cout << "glueTrees: treeList[" << i << "] != treeList[" << j << "]"
				<< std::endl;
			std::cout << "i:    [n = " << (*treeList)[i]->n() << ", (" << (*treeList)[i]->x()
				<< "," << (*treeList)[i]->y() << "," << (*treeList)[i]->z() << ") ]"
				<< std::endl;
			std::cout << "j:    [n = " << (*treeList)[j]->n() << ", (" << (*treeList)[j]->x()
				<< "," << (*treeList)[j]->y() << "," << (*treeList)[j]->z() << ") ]"
				<< std::endl;
		    }
		}

		if ((*treeList)[i]->isDescendant((*treeList)[j].get()))
		{
		    if (debug)
		    {
			std::cout << "glueTrees: treeList[" << i << "] is descendant of treeList[" 
				<< j << "]" << std::endl;
			std::cout << "i:    [n = " << (*treeList)[i]->n() << ", (" << (*treeList)[i]->x()
				<< "," << (*treeList)[i]->y() << "," << (*treeList)[i]->z() << ") ]"
				<< std::endl;
			std::cout << "j:    [n = " << (*treeList)[j]->n() << ", (" << (*treeList)[j]->x()
				<< "," << (*treeList)[j]->y() << "," << (*treeList)[j]->z() << ") ]"
				<< std::endl;
		    }
		    // if we can find treeList[i] on treeList[j] and it's local, then
		    // attach children of treeList[i] to the appropriate node of treeList[j]
		    OctTree<T> *t = new OctTree<T>();
		    t = (*treeList)[j]->findDown((*treeList)[i]->n(), (*treeList)[i]->x(),
			(*treeList)[i]->y(), (*treeList)[i]->z());
		    if (t->equals(((*treeList)[i]).get()))
		    {
		    	if (debug)
		    	{
			    std::cout << "glueTrees: t == treeList[" << i << "]" << std::endl;
			    std::cout << "t:    [n = " << t->n() << ", (" << t->x() << "," << 
			 	t->y() << "," << t->z() << "), rank = " << t->rank() << ", "
				<< "local? " << t->islocal() << " ]" << std::endl;
		    	}
			if (t->islocal())
			{
			    if (debug)
			    {
				std::cout << "glueTrees: t is local" << std::endl;
			    }
			    tmpi = i;
			    FORIJK(
				SharedPtr<OctTree<T> > child = SharedPtr<OctTree<T> >();
				child = (*treeList)[tmpi]->child(i,j,k);
				if (child.get())
				{
			    	    if (debug)
			    	    {
					std::cout << "glueTrees: child(" << i << "," << j << "," << k << ") "
						<< "exists" << std::endl;
			    	    }
				    if (!(t->child(i,j,k).get()) || (t->child(i,j,k)->isremote()))
				    {
					if (debug)
					{
				    	    std::cout << "glueTrees: child(" << i << "," << j << "," << k 
						<< ")" << " either does not exist or is not local," << 
						" so insert new node here" << std::endl;
					}
					t->insert_local_child(child);
				    }
				}
			    );
			    // remove treeList[i]
			    if (debug)
			    {
			    	std::cout << "glueTrees: erase treeList[" << tmpi << "]" << std::endl;
			    }
                	    treeList->erase((treeList->begin()+tmpi));
			    N--;
			    j = i;
			}
			else if ((*treeList)[i]->islocal())
			{
			    if (debug)
			    {
				std::cout << "glueTrees: t is not local" << std::endl;
			    }
			    tmpi = i;
			    OctTree<T> *p = new OctTree<T> ();
			    p = t->parent();
			    if (p)
			    {
			    	if (debug)
			    	{
				    std::cout << "glueTrees: parent exists: n = " << p->n() <<
					", (" << p->x() << "," << p->y() << "," << p->z() << ")" 
					<< std::endl;
			    	}
				p->insert_local_child((*treeList)[i]);
				std::cout << "glueTrees: inserted local child: n = " <<
					(*treeList)[i]->n() << ", (" << (*treeList)[i]->x() << ","
					<< (*treeList)[i]->y() << "," << (*treeList)[i]->z() << ")"
					<< std::endl;
			    	if (debug)
			    	{
				    std::cout << "glueTrees: parent still ok? n = " << p->n() <<
					", (" << p->x() << "," << p->y() << "," << p->z() << ")" 
					<< std::endl;
			    	}
		 		// remove treeList[i]
			    	if (debug)
			    	{
			    	    std::cout << "glueTrees: erase treeList[" << i << "]" << std::endl;
			    	}
                	    	treeList->erase((treeList->begin()+i));
			    	if (debug)
			    	{
				    std::cout << "glueTrees: parent still ok? n = " << p->n() <<
					", (" << p->x() << "," << p->y() << "," << p->z() << ")" 
					<< std::endl;
			    	}
			    	N--;
			    	j = i;
			    }
			}
			else
			{
			    if (debug)
			    {
				std::cout << "glueTrees: this is not the right place for this tree"
					<< std::endl;
			    }
			}
		    }
		}
		else
		{
		    if (debug)
		    {
			std::cout << "glueTrees: treeList[" << i << "] is not descendant of treeList[" 
				<< j << "]" << std::endl;
		    }
		}
	    } while (j < N-1);
	}
	if (debug)
	{
	    std::cout << "glueTrees: at very end of function, about to return" << std::endl;
	    for (int q = 0; q < N; q++)
	    {	
		std::cout << "Subtree " << q << " of " << N << ":" << std::endl;
		(*treeList)[q]->depthFirstTraverseParents();
	    }
	}
    };


    template <class T>
    SharedPtr<OctTree<char> > createGhostTree(std::vector<SharedPtr<OctTree<T> > > *treeList, 
	std::vector<RootList> *ownerList)
    {
	Communicator comm; 
	ProcessID me = comm.rank();
	ProcessID np = comm.nproc();

	SharedPtr<OctTree<char> > ghostTree = SharedPtr<OctTree<char> > (new OctTree<char> ());

//	bool debug = true;
	bool debug = false;

	if (me != 0)
	{
	    archive::MPIOutputArchive arsend(comm, 0);
	    int llen = 0;
	    if (treeList) 
		llen = treeList->size();

	    if (debug)
	    {
		std::cout << "createGhostTree: about to send " << llen << " trees from my list"
			<< std::endl;
	    }

	    arsend & llen;

	    for (int i = 0; i < llen; i++)
	    {
		OctTree<char> *gt = new OctTree<char>();
		gt = (*treeList)[i]->skeletize();
		if (debug)
	 	{
		    std::cout << "createGhostTree: skeletized tree number " << i << ":" << std::endl;
		    gt->depthFirstTraverseAll();
		    std::cout << "createGhostTree: end of tree number " << i << std::endl;
		}
		arsend & *gt;
	    }
	    if (debug)
	    {
		std::cout << "createGhostTree: sent " << llen << " trees from my list"
			<< std::endl;
	    }
	    return ghostTree;
	}

	std::vector<SharedPtr<OctTree<char> > > gtList;
	int llen = 0;
	if (treeList)
	    llen = treeList->size();
	if (debug)
	{
		std::cout << "createGhostTree: about to skeletize my trees" << std::endl;
	}
	for (int i = 0; i < llen; i++)
	{
	    OctTree<char> *gt = new OctTree<char>();
	    gt = (*treeList)[i]->skeletize();
	    gtList.push_back(SharedPtr<OctTree<char> > (gt));
	    if (gt->islocal())
	    {
	        ownerList->push_back(RootList(gt->x(), gt->y(), gt->z(), gt->n(), 0, 0));
	    }
	    else
	    {
		FOREACH_LOCAL_CHILD(SharedPtr<OctTree<char> >, gt,
		    ownerList->push_back(RootList(child->x(), child->y(), child->z(), child->n(), 0, 0));
		);
	    }
	}
	if (debug)
	{
		std::cout << "createGhostTree: about to recv skeletized trees" << std::endl;
	}
	for (ProcessID i = 1; i < np; i++)
	{
	    archive::MPIInputArchive arrecv(comm, i);
	    arrecv & llen;
	    for (int j = 0 ; j < llen; j++)
	    {
		if (debug)
		{
		    std::cout << "createGhostTree: about to receive tree number " << j << " of " << llen
			<< " from process " << i << std::endl;
		}
		OctTree<char> *gt = new OctTree<char>();
		arrecv & *gt;
		if (debug)
		{
		    std::cout << "createGhostTree: just received the following tree:" << std::endl;
		    gt->depthFirstTraverseAll();
		}
		gtList.push_back(SharedPtr<OctTree<char> > (gt));
		if (gt->islocal())
	 	{
		    ownerList->push_back(RootList(gt->x(), gt->y(), gt->z(), gt->n(), i, i));
		    if (debug)
		    {
			std::cout << "createGhostTree: above tree is local " << std::endl;
		    }
		}
		else
		{
		    if (debug)
		    {
			std::cout << "createGhostTree: above tree is remote" << std::endl;
		    }
		    int tmp = i;
		    FOREACH_LOCAL_CHILD(SharedPtr<OctTree<char> >, gt,
		    	ownerList->push_back(RootList(child->x(), child->y(), child->z(), 
				child->n(), tmp, tmp));
			if (debug)
			{
			    std::cout << "createGhostTree: pushed n = " << (*ownerList).back().n << 
				" (" << (*ownerList).back().x << "," << (*ownerList).back().y << "," << 
				(*ownerList).back().z << "), " << (*ownerList).back().current_owner << 
				"->" << (*ownerList).back().future_owner <<  " onto ownerList" << std::endl;
			}
		    );
		}
	    }
	}

	if (debug)
	{
	    llen = gtList.size();
	    std::cout << "createGhostTree: after collecting all the trees from all the procs, we have:"
			<< std::endl;
	    for (int i = 0; i < llen; i++)
	    {
		std::cout << "Subtree:"<< std::endl;
		gtList[i]->depthFirstTraverseAll();
	    }
	    std::cout << std::endl;
	    int olen = ownerList->size();
	    for (int i = 0; i < olen; i++)
	    {
		std::cout << "    ownerList[" << i << "]: n = " << (*ownerList)[i].n << 
		"(" << (*ownerList)[i].x << "," << (*ownerList)[i].y << "," << (*ownerList)[i].z << "), " 
		<< (*ownerList)[i].current_owner << "->" << (*ownerList)[i].future_owner << std::endl;
	    }
	    std::cout << std::endl;
	}
//	sort(gtList.begin(), gtList.end());
	sort(gtList.begin(), gtList.end(), lessPtr<OctTree<char> > ());
	if (debug)
	{
	    std::cout << "createGhostTree: after sorting, before gluing" << std::endl;
	}
	glueTrees(&gtList);
	if (debug) {
	    std::cout << "createGhostTree: after gluing" << std::endl;
	    llen = gtList.size();
	    for (int i = 0; i < llen; i++)
	    {
		std::cout << "Subtree:" << std::endl;
		gtList[i]->depthFirstTraverseAll();
	    }
	}
	sort(ownerList->begin(), ownerList->end());
	if (debug) {
	    std::cout << "createGhostTree: after sorting" << std::endl;
	}
	llen = gtList.size();
	if (llen != 1)
	{
	    std::cout << "createGhostTree: error, there are " << llen << " trees" << std::endl;
	}
	ghostTree = SharedPtr<OctTree<char> > (gtList[0]);
	if (debug) {
	    std::cout << "createGhostTree: after assigning" << std::endl;
	}
	ghostTree->makeAllLocal();
	ghostTree->computeCost();
	if (debug)
	{
	    ghostTree->depthFirstTraverseAll();
	}
	return ghostTree;
    }

    
    std::vector<RootList> createMergeList(std::vector<RootList> globalList, 
	std::vector<RootList> ownerList)
    {
/************************************************************************************************/
/* globalList: list of roots of subtrees and where they belong in the final scheme of things	*/
/*	(future_owner accurate; current_owner inaccurate)					*/
/* ownerList: list of roots of subtrees as currently distributed (current_owner accurate; 	*/
/*	future_owner inaccurate)								*/
/* mergeList: final list of roots of subtrees indicating current and future owners of each	*/
/*	subtree; feed this into exchangeTrees							*/
/************************************************************************************************/
	int olen = ownerList.size(), glen = globalList.size(), mlen, j = 0;
	std::vector<RootList> mergeList;

//	bool debug = true;
	bool debug = false;
	
	for (int i = 0; i < olen; i++)
	{
	    ownerList[i].future_owner = -1;
	    mergeList.push_back(ownerList[i]);
	}
	for (int i = 0; i < glen; i++)
	{
	    globalList[i].current_owner = -1;
	    mergeList.push_back(globalList[i]);
	}

	mlen = mergeList.size();
	if (debug)
	{
	    std::cout << "createMergeList: merge list of length " << mlen << 
			" before sorting" << std::endl;
	    for (int i = 0; i < mlen; i++)
	    {
		std::cout << "MergeList[" << i << "]:" << std::endl;
		std::cout << "    n = " << mergeList[i].n << " (" << mergeList[i].x << "," <<
			mergeList[i].y << "," << mergeList[i].z << "), " << mergeList[i].current_owner
			<< " -> " << mergeList[i].future_owner << std::endl;
	    }
	    std::cout << std::endl;
	}
	sort(mergeList.begin(), mergeList.end());
	if (debug)
	{
	    std::cout << "createMergeList: merge list of length " << mlen << 
			" before processing" << std::endl;
	    for (int i = 0; i < mlen; i++)
	    {
		std::cout << "MergeList[" << i << "]:" << std::endl;
		std::cout << "    n = " << mergeList[i].n << " (" << mergeList[i].x << "," <<
			mergeList[i].y << "," << mergeList[i].z << "), " << mergeList[i].current_owner
			<< " -> " << mergeList[i].future_owner << std::endl;
	    }
	}

	while (j < mlen)
	{
	    for (int i = 1; i < mlen-j; i++)
	    {
		if (mergeList[j].equals(mergeList[j+i]))
		{
		    if (mergeList[j].future_owner==-1)
		    {
			mergeList[j].future_owner = mergeList[j+i].future_owner;
			if (debug)
			{
			    std::cout << "createMergeList: ";
			    std::cout << "mergeList[" << j << "] == mergeList[" << j+i << "]" << std::endl;
			    std::cout << "    n = " << mergeList[j].n << " (" << mergeList[j].x << "," <<
				mergeList[j].y << "," << mergeList[j].z << "), " << 
				mergeList[j].current_owner << " -> " << mergeList[j].future_owner 
				<< std::endl;
			}
		    }
		    else
		    {
			mergeList[j].current_owner = mergeList[j+i].current_owner;
			if (debug)
			{
			    std::cout << "createMergeList: ";
			    std::cout << "mergeList[" << j << "] == mergeList[" << j+i << "]" << std::endl;
			    std::cout << "    n = " << mergeList[j].n << " (" << mergeList[j].x << "," <<
				mergeList[j].y << "," << mergeList[j].z << "), " << 
				mergeList[j].current_owner << " -> " << mergeList[j].future_owner 
				<< std::endl;
			}
		    }
		    mergeList.erase(mergeList.begin()+j+i);
		    if (debug)
		    {
			std::cout << "createMergeList: erased entry " << j+i << std::endl;
		    }
		    mlen--;
		    break;
		}
		else if (mergeList[j].isDescendant(mergeList[j+i]))
		{
		    if ((mergeList[j].future_owner==-1) && (mergeList[j+i].future_owner != -1))
		    {
			mergeList[j].future_owner = mergeList[j+i].future_owner;
			if (debug)
			{
			    std::cout << "createMergeList: changed future owner of entry " << j <<
			    	"    n = " << mergeList[j].n << " (" << mergeList[j].x << "," <<
				mergeList[j].y << "," << mergeList[j].z << "), " << 
				mergeList[j].current_owner << " -> " << mergeList[j].future_owner 
				<< "  to " << mergeList[j].future_owner << std::endl;
			}
			break;
		    }
		    else if ((mergeList[j].current_owner==-1) && (mergeList[j+i].current_owner != -1))
		    {
			mergeList[j].current_owner = mergeList[j+i].current_owner;
			if (debug)
			{
			    std::cout << "createMergeList: changed current owner of entry " << j <<
			    	"    n = " << mergeList[j].n << " (" << mergeList[j].x << "," <<
				mergeList[j].y << "," << mergeList[j].z << "), " << 
				mergeList[j].current_owner << " -> " << mergeList[j].future_owner 
				<< "  to " << mergeList[j].current_owner << std::endl;
			}
			break;
		    }
	 	}
	    }
	    j++;
	}

	if (debug)
	{
	    std::cout << "createMergeList: final list: " << std::endl;
	    for (int i = 0; i < mlen; i++)
	    {
		std::cout << "MergeList[" << i << "]:" << std::endl;
		std::cout << "    n = " << mergeList[i].n << " (" << mergeList[i].x << "," <<
			mergeList[i].y << "," << mergeList[i].z << "), " << mergeList[i].current_owner
			<< " -> " << mergeList[i].future_owner << std::endl;
	    }
	}

	return mergeList;
    }


    template <class T>
    void serialLoadBalance(std::vector<SharedPtr<OctTree<T> > > *treeList)
    {
	Communicator comm;
	ProcessID np, me;
	std::vector<RootList> ownerList, mergeList, globalList;
	SharedPtr<OctTree<char> > ghostTree = SharedPtr<OctTree<char> > ();

//	bool debug = true;
	bool debug = false;

//	bool dumptree = true;
	bool dumptree = false;

	me = comm.rank();
	np = comm.nproc();

        MPI::COMM_WORLD.Barrier();
	if (debug)
	{
	    std::cout << "serialLoadBalance: before creating ghostTree" << std::endl;
	}
	ghostTree = createGhostTree(treeList, &ownerList);
	if (debug)
	{
	    std::cout << "serialLoadBalance: after creating ghostTree" << std::endl;
	}
	ghostTree->treeDiagnostics();
        MPI::COMM_WORLD.Barrier();
	if (me == 0)
	{
	    ghostTree->serialPartition(np, &globalList);
	    std::cout << "    number of broken links = " << globalList.size() << std::endl;
	    if (debug)
	    {
		std::cout << "serialLoadBalance: after serialPartition of ghostTree" << std::endl;
	    }
	    if (dumptree)
	    {
		ghostTree->dumpTree();
	    }
	    mergeList = createMergeList(globalList, ownerList);
	    if (debug)
	    {
		std::cout << "serialLoadBalance: after creating mergeList" << std::endl;
	    }
	    for (ProcessID i = 1; i < np; i++)
	    {
		archive::MPIOutputArchive arsend(comm, i);
		arsend & mergeList;
	    	if (debug)
	    	{
		    std::cout << "serialLoadBalance: after sending mergeList" << std::endl;
	    	}
	    }
            MPI::COMM_WORLD.Barrier();
	}
	else
	{
	    archive::MPIInputArchive arrecv(comm, 0);
	    arrecv & mergeList;
            MPI::COMM_WORLD.Barrier();
	    if (debug)
	    {
		std::cout << "serialLoadBalance: after receiving mergeList" << std::endl;
	    }
	}
        MPI::COMM_WORLD.Barrier();
	if (debug)
	{
	    std::cout << "serialLoadBalance: before exchanging trees" << std::endl;
	}
	exchangeTrees(&mergeList, treeList, true);
	if (debug)
	{
	    std::cout << "serialLoadBalance: after exchanging trees" << std::endl;
	}
	if (me == 0)
	{
	    for (ProcessID i = 1; i < np; i++)
	    {
		archive::MPIOutputArchive arsend(comm, i);
		arsend & globalList;
		if (debug)
		{
		    std::cout << "serialLoadBalance: just finished sending globalList to " << i
			<< std::endl;
		}
	    }
	}
	else
	{
	    if (debug)
	    {
		std::cout << "serialLoadBalance: about to receive globalList from 0" << std::endl;
	    }
	    archive::MPIInputArchive arrecv(comm, 0);
	    arrecv & globalList;
	}

	if (debug)
	{
	    std::cout << "serialLoadBalance: after send/recv globalList" << std::endl;
	}

	int glen = globalList.size();
	int tlen = treeList->size();
	if (debug)
	{
	    std::cout << "serialLoadBalance: length of treeList = " << tlen << std::endl;
	}
	for (int i = 0; i < tlen; i++)
	{
	    RootList rli = RootList(((*treeList)[i]).get(), me, me);
	    for (int j = 0; j < glen; j++)
	    {
		if (rli.isDescendant(globalList[j]))
		{
		    if (debug)
		    {
			std::cout << "serialLoadBalance: about to make tree n = " << rli.n <<
				", (" << rli.x << "," << rli.y << "," << rli.z << ") into child of  n = " <<
				rli.n-1 << ",(" << rli.x/2 << "," << rli.y/2 << "," << rli.z/2 << 
				"), rank =" << globalList[j].future_owner << std::endl;
		    }
		    SharedPtr<OctTree<T> > pptr = SharedPtr<OctTree<T> > ((*treeList)[i]->setParentChild(new 
				OctTree<T> (rli.n-1, rli.x/2, rli.y/2, rli.z/2, true, 0, 
				globalList[j].future_owner, &comm)));
		    (*treeList)[i] = SharedPtr<OctTree<T> >(pptr);
		    break;
		}
		else if (globalList[j].isDescendant(rli))
		{
		    OctTree<T> *u = new OctTree<T>();
		    if (debug)
		    {
		    	std::cout << "serialLoadBalance: about to findDown for tree n = " << 
				globalList[j].n << ", (" << globalList[j].x << "," << 
				globalList[j].y << "," << globalList[j].z << ")" << std::endl;
		    }
		    u = (*treeList)[i]->findDown(globalList[j].n, globalList[j].x, globalList[j].y,
				globalList[j].z);
		    if (debug)
		    {
		    	std::cout << "serialLoadBalance: just finished findDown for tree n = " << 
				globalList[j].n << ", (" << globalList[j].x << "," << 
				globalList[j].y << "," << globalList[j].z << ")" << std::endl;
		    }
		    if ((u) && (globalList[j].n == u->n()) && (globalList[j].x == u->x()) && 
			(globalList[j].y == u->y()) && (globalList[j].z == u->z()))
		    {
			if (debug)
			{
			    std::cout << "serialLoadBalance: setting rank of tree n = " << u->n() <<
				", (" << u->x() << "," << u->y() << "," << u->z() << ") to " <<
				globalList[j].future_owner << std::endl;
			}
			u->setRank(globalList[j].future_owner);
			if (debug)
			{
			    std::cout << "serialLoadBalance: done setting rank of tree n = " << u->n() <<
				", (" << u->x() << "," << u->y() << "," << u->z() << ") to " <<
				globalList[j].future_owner << std::endl;
			}
		    }
		}
	    }
	}

	if (debug)
	{
	    std::cout << "serialLoadBalance: after making all those parents and stuff" << std::endl;
	    int tlen = treeList->size();
	    for (int i = 0; i < tlen; i++)
	    {
		std::cout << "Subtree: " << std::endl;
		(*treeList)[i]->depthFirstTraverseAll();
	    }
	    std::cout << std::endl;
	}

	glueTrees(treeList);
	tlen = treeList->size();
	if (debug)
	{
	    std::cout << "serialLoadBalance: waiting for barrier" << std::endl;
	    MPI::COMM_WORLD.Barrier();
	    std::cout << "serialLoadBalance: after barrier" << std::endl;
	}
    }


    template <class T>
    SharedPtr<OctTree<T> > findPtr(SharedPtr<OctTree<T> > tree, Level n, Translation x, 
	Translation y, Translation z)
    {
	OctTree<T> *tfound = new OctTree<T>();
	tfound = tree->find(n, x, y, z);
	if (tfound)
	{
	    if (tfound->parent())
	    {
		Level nn = n - tfound->parent()->n();
                Translation xx = (x>>nn), yy = (y>>nn), zz = (z>>nn);
		return tfound->parent()->child(xx, yy, zz);
	    }
	    else
		return SharedPtr<OctTree<T> > (tfound);
	}
	else
	{
	    return 0;
	}
    };


}

//template void serialLoadBalance<double>(std::vector< SharedPtr<OctTree<double> > >*);
//#endif
