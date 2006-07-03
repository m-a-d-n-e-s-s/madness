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
#define GLOBAL_COST 2
#define GLOB_PART_SUB_SIGNAL 3

namespace madness {

    template <class T>
    struct less {
        bool operator() (const T p1, const T p2) {
            if (!p1)
                return true;
            if (!p2)
                return false;
            return (!(*p1 < *p2));
        }
    };

    /****************************************************************************************/
    /* class RootList: shorthand, minimal way of communicating that an OctTree<T> needs to  */
    /* be moved, and from current_owner to future_owner.  All public for ease of use.	*/
    /****************************************************************************************/

    class RootList {
    public:
        Translation x, y, z;
        Level n;
        ProcessID current_owner, future_owner;
        RootList(): x(0), y(0), z(0), n(-1), current_owner(0), future_owner(0) {};
        RootList(Translation x0, Translation y0, Translation z0, Level n0, ProcessID current,
                 ProcessID future): x(x0), y(y0), z(z0), n(n0), current_owner(current),
        future_owner(future) {};
        template <class T>
        RootList(OctTree<T> *t, ProcessID current, ProcessID future): x(t->x()), y(t->y()),
        z(t->z()), n(t->n()), current_owner(current), future_owner(future) {};

        template <class Archive>
        inline void serialize(const Archive& ar) {
            ar & x & y & z & n & current_owner & future_owner;
        }

        bool equals(RootList r) {
            if ((n == r.n) && (x == r.x) && (y == r.y) && (z == r.z))
                return true;
            else
                return false;
        }

        bool isDescendant(RootList r) {
            if (r.n >= n)
                return false;

            Level dn = (Level) pow(2, n-r.n);
            if ((x/dn == r.x) && (y/dn == r.y) && (z/dn == r.z))
                return true;
            else
                return false;
        }


        friend bool operator < (const RootList& t1, const RootList& t2) {
            //std::cout << "operator <: beginning" << std::endl;
            if (t1.n > t2.n) {
                //std::cout << "operator <: t1.n > t2.n" << std::endl;
                return true;
            } else if (t1.n < t2.n) {
                //std::cout << "operator <: t1.n < t2.n" << std::endl;
                return false;
            } else {
                //std::cout << "operator <: computing n1 and n2" << std::endl;
                Translation s1, s2, n1 = (Translation) pow(2,t1.n-1), n2 = (Translation) pow(2,t2.n-1);
                //std::cout << "operator <: n1 = " << n1 << ", n2 = " << n2 << std::endl;
                if (n1 == 0) n1 = 1;
                if (n2 == 0) n2 = 1;
                s1 = (t1.x/n1)*4 + (t1.y/n1)*2 + t1.z/n1;
                s2 = (t2.x/n2)*4 + (t2.y/n2)*2 + t2.z/n2;
                //std::cout << "operator <: computed s1 and s2" << std::endl;
                if (s1 < s2)
                    return true;
                else if (s1 > s2)
                    return false;
                else {
                    if (t1.x < t2.x)
                        return true;
                    else if (t1.x > t2.x)
                        return false;
                    else {
                        if (t1.y < t2.y)
                            return true;
                        else if (t1.y > t2.y)
                            return false;
                        else {
                            if (t1.z < t2.z)
                                return true;
                            else
                                return false;
                        }
                    }
                }
            }
        }
    };

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
                      if (debug) {
                              std::cout << "globalPartitionSubtree: don't bother with child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                      }
                      else {
                          keepgoing = true;
                      }
                      if ((child->islocal()) && (keepgoing)) {
                      if (debug) {
                              std::cout << "globalPartitionSubtree: keep going with child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                          Cost cost = child->getLocalSubtreeCost();
                          if (debug) {
                              std::cout << "globalPartitionSubtree: child (" << child->x() << "," <<
                              child->y() << "," << child->z() << "), with cost " << cost <<
                              " is local" << std::endl;
                          }
                          if (cost <= *spaceleft) {
                              *spaceleft -= cost;
                              tmpList.push_back(RootList(child, me, p));
                              child->setSendto(p);
                              if (debug) {
                                  std::cout << "globalPartitonSubtree: added child to list" << std::endl;
                              }
                              std::vector <RootList> remoteList;
                              remoteList = findRemoteChildren(child, p);
                              int rsize = remoteList.size();
                              if (debug) {
                                  std::cout << "globalPartitionSubtree: found " << rsize << "  new remote children"
                                  << std::endl;
                              }
                              for (int i = 0; i < rsize; i++) {
                                  tmpList.push_back(remoteList[i]);
                              }
                          } else {
                              if (debug) {
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
                      if (debug) {
                              std::cout << "globalPartitionSubtree: keep going with remote child " <<
                              "(" << child->x() << "," << child->y() << "," << child->z() << ")"
                              << std::endl;
                          }
                          archive::MPIOutputArchive arsend(comm, child->rank());
                          archive::MPIInputArchive arrecv(comm, child->rank());
                          if (debug) {
                              std::cout << "globalPartitionSubtree: about to send GLOB_PART_SUB_SIGNAL "
                              << GLOB_PART_SUB_SIGNAL << " to processor " << child->rank() <<
                              " concerning partition " << p << " and space left " << *spaceleft
                              << std::endl;
                          }
                          arsend & GLOB_PART_SUB_SIGNAL & RootList(child,child->rank(), p) & *spaceleft;
                          arrecv & success;
                          if (debug) {
                              std::cout << "globalPartitionSubtree:  about to add " << success <<
                              " new roots to list" << std::endl;
                          }
                          if (success > 0) {
                              for (int i = 0; i < success; i++) {
                                  arrecv & tmproot & costofroot;
                                  *spaceleft -= costofroot;
                                  tmproot.future_owner = p;
                                  localList.push_back(tmproot);
                                  if (debug) {
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
                                  if (debug) {
                                      std::cout << "globalPartitionSubtree: put remote subtree n = "
                                      << remoteList[i].n << ", (" << remoteList[i].x << "," <<
                                      remoteList[i].y << "," << remoteList[i].z << "), " <<
                                      remoteList[i].current_owner << " -> " <<
                                      remoteList[i].future_owner << " on globalList" << std::endl;
                                  }
                              } else {
                                  if (debug) {
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
                          if (debug) {
                          std::cout << "partitionWorkerSub: at beginning of FOREACH_CHILD"
                          << std::endl;
                      }
                      partitionWorkerSub(child, localList, spaceleft, costList, remoteList, p);
                      if (*spaceleft == 0) {
                          if (debug) {
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
        madness::redirectio(comm);
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
                std::cout << "exchangeTrees: root " << i << ": n = " << root.n << ", " <<
                "(" << root.x << "," << root.y << "," << root.z << ")" << "current = "
                << root.current_owner << ", future = " << root.future_owner << std::endl;
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
                        (*treeList)[j]->depthFirstTraverseAll();
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
                if (t != (*treeList)[j]) {
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
                    SharedPtr<OctTree<T> > tptr = p->childPtr(x,y,z);
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
                    t->setParentChild(pprime);
                    if (debug) {
                        std::cout << "exchangeTrees: just set t's parent" << std::endl;
                    }
                    OctTree<T> *q = new OctTree<T>();
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
                        (*treeList)[j]->depthFirstTraverseAll();
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
                else {
                    Translation x = (*globalList)[i].x, y = (*globalList)[i].y, z = (*globalList)[i].z;
                    x -= 2*(x/2);
                    y -= 2*(y/2);
                    z -= 2*(z/2);
                    p->insert_remote_child(x, y, z, (*globalList)[i].future_owner);
		    bool deleteit = true;
		    FOREACH_CHILD(OctTree<T>, p,
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
                    (*treeList)[i]->depthFirstTraverseAll();
                }
            }
        }
        if (debug) {
            int tlength = treeList->size();
            std::cout << "exchangeTrees: before sorting and gluing together, list of length " <<
            tlength << ":" << std::endl;
            for (int i = 0; i < tlength; i++) {
                std::cout << "tree " << i << " of " << tlength << ":" << std::endl;
                (*treeList)[i]->depthFirstTraverseAll();
            }
        }

//	sort((*treeList).begin(), (*treeList).end(), less<OctTree<T>* > ());
	sort((*treeList).begin(), (*treeList).end());
        if (glue)
            glueTrees(treeList);
	if (debug)
	{
	    std::cout << "exchangeTrees: about to return" << std::endl;
	}
    }

    template <class T>
    void sendSubtree(OctTree<T> *tree, ProcessID me, ProcessID dest) {
        Communicator comm;
        archive::MPIOutputArchive ardest(comm, dest);

        bool debug = false;
//	bool debug = true;

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

        bool debug = false;
//	bool debug = true;

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
        madness::redirectio(comm);
        comm.print();
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
        madness::redirectio(comm);
        comm.print();
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
        madness::redirectio(comm);
        archive::MPIOutputArchive ardest(comm, dest);

        ardest & msg;
    }

    template <class T>
    void recvMsg(T *msg, ProcessID me, ProcessID source) {
        Communicator comm;
        madness::redirectio(comm);
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
        int size = treeList->size();
        int nconsidered = size, nunchanged = 0;
        bool flag = false;

//        bool debug = true;
	bool debug = false;

        for (int i = 0; i < size; i++) {
            if (debug) {
                std::cout << "glueTrees: at beginning of i loop, i = " << i << ", nunchanged = " <<
                nunchanged << ", nconsidered = " << nconsidered << std::endl;
            }
//            OctTree<T> *t = new OctTree<T>();
            SharedPtr<OctTree<T> > t = SharedPtr<OctTree<T> >();
            t = (*treeList)[nunchanged];
            if (debug) {
                std::cout << "glueTrees: after assigning t to be (*treeList)[" << nunchanged << "] = "
                << "(" << t->x() << "," << t->y() << "," << t->z() << ")" << std::endl;
            }

            for (int j = 1; j < nconsidered; j++) {
                if (debug) {
                    std::cout << "glueTrees: at beginning of j loop, j = " << j << std::endl;
                }
		if (t->n() == (*treeList)[nunchanged+j]->n()) {
		// if this tree is the same level, it could be the same tree, in which case
		// we need to replace one tree with another, inserting its children as appropriate
		    if (debug)
		    {
			std::cout << "glueTrees: the same node is duplicated" << std::endl;
		    }
		    if (t->equals((*treeList)[nunchanged+j])) {
			int tmp = nunchanged+j;
			if (t->isremote())
			{
			    if (debug)
			    {
				std::cout << "glueTrees: the first one is remote" << std::endl;
			    }
			    FORIJK(
//			        if (((child->isremote())||(!child)) && (t->childPtr(i,j,k)->islocal()))
			        if (t->child(i,j,k) && (t->childPtr(i,j,k)->islocal()))
			        {
				    (*treeList)[tmp]->setChild(i,j,k,t->childPtr(i,j,k)); 
					// set the child ptr to this local node
				    if (debug)
				    {
					std::cout << "glueTrees: set child ptr to local ptr (" << i << ","
						<< j << "," << k << ")" << std::endl;
				    }
			        }
				else
				{
				    if (debug)
				    {
					std::cout << "glueTrees: child (" << i << "," << j << "," << k << 
						") does not exist or is not local" << std::endl;
				    }
				}
			    );
                            treeList->erase((treeList->begin()+nunchanged));
                            flag = true;
			}
			else if ((*treeList)[tmp]->isremote())
			{
			    if (debug)
			    {
				std::cout << "glueTrees: the second one is remote" << std::endl;
			    }
			    SharedPtr<OctTree<T> > s = SharedPtr<OctTree<T> > ((*treeList)[tmp]);
			    (*treeList)[tmp] = t;
			    FORIJK(
			        if (s->child(i,j,k) && (s->childPtr(i,j,k)->islocal()))
			        {
				    t->setChild(i,j,k,s->childPtr(i,j,k)); 
					// set the child ptr to this local node
				    if (debug)
				    {
					std::cout << "glueTrees: set child ptr to local ptr (" << i << ","
						<< j << "," << k << ")" << std::endl;
				    }
			        }
				else
				{
				    if (debug)
				    {
					std::cout << "glueTrees: child (" << i << "," << j << "," << k << 
						") does not exist or is not local" << std::endl;
				    }
				}
			    );
//                            treeList->erase((treeList->begin()+tmp));
//			    (*treeList)[tmp] = t;
			    treeList->erase((treeList->begin()+nunchanged));
			    if (debug)
			    {
				std::cout << "glueTrees: erased entry number " << tmp << " in treeList"
					<< std::endl;
			    }
                            flag = true;
			}
			else
			{
			    std::cout << "glueTrees: Huh?  This should not happen" << std::endl;
			}
			break;
		    }
                } else if (t->n() > (*treeList)[nunchanged+j]->n()) {
                    if (t->isDescendant((*treeList)[nunchanged+j])) {
                        if (debug) {
                            std::cout << "glueTrees: t (" << t->x() << "," << t->y() << "," << t->z() << ")"
                            << " is a descendant of (" << (*treeList)[nunchanged+j]->x() << "," <<
                            (*treeList)[nunchanged+j]->y() << "," << (*treeList)[nunchanged+j]->z()
                            << ")" << std::endl;
                        }
                        // find t on this tree
                        OctTree<T> *u = new OctTree<T>();
			u = (*treeList)[nunchanged+j]->findDown(t->n(), t->x(), t->y(), t->z());
                        if (u) {
                            if (debug) {
                                std::cout << "glueTrees: found the tree on the list" << std::endl;
                            }
                            OctTree<T> *p = new OctTree<T>();
                            p = u->parent();
                            if (debug) {
                                std::cout << "glueTrees: made parent to u" << std::endl;
                            }
                            p->insert_local_child(t);
                            if (debug) {
                                std::cout << "glueTrees: inserted local child" << std::endl;
                            }
                            treeList->erase((treeList->begin()+nunchanged));
                            if (debug) {
                                std::cout << "glueTrees: erased tree from treeList" << std::endl;
                            }
                            flag = true;
                        }
                        break;
                    }
                }
            }
	    if (debug)
	    {
		std::cout << "glueTrees: outside of loop, about to change nunchanged and nconsidered"
			<< std::endl;
	    }
            if (!(flag))
                nunchanged++;
            nconsidered--;
            flag = false;
        }
	if (debug)
	{
	    std::cout << "glueTrees: does deletion occur before or after this point?" << std::endl;
	}
    }


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
		FOREACH_LOCAL_CHILD(OctTree<char>, gt,
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
		    FOREACH_LOCAL_CHILD(OctTree<char>, gt,
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
	sort(gtList.begin(), gtList.end());
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

	bool debug = false;
//	bool debug = true;

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
        MPI::COMM_WORLD.Barrier();
	if (me == 0)
	{
	    ghostTree->serialPartition(np, &globalList);
	    if (debug)
	    {
		std::cout << "serialLoadBalance: after serialPartition of ghostTree" << std::endl;
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
	    }
	}
	else
	{
	    archive::MPIInputArchive arrecv(comm, 0);
	    arrecv & globalList;
	}

	int glen = globalList.size();
	int tlen = treeList->size();
	for (int i = 0; i < tlen; i++)
	{
	    RootList rli = RootList(((*treeList)[i]).get(), me, me);
	    for (int j = 0; j < glen; j++)
	    {
		if (rli.isDescendant(globalList[j]))
		{
		    SharedPtr<OctTree<T> > pptr = SharedPtr<OctTree<T> > ((*treeList)[i]->setParentChild(new 
				OctTree<T> (globalList[j].n, globalList[j].x, globalList[j].y, 
				globalList[j].z, true, 0, globalList[j].future_owner, &comm)));
		    (*treeList)[i] = pptr;
		    break;
		}
	    }
	}
	if (debug)
	{
	    std::cout << "serialLoadBalance: after making all those parents and stuff" << std::endl;
	    tlen = treeList->size();
	    for (int i = 0; i < tlen; i++)
	    {
		std::cout << "Subtree: " << std::endl;
		(*treeList)[i]->depthFirstTraverseAll();
	    }
	    std::cout << std::endl;
	}

	glueTrees(treeList);
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
    }

}

#endif
