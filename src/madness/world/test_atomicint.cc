/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#include <madness/world/world.h>
#include <madness/world/thread.h>
#include <iostream>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>

#include <madness/world/atomicint.h>

#define MAX_NUM_THREADS 1000

using namespace madness;
using namespace std;

AtomicInt sum;
AtomicInt ndone;

void* doit(void *args) {
  // int* a = new int; // work around libtcmalloc bug?
    for (int j=0; j<1000; ++j) {
        for (int i=0; i<100000; ++i) {
            sum++;
        }
        sched_yield();
    }
    ndone++;

    // delete a; // NAR 3/11/2013, probably caused from missing madness::initialize

    return 0;
}


class Greet : public PoolTaskInterface {
public:
    void run(const TaskThreadEnv& env) {
        std::cout << "HI\n";
    }


};

const int NDO = 10000000;

class Adder : public PoolTaskInterface {
public:
    void run(const TaskThreadEnv& env) {
        if (sum >= NDO) {
            ndone++;
        }
        else {
            sum++;
        }
    }
};

int get_nthread() {
    int nthread;
    char *cnthread = getenv("MAD_NUM_THREADS");
    if (cnthread) {
        int result = sscanf(cnthread, "%d", &nthread);
        if (result != 1)
            MADNESS_EXCEPTION("MAD_NUM_THREADS is not an integer", result);
    }
    else {
        nthread = 3; /* this was the hard-coded value in the test previously */
    }
    return nthread;
}

int main(int argc, char** argv) {
    bool smalltest = false;
    if (getenv("MAD_SMALL_TESTS")) smalltest=true;
    for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
    std::cout << "small test : " << smalltest << std::endl;
    if (smalltest) return 0;
    
    madness::initialize(argc,argv);
    const int nthread = get_nthread();
    if (nthread > MAX_NUM_THREADS)
        MADNESS_EXCEPTION("MAD_NUM_THREADS is too large", nthread);
    Thread* threads = new Thread[nthread];

    try {
        sum = ndone = 0;
        for (int i=0; i<nthread; ++i) threads[i].start(doit,0);
        while (ndone != nthread) sleep(1);
        cout << "SUM " << sum << endl;

        sum = ndone = 0;
        ThreadPool::add(new Greet());
        for (int i=0; i<(NDO+1000); ++i) ThreadPool::add(new Adder());
        while (!ndone) {
            sleep(1);
            cout << int(ndone) << " " << int(sum) << endl;
        }
        sleep(1);
        cout << "SUM " << int(sum) << endl;
    }
    catch (const char * e) {
        cout << "string exception: " << e << endl;
    }

    delete [] threads;

    madness::finalize();
    return 0;
}
