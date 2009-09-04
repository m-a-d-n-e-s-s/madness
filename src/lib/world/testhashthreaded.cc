#include <world/worldhash.h>
#include <world/worldhashmap.h>
#include <world/worldtime.h>
#include <world/atomicint.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>

/// \file testhashthreaded.cc
/// \brief Test code for parallel hash

using namespace std;
using namespace madness;

void errmsg(const char *msg, int status) {
    cerr << msg << " " << status << std::endl;
    exit(1);
}

double drand() {
    return random()*(1.0/RAND_MAX);
}

void test_coverage() {
    // This test aims for complete code coverage for whatever that
    // is worth, and tests for basic sequential correctness.
    ConcurrentHashMap<int,int> a;
    typedef ConcurrentHashMap<int,int>::datumT datumT;
    typedef ConcurrentHashMap<int,int>::iterator iteratorT;
    typedef ConcurrentHashMap<int,int>::const_iterator const_iteratorT;


    a[-1] = -99;
    if (a[-1] != -99) cout << "was expecting -99 " << a[-1] << endl;
    if (a.size() != 1) cout << "was expecting size to be 1" << endl;

    for (int i=0; i<10000; i++) a.insert(datumT(i,i*99));

    for (int i=0; i<10000; i++) {
        pair<iteratorT,bool> r = a.insert(datumT(i,i*99));
        if (r.second) cout << "expected second insert to fail " << i << endl;
        if (r.first->first != i) cout << "key mismatch " << i << " " << r.first->first << endl;
        if (r.first->second != i*99) cout << "value mismatch " << i << " " << r.first->second << endl;
    }
    if (a.size() != 10001) cout << "was expecting size to be 10001 " << a.size() << endl;

    for (int i=0; i<10000; i++) {
        iteratorT it = a.find(i);
        if (it == a.end()) cout << "expected to find this element " << i << endl;
        if (it->first != i) cout << "key mismatch on find" << i << " " << it->first << endl;
        if (it->second != 99*i) cout << "value mismatch on find" << i << " " << it->second << endl;
    }

    const ConcurrentHashMap<int,int>* ca = &a;
    for (int i=0; i<10000; i++) {
        const_iteratorT it = ca->find(i);
        if (it == ca->end()) cout << "expected to find this element " << i << endl;
        if (it->first != i) cout << "key mismatch on find" << i << " " << it->first << endl;
        if (it->second != 99*i) cout << "value mismatch on find" << i << " " << it->second << endl;
    }

    size_t count = 0;
    for (iteratorT it=a.begin(); it!=a.end(); ++it) {
        count++;
        if (it->second != 99*it->first) cout << "key/value mismatch" << it->first << " " << it->second << endl;
    }
    if (count != 10001) cout << "a. count should have been 10001 " << count << endl;

    count= 0;
    for (const_iteratorT it=ca->begin(); it!=ca->end(); ++it) {
        count++;
        if (it->second != 99*it->first) cout << "key/value mismatch" << it->first << " " << it->second << endl;
    }
    if (count != 10001) cout << "b. count should have been 10001 " << count << endl;

    count = 10001;
    for (int i=0; i<2000; i+=2) {
        iteratorT it = a.find(i);
        a.erase(i);
        count--;
        if (a.size() != count) cout << "size should have been " << count << " " << a.size() << endl;
        it = a.find(i);
        if (it != a.end()) cout << "this was just deleted but was found " << i << endl;
    }

    for (int i=2001; i<4000; i+=2) {
        iteratorT it = a.find(i);
        if (a.erase(i) != 1) cout << "expected to have deleted one element " << i << endl;
        count--;
        if (a.size() != count) cout << "c. size should have been " << count << " " << a.size() << endl;
        it = a.find(i);
        if (it != a.end()) cout << "this was just deleted but was found " << i << endl;
    }

    for (iteratorT it=a.begin(); it!=a.end();) {
        int i = it->first;
        iteratorT prev = it++;
        a.erase(prev);
        count--;
        if (a.size() != count) cout << "d. size should have been " << count << " " << a.size() << endl;
        iteratorT fit = a.find(i);
        if (fit != a.end()) cout << "this was just deleted but was found " << i << endl;

    }

    if (a.size() != 0) cout << "e. size should have been 0 " << a.size() << endl;

    for (int i=0; i<10000; i++) a.insert(datumT(i,i*99));
}


vector<int> random_perm(int n) {
    vector<int> v(n);
    for (int i=0; i<n; i++) v[i] = i;
    for (int i=0; i<n; i++) swap(v[i],v[int(drand()*n)]);
    return v;
}


void test_time() {
    // Examine interaction between nbins and nentries by looping thru
    // bin sizes and measuring time to insert and then delete varying
    // number of keys in random order
    typedef ConcurrentHashMap<int,int>::datumT datumT;
    for (int nbins=100; nbins<=10000; nbins*=10) {
        for (int nentries=nbins; nentries<=nbins*100; nentries*=10) {
            ConcurrentHashMap<int,double> a(nbins);
            vector<int> v = random_perm(nentries);
            double insert_used = madness::cpu_time();
            for (int i=0; i<nentries; i++) {
                a.insert(datumT(i,i));
            }
            insert_used = madness::cpu_time()-insert_used;
            v = random_perm(nentries);
            double del_used = madness::cpu_time();
            for (int i=0; i<nentries; i++) {
                a.erase(i);
            }
            del_used = madness::cpu_time()-del_used;
            printf("nbin=%8d   nent=%8d   insert=%.1es/call   del=%.1es/call\n",
                   nbins, nentries, insert_used/nentries, del_used/nentries);
        }
    }
}

void do_test_random(ConcurrentHashMap<int,double>& a, size_t& count, double& sum) {
    typedef ConcurrentHashMap<int,double>::datumT datumT;
    typedef ConcurrentHashMap<int,double>::iterator iteratorT;
    // Randomly generate keys in range 4*nbin and randomly insert or
    // delete that entry.  Maintain expected sum and count of values
    // and verify at end.
    const int nbin = 131; // A small # bins means more chance of bad thread interaction
    count = 0;
    sum = 0.0;
    for (int i=0; i<100000000; i++) {
        int key = int(drand()*4*nbin);
        double value = key;
        bool choice = (drand() < 0.5);
        if (choice) {
            pair<iteratorT,bool> it = a.insert(datumT(key,value));
            if (it.second) {
                sum += value;
                count++;
            }
        }
        else {
            int ndel = a.erase(key);
            if (ndel == 1) {
                sum -= value;
                count --;
            }
        }
    }
}

void test_random() {
    ConcurrentHashMap<int,double> a(131);
    typedef ConcurrentHashMap<int,double>::iterator iteratorT;

    size_t count;
    double sum;
    do_test_random(a, count, sum);

    double end_sum = 0.0;
    size_t end_count = 0;
    iteratorT ppp=a.begin();
    for (iteratorT it=a.begin(); it!=a.end(); ++it) {
        end_count++;
        end_sum += it->second;
    }

    if (sum != end_sum || count != end_count || count != a.size()) {
        cout << "expected sum and count " << sum << " " << count << endl;
        cout << "  actual sum and count " << end_sum << " " << end_count << " " << a.size() << " " << endl;
    }
}

madness::AtomicInt ndone;

class Worker : public madness::ThreadBase {
private:
    ConcurrentHashMap<int,double>& a; // Better would be a shared pointer
    size_t& count;
    double& sum;

public:
    Worker(ConcurrentHashMap<int,double>& a, size_t& count, double& sum)
            : ThreadBase(), a(a), count(count), sum(sum) {
        start();
    }

    void run() {
        do_test_random(a, count, sum);

        ndone++;
    }
};



void test_thread() {
    ConcurrentHashMap<int,double> a(131);
    typedef ConcurrentHashMap<int,double>::datumT datumT;
    typedef ConcurrentHashMap<int,double>::iterator iteratorT;
    typedef ConcurrentHashMap<int,double>::const_iterator const_iteratorT;
    const int nthread = 2;
    size_t counts[nthread];
    double sums[nthread];

    ndone = 0;

    Worker worker1(a,counts[0],sums[0]);
    Worker worker2(a,counts[1],sums[1]);
    while (ndone != 2) sched_yield();

    size_t count = 0;
    double sum = 0;
    for (int i=0; i<nthread; i++) {
        sum += sums[i];
        count += counts[i];
    }

    double end_sum = 0.0;
    size_t end_count = 0;
    iteratorT ppp=a.begin();
    for (iteratorT it=a.begin(); it!=a.end(); ++it) {
        end_count++;
        end_sum += it->second;
    }

    if (sum != end_sum || count != end_count || count != a.size()) {
        cout << "threaded: expected sum and count " << sum << " " << count << endl;
        cout << "            actual sum and count " << end_sum << " " << end_count << " " << a.size() << " " << endl;
    }
}


class Peasant : public madness::ThreadBase {
private:
    ConcurrentHashMap<int,double>& a; // Better would be a shared pointer

public:
    Peasant(ConcurrentHashMap<int,double>& a)
            : ThreadBase(), a(a) {
        start();
    }

    void run() {
        for (int i=0; i<10000000; i++) {
            ConcurrentHashMap<int,double>::accessor r;
            if (!a.find(r, 1)) throw "OK ... where is it?";
            r->second++;
        }

        ndone++;
    }
};


void test_accessors() {
    ConcurrentHashMap<int,double> a(131);
    typedef ConcurrentHashMap<int,double>::datumT datumT;
    typedef ConcurrentHashMap<int,double>::accessor accessorT;

    ndone = 0;

    accessorT result;
    if (a.find(result,1)) throw "should not have found this";

    a[1] = 0.0;

    if (!a.find(result,1)) throw "should have found this";
    if (result->second != 0.0) throw "should have been zero";

    if (!a.find(result,1)) throw "should have found this";
    if (result->second != 0.0) throw "should have been zero";


    Peasant a1(a),a2(a);
    result.release();
    while (ndone != 2) sched_yield();

    if (a[1] != 20000000) throw "Ooops";
}

int main() {

    try {
        test_coverage();
        test_random();
        test_time();
        test_thread();
        test_accessors();

        cout << "Things seem to be working!\n";
    }
    catch (const char* s) {
        cout << "STRING EXCEPTION: " << s << endl;
    }
    catch (...) {
        cout << "UNKNOWN EXCEPTION: " << endl;
    }

    return 0;
}
