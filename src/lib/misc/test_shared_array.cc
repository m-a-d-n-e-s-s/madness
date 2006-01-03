#ifndef SHARED_ARRAY_H
#define SHARED_ARRAY_H

/// \file test_shared_array.cc
/// \brief Minimal shared_array solely for machines on which we do not have BOOST

class shared_counter {
private:
  int count;
  shared_counter(const shared_counter& x) {}; // verboten
  void operator=(const shared_counter& x) {}; // verboten

public:
  shared_counter() : count(1) {};

  inline int get() {return count;};

  inline int inc() {return ++count;};

  inline int dec() {return --count;};
};


// Easy to extend this to general shared object by templating
// appropriate destructor

template <class T>
class shared_array {
private:
  T* p;
  shared_counter *count;

  void dec() {
    if (count) {
      if (count->dec() == 0) {
	if (p) {
	  //cout << "freeing array\n";
	  delete [] p;
	}
	delete count;
      }
    }
  };

public:

  shared_array(T* ptr = 0) : p(ptr), count(0) {
    if (p) count = new shared_counter;
  };

  shared_array(const shared_array& s) : p(s.p), count(s.count) {
    if (count) count->inc();
  };

  ~shared_array() {dec();};

  shared_array& operator=(const shared_array& s) {
    if (this != &s) {
      dec();
      p = s.p;
      count = s.count;
      if (count) count->inc();
    }
    return *this;
  };

  inline int use_count() const {
    if (count) return count->get();
    else return 1;
  };

  inline T* get() const {return p;};

};

/*

#include <iostream>

using namespace std;

int main() {

  shared_array<int> a(new int(100));
  cout << "a should be 1 " << a.use_count() << " " << (void *) a.get() << endl;

  shared_array<int> b(new int(100));
  cout << "b should be 1 " << b.use_count() << " " << (void *) b.get() << endl;

  cout << "assignment ... should see destructor if printing there\n";
  b = a;
  cout << "a should be 2 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 2 " << b.use_count() << " " << (void *) b.get() << endl;

  cout << "copy constructor\n";
  shared_array<int> c(a);
  cout << "a should be 3 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 3 " << b.use_count() << " " << (void *) b.get() << endl;
  cout << "c should be 3 " << c.use_count() << " " << (void *) c.get() << endl;

  cout << "default constructor\n";
  shared_array<int> d;
  cout << "a should be 3 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 3 " << b.use_count() << " " << (void *) b.get() << endl;
  cout << "c should be 3 " << c.use_count() << " " << (void *) c.get() << endl;
  cout << "d should be 1 " << d.use_count() << " " << (void *) d.get() << endl;

  cout << "a=d\n";
  a = d;
  cout << "a should be 1 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 2 " << b.use_count() << " " << (void *) b.get() << endl;
  cout << "c should be 2 " << c.use_count() << " " << (void *) c.get() << endl;
  cout << "d should be 1 " << d.use_count() << " " << (void *) d.get() << endl;

  cout << "b=d\n";
  b = d;
  cout << "a should be 1 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 1 " << b.use_count() << " " << (void *) b.get() << endl;
  cout << "c should be 1 " << c.use_count() << " " << (void *) c.get() << endl;
  cout << "d should be 1 " << d.use_count() << " " << (void *) d.get() << endl;

  cout << "c=d ... should see destructor of original a if printing there\n";
  c = d;
  cout << "a should be 1 " << a.use_count() << " " << (void *) a.get() << endl;
  cout << "b should be 1 " << b.use_count() << " " << (void *) b.get() << endl;
  cout << "c should be 1 " << c.use_count() << " " << (void *) c.get() << endl;
  cout << "d should be 1 " << d.use_count() << " " << (void *) d.get() << endl;

  return 0;
}

*/


#endif
