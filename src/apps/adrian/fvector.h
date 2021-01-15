#include <cassert>
#include <cstddef>

// grabs the underlying type of a type...for example int
template <typename T> struct underlying_type_traits {
  typedef T underlying_type;
};

// macro
#define UNDERLYING_TYPE(T) typename underlying_type_traits<T>::underlying_type

template <typename T> inline void swap(T &x, T &y) {
  // This underlying type does not call a destructor
  UNDERLYING_TYPE(T) tmp;
  move_raw(x, tmp);
  move_raw(y, x);
  // at the end of this move tmp is in an undefined state
  // this is the purpose of underlying_type
  move_raw(tmp, y);
}

template <typename T> inline void cycle_left(T &x1, T &x2, T &x3) {
  // x1->x2-x3-x1
  UNDERLYING_TYPE(T) tmp;
  move_raw(x1, tmp); // tmp=x1
  move_raw(x2, x1);  // x1=x2
  move_raw(x3, x2);  // x2=x3
  move_raw(tmp, x3); // x3=x1
}
template <typename T> inline void cycle_right(T &x1, T &x2, T &x3) {
  // x1->x2-x3-x1
  cycle_left(x3, x2, x1); // just call cycle in reverse order
}
class fvector_int {
private:
  size_t length; // size of allocated area
  int *v;        // v points to the allocated area
public:
  // default constructor
  fvector_int() : length(size_t(0)), v(NULL) {}
  // copy constructor
  fvector_int(const fvector_int &x);
  explicit fvector_int(size_t n) : v(new int[n]), length(n) {}
  ~fvector_int() { delete[] v; }
  fvector_int &operator=(const fvector_int &x);
  friend size_t size(const fvector_int &x) { return x.length; }
  int &operator[](size_t n) {
    assert(n < length);
    return v[n];
  }
  const int &operator[](size_t n) const {
    assert(n < length);
    return v[n];
  }
  struct underlying_type {
    size_t length;
    int *v;
  };
  friend void move_raw(fvector_int &x, underlying_type &y) {
    y.length = x.length;
    y.v = x.v;
  }

  friend void move_raw(underlying_type &x, fvector_int &y) {
    y.length = x.length;
    y.v = x.v;
  }
  friend void move_raw(fvector_int &x, fvector_int &y) {
    y.length = x.length;
    y.v = x.v;
  }
};
template <> struct underlying_type_traits<fvector_int> {
  typedef fvector_int::underlying_type underlying_type;
};

fvector_int::fvector_int(const fvector_int &x)
    : length(size(x)), v(new int[size(x)]) {
  for (std::size_t i = 0; i < length; ++i)
    (*this)[i] = x[i];
}

fvector_int &fvector_int ::operator=(const fvector_int &x) {
  if (this != &x) {
    if (size(*this) == size(x)) {
      for (size_t i = 0; i < size(*this); ++i) {
        int tmp(x[i]); // uses more mem
        (*this)[i] = x[i];
      }
    } else {
      fvector_int tmp(x);
      swap(*this, tmp);
    }
  }
  return *this;
} // problem...if there is an exception during the construction

inline void move(fvector_int &x, fvector_int &y) { swap(x, y); }

bool operator==(const fvector_int &x, const fvector_int &y) {
  if (size(x) != size(y))
    return false;
  for (size_t i = 0; i < size(x); ++i)
    if (x[i] != y[i])
      return false;
  return true;
}
bool operator!=(const fvector_int &x, const fvector_int &y) {
  return !(x == y);
}
bool operator<(const fvector_int &x, const fvector_int &y) {
  for (size_t i(0);; ++i) {
    if (i >= size(y))
      return false;
    if (i >= size(x))
      return true;
    if (y[i] < x[i])
      return false;
    if (x[i] < y[i])
      return true;
  }
}
inline bool operator>(const fvector_int &x, const fvector_int &y) {
  return y < x;
}

inline bool operator<=(const fvector_int &x, const fvector_int &y) {
  return !(y < x);
}

bool operator>=(const fvector_int &x, const fvector_int &y) { return !(x < y); }

size_t areaof(const fvector_int &x) {
  return size(x) * sizeof(int) + sizeof(fvector_int);
}

double memory_utilization(const fvector_int &x) {
  double useful(size(x) * sizeof(int));
  double total(areaof(x));
  return useful / total;
}