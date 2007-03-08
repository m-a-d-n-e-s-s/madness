#ifndef WORLD_ARRAY
#define WORLD_ARRAY

#include <vector>

namespace madness {

    /// A simple, fixed dimension array 

    /// Like a vector, but better!  More functionality, eliminates
    /// memory allocation cost, is just POD so can be copied easily
    /// and allocated on the stack, and the known dimension permits
    /// agressive compiler optimizations.
    ///
    /// Presently assumes contents are number like things but
    /// this would best be moved out to be a specialization.
    template <typename T, std::size_t N>
    class Array {
    private:
        T v[N];
    public:
        typedef T* iterator;
        typedef const T* const_iterator;

        /// Default constructor does not initialize anything
        Array() {};

        /// Initialize all elements to value t
        Array(T t) {
            for (std::size_t i=0; i<N; i++) v[i] = t;
        };

        /// Construct from a C++ array of the same dimension
        Array(const T (&t)[N]) {
            for (std::size_t i=0; i<N; i++) v[i] = t[i];
        };

        /// Construct from an STL vector of the same dimension
        Array(const std::vector<T> t) {
            MADNESS_ASSERT(t.size() == N);
            for (std::size_t i=0; i<N; i++) v[i] = t[i];
        };
        
        /// Copy constructor is deep
        Array(const Array<T,N>& other) {
            for (std::size_t i=0; i<N; i++) v[i] = other.v[i];
        };

        /// Assignment is deep
        Array& operator=(const Array<T,N>& other) {
            for (std::size_t i=0; i<N; i++) v[i] = other.v[i];
        };

        /// Assignment is deep
        Array& operator=(const std::vector<T>& other) {
            for (std::size_t i=0; i<N; i++) v[i] = other[i];
        };

        /// Fill from scalar value
        Array& operator=(T t) {
            for (std::size_t i=0; i<N; i++) v[i] = t;
        };

        /// Test for element-wise equality
        bool operator==(const Array<T,N>& other) const {
            for (std::size_t i=0; i<N; i++) 
                if (v[i] != other.v[i]) return false;
            return true;
        };

        /// Same as !(*this==other)
        bool operator!=(const Array<T,N>& other) const {
            return !(*this==other);
        };

        /// Tests a<b in sense of lexically ordered index

        /// Can be used to impose a standard ordering on arrays.
        bool operator<(const Array<T,N>& other) const {
            for (std::size_t i=0; i<N; i++) {
                if (v[i] < other.v[i]) return true;
                else if (v[i] > other.v[i]) return false;
            }
            if (v[N-1] == other.v[N-1]) return false; // equality
            return true;
        };

        /// Indexing
        T& operator[](std::size_t i) {
            return v[i];
        };

        /// Indexing
        const T& operator[](std::size_t i) const {
            return v[i];
        };

        /// Element-wise multiplcation by a scalar
        template <typename Q>
        Array<T,N> operator*(Q q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i] * q;
            return r;
        };

        /// In-place element-wise multiplcation by a scalar
        template <typename Q>
        Array<T,N>& operator*=(Q q) {
            for (std::size_t i=0; i<N; i++) v[i] *= q;
            return *this;
        };

        /// Element-wise multiplcation by another array
        template <typename Q>
        Array<T,N> operator*(const Array<Q,N>& q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i]*q[i];
            return r;
        };

        /// Element-wise addition of a scalar
        template <typename Q>
        Array<T,N> operator+(Q q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i] + q;
            return r;
        };

        /// In-place element-wise addition of a scalar
        template <typename Q>
        Array<T,N>& operator+=(Q q) {
            for (std::size_t i=0; i<N; i++) v[i] += q;
            return *this;
        };

        /// Element-wise addition of another array
        template <typename Q>
        Array<T,N> operator+(const Array<Q,N>& q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i] + q[i];
            return r;
        };

        /// Element-wise subtraction of a scalar
        template <typename Q>
        Array<T,N> operator-(Q q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i] - q;
            return r;
        };

        /// Element-wise subtraction of another array
        template <typename Q>
        Array<T,N> operator-(const Array<Q,N>& q) const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = v[i] - q[i];
            return r;
        };

        /// Unary negation
        template <typename Q>
        Array<T,N> operator-() const {
            Array<T,N> r;
            for (std::size_t i=0; i<N; i++) r[i] = -v[i];
            return r;
        };


        /// STL iterator support
        iterator begin() {
            return v;
        };

        /// STL iterator support
        const_iterator begin() const {
            return v;
        };

        /// STL iterator support
        iterator end() {
            return v+N;
        };

        /// STL iterator support
        const_iterator end() const {
            return v+N;
        };

        /// Length of the array
        std::size_t size() const {
            return N;
        };
        
        /// Support for MADNESS serialization
        template <typename Archive>
        void serialize(Archive& ar) {
            ar & v;
        };

        /// Support for MADNESS hashing
        hashT hash() const {
            return madness::hash(v);
        };
    };

    /// Output array to stream for human consumtion
    template <std::size_t N, typename T>
    std::ostream& operator<<(std::ostream& s, const Array<T,N>& a) {
        s << "[";
        for (std::size_t i=0; i<N; i++) {
            s << a[i];
            if (i != (N-1)) s << ",";
        }
        s << "]";
        return s;
    }
}

#endif
