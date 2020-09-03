#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>

#include <madness/misc/ran.h>
#include <madness/world/print.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>

#include "mc.h"

using namespace madness;

namespace std {
    template <> struct hash<std::pair<size_t,size_t>> {
        size_t operator()(const std::pair<size_t,size_t>& ij) const {
            return (ij.first<<32) | ij.second;
        }
    };
}

namespace madness {
    
    // Simple sparse vector class 
    template <typename T=double>
    class sparse_vector {
        size_t M;
        std::vector<size_t> indices;
        std::vector<T> data;
    public:

        // Make an empty sparse vector of length M
        sparse_vector(size_t M) : M(M) {}
        
        // Make a sparse vector from a dense vector with screening
        template <typename vecT>
        sparse_vector(const vecT& x, T tol=0) : M(x.size()) {
            for (size_t i=0; i<M; i++) {
                if (std::abs(x[i]) > tol) {
                    indices.push_back(i);
                    data.push_back(x[i]);
                }
            }
        }

        // Make a dense tensor 
        Tensor<T> dense() const {
            Tensor<T> v(M);
            for (auto [i,x] : *this) v[i] = x;
            return v;
        }
        
        // Return the indices of the non-zero elements
        const std::vector<size_t>& get_indices() const {return indices;}
        std::vector<size_t>& get_indices() {return indices;}

        // Return the values of the non-zero elements
        const std::vector<T>& get_data() const {return data;}
        std::vector<T>& get_data() {return data;}

        // Returns the length (if dense) of the vector
        size_t size() const {return M;}

        // Returns the number of nonzero elements of the vector
        size_t nnz() const {return indices.size();}
        
        // Push (i,value) onto the end of the indices and data
        void push_back(size_t i, T value) {
            MADNESS_ASSERT(i>=0 && i<M);
            indices.push_back(i);
            data.push_back(value);
        }
        
        // Sort the index vector into increasing order
        void sort_indices() {
            using pairT = std::pair<size_t,T>;
            std::vector<pairT> x; x.reserve(nnz());
            for (size_t i=0; i<data.size(); i++) x.push_back({indices[i],data[i]});
            std::sort(x.begin(), x.end(), [](const pairT& a, const pairT& b){return a.first<b.first;});
            indices = std::vector<size_t>(indices.size());
            data = std::vector<T>(data.size());
            for (size_t i=0; i<data.size(); i++) {indices[i]=x[i].first,data[i]=x[i].second;}
        }

        // Compress data keeping entries with absolute value greater than tol
        void compress(T tol=0.0) {
            std::vector<size_t> ind; ind.reserve(nnz());
            std::vector<T> dat; dat.reserve(nnz());
            for (auto [i,x] : *this) {
                if (std::abs(x)>tol) {
                    ind.push_back(i);
                    dat.push_back(x);
                }
            }
            indices = ind; // Assign/copy should shrink capacity to match size but perhaps move is done instead?
            data = dat;
        }

        // return (true,value) if index i exists, (false,0.0) otherwise --- this is presently implemented with dumb linear search and so is slooooow
        std::pair<bool,T> find(size_t i) const {
            for (auto [j,x] : *this) if (j == i) return std::make_pair(true,x);
            return {false,0.0};
        }

        // Non-const iterator dereferences into pair(index&,value&)
        struct iterator {
            std::vector<size_t>& indices;
            std::vector<T>& data;
            size_t i;
            iterator(const iterator& other) : indices(other.indices), data(other.data), i(other.i) {}
            iterator(std::vector<size_t>& indices, std::vector<T>& data, size_t start=0) : indices(indices), data(data), i(start) {}
            iterator& operator++() {++i; return *this;}
            bool operator!=(const iterator& other) const {return i!=other.i;}
            std::pair<size_t&,T&> operator*() {return std::pair<size_t&,T&>(indices[i],data[i]);}
        };

        // Const iterator dereferences into pair(index,value)
        struct const_iterator {
            const std::vector<size_t>& indices;
            const std::vector<T>& data;
            size_t i;
            const_iterator(const iterator& other) : indices(other.indices), data(other.data), i(other.i) {}
            const_iterator(const std::vector<size_t>& indices, const std::vector<T>& data, size_t start=0) : indices(indices), data(data), i(start) {}
            const_iterator& operator++() {++i; return *this;}
            bool operator!=(const const_iterator& other) const {return i!=other.i;}
            std::pair<size_t,T> operator*() {return std::pair<size_t,T>(indices[i],data[i]);}
        };

        // Inplace renumber the indices so that inew = map(iold)
        // For convenience we also sort the new vector of indices/values, which also minimizes the storage
        // map should be a permutation
        void renumber(const std::vector<size_t>& map) {
          sparse_vector<T> v(M);
          for (auto [i,x] : *this) {v.push_back(map[i],x);}
          v.sort_indices();
          *this = v;
        }

        const_iterator begin() const {return const_iterator(indices,data,0);}
        const_iterator end() const {return const_iterator(indices,data,indices.size());}
        iterator begin() {return iterator(indices,data,0);}
        iterator end() {return iterator(indices,data,indices.size());}
    };
    
    // Print sparse vector
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const sparse_vector<T>& v) {
        using madness::operators::operator<<;        
        s<<"(";
        s<<v.get_indices();
        s<<", ";
        s<<v.get_data();
        s<<")";
        return s;
    }
    
    // Simple class storing sparse matrix as map {i,j}->value --- inherits most of interface from std::unordered_map
    template <typename T=double>
    class sparse_matrix_coo : public  std::unordered_map<std::pair<size_t,size_t>,T> {
        size_t N,M;
    public:
        sparse_matrix_coo(size_t N, size_t M) : N(N), M(M) {}

        // Return the number of rows or columns in the matrix
        size_t nrow() const {return N;}
        size_t ncol() const {return M;}    
    };
    
    // Simple class storing sparse matrix as compressed sparse row
    template <typename T=double>
    class sparse_matrix_csr {
        size_t N,M;
        std::vector<sparse_vector<T>> data;
    public:
        typedef sparse_vector<T>* iterator;
        typedef const sparse_vector<T>* const_iterator;
        
        sparse_matrix_csr(size_t N, size_t M) : N(N), M(M), data(N,M) {}

      // Make a sparse matrix from a dense matrix with optional thresholding
      sparse_matrix_csr(const Tensor<T>& A, T tol = 0.0) 
      	: N(A.dim(0)), M(A.dim(1)), data(A.dim(0),A.dim(1))
      {
        for (size_t i=0; i<N; i++) {
            sparse_vector<T>& row = get_row(i);
            for (size_t j=0; j<M; j++) {
                T value = A(i,j);
                if (std::abs(value) > tol) row.push_back(j,value);
            }
            row.get_indices().shrink_to_fit();
            row.get_data().shrink_to_fit();
        }
      }
        
        sparse_matrix_csr(const sparse_matrix_coo<T>& A) : N(A.nrow()), M(A.ncol()), data(N,M) { 
            for (const auto& [ij,value] : A) {
                const auto [i,j] = ij;
                if (i<0 || i>=N) print("bad", i, j, value);
                get_row(i).push_back(j,value);
            }
            for (auto& row : data) row.sort_indices();
        }

      Tensor<T> dense() {
        Tensor<T> B(nrow(),ncol());
        size_t i=0;
        for (const auto& row: *this) {
 	  for (const auto& [j,value] : row) B(i,j) = value;
	  i++;
        }
        return B;
      }

      // Returns the sparse vector representing row i
      const sparse_vector<T>& get_row(size_t i) const {return data.at(i);}
      sparse_vector<T>& get_row(size_t i) {return data.at(i);}

        const_iterator begin() const {return &data[0];}
        const_iterator end() const {return (&data[0])+N;}
        iterator begin() {return &data[0];}
        iterator end() {return (&data[0])+N;}
        
        // Return the number of rows or columns in the matrix
        size_t nrow() const {return N;}
        size_t ncol() const {return M;}

        // Compress out zero/small elements
        void compress() {
            for (auto& row : *this) row.compress();
        }

        // Inplace renumber indices so that inew = imap(iold), jnew = jmap(jold)
        // For convenience we also sort the new vectors of indices/values, which also minimizes the storage
        // imap/jmap should be a permutation
        void renumber(const std::vector<size_t>& imap, const std::vector<size_t>& jmap) {
          for (auto& row : *this) row.renumber(jmap);
          std::vector<sparse_vector<T>> rows(N,M);
          for (size_t i=0; i<N; i++) rows[imap[i]] = std::move(data[i]);
          data = std::move(rows);
        }

        // Return a new matrix with elements in the specified slices.
        // Note a matrix is always returned and zero size slice will
        // throw.  Also, MADNESS slices are inclusive ranges whereas
        // Python slice have an exclusive endpoint.  Look at slice.h
        // and tensor.h for more info.
        sparse_matrix_csr<T> operator()(const Slice& rows, const Slice& cols) const {
            long rowstart = rows.start;  if (rowstart < 0) rowstart += N;
            long rowend = rows.end;  if (rowend < 0) rowend += N;
            long rowstep = rows.step;
            long colstart = cols.start;  if (colstart < 0) colstart += M;
            long colend = cols.end;  if (colend < 0) colend += M;
            long colstep = cols.step;

            if (rowstart == rowend && rowstep==0) throw "empty row slice not permitted";
            if (colstart == colend && colstep==0) throw "empty col slice not permitted";
            if (colstep!=1 || rowstep!=1) throw "non-unit step not supported yet";
            MADNESS_ASSERT(rowstart <= rowend);
            MADNESS_ASSERT(colstart <= colend);

            sparse_matrix_csr<T> A(rowend-rowstart+1, colend-colstart+1);
            for (long i=rowstart; i<=rowend; i++) {
                const auto& oldrow = get_row(i);
                auto& newrow = A.get_row(i-rowstart);
                for (auto [j,x] : oldrow) {
                    if (long(j)>=colstart && long(j)<=colend)
                        newrow.push_back(j-colstart,x);
                }
            }
            return A;
        }
    };

    // Simple class storing sparse matrix as compressed sparse row with tiling
    template <typename T=double>
    class sparse_matrix_csr_tiled {
    public:
        struct triple {
            size_t i;
            size_t j;
            T value;
        };
        struct lookup {
            size_t nu;
            size_t size;
            size_t offset;
        };
        typedef std::vector<triple> rowtype;

        size_t N,M,tilesize;
        std::vector<rowtype> data;
        std::vector<std::vector<lookup>> index;

    public:
        sparse_matrix_csr_tiled(size_t N, size_t M, size_t tilesize) : N(N), M(M), tilesize(tilesize), data(), index((M-1)/tilesize + 1) {}

        // Returns the sparse vector representing row i
        const rowtype& get_row(size_t i) const {return data.at(i);}
        rowtype& get_row(size_t i) {return data.at(i);}

        // Once data is loaded we need to complete the data structure
        void complete() {
            for (auto& row : data) {
                if (row.size() > 0) {
                    // 1) sort the rows so data in a tile is contiguous
                    std::sort(row.begin(), row.end(),
                              [](const triple& left, const triple& right) {return left.j < right.j;});

                    size_t itile = row[0].i / tilesize;
                    auto& rowindex = index[itile];

                    // 2) Build the index vector
                    lookup current = {size_t(-1),0,0};
                    size_t offset = 0;
                    for (const triple& entry : row) {
                        size_t jtile = entry.j/tilesize;
                        if (jtile != current.nu) {
                            if (current.size > 0) rowindex.push_back(current);
                            current = {jtile,0,offset};
                        }
                        current.size++;
                        offset++;
                    }
                    if (current.size > 0) rowindex.push_back(current);
                }
            }
        }
    };


    // Print sparse matrix CSR
    template <typename T>
    std::ostream& operator<<(std::ostream& s, const sparse_matrix_csr<T>& A) {
        using madness::operators::operator<<;        
        s<<"[";
	size_t i=0;
	for (auto& row : A) {
	  s << i << ":" << row;
	  if (i != A.nrow()-1) s << ",";
	  else s << "]";
	  i++;
	}
        return s;
    }


    // Read A test matrix from the configuration interaction (CI) code
    sparse_matrix_csr<double> load_ci_matrix(const std::string& filename) {
        std::fstream file(filename,std::ios::in);
        size_t N;
        file >> N;
	print("read ", N);
        //sparse_matrix_coo<double> A(N,N);
        sparse_matrix_csr<double> A(N,N);
        int i, j; // MUST BE SIGNED
        double value;
        size_t count=0;
        while (file >> i >> j >> value) {
	  if (i<0) {print("breaking with ", count); break;}
	  i--; j--; // fortran indexes from 1, c++ from 0
	  //MADNESS_ASSERT(i>=0 && i<N && j>=0 && j<N);
	  if (i<0 || j<0 || i>=long(N) || j>=long(N)) {
	    print("bad", count, ":", i, j, N, value);
	    throw "bad";
	  }
            //A[{i,j}] = value; A[{j,i}] = value; // read full square
	  A.get_row(i).push_back(j,value);
	  if (i != j) A.get_row(j).push_back(i,value);
	  count++;
        }
        print("load_ci_matrix:",N,count);
        
        //sparse_matrix_csr<double> B(A);
        //print(B.nrow());
        for (auto& row : A) row.sort_indices();
        return A;
    }


    // Read A test matrix from the configuration interaction (CI) code
    sparse_matrix_csr_tiled<double> load_ci_matrix_tiled(const std::string& filename, size_t tilesize) {
        std::fstream file(filename,std::ios::in);
        size_t N;
        file >> N;
	print("read ", N);

        sparse_matrix_csr_tiled<double> A(N,N,tilesize);
        int ii, jj; // MUST BE SIGNED
        double value;
        size_t count=0;
        while (file >> ii >> jj >> value) {
	  if (ii<0) {print("breaking with ", count); break;}
	  size_t i = ii-1;
          size_t j = jj-1;
          size_t itile = i/tilesize;
          size_t jtile = j/tilesize;
	  A.get_row(itile).push_back({i, j, value});
	  if (i != j) A.get_row(jtile).push_back({j, i, value});
	  count++;
        }
        print("load_ci_matrix:",N,count);
        
        A.complete();

        return A;
    }


    // // Converts a sparse_matrix_csr into a MADNESS tensor
    // template <typename T>
    // Tensor<T> sparse_matrix_to_tensor(const sparse_matrix_csr<T>& A) {
    //     Tensor<T> B(A.nrow(),A.ncol());
    //     size_t i=0;
    //     for (const auto& row: A) {
    //         for (const auto& [j,value] : row) B(i,j) = value;
    //         i++;
    //     }
    //     return B;
    // }

    // // Converts a MADNESS tensor into a sparse_matrix_csr
    // template <typename T>
    // sparse_matrix_csr<T> tensor_to_sparse_matrix(const Tensor<T>& A, T tol=0.0) {
    //     size_t N=A.dim(0), M=A.dim(1);
    //     sparse_matrix_csr<T> B(N,M);

    //     for (size_t i=0; i<N; i++) {
    //         sparse_vector<T>& row = B.get_row(i);
    //         for (size_t j=0; j<M; j++) {
    //             T value = A(i,j);
    //             if (std::abs(value) > tol) row.push_back(j,value);
    //         }
    //         row.get_indices().shrink_to_fit();
    //         row.get_data().shrink_to_fit();
    //     }
    //     return B;
    // }
    
    // Multiplies a sparse matrix (CSR) onto a vector (that can be either std::vector or madness::Tensor)
    template <typename T, typename vecT>
    vecT sparse_matrix_times_vector(const sparse_matrix_csr<T>& A, const vecT& v) {
        MADNESS_ASSERT(A.ncol() == size_t(v.size()));
        vecT Av(v.size());
        size_t i = 0;
        for (const auto& row: A) {
            T sum = 0.0;
            for (const auto& [j,aij] : row) sum += aij * v[j];
            Av[i++] = sum;
        }
        return Av;
    }
    
    // Power iteration with shift intended to crudely estimate extreme eigenvalues
    template <typename T>
    T simple_power_iteration(const sparse_matrix_csr<T>& A, T shift, double tol=0.001, bool ranguess=true) {
        Tensor<double> x(A.nrow());
        if (ranguess) x.fillrandom();
        else x[0] = 1.0;
        
        T eprev = 0.0;
        for (size_t iter=0; iter<1000; iter++) {
            x.scale(1.0/x.normf());
            Tensor<T> Ax = sparse_matrix_times_vector(A, x);
            T eval = x.trace(Ax);
            //print(iter,eval);
            if ((iter%10)==1) {
                if (std::abs((eval-eprev)/eval) < tol) return eval;
                eprev = eval;
            }
            x = Ax - shift*x;
        }
        return eprev;
    }
    
    // More efficient power iteration for lowest eigenvalue using Chebyshev nodes
    template <typename T>
    T power_iteration(const sparse_matrix_csr<T>& A, T e0, T e1, T ehi, double tol=0.0001) {
        Tensor<double> x(A.nrow());
        x[0] = 1.0;
        
        T xlo = (e1-e0)/(ehi-e1);
        T xhi = 1.0;
        //T pts[] = {0.0};
        //T pts[] = {std::sqrt(2.0)/2.0,-std::sqrt(2.0)/2.0};
        //T pts[] = {std::sqrt(3.0)/2.0,0.0,-std::sqrt(3.0)/2.0};
        T pts[] = {-std::sqrt(2.0+std::sqrt(2.0))/2.0,-std::sqrt(2.0-std::sqrt(2.0))/2.0,std::sqrt(2.0-std::sqrt(2.0))/2.0,std::sqrt(2.0+std::sqrt(2.0))/2.0};
        T eprev = 0.0;
        for (size_t iter=0; iter<1000; iter++) {
            T eval;
            for (T pt : pts) {
                T omega = 1.0/((0.5*(pt+1.0)*(xhi-xlo)+ xlo));
                //omega = 1.0;
                x.scale(1.0/x.normf());
                Tensor<T> Ax = sparse_matrix_times_vector(A, x);
                eval = x.trace(Ax);
                x = x + omega*((-1.0/Ax.normf())*Ax - x);
            }
            print(iter,eval);
            if ((iter%10)==1) {
                if (std::abs((eval-eprev)/eval) < tol) return eval;
                eprev = eval;
            }
        }
        return eprev;
    }
  
  template <typename T>
  class vector_sampler {
    std::vector<size_t> indices;
    std::vector<T> P;
    std::vector<T> signs;
    size_t N;
    T X;
    
  public:
    vector_sampler() {} // need this to be able to store in a vector

    vector_sampler(const std::vector<T>& x, const std::vector<size_t>& indices = std::vector<size_t>())
      : indices(indices)
      , P(x.size())
      , signs(x.size())
      , N(x.size())
    {
      if (indices.size() > 0) MADNESS_ASSERT(indices.size() == N);
      T sum = 0.0;
      for (size_t i=0; i<N; i++) {
	sum += std::abs(x[i]);
	P[i] = sum;
	signs[i] = x[i]>=0.0 ? 1.0 : -1.0;
      }
      X = sum;
      sum = 1.0/sum;
      for (size_t i=0; i<N; i++) P[i] *= sum;
    }

    vector_sampler(const sparse_vector<T>& x) 
      : vector_sampler(x.get_data(), x.get_indices())
    {}
    
    std::pair<size_t,T> sample() const {
      T xi = RandomValue<T>();
      size_t a=0, b=N-1;
      while ((b-a) > 1) { // binary search for i
	size_t m = (a+b)/2;
	if (P[m]>=xi) b = m;
	else a = m;
      }
      size_t i;
      if (P[a]>=xi) i = a;
      else i = b;
      T isgn = signs[i];
      if (indices.size()) i = indices[i];
      return {i,isgn};
    }
    
    T norm() const {
      return X;
    }
  };

  template <typename T>
  class sparse_matrix_sampler {
    std::vector<vector_sampler<T>> data;
  public:
    sparse_matrix_sampler(const sparse_matrix_csr<T>& a) {
      data.reserve(a.nrow());
      for (const auto& row: a) data.push_back(vector_sampler<T>(row));
    }

    // sample an element from row i of the matrix
    std::pair<size_t, T> sample(size_t i) const {return data[i].sample();}

    // return the norm (sum of absolute values) of row i of the matrix
    T norm(size_t i) const {return data[i].norm();}
  };
}

// Computes the Shur energy expression using dense matrix algebra
template <typename T>
T ShurDense(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
    }

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j);
        }
    }

    // Tensor<double> evals,evecs;
    // syev(A, evecs, evals);
    // print("dense Aevals",evals[0],evals[A.dim(0)-1]);

    // Solve the linear equation
    Tensor<T> y;
    gesv(A, x, y);

    // // Project x and y onto Aevecs out of curiosity
    // Tensor<T> xA = inner(x,evecs);
    // Tensor<T> yA = inner(y,evecs);
    // print("");
    // print("xA",xA);
    // print("yA",yA);
    // std::exit(0);

    // Compute new energy
    return E0 - x.trace(y);
}

// Computes the Shur energy expression using (mostly) sparse matrix algebra
// (we need to extend the sparse matrix class a little to avoid the dense operations)
template <typename T>
T ShurSparse(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    Tensor<double> s(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
        if (x(i)<0) {
            s(i) = -1.0;
            x(i) = -x(i);
        }
        else {
            s(i) = 1.0;
        }
    }

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j)*s(i)*s(j);
        }
    }

    // Convert to sparse
    sparse_matrix_csr<T> Asp(A,1e-6);

    // Iteratively solve the linear equation --- for now just bute force iteration --- really need to estimate the range of the spectrum using power method
    // then switch to faster iteration to solve.

    Tensor<T> y(N-1);
    double omega = 0.773;
    for (int iter=0; iter<1000; iter++) {
        Tensor<T> r = x - sparse_matrix_times_vector(Asp,y);
        double err = r.normf();
        //print(iter,err);
        y += omega*r;
        if (err < 1e-6) break;
    }

    // Compute new energy
    return E0 - x.trace(y);
}

// Computes the Shur energy expression using fully sparse matrix algebra
template <typename T>
T ShurSparse(const sparse_matrix_csr<T>& H, T E) {
    size_t N = H.nrow();

    // Reference state energy
    T E0 = H.get_row(0).get_data()[0];

    // Extract Shur complement matrix and vector
    sparse_matrix_csr<T> A = H(Slice(1,-1),Slice(1,-1));
    Tensor<T> x = copy(H.get_row(0).dense()(Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        auto [found, d] = A.get_row(i).find(i);
        if (!found) throw "diagonal element is missing?";
        d = 1.0/std::sqrt(d-E);
        D(i) = d;
        x(i) *= d;
    }
    //print("xxxxxxxx", x);

    // Scale A and shift its diagonal (which becomes 1)
    for (size_t i=0; i<N-1; i++) {
        auto& rowi = A.get_row(i);
        auto& ind = rowi.get_indices();
        auto& dat = rowi.get_data();
        for (size_t jj=0; jj<ind.size(); jj++) {
            size_t j = ind[jj];
            if (i == j) {
                dat[jj] = 1.0;
            }
            else {
                dat[jj] *= D(i)*D(j);
            }
        }
    }
    //print("AAAAAAAAAA", A);

    // Iteratively solve the linear equation --- for now just bute force iteration --- really need to estimate the range of the spectrum using power method
    // then switch to faster iteration to solve.

    Tensor<T> y(N-1);
    double omega = 0.773;
    for (size_t iter=0; iter<1000; iter++) {
        Tensor<T> r = x - sparse_matrix_times_vector(A,y);
        double err = r.normf();
        //print(iter,err);
        y += omega*r;
        if (err < 1e-6) break;
    }

    // Compute new energy
    return E0 - x.trace(y);
}


// Computes the Shur energy expression using (mostly) sparse matrix algebra and Monte Carlo
// (we need to extend the sparse matrix class a little to avoid the dense operations)
template <typename T>
T ShurSparseMC(const Tensor<T>& H, T E) {
    size_t N = H.dim(0);

    // Reference state energy
    T E0 = H(0,0);

    // Extract Shur complement matrix and vector
    Tensor<T> A = copy(H(Slice(1,-1),Slice(1,-1)));
    Tensor<T> x = copy(H(0,Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<T> D(N-1);
    Tensor<T> s(N-1);
    for (size_t i=0; i<N-1; i++) {
        T d = A(i,i) - E;
        A(i,i) = d;
        d = 1.0/std::sqrt(d);
        D(i) = d;
        x(i) *= d;
        // if (x(i)<0) {  // experiment with using sign information of x
        //     s(i) = -1.0;
        //     x(i) = -x(i);
        // }
        // else {
           s(i) = 1.0;
        // }
    }

    // Find highest index non-zero element of x so we can optimize truncation etc.
    int ximax = N-2;
    while (x[ximax] == 0.0) ximax--;
    print("ximax", ximax);

    // Scale A
    for (size_t i=0; i<N-1; i++) {
        for (size_t j=0; j<N-1; j++) {
            A(i,j) *= D(i)*D(j);
        }
    }

    // Tensor<double> evals,evecs;
    // syev(A, evecs, evals);
    // print("Aevals",evals[0],evals[1],evals[2],evals[3],evals[A.dim(0)-1]);


    // For MC algorithm zero the diagonal entries of A which we know are all 1
    for (size_t i=0; i<N-1; i++) {
        if (std::abs(A(i,i)-1.0) > 1e-12) print("bad diag", i, A(i,i));
        A(i,i) = 0.0;
    }

    // // Examine spectra of +ve and -ve parts of A
    // Tensor<T> Apos=copy(A), Aneg=copy(A);
    // for (size_t i=0; i<N-1; i++) {
    //     for (size_t j=0; j<N-1; j++) {
    //         if (A(i,j)<0) {
    //             Apos(i,j) = 0.0;
    //         }
    //         else {
    //             Aneg(i,j) = 0.0;
    //         }
    //     }
    // }
    // Tensor<double> evals,evecs;
    // {syev(Aneg, evecs, evals); print("Aneg",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos, evecs, evals); print("Apos",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos+Aneg, evecs, evals); print("Apos+Aneg",evals[0],evals[A.dim(0)-1]);}
    // {syev(Apos-Aneg, evecs, evals); print("Apos-Aneg",evals[0],evals[A.dim(0)-1]);}
    
    // Convert to sparse
    sparse_matrix_csr<T> Asp(A,1e-6);

    double omega = 0.1; //0.773;
    Tensor<T> y = omega*copy(x); // can generate a better initial guess

    double sum = 0.0;
    double sumsq = 0.0;
    double sum2 = 0.0;
    double sumsq2 = 0.0;
    size_t count = 0;

    const double cut = 0.01;
    const double pkill = 0.75;

    {
        T eval0 = simple_power_iteration(Asp, -1.0);
        T eval1 = simple_power_iteration(Asp, eval0);
        print("Aspevals", eval0, eval1);
    }

    for (size_t iter=0; iter<3500; iter++) {
        Tensor<T> ynew(N-1);
        sparse_matrix_coo<double> Asample(N,N);
        for (size_t i=0; i<N-1; i++) {
            if (y[i]) {
                int nrep = 3;
                for (int rep=0; rep<nrep; rep++) {
                    // Randomly sample non-zero element in row A[i,*]
                    auto row = Asp.get_row(i);
                    size_t n = row.get_indices().size();
                    size_t j = n*RandomValue<double>();
                    size_t jj = row.get_indices()[j];
                    double aji = s(i)*s(jj)*row.get_data()[j];
                    Asample[{i,jj}] = aji;
                    Asample[{jj,i}] = aji;
                    double gji = double(nrep)/double(n);
                    ynew[jj] -= 0.5*omega*aji*y[i] / gji; 
                    ynew[i]  -= 0.5*omega*aji*y[jj] / gji; // Enforcing symmetry is important for stability (part of detailed balance?)
                }
                //for (auto [j,aji] : Asp.get_row(i)) ynew[j] -= omega*aji*y[i];
            }
            //ynew[i] += (1.0-omega)*y[i] + omega*x[i];
        }
        sparse_matrix_csr<double> B(Asample);
        T eval0 = simple_power_iteration(B, -1.0);
        T eval1 = simple_power_iteration(B, eval0);
        print("evals", eval0, eval1);

        //ynew -= omega*sparse_matrix_times_vector(Asp,y);
        double xzy = x.trace(x) + x.trace(ynew)/(omega); // using (1-z)^-1 = 1 + z (1-z)^-1
        ynew += (1.0-omega)*y + omega*x;
        double ecor = ynew.trace(x);
        //print(iter,ecor,(ynew-y).normf());
        y = ynew;
        if (iter >500) {
            sum += ecor;
            sumsq += ecor*ecor;
            sum2 += xzy;
            sumsq2 += xzy*xzy;
            count += 1;
        }

        size_t nkill = 0.0;
        for (size_t i=0; i<N-1; i++) {
            if (std::abs(y[i]) < cut) {
                if (RandomValue<double>() > pkill) {
                    y[i] *= 1.0/(1.0-pkill);
                }
                else {
                    nkill ++;
                    y[i] = 0.0;
                }
            }
        }
        //print("nkill", nkill);
    }
    sum /= count;
    sumsq /= count;
    double err = std::sqrt((sumsq - sum*sum)/count);
    sum2 /= count;
    sumsq2 /= count;
    double err2 = std::sqrt((sumsq2 - sum2*sum2)/count);
    print("Correlation energy ",-sum,"stderr",err,"total energy", E0 - sum, "F(E)", E0-sum-E);
    print("Correlation energy2",-sum2,"stderr",err2,"total energy", E0 - sum2, "F(E)", E0-sum2-E);
    //print(y);

    // Compute new energy
    return E0 - sum;
}

// Computes the Shur energy expression using fully sparse matrix algebra and Monte Carlo
template <typename T>
T ShurSparseMC(const sparse_matrix_csr<T>& H, T E) {
    size_t N = H.nrow();

    // Reference state energy
    T E0 = H.get_row(0).get_data()[0];

    // Extract Shur complement matrix and vector
    sparse_matrix_csr<T> A = H(Slice(1,-1),Slice(1,-1));
    Tensor<T> x = copy(H.get_row(0).dense()(Slice(1,-1)));

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        auto [found, d] = A.get_row(i).find(i);
        if (!found) throw "diagonal element is missing?";
        d = 1.0/std::sqrt(d-E);
        D(i) = d;
        x(i) *= d;
    }

    // Scale A and shift its diagonal (which becomes 1 and here we set to zero since it is handled exactly)
    for (size_t i=0; i<N-1; i++) {
        auto& rowi = A.get_row(i);
        auto& ind = rowi.get_indices();
        auto& dat = rowi.get_data();
        for (size_t jj=0; jj<ind.size(); jj++) {
            size_t j = ind[jj];
            if (i == j) {
                dat[jj] = 0.0;
            }
            else {
                dat[jj] *= D(i)*D(j);
            }
        }
    }

    // Compress out the zeros (diagonal elements)
    A.compress();

    double omega = 0.1; //0.773;
    Tensor<T> y = copy(x); // can generate a better initial guess

    double sum = 0.0;
    double sumsq = 0.0;
    double sum2 = 0.0;
    double sumsq2 = 0.0;
    size_t count = 0;

    const double cut = 0.001;
    const double pkill = 0.75;

    {
        T eval0 = simple_power_iteration(A, -1.0);
        T eval1 = simple_power_iteration(A, eval0);
        print("Aevals", eval0, eval1);
    }

    for (size_t iter=0; iter<3800; iter++) {
        Tensor<T> ynew(N-1);
        sparse_matrix_coo<double> Asample(N,N);
        for (size_t i=0; i<N-1; i++) {
            if (y[i]) {
                int nrep = 10;
                for (int rep=0; rep<nrep; rep++) {
                    // Randomly sample non-zero element in row A[i,*]
                    auto row = A.get_row(i);
                    size_t n = row.get_indices().size();
                    size_t j = n*RandomValue<double>();
                    size_t jj = row.get_indices()[j];
                    double aji = row.get_data()[j];
                    Asample[{i,jj}] = aji;
                    Asample[{jj,i}] = aji;
                    double gji = double(nrep)/double(n);
                    ynew[jj] -= 0.5*omega*aji*y[i] / gji; 
                    ynew[i]  -= 0.5*omega*aji*y[jj] / gji; // Enforcing symmetry is important for stability (part of detailed balance?)
                }
                //for (auto [j,aji] : A.get_row(i)) ynew[j] -= omega*aji*y[i];
            }
            //ynew[i] += (1.0-omega)*y[i] + omega*x[i];
        }
        // sparse_matrix_csr<double> B(Asample);
        // T eval0 = simple_power_iteration(B, -1.0);
        // T eval1 = simple_power_iteration(B, eval0);
        // print("evals", eval0, eval1);

        //ynew -= omega*sparse_matrix_times_vector(A,y);
        double xzy = x.trace(x) + x.trace(ynew)/(omega); // using (1-z)^-1 = 1 + z (1-z)^-1
        ynew += (1.0-omega)*y + omega*x;
        double ecor = ynew.trace(x);
        y = ynew;
        if (iter > 800) {
            sum += ecor;
            sumsq += ecor*ecor;
            sum2 += xzy;
            sumsq2 += xzy*xzy;
            count += 1;
        }

        size_t nnz = 0;
        size_t nkill = 0;
        for (size_t i=0; i<N-1; i++) {
            if (std::abs(y[i]) < cut) {
                if (RandomValue<double>() > pkill) {
                    y[i] *= 1.0/(1.0-pkill);
                }
                else {
                    nkill ++;
                    y[i] = 0.0;
                }
            }
            if (std::abs(y[i])) nnz++;
        }
        if ((iter%10) == 1) print(iter,ecor,xzy,(ynew-y).normf(),nkill,nnz,nnz/double(N-1));
    }
    sum /= count;
    sumsq /= count;
    double err = std::sqrt((sumsq - sum*sum)/count);
    sum2 /= count;
    sumsq2 /= count;
    double err2 = std::sqrt((sumsq2 - sum2*sum2)/count);
    print("Correlation energy ",-sum,"stderr",err,"total energy", E0 - sum, "F(E)", E0-sum-E);
    print("Correlation energy2",-sum2,"stderr",err2,"total energy", E0 - sum2, "F(E)", E0-sum2-E);
    //print(y);

    // Compute new energy
    return E0 - sum;
}

/*

  y = y + (x - A*y)

  y = A^-1 x

    = (1 - Z)^-1 x = (1 + Z + Z^2 + ...) x

 */

// Computes the Shur energy expression using fully sparse matrix
// algebra and Monte Carlo This routine has code that demonstrates the
// sign problem by computing with combinations of A+ and A-
template <typename T>
T ShurSparseMCX(const sparse_matrix_csr<T>& H, T E) {
    size_t N = H.nrow();

    print("HHHH", N);
    //print(H);

    // Reference state energy
    T E0 = H.get_row(0).get_data()[0];

    // Extract Shur complement matrix and vector
    sparse_matrix_csr<T> A = H(Slice(1,-1),Slice(1,-1));
    Tensor<T> x = copy(H.get_row(0).dense()(Slice(1,-1)));

    // For convenience of testing algorithms, reorder the rows/columns
    // so that the non-zero elements of x come first in the vector.
    std::vector<size_t> map(N-1);
    size_t countx = 0;
    for (size_t i=0; i<N-1; i++) if (std::abs(x[i]) > 0) countx++;
    print("countx", countx);
    size_t inew=0, jnew=countx;
    for (size_t i=0; i<N-1; i++) {
      if (std::abs(x[i]) > 0) map[i] = inew++;
      else map[i] = jnew++;
    }
    MADNESS_ASSERT(inew==countx && jnew==N-1);
    //print(map);
    A.renumber(map,map);
    {
      Tensor<T> y = copy(x);
      for (size_t i=0; i<N-1; i++) x[map[i]] = y[i];
    }

    T DELTA = E; // (A-E) = (D-DELTA) + (A-D-E+DELTA)

    T lambda = 1.0;

    // Shift and extract diagonal, and scale x
    Tensor<double> D(N-1);
    for (size_t i=0; i<N-1; i++) {
        const auto [found, d] = A.get_row(i).find(i);
        if (!found) throw "diagonal element is missing?";
        T rsqd = 1.0/std::sqrt(d-DELTA);
        D(i) = rsqd;
        x(i) *= rsqd;
    }

    double omega = 0.733;

    // Scale and shift A to make Zomega = I - omega * A
    for (size_t i=0; i<N-1; i++) {
        auto& rowi = A.get_row(i);
        auto& ind = rowi.get_indices();
        auto& dat = rowi.get_data();
        for (size_t jj=0; jj<ind.size(); jj++) {
            size_t j = ind[jj];
	    if (i == j) dat[jj] = DELTA - E;
	    dat[jj] *= -lambda*omega*D(i)*D(j);
	    if (i == j) dat[jj] += lambda*(1.0 - omega);
        }
    }

    // Compress out any zeros
    A.compress();

    // {
    //   // Compute the exact solution
    //   // print("Z", A);
    //   Tensor<T> tA = A.dense();
    //   tA *= -1.0;
    //   for (size_t i=0; i<N-1; i++) tA(i,i) += 1.0;
    //   Tensor<T> ty;
    //   gesv(tA, x, ty);
    //   print("exact soln", ty.trace(x)*omega, E0 - ty.trace(x)*omega);
    // }
    
    // // Test with the absolute value of the matrix elements
    // print("USING ABSOLUTE VALUE OF THE MATRIX");
    // for (auto& row: A) {
    //   for (auto [j,value] : row) {
    // 	value = std::abs(value);//value<0 ? value : 0.0;
    //   }
    // }

    {
        T eval0 = simple_power_iteration(A, 0.0);
        T eval1 = simple_power_iteration(A, eval0);
        print("Zevals", eval0, eval1);
    }

    // if (countx != N-1) {
    //     sparse_matrix_csr<T> AA = A(Slice(countx,-1),Slice(countx,-1));
    //     T eval0 = simple_power_iteration(AA, 0.0);
    //     T eval1 = simple_power_iteration(AA, eval0);
    //     print("TZevals", eval0, eval1);
    // }

    // {
    //   Tensor<double> evals,evecs;
    //   syev(A.dense(), evecs, evals);
    //   print("dense Zevals",evals);
    // }

    size_t nstep = 40;
    {
      // test with simple von Neumann iteration
      Tensor<T> y = copy(x);
      Tensor<T> xk = copy(x);
      for (size_t iter=0; iter<nstep; iter++) {
	print("vN", iter, x.trace(xk)*omega, y.trace(x)*omega);
	xk = sparse_matrix_times_vector(A,xk);
	y += xk;
      }
    }

    {
      // test with iterative solver iteration
      Tensor<T> y = copy(x);
      for (size_t iter=0; iter<nstep; iter++) {
	y = x + sparse_matrix_times_vector(A,y);
	print("it", iter, y.trace(x)*omega);
      }
    }

    // Now do the MC version of von Neumann iteration
    size_t nsample = 1000000;
    std::vector<double> sum(nstep,0.0);
    std::vector<double> sumsq(nstep,0.0);
    std::vector<double> stderr(nstep,0.0);

    vector_sampler<T> xsampler(x);
    sparse_matrix_sampler<T> Zsampler(A);

    //for (size_t i=0; i<N-1; i++) print(i,Zsampler.norm(i));

    for (size_t sample=0; sample<nsample; sample++) {
      auto [i,wi] = xsampler.sample();
      wi *= xsampler.norm();
      for (size_t step=0; step<nstep; step++) {
	sum[step] += x[i]*wi;
        sumsq[step] += (x[i]*wi)*(x[i]*wi);
	auto [j,s] = Zsampler.sample(i);
	wi *= s*Zsampler.norm(i);
	i = j;
      }
    }

    for (auto& s : sum) s *= omega/nsample;
    for (auto& s : sumsq) s *= omega*omega/nsample;
    for (size_t k=0; k<nstep; k++) stderr[k] = std::sqrt((sumsq[k] - sum[k]*sum[k])/nsample);
    for (size_t k=0; k<nstep; k++) print(k, sum[k], stderr[k]);

    T s = 0.0, ss = 0.0, lamk = 1.0;
    for (auto v : sum) {
      s += v; 
      ss += lamk*v; 
      print(v, lamk, s, ss);
      lamk /= lambda;
    }

    print("x.x", x.trace(x)*omega);
    print("sss", s, ss, E0 - ss );

    return 1.0;
}

void test_sampler() {
  size_t N = 1000;
  std::vector<double> x(N), y(N);
  std::vector<size_t> indices(N);
  for (size_t i=0; i<N; i++) {
      // size_t j = RandomValue<double>()*i;
      // indices[i] = indices[j];
      // indices[j] = i;
      indices[i] = N-i-1;
      x[i] = RandomValue<double>()-0.8;
      y[i] = RandomValue<double>()-0.3;
    }
    double exact = 0.0, exact2 = 0.0;
    for (size_t i=0; i<N; i++) {
      exact += x[i]*y[i];
      exact2 += x[i]*y[indices[i]];
    }
    vector_sampler<double> X(x), XI(x,indices);
    double sum = 0.0, sum2 = 0.0;

    size_t nsample = 10000000;
    for (size_t sample=0; sample<nsample; sample++) {
      {
	auto [i,s] = X.sample();
	sum += y[i]*s;
      }
      {
	auto [i,s] = XI.sample();
	sum2 += y[i]*s;
      }
    }
    print(X.norm()*sum/nsample, exact);
    print(XI.norm()*sum2/nsample, exact2);
  }
    

int main() {
    std::cout.precision(10);    
    sparse_matrix_csr A = load_ci_matrix("h2o-4037.txt");

    auto AA = load_ci_matrix_tiled("h2o-4037.txt", 2);

    size_t Nuse = A.nrow()-1;
    //size_t Nuse = 500;
    A = A(Slice(0,Nuse),Slice(0,Nuse)); 
    print(A.nrow(),A.nrow());

    //sparse_matrix_csr A = load_ci_matrix("sparse.txt");
    //sparse_matrix_csr A = load_ci_matrix("hooh-15252.txt");
    //sparse_matrix_csr A = load_ci_matrix("hooh-44034.txt");
    //print(A.get_row(0));
    
    //Tensor<double> B = A.dense();
    // Tensor<double> evals,evecs;
    // syev(B, evecs, evals);
    // print("Hevals",evals[0],evals[B.dim(0)-1]);

    // Tensor<double> x(A.nrow());
    // x.fillrandom();
    // Tensor<double> Bx = inner(B,x);
    // Tensor<double> Ax = sparse_matrix_times_vector(A, x);
    // print("error in sparse mxv", (Ax - Bx).normf());

    // double elo = simple_power_iteration(A,-65.0,1e-4,false);
    // double ehi = simple_power_iteration(A,elo,1e-4);
    // print(elo,ehi);
    // elo -= 0.1;
    // ehi += 0.1;
    // print(elo,ehi);

    // power_iteration(A, elo, elo+0.2, ehi, 1e-7);

    //    print("sparse1 ", -76.1, ShurSparse(B,-76.1));
    double E0 = -76.0; //-150.807 - 0.4;
    print("sparse2 ", E0, ShurSparse(A,E0));
    //    //print(" dense", -76.1,ShurDense(B,-76.1));

    //test_sampler();

    ShurSparseMCX(A,E0);

    
    
    // // Use secant method to find the root ... better would be to use all points computed and fit a low-order polyn to manage noise
    // //double E0 = -76.1; // Initial guess
    // double F0 = ShurSparse(A,E0) - E0;
    // double E1 = E0 + 0.1;
    // while (std::abs(E0-E1)>1e-4) {
    //     double F1 = ShurSparse(A,E1) - E1;
    //     print(E0,F0,E1,F1);
    //     double E2 = E1 - F1*(E1-E0)/(F1-F0);
    //     // Save the old point closest to the new point
    //     if (std::abs(E2-E1) < std::abs(E2-E0)) {
    //         E0 = E1; F0 = F1; E1 = E2;
    //     }
    //     else {
    //         E1 = E2;
    //     }
    // }
    
    return 0;
}
