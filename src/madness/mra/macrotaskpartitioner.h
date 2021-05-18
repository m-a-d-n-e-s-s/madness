//
// Created by Florian Bischoff on 12/10/20.
//

#ifndef MADNESS_MACROTASKPARTITIONER_H
#define MADNESS_MACROTASKPARTITIONER_H

#include<vector>
#include<list>
#include<string>
#include<iomanip>
#include<sstream>
#include<madness/world/madness_exception.h>


namespace madness {
template<typename ... Ts>
constexpr auto decay_types(std::tuple<Ts...> const &)
-> std::tuple<std::remove_cv_t<std::remove_reference_t<Ts>>...>;

template<typename T>
using decay_tuple = decltype(decay_types(std::declval<T>()));

template<typename>
struct is_madness_function_vector : std::false_type {
};

template<typename T, std::size_t NDIM>
struct is_madness_function_vector<std::vector<typename madness::Function<T, NDIM>>> : std::true_type {
};

template<typename tupleT, std::size_t I>
constexpr std::size_t get_index_of_first_vector_argument() {

    typedef decay_tuple <tupleT> argtupleT;   // removes const, &, etc

    if constexpr(I >= std::tuple_size_v<tupleT>) {
        // Last case, if nothing is left to iterate, then exit the function
//        MADNESS_EXCEPTION("there is no madness function vector argument in the list, cannot partition the tasks", 1);
        return I;
    } else {
        using typeT = typename std::tuple_element<I, argtupleT>::type;       // use decay types for determining a vector
        if constexpr (is_madness_function_vector<typeT>::value) {
            return I;
        } else {
            // Going for next element.
            return get_index_of_first_vector_argument<tupleT,I+1>();
        }
    }
}

template<typename tupleT, std::size_t I>
constexpr std::size_t get_index_of_second_vector_argument() {
    constexpr std::size_t index0=get_index_of_first_vector_argument<tupleT,0>();
    return get_index_of_first_vector_argument<tupleT,index0+1>();
}

class Batch_1D {
    friend class MacroTaskPartitioner;

public:
    std::size_t begin=0, end=-1; ///< first and first past last index [begin,end)

    Batch_1D() {}
    Batch_1D(const Slice& s) : begin(s.start), end(s.end) {
        MADNESS_CHECK(s.step==1);
    }
    Batch_1D(const std::size_t& begin, const std::size_t& end) : begin(begin), end(end) {}

    bool operator==(const Batch_1D& other) const {
        return (end==other.end && begin==other.begin);
    }

    std::size_t size() const {
        return is_full_size() ? SIZE_MAX : end - begin;
    }

    bool is_full_size() const {
        return (begin==0 and end==-1);
    }

    /// select the relevant vector elements from the argument tuple
    template<typename tupleT, std::size_t appearance>
    tupleT copy_batch(const tupleT& arg) const {

        constexpr std::size_t index = (appearance==0) ?  get_index_of_first_vector_argument<tupleT,0>() :
                                      get_index_of_second_vector_argument<tupleT,0>();
        if constexpr (index>=std::tuple_size<tupleT>::value) {
            MADNESS_CHECK(is_full_size());  // better: not set by user..
            return arg;
        } else {
            auto v = std::get<index>(arg);
            decltype(v) v_batch;
            std::size_t end1 = (end > 0) ? end : v.end() - v.begin();
            std::copy(v.begin() + begin, v.begin() + end1, std::back_inserter(v_batch));
            tupleT batched_arg = arg;
            std::get<index>(batched_arg) = v_batch;
            return batched_arg;
        }
    }

    template<typename vecT>
    vecT insert_batch(vecT v, const vecT& v_batch) const {
        MADNESS_CHECK(v_batch.size()==this->size() or this->is_full_size());
        std::copy(v_batch.begin(), v_batch.end(), v.begin()+begin);
        return v;
    }

    friend std::ostream& operator<<(std::ostream& os, const Batch_1D& batch) {
        std::stringstream ss;
        if (batch.is_full_size()) ss << "[  ----  )";
        else ss << "[" << std::setw(3) << batch.begin << ", " << std::setw(3) << batch.end << ")";
        os << ss.str();
        return os;
    }
};

/// a batch consists of a 2D-input batch and a 1D-output batch: K-batch <- (I-batch, J-batch)
class Batch {
public:
    friend class MacroTaskPartitioner;
    Batch_1D result;
    std::vector<Batch_1D> input;

    Batch() {}
    Batch(const Batch& other) {
        *this=other;
    }
    Batch& operator=(const Batch& other) {
        if (this==&other) return *this;
        result=other.result;
        input=other.input;
        return *this;
    }

    Batch(Batch_1D input1,  Batch_1D result) : input{input1}, result(result) {}
    Batch(Batch_1D input1, Batch_1D input2, Batch_1D result)
            : input{input1,input2}, result(result) {}

    /// select the relevant vector elements from the argument tuple
    template<typename tupleT>
    tupleT copy_input_batch(const tupleT& arg) const {
        if (input.size()==0) return arg;
        tupleT arg1=input[0].template copy_batch<tupleT,0>(arg);
        if (input.size()>1) arg1=input[1].template copy_batch<tupleT,1>(arg1);
        MADNESS_CHECK(input.size()<=2);
        return arg1;
    }

    template<typename vecT>
    vecT insert_result_batch(vecT v, const vecT& v_batch) const {
        return result.template insert_batch(v,v_batch);
    }

    friend std::ostream& operator<<(std::ostream& os, const Batch& batch) {
        os << batch.result<< " <-- [ ";
        if (batch.input.size()>0)  os << batch.input[0];
        if (batch.input.size()>1) os << ", " << batch.input[1];
        os << " ]";
        return os;
    }
};

class MacroTaskPartitioner {
    friend class Batch;

public:
    typedef std::list<Batch> partitionT;
    std::size_t min_batch_size=5;
    std::size_t max_batch_size = 10;
    std::size_t nsubworld=1;
    std::string policy = "guided";
    std::size_t dimension = 1;

    MacroTaskPartitioner() {
    }

    MacroTaskPartitioner& set_nsubworld(const long& n) {
        nsubworld=n;
        return *this;
    }
    MacroTaskPartitioner& set_policy(const std::string& n) {
        policy=n;
        return *this;
    }
    MacroTaskPartitioner& set_dimension(const std::size_t& n) {
        dimension=n;
        return *this;
    }
    MacroTaskPartitioner& set_min_batch_size(const long& n) {
        min_batch_size=n;
        return *this;
    }
    MacroTaskPartitioner& set_max_batch_size(const long& n) {
        max_batch_size=n;
        return *this;
    }

    template<typename tupleT>
    partitionT partition_tasks(const tupleT& argtuple) const {

        constexpr std::size_t I1 = get_index_of_first_vector_argument<tupleT, 0>();
        constexpr std::size_t I2 = get_index_of_second_vector_argument<tupleT, 1>();
        if constexpr (I2 < std::tuple_size_v<tupleT>) {     // found at least 2 vectors of madness functions
            if (dimension == 1) {
                constexpr std::size_t I = get_index_of_first_vector_argument<tupleT, 0>();
                auto v = std::get<I>(argtuple);
                return do_1d_partition(v.size(), "guided");
            } else if (dimension == 2) {
                constexpr std::size_t I = get_index_of_first_vector_argument<tupleT, 0>();
                auto v = std::get<I>(argtuple);
                constexpr std::size_t I2 = get_index_of_second_vector_argument<tupleT, 0>();
                auto v2 = std::get<I2>(argtuple);
                return do_2d_partition(v.size(), v2.size(), "guided");
            }
        }
        else if constexpr (I1 < std::tuple_size_v<tupleT>) { // found 1 vector
            MADNESS_CHECK(dimension==1);
            auto v = std::get<I1>(argtuple);
            return do_1d_partition(v.size(),"guided");
        }
        else {
            std::string msg="confused partitioning dimension: "+std::to_string(dimension) +" typeid " + typeid(tupleT).name();
            MADNESS_EXCEPTION(msg.c_str(),1);
        }
        return partitionT();
    }

    partitionT do_1d_partition(const std::size_t vsize, const std::string policy) const {
        partitionT result;
        if (policy == "guided") {
            std::size_t begin = 0;
            std::size_t end = 0;
            while (end < vsize) {
                end += std::min(max_batch_size, std::max(min_batch_size, ((vsize - end) / nsubworld)));
                end = std::min(end, vsize);
                result.push_back(Batch(Batch_1D(begin, end),Batch_1D(begin,end)));
                begin = end;
            }
        } else {
            std::string msg = "unknown partitioning policy: " + policy;
            MADNESS_EXCEPTION(msg.c_str(), 1);
        }
        return result;
    }

    partitionT do_2d_partition(const std::size_t vsize, const std::size_t v2size, const std::string policy) const {
        partitionT partition1=do_1d_partition(vsize,policy);
        partitionT partition2=do_1d_partition(v2size,policy);
        partitionT result;
        for (auto p1 : partition1) {
            for (auto p2 : partition2) {
                result.push_back(Batch(p1.input[0],p2.input[0],p1.result));
            }
        }
        return result;

    }
};

}

#endif //MADNESS_MACROTASKPARTITIONER_H
