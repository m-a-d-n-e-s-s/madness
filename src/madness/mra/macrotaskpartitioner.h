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

template<typename Q> struct is_vector : std::false_type { };
template<typename Q> struct is_vector<std::vector<Q>> : std::true_type { };

/// given a tuple return the index of the first argument that is a vector of Function<T,NDIM>
template<typename tupleT, std::size_t I>
constexpr std::size_t get_index_of_first_vector_argument() {

    typedef decay_tuple <tupleT> argtupleT;   // removes const, &, etc

    if constexpr(I >= std::tuple_size_v<tupleT>) {
        // Last case, if nothing is left to iterate, then exit the function
//        MADNESS_EXCEPTION("there is no madness function vector argument in the list, cannot partition the tasks", 1);
        return I;
    } else {
        using typeT = typename std::tuple_element<I, argtupleT>::type;// use decay types for determining a vector
        if constexpr (is_vector<typeT>::value) {
            return I;
        } else {
            // Going for next element.
            return get_index_of_first_vector_argument<tupleT,I+1>();
        }
    }
}

/// given a tuple return the index of the second argument that is a vector of Function<T,NDIM>
template<typename tupleT, std::size_t I>
constexpr std::size_t get_index_of_second_vector_argument() {
    constexpr std::size_t index0=get_index_of_first_vector_argument<tupleT,0>();
    return get_index_of_first_vector_argument<tupleT,index0+1>();
}

class Batch_1D {
    friend class MacroTaskPartitioner;

public:
    long begin=0, end=-1; ///< first and first past last index [begin,end)

    Batch_1D() {}
    Batch_1D(const Slice& s) : begin(s.start), end(s.end) {
        MADNESS_CHECK(s.step==1);
    }
    Batch_1D(const long& begin, const long& end) : begin(begin), end(end) {}

    bool operator==(const Batch_1D& other) const {
        return (end==other.end && begin==other.begin);
    }

    long size() const {
        return end - begin;
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
            auto v_batch=copy_batch(v);
            tupleT batched_arg = arg;
            std::get<index>(batched_arg) = v_batch;
            return batched_arg;
        }
    }

    /// given vector v, copy vector elements of v_batch into vector
    template<typename vecT>
    vecT insert_batch(vecT v, const vecT& v_batch) const {
      MADNESS_CHECK(v_batch.size()==size_t(this->size()) or this->is_full_size());
        std::copy(v_batch.begin(), v_batch.end(), v.begin()+begin);
        return v;
    }

    /// given vector v, return those vector elements inside this batch
    template<typename vecT>
    vecT copy_batch(const vecT v) const {
        vecT result_batch;
        long end1 = (end > 0) ? end : v.end() - v.begin();
        std::copy(v.begin() + begin, v.begin() + end1, std::back_inserter(result_batch));
        return result_batch;
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
    std::vector<Batch_1D> input;
    Batch_1D result;

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

    std::size_t size_of_input() const {
        std::size_t result=1;
        for (auto i: input) result*=i.size();
        return  result;
    }

    /// select the relevant vector elements from the argument tuple
    template<typename tupleT>
    tupleT copy_input_batch(const tupleT& arg) const {
        if (input.size()==0) return arg;
        tupleT arg1=input[0].template copy_batch<tupleT,0>(arg);
        if (input.size()>1) arg1=input[1].template copy_batch<tupleT,1>(arg1);
        MADNESS_CHECK(input.size()<=2);
        return arg1;
    }

    /// copy v_batch into the result vector
    template<typename vecT>
    vecT insert_result_batch(vecT v, const vecT& v_batch) const {
        return result.template insert_batch(v,v_batch);
    }

    /// pretty print this batch
    friend std::ostream& operator<<(std::ostream& os, const Batch& batch) {
        std::stringstream ss;
        ss << batch.result<< " <-- ";
        if (batch.input.size()>0)  ss << batch.input[0];
        if (batch.input.size()>1) ss << ", " << batch.input[1];
        os << ss.str();
        return os;
    }
};

/// partition one (two) vectors into 1D (2D) batches.

/// derive from this class and override the \it{do_partitioning} method if you want to implement
/// your custom partitioner
class MacroTaskPartitioner {
    friend class Batch;

public:
    typedef std::list<std::pair<Batch,double>> partitionT;
    std::size_t min_batch_size=5;           ///< minimum batch size
    std::size_t max_batch_size = 10;        ///< maximum batch size (for memory management)
    std::size_t nsubworld=1;                ///< number of worlds (try to have enough batches for all worlds)
    std::string policy = "guided";          ///< how to partition the batches
    std::size_t dimension = 1;              ///< partition one or two vectors

    MacroTaskPartitioner() {}

    virtual ~MacroTaskPartitioner() {}

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

    /// this will be called by MacroTask, it will *always* partition first (and possibly second) vector of arguments
    template<typename tupleT>
    partitionT partition_tasks(const tupleT& argtuple) const {

        constexpr std::size_t I1 = get_index_of_first_vector_argument<tupleT, 0>();
        constexpr std::size_t I2 = get_index_of_second_vector_argument<tupleT, 1>();
        std::size_t vsize1=1,vsize2=1;
        if constexpr (I2 < std::tuple_size_v<tupleT>) {     // found at least 2 vectors of madness functions
            constexpr std::size_t I2 = get_index_of_second_vector_argument<tupleT, 0>();
            vsize2 = std::get<I2>(argtuple).size();
        }
        if constexpr (I1 < std::tuple_size_v<tupleT>) { // found at least 1 vector
            constexpr std::size_t I1 = get_index_of_first_vector_argument<tupleT, 0>();
            vsize1 = std::get<I1>(argtuple).size();
        } else {
            std::string msg="confused partitioning dimension: "+std::to_string(dimension) +" typeid " + typeid(tupleT).name();
            MADNESS_EXCEPTION(msg.c_str(),1);
        }

        return do_partitioning(vsize1,vsize2,policy);
    }

    /// override this if you want your own partitioning
    virtual partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                               const std::string policy) const {
        if (dimension == 1) {
            return do_1d_partition(vsize1, policy);
        } else if (dimension == 2) {
            return do_2d_partition(vsize1, vsize2, policy);
        }
        return partitionT();
    }

    partitionT do_1d_partition(const std::size_t vsize, const std::string policy) const {
        partitionT result;
        if (policy == "guided") {
            long begin = 0;
            long end = 0;
            while (end < long(vsize)) {
                end += std::min(max_batch_size, std::max(min_batch_size, ((vsize - end) / nsubworld)));
                end = std::min(end, long(vsize));
                Batch batch(Batch_1D(begin, end),Batch_1D(begin,end));
                double priority=MacroTaskPartitioner::compute_priority(batch);
                result.push_back(std::make_pair(batch,priority));
                begin = end;
            }
        } else {
            std::string msg = "unknown partitioning policy: " + policy;
            MADNESS_EXCEPTION(msg.c_str(), 1);
        }
        return result;
    }

    /// outer product of 2 1d-partitionings -- result batches correspond to first input batches

    /// [begin1,end1)   <-- [begin1,end1) [begin2,end2)
    partitionT do_2d_partition(const std::size_t vsize, const std::size_t v2size, const std::string policy) const {
        partitionT partition1=do_1d_partition(vsize,policy);
        partitionT partition2=do_1d_partition(v2size,policy);
        partitionT result;
        for (auto p1 : partition1) {
            for (auto p2 : partition2) {
                Batch batch(p1.first.input[0],p2.first.input[0],p1.first.result);
                double priority=compute_priority(batch);
                result.push_back(std::make_pair(batch,priority));
            }
        }
        return result;
    }

    virtual double compute_priority(const Batch& batch) const {
        return batch.size_of_input();
    }

};

}

#endif //MADNESS_MACROTASKPARTITIONER_H
