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

template<size_t I = 0, typename tupleT>
constexpr std::size_t get_index_of_first_vector_argument() {

    typedef decay_tuple <tupleT> argtupleT;   // removes const, &, etc

    if constexpr(I >= std::tuple_size_v<tupleT>) {
        // Last case, if nothing is left to iterate, then exit the function
        MADNESS_EXCEPTION("there is no madness function vector argument in the list, cannot partition the tasks", 1);
        return I;
    } else {
        // Print the tuple and go to next element
        using typeT = typename std::tuple_element<I, argtupleT>::type;       // use decay types for determining a vector
        if constexpr (is_madness_function_vector<typeT>::value) {
            return I;
        } else {
            // Going for next element.
            return get_index_of_first_vector_argument<I + 1, tupleT>();
        }
    }
}


class Batch {
    friend class MacroTaskPartitioner;

public:
    std::size_t begin, end; ///< first and first past last index [begin,end)

    Batch(const std::size_t& begin, const std::size_t& end) : begin(begin), end(end) {}

    /// select the relevant vector elements from the argument tuple
    template<typename tupleT>
    tupleT operator()(const tupleT& arg) const {

        constexpr std::size_t index = get_index_of_first_vector_argument<0, tupleT>();
        auto v = std::get<index>(arg);
        decltype(v) v_batch;
        std::copy(v.begin() + begin, v.begin() + end, std::back_inserter(v_batch));
        tupleT batched_arg = arg;
        std::get<index>(batched_arg) = v_batch;
        return batched_arg;
    }

    friend std::ostream& operator<<(std::ostream& os, const Batch& batch) {
        std::stringstream ss;
        ss << "[" << std::setw(3) << batch.begin << ", " << batch.end << ")";
        os << ss.str();
        return os;
    }

};

class MacroTaskPartitioner {
    friend class Batch;

public:
    typedef std::list<Batch> partitionT;
    std::size_t min_batch_size;
    std::size_t max_batch_size = 10;
    std::size_t nsubworld;
    std::string policy = "guided";

    MacroTaskPartitioner(long nsubworld, long min_batch_size = 5, std::string policy = "guided")
            : min_batch_size(min_batch_size), nsubworld(nsubworld), policy(policy) {}

    template<typename tupleT>
    partitionT partition_tasks(const tupleT& argtuple) const {

        constexpr std::size_t I= get_index_of_first_vector_argument<0, tupleT>();
        auto v=std::get<I>(argtuple);
        partitionT result;

        if (policy == "guided") {
            std::size_t begin = 0;
            std::size_t end = 0;
            while (end < v.size()) {
                end += std::min(max_batch_size, std::max(min_batch_size, ((v.size() - end) / nsubworld)));
                end = std::min(end, v.size());
                result.push_back(Batch(begin, end));
                begin = end;
            }
        } else {
            std::string msg = "unknown partitioning policy: " + policy;
            MADNESS_EXCEPTION(msg.c_str(), 1);
        }
        return result;
    }
};



///// partition orbitals into sets for the multiworld/macrotask algorithms
//class MacroTaskPartitioner {
//
//public:
//    static std::vector<std::pair<long,long> > partition_for_exchange(
//            long min_batch_size, long nsubworld, long nocc, std::string policy="guided");
//
//
//
//};

}

#endif //MADNESS_MACROTASKPARTITIONER_H
