//
// Created by Florian Bischoff on 2/15/24.
//

#ifndef MP3_H
#define MP3_H


#include <madness/mra/mra.h>
#include<madness/mra/commandlineparser.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/CCPotentials.h>
#include<madness/mra/QCCalculationParametersBase.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <madness/mra/macrotaskq.h>

namespace madness {
class MP3 : public CCPotentials {
public:

    MP3(World& world, const std::shared_ptr<Nemo> nemo, const CCParameters& param)
        : CCPotentials(world,nemo,param) {}
    MP3(const CCPotentials& ops) : CCPotentials(ops) {}

    double mp3_energy_contribution(const Pairs<CCPair>& mp2pairs) const;

    /// compute the MP3 energy contribution, macrotask version
    double mp3_energy_contribution_macrotask_driver(const Pairs<CCPair>& mp2pairs) const;

private:
    /// helper class for calculating the MP3 energy contributions
    class MacroTaskMP3 : public MacroTaskOperationBase {

        class ConstantPartPartitioner : public MacroTaskPartitioner {
        public:
            ConstantPartPartitioner() {};

            partitionT do_partitioning(const std::size_t& vsize1, const std::size_t& vsize2,
                                       const std::string policy) const override {
                partitionT p;
                for (size_t i = 0; i < vsize1; i++) {
                    Batch batch(Batch_1D(i,i+1), Batch_1D(i,i+1));
                    p.push_back(std::make_pair(batch,1.0));
                }
                return p;
            }
        };

    public:
        MacroTaskMP3(){partitioner.reset(new ConstantPartPartitioner());}

        typedef std::tuple<
                const std::string&,
                const std::vector<int>&,
                const std::vector<std::vector<CCPairFunction<double,6>>>& ,                 // all pairs ij
                const std::vector<Function<double,3>>&,
                const std::vector<Function<double,3>>&,
                const CCParameters&,
                const Molecule&,
                const Function<double,3>&,
                const std::vector<std::string>& > argtupleT;

        using resultT =std::shared_ptr<ScalarResult<double>>;

        resultT allocator(World& world, const argtupleT& argtuple) const {
            return std::shared_ptr<ScalarResult<double>>(new ScalarResult<double>(world));
        }

        resultT operator() (const std::string& diagram,                             // which diagram to calculate
                            const std::vector<int>& ij_vec,                         // dummy vector of size npair
                            const std::vector<std::vector<CCPairFunction<double,6>>>& pair_square,                 // all pairs ij
                            const std::vector<Function<double,3>>& mo_ket,          // the orbitals
                            const std::vector<Function<double,3>>& mo_bra,          // the orbitals*R2
                            const CCParameters& parameters,
                            const Molecule& molecule,
                            const Function<double,3>& Rsquare,
                            const std::vector<std::string>& argument) const {

            // the partitioner will break the input vector of pairs into single pairs
            MADNESS_CHECK(ij_vec.size()==1);

            // determine the orbital indices i and j for the pair
            // active occupied orbitals, the total length of pair_triangular is nact*(nact+1)
            const long nact=mo_ket.size()-parameters.freeze();
            MADNESS_CHECK(pair_square.size()==nact*nact);

            // the batch index is the ij composite index [0,nact*(nact+1)-1]
            const long ij=batch.result.begin;
            MADNESS_CHECK(batch.result.size()==1);

            // turn composite index ij into i and j, taking care of frozen orbitals
            PairVectorMap tri_map=PairVectorMap::triangular_map(parameters.freeze(),mo_ket.size());
            auto ij_to_i_and_j = [&tri_map](const int ij) { return tri_map.map[ij]; };
            auto [i,j]=ij_to_i_and_j(ij);

            PairVectorMap square_map=PairVectorMap::quadratic_map(parameters.freeze(),mo_ket.size());
            auto clusterfunctions=Pairs<std::vector<CCPairFunction<double,6>>>::vector2pairs(pair_square,square_map);

            double result=0.0;
            World& world=Rsquare.world();
            if (diagram=="cd")
                result= MP3::compute_mp3_cd(world,i,j,clusterfunctions,mo_ket,mo_bra,parameters,molecule,Rsquare,argument);
//            else if (diagram=="ef")
//                result= MP3::compute_mp3_ef(pair_triangular,pair_square,mo_ket,mo_bra,parameters,molecule,Rsquare,argument);
//            else if (diagram=="ghij")
//                result= MP3::compute_mp3_ghij(pair_triangular,pair_square,mo_ket,mo_bra,parameters,molecule,Rsquare,argument);
//            else if (diagram=="klmn")
//                result= MP3::compute_mp3_klmn(pair_triangular,pair_square,mo_ket,mo_bra,parameters,molecule,Rsquare,argument);
            else {
                std::string msg = "Unknown MP3 diagram: " + diagram;
                MADNESS_EXCEPTION(msg.c_str(), 1);
            }
            auto result1=std::shared_ptr<ScalarResult<double>>(new ScalarResult<double>(world));
            *result1=result;
            return result1;

        };


    };


    double compute_mp3_cd(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef_with_permutational_symmetry(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef_low_scaling(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;
    double compute_mp3_ef_as_overlap(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;
    double compute_mp3_ghij(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;
    double compute_mp3_ghij_fast(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;
    double compute_mp3_klmn(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_klmn_fast(const Pairs<CCPair>& mp2pairs) const;
    double mp3_test(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;

    static double compute_mp3_cd(World& world,
                                 const long i, const long j,
                                 const Pairs<std::vector<CCPairFunction<double,6>>>& pair_square,
                                 const std::vector<Function<double,3>>& mo_ket,
                                 const std::vector<Function<double,3>>& mo_bra,
                                 const CCParameters& parameters,
                                 const Molecule& molecule,
                                 const Function<double,3>& Rsquare,
                                 const std::vector<std::string>& argument);
    static double compute_mp3_ef(World& world,
                                 const std::vector<CCPair>& pair_triangular,
                                 const std::vector<CCPair>& pair_square,
                                 const std::vector<Function<double,3>>& mo_ket,
                                 const std::vector<Function<double,3>>& mo_bra,
                                 const CCParameters& parameters,
                                 const Molecule& molecule,
                                 const Function<double,3>& Rsquare,
                                 const std::vector<std::string>& argument);
    static double compute_mp3_ghij(World& world,
                                   const std::vector<CCPair>& pair_triangular,
                                   const std::vector<CCPair>& pair_square,
                                   const std::vector<Function<double,3>>& mo_ket,
                                   const std::vector<Function<double,3>>& mo_bra,
                                   const CCParameters& parameters,
                                   const Molecule& molecule,
                                   const Function<double,3>& Rsquare,
                                   const std::vector<std::string>& argument);
    static double compute_mp3_klmn(World& world,
                                   const std::vector<CCPair>& pair_triangular,
                                   const std::vector<CCPair>& pair_square,
                                   const std::vector<Function<double,3>>& mo_ket,
                                   const std::vector<Function<double,3>>& mo_bra,
                                   const CCParameters& parameters,
                                   const Molecule& molecule,
                                   const Function<double,3>& Rsquare,
                                   const std::vector<std::string>& argument);

};
}


#endif //MP3_H
