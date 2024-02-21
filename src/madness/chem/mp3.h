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


    double compute_mp3_cd(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef_with_permutational_symmetry(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_ef_low_scaling(const Pairs<CCPair>& mp2pairs, const Pairs<std::vector<CCPairFunction<double,6>>> clusterfunctions) const;
    double compute_mp3_ghij(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_klmn(const Pairs<CCPair>& mp2pairs) const;
    double compute_mp3_klmn_fast(const Pairs<CCPair>& mp2pairs) const;
    double mp3_energy_contribution(const Pairs<CCPair>& mp2pairs) const;

};
}


#endif //MP3_H
