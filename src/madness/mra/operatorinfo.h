//
// Created by Florian Bischoff on 11/2/23.
//

#ifndef MADNESS_OPERATORINFO_H
#define MADNESS_OPERATORINFO_H

namespace madness {

/// operator types
enum OpType {
    OT_UNDEFINED,
    OT_ONE,         /// indicates the identity
    OT_G12,         /// 1/r
    OT_SLATER,      /// exp(-r)
    OT_GAUSS,       /// exp(-r2)
    OT_F12,         /// 1-exp(-r)
    OT_FG12,        /// (1-exp(-r))/r
    OT_F212,        /// (1-exp(-r))^2
    OT_F2G12,       /// (1-exp(-r))^2/r = 1/r + exp(-2r)/r - 2 exp(-r)/r
    OT_BSH,         /// exp(-r)/r
    OT_SIZE         /// for ending loops
};

/// operator type to string
template<std::size_t N=1>   // dummy template argument to avoid ambiguity with the other operator<<
std::ostream& operator<<(std::ostream& os, const OpType type) {
    auto name = [](OpType type) {
        switch (type) {
            case OpType::OT_UNDEFINED:
                return "undefined";
            case OpType::OT_ONE:
                return "identity";
            case OpType::OT_G12:
                return "g12";
            case OpType::OT_SLATER:
                return "slater";
            case OpType::OT_GAUSS:
                return "gauss";
            case OpType::OT_F12:
                return "f12";
            case OpType::OT_FG12:
                return "fg12";
            case OpType::OT_F212:
                return "f12^2";
            case OpType::OT_F2G12:
                return "f12^2g";
            case OpType::OT_BSH:
                return "bsh";
            default:
                return "undefined";
        }
    };
    os << name(type);
    return os;
}

struct OperatorInfo {
    OperatorInfo() = default;
    OperatorInfo(double mu, double lo, double thresh, OpType type, std::optional<bool> truncate = {}) : mu(mu), lo(lo), thresh(thresh), type(type), truncate_lowexp_gaussians(truncate) { }
    double mu=0.0;     ///< some introspection
    double lo=1.e-5;
    double thresh=1.e-4;
    OpType type=OT_UNDEFINED;    ///< introspection
    double hi=-1.0;
    bool debug=false;
    std::optional<bool> truncate_lowexp_gaussians;  // if given, overrides the default for whether to truncate low-exponent gaussians
};



}
#endif //MADNESS_OPERATORINFO_H
