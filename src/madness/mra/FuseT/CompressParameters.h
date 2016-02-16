//
// Ghaly
//
#ifndef __MADNESS_MRA_FUSET_COMPRESSPARAMETERS_H__INCLUDED__
#define __MADNESS_MRA_FUSET_COMPRESSPARAMETERS_H__INCLUDED__

#include <madness/world/MADworld.h>
#include "../../tensor/gentensor.h"
#include <madness/world/archive.h>

namespace madness
{
	template <typename T>
	class CompressParameters 
	{
	public:
		typedef GenTensor<T> coeffT;

	public:	
		coeffT				_coeff;	// for CompressOp, ReconstructOp

		bool				_isVector;
		std::vector<coeffT> _vCoeff;

	public:
		CompressParameters() 
		{
			_isVector = false;
		}

		CompressParameters(coeffT coeff)
		{
			_coeff = coeff;
		}
	
		~CompressParameters() { }

		template <typename Archive> 
		void serialize(const Archive& ar) 
		{
			ar & _coeff;
			ar & _isVector;
			ar & _vCoeff;
		}
	};
}
#endif	// 
