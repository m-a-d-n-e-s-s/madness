//
// Ghaly
//
#ifndef __MADNESS_MRA_FUSET_BASEPARAMETERS_H__INCLUDED__
#define __MADNESS_MRA_FUSET_BASEPARAMETERS_H__INCLUDED__

#include <madness/world/MADworld.h>
#include "../../tensor/gentensor.h"
#include <madness/world/archive.h>

namespace madness
{
	template <typename T, std::size_t NDIM>
	class BaseParameters 
	{
	public:
		typedef GenTensor<T> coeffT;

	public:	
		coeffT				_coeff;	// for CompressOp, ReconstructOp

		bool				_isVector;
		std::vector<coeffT> _vCoeff;

	public:
		BaseParameters() 
		{
			_isVector = false;
		}

			//d(child_patch(kit.key())) += v[i].get().full_tensor_copy();	
			//
		BaseParameters(coeffT coeff)
		{
			_coeff = coeff;
		}
	
		~BaseParameters() { }

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
