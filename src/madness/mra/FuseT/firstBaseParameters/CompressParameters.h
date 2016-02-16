//
// Ghaly
//
#include "BaseParameters.h"

#include "../../tensor/gentensor.h"
#include "madness/world/archive.h"


namespace madness
{

	template <typename T, std::size_t NDIM>
	class CompressParameters:BaseParameters<T,NDIM>
	{
	public:
		typedef GenTensor<T> coeffT;

	private:
		coeffT	_coeff;

	public:
					CompressParameters() : _coeff() {}
		coeffT&		getCoeff()
		{
			return const_cast<coeffT&>(this->_coeff);
		}

		// Future<coeffT > result(node.coeff());
		//
		//	node <- nodeT <- FunctionNode
		//	node.coeff() // coeffT& coeff() { return const_cast<coeffT&>(_coeffs);
		//
		void		setCoeff(const coeffT& coeffs)
		{
			this->_coeff = coeffs;
		}
	};

}
