//
//
//
//
#ifndef __MADNESS_MRA_FUSET_FUSETCONTAINER__INCLUDED__
#define __MADNESS_MRA_FUSET_FUSETCONTAINER__INCLUDED__

#include "BaseParameters.h"

namespace madness
{
	template<typename classT>
	class FuseTContainer
	{
	public:
		classT	parameter;		

	public:
		FuseTContainer() {}

		void setParameter(classT& a);
	};

	template<typename classT>
	void 
	FuseTContainer<classT>::setParameter(classT& inputT)
	{
		this->parameter = inputT;
	}

};

#endif //
