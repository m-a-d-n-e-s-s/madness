//
// Ghaly
//
#ifndef __MADNESS_MRA_FUSET_BASEPARAMETERS_H__INCLUDED__
#define __MADNESS_MRA_FUSET_BASEPARAMETERS_H__INCLUDED__

#include <madness/world/MADworld.h>

namespace madness
{
	template <typename T, std::size_t NDIM>
	class BaseParameters 
	{
	private:
		int temp;

	public:
		BaseParameters() { }

		virtual ~BaseParameters() {}

		void setTemp(int value)
		{
			this->temp = value;
		}
		
		int getTemp()
		{
			return this->temp;
		}

		template<typename Archive> 
		void serialize(const Archive& ar) {}
	};
}
#endif	// 
