#ifndef AD_SCALAR_H
#define AD_SCALAR_H

namespace ad
{
	typedef double Scalar;

	template <typename T>
	inline T max(T a, T b)
	{
	  return a > b ? a : b;
	}

	template <typename T>
	inline T min(T a, T b)
	{
	  return a < b ? a : b;
	}

};

#endif  // AD_SCALAR_H