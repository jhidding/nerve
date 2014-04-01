#pragma once
#include "../base/mvector.hh"
#include "../base/array.hh"

namespace Misc
{
	using System::mVector;

	template <typename T>
	inline double det2(
		mVector<T, 2> const &a,
		mVector<T, 2> const &b)
	{
		return (a[0] * b[1] - b[0] * a[1]);
	}

	template <typename T>
	inline double det3(
		mVector<T, 3> const &a,
		mVector<T, 3> const &b,
		mVector<T, 3> const &c)
	{
		return (  a[0] * (b[1] * c[2] - c[1] * b[2])
			- b[0] * (a[1] * c[2] - c[1] * a[2])
			+ c[0] * (a[1] * b[2] - b[1] * a[2]));
	}

	template <int R> 
	double vecdet(Array<mVector<double, R>> A, mVector<double, R> const &b);

	template <>
	double vecdet<2>(Array<mVector<double, 2>> A, mVector<double, 2> const &b)
	{
		return det2(A[0] - b, A[1] - b);
	}

	template <>
	double vecdet<3>(Array<mVector<double, 3>> A, mVector<double, 3> const &b)
	{
		return det3(A[0] - b, A[1] - b, A[2] - b);
	}
}

