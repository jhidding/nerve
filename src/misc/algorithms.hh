#pragma once
#include "../base/array.hh"
#include "../base/mvector.hh"

namespace Misc
{
	template <typename T, typename Func>
	inline void for_each_decrement(
		Array<T> A,
		Func f)
	{
		for (unsigned char i = 0; i < A.size(); ++i)
		{
			Array<T> B = A.copy(); 
			B->erase(B->begin() + i);
			f(B, A[i]);
		}
	}

	template <typename Func>
	void for_each_integer_in_interval(double a, double b, Func f)
	{
		double v = b - a;
		a = (v < 0.0 ? round_down(a) : round_up(a));

		for (double t = 1.0 - (b - a) / v; t < 1.0; t += fabs(1./v))
			f(t);
	}

	template <typename T, int R, typename Func>
	void for_each_line_intersection(
		mVector<T, R> const &a, mVector<T, R> const &b,
		Func f)
	{
		mVector<T, R> v = b - a;
		for (unsigned k = 0; k < R; ++k)
		{
			for_each_integer_in_interval(a[k], b[k],
				[&] (double t)
			{
				f(a + v * t, k);
			});
		}
	}
}

