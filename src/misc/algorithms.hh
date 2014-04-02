#pragma once
#include "../base/array.hh"
#include "../base/mvector.hh"
#include "../base/misc.hh"

#include <queue>

namespace Misc
{
	using System::Array;
	using System::mVector;
	using System::round_up;
	using System::round_down;

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

	template <typename Point, typename Func>
	void for_each_line_intersection(
		Point const &a, Point const &b,
		Func f)
	{
		Point v = b - a;
		for (unsigned k = 0; k < Point::size(); ++k)
		{
			for_each_integer_in_interval(a[k], b[k],
				[&] (double t)
			{
				f(a + v * t, k);
			});
		}
	}

	template <typename Coord, typename Next, typename Pred, typename Actor>
	bool flood_fill(Coord start, Next next, Pred pred, Actor act)
	{
		std::queue<Coord> Q;

		if (pred(start))
		{
			act(start);
			Q.push(start);
		}
		else
		{
			return false;
		}

		while (!Q.empty())
		{
			for (Coord const &n : next(Q.front()))
			{
				if (pred(n))
				{
					act(n);
					Q.push(n);
				}
			}

			Q.pop();
		}

		return true;
	}
}

