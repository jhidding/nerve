#pragma once
#include "../base/system.hh"
#include "../base/box.hh"

#include "../misc/algorithms.hh"
#include "../misc/linalg.hh"

#include "cell.hh"

namespace Nerve
{
	using System::ptr;
	using System::Box;
	using System::Array;
	using System::mVector;

	template <unsigned R>
	class Remesh
	{
		public:
			typedef mVector<double, R> Point;
			typedef Array<Point> PointSet;

			static Array<double> from_cells(ptr<Box<R>> box, Array<Cell<R>> data);
	};

	template <unsigned R>
	static Array<double> Remesh<R>::from_cells(ptr<Box<R>> box, Array<Cell<R>> data)
	{
		Array<double> A(box->size());

		for (Cell<R> c : data)
		{
		}
	}
}
