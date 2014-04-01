#pragma once
#include "../base/system.hh"
#include "../base/box.hh"

#include "../misc/algorithms.hh"
#include "../misc/linalg.hh"

namespace Nerve
{
	using System::ptr;
	using System::Box;
	using System::Array;
	using System::mVector;

	template <unsigned R>
	void from_txt(ptr<Box<R>> box, std::istream &fi, std::ostream &fo)
	{
		typedef mVector<double, R> Point;
		typedef Array<Point> PointSet;

		std::cerr << "Reading triangles from file ... ";

		Array<Cell<R>> data;

		/* read triangles */
		unsigned M; fi >> M;
		for (unsigned i = 0; i < M; ++i)
		{
			PointSet pts;

			for (unsigned j = 0; j < R+1; ++j)
			{
				Point p;
				for (unsigned k = 0; k < R; ++k)
					fi >> p[k];
				pts->push_back(p / box->scale());
			}

			data.push_back(Cell<R>(data));
		}

		std::cerr << "read " << data.size() << " cells. [done]\n";

		Array<double> rho = Remesh<R>::from_cells(box, data);

		if (R == 2)
		{
			size_t o = 0;
			for (unsigned j = 0; j < box->N(); ++j)
			{
				for (unsigned i = 0; i < box->N(); ++i, ++o)
				{
					fo << rho[o] << " ";
				}

				fo << "\n";
			}
		}
	}
}

