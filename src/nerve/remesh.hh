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
	using System::MdRange;

	template <unsigned R>
	class Remesh
	{
		public:
			typedef mVector<int, R> iPoint;
			typedef mVector<double, R> Point;
			typedef Array<Point> PointSet;

			static Array<double> from_cells(ptr<Box<R>> box, Array<Cell<R>> data);
	};

	template <unsigned R>
	Array<double> Remesh<R>::from_cells(ptr<Box<R>> box, Array<Cell<R>> data)
	{
		Array<double> A(box->size(), 0.0);
		Array<size_t> bread_crumbs(box->size(), 0);

		for (size_t i = 0; i < data.size(); ++i)
		{
			double f = 1./fabs(data[i].volume());

			// draw the contours
			data[i].for_each_vefa(box, [&] (iPoint const &p, double a)
			{
				size_t idx = box->idx(p);
				bread_crumbs[idx] = i*3+1;
				A[idx] += f * a;
			});

			Misc::flood_fill(
				iPoint(data[i].centre_of_mass()),
				[&] (iPoint const &x)
			{
				size_t idx = box->idx(x);
				Array<iPoint> n;

				// if on boundary, only continue on boundary
				if (bread_crumbs[idx] == i*3 + 3)
				{
					MdRange<R> DX(iPoint(3)); iPoint one(1);
					for (auto &dx : DX)
					{
						iPoint p = x + dx - one;
						size_t idx = box->idx(p);
						size_t bc = bread_crumbs[idx];
						if ((bc == i*3+1) or 
						    (bc <  i*3+1) and data[i].contains(p))
							n->push_back(p);
					}
				}
				else {
					MdRange<R> DX(iPoint(3)); iPoint one(1);
					for (auto &dx : DX)
					{
						n->push_back(x + dx - one);
					}
				}

				return n;
			},
				[&] (iPoint const &x)
			{
				size_t idx = box->idx(x);
				if (bread_crumbs[idx] > i*3 + 1) 
				{
					return false;
				}

				return true;
			},
				[&] (iPoint const &x)
			{
				size_t idx = box->idx(x);
				if (bread_crumbs[idx] == i*3 + 1)
				{
					bread_crumbs[idx] = i*3 + 3;
					if (not data[i].contains(x)) return;
				}
				else
					bread_crumbs[idx] = i*3 + 2;

				// the value added is +/-(x[0]*x[1]*...)
				for (iPoint const &dx : box->block)
				{
					iPoint p = (x - dx) % box->N();
					iPoint x = p + dx;
					double a = x.prod();

					size_t idx = box->idx(p);
					int sign = (dx.sum() % 2 == 0 ? 1 : -1);
					A[idx] += sign * f * a;
				}
			});

			/*
			// flood-fill within breadcrumbs
			Misc::flood_fill(
				// start position
				iPoint(data[i].centre_of_mass()),

				// give back neighbours
				[&] (iPoint const &x)
			{
				size_t idx = box->idx(x);
				// if we are on the boundary return no neighbours
				if (bread_crumbs[idx] == i*2+1)
					return Array<iPoint>(0);

				Array<iPoint> n;
				for (unsigned k = 0; k < R; ++k)
				{
					n->push_back(x - box->dx[k]);
					n->push_back(x + box->dx[k]);
				}
				return n;
			},

				// predicate
				[&] (iPoint const &x)
			{
				size_t idx = box->idx(x);

				// if on boundary, check for insideness
				if (bread_crumbs[idx] == i*2 + 1)
					return data[i].contains(x);
				// otherwise, see if already done
				else
					return bread_crumbs[idx] != i*2 + 2;
			},
				
				// action
				[&] (iPoint const &x)
			{
				std::cout << x << std::endl;

				size_t idx = box->idx(x);
				// never overdraw the boundary
				if (bread_crumbs[idx] != i*2+1)
					bread_crumbs[idx] = i*2+2;

				// the value added is +/-(x[0]*x[1]*...)
				double a = x.prod();
				for (iPoint const &dx : box->block)
				{
					size_t idx = box->idx(x - dx);
					int sign = (dx.sum() % 2 == 0 ? 1 : -1);
					A[idx] += -f * a;
				}
			});*/
		}

		return A;
	}
}

