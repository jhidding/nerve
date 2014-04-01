#pragma once
#include "../base/system.hh"
#include "../base/box.hh"

#include "../misc/algorithms.hh"
#include "../misc/linalg.hh"

namespace Nerve
{
	using Misc::vecdet;
	using System::Array;
	using System::mVector;
	using System::fsign;
	using System::round_down;

	using Misc::for_each_decrement;
	using Misc::for_each_line_intersection;

	template <unsigned R>
	class Cell
	{
		typedef mVector<double, R> Point;
		typedef Array<Point> PointSet;

		PointSet m_points;

		public:
			Cell(PointSet p):
				m_points(p)
			{}

			PointSet points() const { return m_points; }

			/* contains?
			 * Checks wether a point is inside or outside the cell
			 */
			bool contains(Point const &q) const
			{
				bool l = true;
				for_each_decrement(m_points, 
					[&] (PointSet A, Point const &b)
				{
					if (vecdet(A, b) * vecdet(A, q) < 0)
						l = false;
				});

				return l;
			}

			/* for-each-vertex-edge-face-adjacency
			 * Assuming the coordinates are defined on a grid, we may compute
			 * each /vefa/ by taking each intersection of a (R-n)-cell with the
			 * grid. We may pass this operator a function that is then executed
			 * on every occurence of a /vefa/. This means 'blub' needs to be some
			 * kind of accumulator function.
			 */
			template <typename Func>
			void for_each_vefa(Func blub) const;
	};

	template <>
	template <typename Func>
	Cell<2>::for_each_vefa<Func>(Func blub) const
	{
		mVector<int, 2> dx[2];
		for (unsigned k = 0; k < 2; ++k)
			dx[k][k] = 1;

		/* Each time we take out a vertex from the cell, we obtain a sub-cell
		 * by descending into the hierarchy of cells, we pass every combination
		 * of vertex and face.
		 */
		for_each_decrement(m_points, [&]
			(PointSet A, Point const &a)
		{
			// we have split the cell into a vertex and its opposite (N-1)-cell.
			// the vertex generates two adjacencies for itself.
			for_each_decrement(A, [&]
				(PointSet B, Point const &b)
			{
				Point c = B[0], l = b - a, p;

				l /= l.norm(); 			   // normalise l
				p[0] = l[1]; p[1] = -l[0];	   // p perp. to l
				p *= (p.dot(c - a) > 0 ? 1 : -1);  // pointing into the cell

				blub(a.floor(), a.dot(l) * a.dot(p) / 2.); 
				                                   // accumulate the result of 
								   // the vertex on the grid cell
			});

			// the oposite edge needs to be checked against intersections with
			// lines on the grid
			Point l = A[1] - A[0], p;
			l /= l.norm(); 				   // normalise l
			p[0] = l[1]; p[1] = -l[0]; 		   // p perp. to l
			// if sign is positive, triangle is on the right of the line A[0]->A[1]
			int sign = (p.dot(a - A[0]) > 0 ? 1 : -1);
			p *= sign; 				   // points into the cell

			// the edges of the grid are in unit directions, so the dot-product
			// is equal to just indexing a vector.
			//                    |      /|
			//         k+  -------+-----*-+------  
			//                    |    /  |
			//                    | X /   |
			//                    |  /    |
			//         k0   ------+-#-----+-------
			//                    |/      |
			//                    *X-dx[k]|
			//                   /|       |
			//         k-   ----*-+-------+-------
			//                 /  |       |
			//                    m-      m+
			for_each_line_intersection(A[0], A[1],
				[&] (Point const &x, unsigned k)
			{
				unsigned m = (k + 1) % 2;
				mVector<int, 2> X; X[k] = round(x[k]); X[m] = round_down(x[m]);

				blub(X,           x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X - dx[k], - x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X,           x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
				blub(X - dx[k], - x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
			});
		});
	}

	template <>
	template <typename Func>
	Cell<3>::for_each_vefa<Func>(Func blub) const
	{
	}
}
