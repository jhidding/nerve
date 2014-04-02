#pragma once
#include "../base/system.hh"
#include "../base/box.hh"

#include "../misc/algorithms.hh"
#include "../misc/linalg.hh"

namespace Nerve
{
	using Misc::vecdet;
	using Misc::det2;
	using Misc::det3;

	using System::Array;
	using System::mVector;
	using System::fsign;
	using System::round_down;
	using System::Box;
	using System::ptr;

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
					if (vecdet<R>(A, b) * vecdet<R>(A, q) < 0)
						l = false;
				});

				return l;
			}

			Point centre_of_mass() const
			{
				Point p = std::accumulate(m_points.begin(), m_points.end(), Point(0));
				return p / m_points.size();
			}

			double volume() const;

			/* for-each-vertex-edge-face-adjacency
			 * Assuming the coordinates are defined on a grid, we may compute
			 * each /vefa/ by taking each intersection of a (R-n)-cell with the
			 * grid. We may pass this operator a function that is then executed
			 * on every occurence of a /vefa/. This means 'blub' needs to be some
			 * kind of accumulator function.
			 */
			template <typename Func>
			void for_each_vefa(ptr<Box<R>> box, Func blub) const;
	};

	template <>
	double Cell<2>::volume() const
	{
		return det2(m_points[1] - m_points[0], 
		            m_points[2] - m_points[0]) / 2;
	}

	template <>
	double Cell<3>::volume() const
	{
		return det3(m_points[1] - m_points[0], 
			    m_points[2] - m_points[0],
			    m_points[3] - m_points[0]) / 6;
	}
	
	template <>
	template <typename Func>
	void Cell<2>::for_each_vefa(ptr<Box<2>> box, Func blub) const
	{
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

				blub(X,                x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X - box->dx[k], - x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X,                x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
				blub(X - box->dx[k], - x.dot(l) * x.dot(p) * fsign(l[k]) / 2.0);
			});
		});
	}

	template <>
	template <typename Func>
	void Cell<3>::for_each_vefa(ptr<Box<3>> box, Func blub) const
	{
	}
}
