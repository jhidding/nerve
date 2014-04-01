#pragma once

#include <algorithm>
#include "../conan/mvector.h"

namespace Conan
{
	for (auto s = D.S.begin(); s != D.S.end(); ++s)
	{
		mVector<T, 3> l_ = d[s->vertices[1]] - d[s->vertices[0]];
		l_ /= Conan::norm(l_);

		for_all_line_intersections(
			d[s->vertices[0]],
			d[s->vertices[1]],
			[&] (mVector<T, 3> const &x, int k)
		{
			// for all adjacent faces
			for (auto f = s->faces.begin(); f != s->faces.end(); ++f)
			{
				mVector<int, 3> X(x);
				mVector<T, 3>   const &n = d_normals[*f];
						p, q, u, l = l_ * fsign(l_[k]);

				int i, j;
				switch (k)
				{
					case 0: i = 1; j = 2; break;
					case 1: i = 0; j = 2; break;
					case 2: i = 0; j = 1; break;
				}

				u = Conan::cross(l, n); // TODO get u on good side of line
				on_intersect(X, Conan::dot(x, l) * Conan::dot(x, u) * Conan::dot(x, n);
				on_intersect(X, Conan::dot(x, l) * Conan::dot(x, u) * Conan::dot(x, -n);
				on_intersect(X - dx[k], Conan::dot(x, -l) * Conan::dot(x, u) * Conan::dot(x, n);
				on_intersect(X - dx[k], Conan::dot(x, -l) * Conan::dot(x, u) * Conan::dot(x, -n);

				q = calc_q(n, k, i, j); // TODO get q on good side of line
				p = calc_p(n, k, i, j);
				p *= fsign(p[k]);

				on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n);
				on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, -n);
				on_intersect(X - dx[k], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, n);
				on_intersect(X - dx[k], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, -n);
			}
		});
	}

	// for all faces
	for (byte f = 0; f < D.F.count(); ++f)
	{
		for_all_face_intersections(
			d[D.F[i].vertices[0]], 
			d[D.F[i].vertices[1]], 
			d[D.F[i].vertices[2]],
			[&] (mVector<T, 3> const &p, int k)
		{
			mVector<int, 3> X(x);
			mVector<T, 3>   const &n = d_normals[*f];
					p, q;

			int i, j;
			switch (k)
			{
				case 0: i = 1; j = 2; break;
				case 1: i = 0; j = 2; break;
				case 2: i = 0; j = 1; break;
			}
		
			q = calc_q(n, i, j, k);
			q *= fsign(q[j]);
			p = calc_p(n, i, j, k);
			p *= fsign(p[i]);

			on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n));
			on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, -n));
			on_intersect(X - dx[i], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, n));
			on_intersect(X - dx[i], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, -n));
			on_intersect(X - dx[j], Conan::dot(x, -q) * Conan::dot(x, p) * Conan::dot(x, n));
			on_intersect(X - dx[j], Conan::dot(x, -q) * Conan::dot(x, p) * Conan::dot(x, -n));
			on_intersect(X - dx[j] - dx[i], Conan::dot(x, -q) * Conan::dot(x, -p) * Conan::dot(x, n));
			on_intersect(X - dx[j] - dx[i], Conan::dot(x, -q) * Conan::dot(x, -p) * Conan::dot(x, -n));

			q = calc_q(n, j, i, k);
			q *= fsign(q[i]);
			p = calc_p(n, j, i, k);
			p *= fsign(p[j]);

			on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n));
			on_intersect(X, Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, -n));
			on_intersect(X - dx[j], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, n));
			on_intersect(X - dx[j], Conan::dot(x, q) * Conan::dot(x, -p) * Conan::dot(x, -n));
			on_intersect(X - dx[i], Conan::dot(x, -q) * Conan::dot(x, p) * Conan::dot(x, n));
			on_intersect(X - dx[i], Conan::dot(x, -q) * Conan::dot(x, p) * Conan::dot(x, -n));
			on_intersect(X - dx[i] - dx[j], Conan::dot(x, -q) * Conan::dot(x, -p) * Conan::dot(x, n));
			on_intersect(X - dx[i] - dx[j], Conan::dot(x, -q) * Conan::dot(x, -p) * Conan::dot(x, -n));
		});
	}

namespace Hex
{
	template <typename T, int Rank>
	inline double det(
		mVector<T, Rank> a,
		mVector<T, Rank> b,
		mVector<T, Rank> c)
	{
		return (  a[0] * (b[1] * c[2] - c[1] * b[2])
			- b[0] * (a[1] * c[2] - c[1] * a[2])
			+ c[0] * (a[1] * b[2] - b[1] * a[2]));
	}

	class Cell
	{
		unsigned char v[4];
		public:
			Cell(unsigned v1, unsigned v2, unsigned v3, unsigned v4)
			{
				v[0] = v1;
				v[1] = v2;
				v[2] = v3;
				v[3] = v4;
			}

			template <typename T>
			double volume(mVector<T, 3> cube[8])
			{
				return 1./6 * det(
						cube[v[1]] - cube[v[0]],
						cube[v[2]] - cube[v[0]],
						cube[v[3]] - cube[v[0]]);
			}
	};

	struct Decomposition
	{
		unsigned 	id;
		unsigned 	n;
		Cell		cells[6];

		template <typename T>
		double volume(mVector<T, 3> cube[8])
		{
			double v = 0.0;
			for (unsigned i = 0; i < n; ++i)
			{
				v += cells[i].volume(cube);
			}
			return v;
		}
	};

	extern unsigned n_decos;
	extern Decomposition decos[];
} // namespace Hex

} // namespace Conan

