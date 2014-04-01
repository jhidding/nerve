#pragma once
#include "../conan/mvector.h"
#include "misc.h"
#include <vector>

using Conan::mVector;

template <int R> 
double vecdet(std::vector<mVector<double, R>> const &A, mVector<double, R> const &b);

template <>
double vecdet<2>(std::vector<mVector<double, 2>> const &A, mVector<double, 2> const &b)
{
	return det2(A[0] - b, A[1] - b);
}

template <>
double vecdet<3>(std::vector<mVector<double, 3>> const &A, mVector<double, 3> const &b)
{
	return det3(A[0] - b, A[1] - b, A[2] - b);
}

template <int R>
class Standalone_cell
{
	std::vector<mVector<double, R>> p;

	public:
		typedef mVector<double, R> Point;
		typedef std::vector<Point> Vec;

		Standalone_cell(std::vector<mVector<double, R>> const &s):
			p(s) {}

		Vec const &points() const { return p; }

		bool contains(mVector<double, R> const &q) const
		{
			bool l = true;
			for_each_decrement(p, 
				[&] (std::vector<mVector<double,R>> const &A, mVector<double, R> const &b)
			{
				if (vecdet(A, b) * vecdet(A, q) < 0)
					l = false;
			});

			return l;
		}

		template <typename Func>
		void for_each_vefa(Func blub) const;
};

template <>
template <typename Func>
void Standalone_cell<3>::for_each_vefa(Func blub) const
{
	mVector<int, 3> dx[3];
	for (unsigned k = 0; k < 3; ++k)
		dx[k][k] = 1;

	// for each cell
	for_each_decrement(p,
		[&] (Vec const &A, Point const &a)
	{
		// add vertex [a], edge [ab], face [abc], adjacency
		// one time for each permutation of [abcd]
		for_each_decrement(A,
			[&] (Vec const &B, Point const &b)
		{
			Point l = b - a;
			l /= Conan::norm(l);

			for_each_decrement(B,
				[&] (Vec const &C, Point const &c)
			{
				Point d = C[0], n = Conan::cross(b - a, c - a);
				n /= Conan::norm(n);
				Point q = Conan::cross(l, n);
				q *= (Conan::dot(q, c - a) > 0 ? 1 : -1); // pointing into face
				n *= (Conan::dot(n, d - a) > 0 ? 1 : -1); // pointing into cell

				blub(a.floor(), - Conan::dot(a, l) * Conan::dot(a, q) * Conan::dot(a, n) / 6.);
			});
		});

		// determine normal of face;
		Point n = Conan::cross(A[1] - A[0], A[2] - A[0]);
		n /= Conan::norm(n);
		n *= (Conan::dot(n, a - A[0]) > 0 ? 1 : -1); // pointing into cell

		for_all_ef_intersections(A,
			[&] (Point const &x, Point const &l, Point const &b, int k) // line intersection
		{
			int i = (k + 1) % 3, j = (k + 2) % 3;
			mVector<int, 3> X; X[k] = round_down(x[k] + 0.5); 
			X[i] = round_down(x[i]); X[j] = round_down(x[j]);

			Point u = Conan::cross(n, l);
			u *= (Conan::dot(Conan::cross(l, b - x), n) > 0 ? 1 : -1); // pointing into face

			blub(X, - Conan::dot(x, u) * Conan::dot(x, l) * Conan::dot(x, n) * fsign(l[k]) / 6.);
			blub(X - dx[k], Conan::dot(x, u) * Conan::dot(x, l) * Conan::dot(x, n) * fsign(l[k]) / 6.);

			Point q = calc_q(n, k, i, j), p = Conan::cross(n, q);
			q *= (Conan::dot(q, u) > 0 ? 1 : -1);

			blub(X, - Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) * fsign(p[k]) / 6.);
			blub(X - dx[k], Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) * fsign(p[k]) / 6.);

			p = Conan::cross(q, Point(dx[k]));
			p *= (Conan::dot(p, n) > 0 ? 1 : -1);
			blub(X, - Conan::dot(x, q) * Conan::dot(x, p) * x[k] / 6.);
			blub(X - dx[k], Conan::dot(x, q) * Conan::dot(x, p) * x[k] / 6.);
		},
			[&] (Point const &x, int k) // face intersection
		{
			int i = (k + 1) % 3, j = (k + 2) % 3;
			mVector<int, 3> X; X[k] = round_down(x[k]); 
			X[i] = round_down(x[i] + 0.5); X[j] = round_down(x[j] + 0.5);
			Point p, q;

			q = calc_q(n, i, j, k);
			q *= fsign(q[j]);
			p = Conan::cross(n, q);
			p *= fsign(p[i]);

			double z = - Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) / 6.;
			blub(X, z);
			blub(X - dx[i], -z);
			blub(X - dx[j], -z);
			blub(X - dx[j] - dx[i], z);

			p = Conan::cross(q, Point(dx[i]));
			p *= (Conan::dot(p, n) > 0 ? 1 : -1);
			z = - Conan::dot(x, q) * Conan::dot(x, p) * x[i] / 6.;
			blub(X, z);
			blub(X - dx[i], -z);
			blub(X - dx[j], -z);
			blub(X - dx[j] - dx[i], z);

			q = calc_q(n, j, i, k);
			q *= fsign(q[i]);
			p = Conan::cross(n, q);
			p *= fsign(p[j]);

			double w = - Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) / 6.;
			blub(X, w);
			blub(X - dx[j], -w);
			blub(X - dx[i], -w);
			blub(X - dx[i] - dx[j], w);

			p = Conan::cross(q, Point(dx[j]));
			p *= (Conan::dot(p, n) > 0 ? 1 : -1);
			w = - Conan::dot(x, q) * Conan::dot(x, p) * x[j] / 6.;
			blub(X, w);
			blub(X - dx[j], -w);
			blub(X - dx[i], -w);
			blub(X - dx[i] - dx[j], w);

			w = - x[k] * fsign(n[k]) * x[i] * x[j] / 3.;
			blub(X, w);
			blub(X - dx[j], -w);
			blub(X - dx[i], -w);
			blub(X - dx[i] - dx[j], w);
		});
	});
}

template <>
template <typename Func>
void Standalone_cell<2>::for_each_vefa(Func blub) const
{
	mVector<int, 2> dx[2];
	for (unsigned k = 0; k < 2; ++k)
		dx[k][k] = 1;

	for_each_decrement(p, [&]
		(std::vector<mVector<double, 2>> const &A, mVector<double, 2> const &a)
	{
		for_each_decrement(A, [&]
			(std::vector<mVector<double, 2>> const &B, mVector<double, 2> const &b)
		{
			mVector<double, 2> 
				c = B[0],
				l = b - a,
				p;

			l /= Conan::norm(l);
			p[0] = l[1]; p[1] = -l[0];
			p *= (Conan::dot(p, c - a) > 0 ? 1 : -1);
			blub(a.floor(), Conan::dot(a, l) * Conan::dot(a, p) / 2.);
		});

		mVector<double, 2>
			l = A[1] - A[0],
			p;

		l /= Conan::norm(l);
		p[0] = l[1]; p[1] = -l[0];
		// if sign is positive, triangle is on the right of the line A[0]->A[1]
		int sign = (Conan::dot(p, a - A[0]) > 0 ? 1 : -1);

		// p points into the triangle
		p *= sign;

		for_all_line_intersections(A[0], A[1],
			[&] (mVector<double, 2> const &x, unsigned k)
		{
			unsigned m = (k + 1) % 2;
			mVector<int, 2> X; X[k] = round(x[k]); X[m] = floor(x[m]);
			//std::cerr << x << std::endl << X << "\n\n\n";
			//std::cerr << x << " " << l << " " << 1 << std::endl;
			//std::cerr << x << " " << p << " " << 2 << std::endl;

			blub(X, x[k] * x[m] * fsign(p[m]) / 2.0);
			blub(X - dx[k], - x[k] * x[m] * fsign(p[m]) / 2.0);
			blub(X, Conan::dot(x, l) * Conan::dot(x, p) * fsign(l[k]) / 2.0);
			blub(X - dx[k], - Conan::dot(x, l) * Conan::dot(x, p) * fsign(l[k]) / 2.0);
		});
	});
}

