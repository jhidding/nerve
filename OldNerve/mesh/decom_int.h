#pragma once
#include "misc.h"

template <typename T, int R>
int arg_min(mVector<T, R> const &v);

template <typename T>
int arg_min(mVector<T, 2> const &v)
{ return (v[0] < v[1] ? 0 : 1); }

template <typename T>
int arg_min(mVector<T, 3> const &v)
{ return (v[0] < v[1] ? (v[0] < v[2] ? 0 : 2) : (v[1] < v[2] ? 1 : 2)); }



inline int fsign(double a)
{
	if (a > 0.0) return 1;
	if (a < 0.0) return -1;
	return 0;
}

template <typename Func>
void Decomposition<2>::for_each_intersection(Func blub) const
{
	mVector<int, 2> dx[2];
	for (unsigned k = 0; k < 2; ++k)
		dx[k][k] = 1;

	// for each cell
	for (byte c_ = 0;  c_ < 2; ++c_)
	{
		double density = D.C[c_].unit_volume / d_cell_volume[c_];
		Decom_data_2::Cell const &c = D.C[c_];
		std::cerr << (d[c.vertices[0]] + d[c.vertices[1]] + d[c.vertices[2]])/3. << std::endl;

		for_each_decrement(c.vertices, [&]
			(Vec const &A, byte a)
		{
			for_each_decrement(A, [&]
				(Vec const &B, byte b)
			{
				mVector<double, 2> 
					c = d[B[0]],
					l = d[b] - d[a],
					p;

				l /= Conan::norm(l);
				p[0] = l[1]; p[1] = -l[0];
				p *= (Conan::dot(p, c - d[a]) > 0 ? 1 : -1);
				blub(d[a].floor(), density * Conan::dot(d[a], l) * Conan::dot(d[a], p) / 2.);
			});

			mVector<double, 2>
				l = d[A[1]] - d[A[0]],
				p;

			l /= Conan::norm(l);
			p[0] = l[1]; p[1] = -l[0];
			// if sign is positive, triangle is on the right of the line A[0]->A[1]
			int sign = (Conan::dot(p, d[a] - d[A[0]]) > 0 ? 1 : -1);

			// p points into the triangle
			p *= sign;

			for_all_line_intersections(d[A[0]], d[A[1]],
				[&] (mVector<double, 2> const &x, unsigned k)
			{
				unsigned m = (k + 1) % 2;
				mVector<int, 2> X; X[k] = round(x[k]); X[m] = floor(x[m]);
				//std::cerr << x << std::endl << X << "\n\n\n";
				//std::cerr << x << " " << l << " " << 1 << std::endl;
				//std::cerr << x << " " << p << " " << 2 << std::endl;

				blub(X, density * x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X - dx[k], - density * x[k] * x[m] * fsign(p[m]) / 2.0);
				blub(X, density * Conan::dot(x, l) * Conan::dot(x, p) * fsign(l[k]) / 2.0);
				blub(X - dx[k], - density * Conan::dot(x, l) * Conan::dot(x, p) * fsign(l[k]) / 2.0);
			});
		});
	}
}

template <typename Func>
void Decomposition<2>::loop_vefas_for_vertex(byte v_, Func blub) const
{
	Decom_data_2::Vertex const &v = D.V[v_];
	// for all incident edges
	for (auto e_ = v.begin(); e_ != v.end(); ++e_)
	{
		// e_ points to pair(target vertex idx, resulting edge)
		Decom_data_2::Edge const &e = D.E[e_->second];
		mVector<double, 2> 
			l = d[e_->first] - d[v_],
			p;

		l /= Conan::norm(l);
		// p perpendicular to l -> quarter turn
		p[0] = l[1];
		p[1] = -l[0];

		double z = Conan::dot(d[v_], l) * Conan::dot(d[v_], p);
		// both adjacent faces/cells (nomenclature compromise?)
		for (auto c_ = e.begin(); c_ != e.end(); ++c_)
		{
			Decom_data_2::Cell const &c = D.C[c_->second];
			double density = c.unit_volume / d_cell_volume[c_->second];
			blub(density * z * (Conan::dot(d[c_->first] - d[v_], p) > 0 ? 1 : -1) / 2.0);
		}
	}
}

template <typename Func>
void Decomposition<3>::loop_vefas_for_vertex(byte v_, Func blub) const
{
	Decom_data_3::Vertex const &v = D.V[v_];
	for (auto e_ = v.begin(); e_ != v.end(); ++e_)
	{
		Decom_data_3::Edge const &e = D.E[e_->second];
		mVector<double, 3> 
			l = d[e_->first] - d[v_];

		for (auto f_ = e.begin(); f_ != e.end(); ++f_)
		{
			Decom_data_3::Face const &f = D.F[f_->second];
			mVector<double, 3>
				n = d_normals[f_->second],
				p = Conan::cross(n, l);

			p *= (Conan::dot(d[f_->first] - d[v_], p) > 0 ? 1 : -1);

			for (auto c_ = f.begin(); c_ != f.end(); ++c_)
			{
				Decom_data_3::Cell const &c = D.C[c_->second];
				double density = c.unit_volume / d_cell_volume[c_->second];
				n *= (Conan::dot(d[c_->first] - d[v_], n) > 0 ? 1 : -1);

				blub(density * Conan::dot(d[v_], l) * Conan::dot(d[v_], p) * Conan::dot(d[v_], n));
			}
		}
	}
}

template <typename Func>
void Decomposition<3>::for_each_intersection(Func blub) const
{
	mVector<int, 3> dx[3];
	for (unsigned k = 0; k < 3; ++k)
		dx[k][k] = 1;

	// several problems force me to loop first over the cells
	// for all segments
	for (byte c_ = 0; c_ < 5; ++c_)
	{
		double density = D.C[c_].unit_volume / d_cell_volume[c_];

		// for all faces
		// for each cell
		for_each_decrement(D.C[c_].vertices,
			[&] (Vec const &A, byte a)
		{
			// add vertex [a], edge [ab], face [abc], adjacency
			// one time for each permutation of [abcd]
			for_each_decrement(A,
				[&] (Vec const &B, byte b)
			{
				Point l = d[b] - d[a];
				l /= Conan::norm(l);

				for_each_decrement(B,
					[&] (Vec const &C, byte c)
				{
					byte dd = C[0]; Point n = Conan::cross(d[b] - d[a], d[c] - d[a]);
					n /= Conan::norm(n);
					Point q = Conan::cross(l, n);
					q *= (Conan::dot(q, d[c] - d[a]) > 0 ? 1 : -1); // pointing into face
					n *= (Conan::dot(n, d[dd] - d[a]) > 0 ? 1 : -1); // pointing into cell

					blub(d[a].floor(), - density * Conan::dot(d[a], l) * Conan::dot(d[a], q) * Conan::dot(d[a], n) / 6.);
				});
			});

			// determine normal of face;
			Point n = Conan::cross(d[A[1]] - d[A[0]], d[A[2]] - d[A[0]]);
			n /= Conan::norm(n);
			n *= (Conan::dot(n, d[a] - d[A[0]]) > 0 ? 1 : -1); // pointing into cell

			std::vector<Point> V(3);
			V[0] = d[A[0]]; V[1] = d[A[1]]; V[2] = d[A[2]];
			for_all_ef_intersections(V,
				[&] (Point const &x, Point const &l, Point const &b, int k) // line intersection
			{
				int i = (k + 1) % 3, j = (k + 2) % 3;
				mVector<int, 3> X; X[k] = round_down(x[k] + 0.5); 
				X[i] = round_down(x[i]); X[j] = round_down(x[j]);

				Point u = Conan::cross(n, l);
				u *= (Conan::dot(Conan::cross(l, b - x), n) > 0 ? 1 : -1); // pointing into face

				double z = -density * Conan::dot(x, u) * Conan::dot(x, l) * Conan::dot(x, n) * fsign(l[k]) / 6.;
				blub(X, z);
				blub(X - dx[k], -z);

				Point q = calc_q(n, k, i, j), p = Conan::cross(n, q);
				q *= (Conan::dot(q, u) > 0 ? 1 : -1);
				z = -density * Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) * fsign(p[k]) / 6.;
				blub(X, z);
				blub(X - dx[k], -z);

				p = Conan::cross(q, Point(dx[k]));
				p *= (Conan::dot(p, n) > 0 ? 1 : -1);
				z = -density * Conan::dot(x, q) * Conan::dot(x, p) * x[k] / 6.; 
				blub(X, z);
				blub(X - dx[k], -z);
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

				double z = - density * Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) / 6.;
				blub(X, z);
				blub(X - dx[i], -z);
				blub(X - dx[j], -z);
				blub(X - dx[j] - dx[i], z);

				p = Conan::cross(q, Point(dx[i]));
				p *= (Conan::dot(p, n) > 0 ? 1 : -1);
				z = - density * Conan::dot(x, q) * Conan::dot(x, p) * x[i] / 6.;
				blub(X, z);
				blub(X - dx[i], -z);
				blub(X - dx[j], -z);
				blub(X - dx[j] - dx[i], z);

				q = calc_q(n, j, i, k);
				q *= fsign(q[i]);
				p = Conan::cross(n, q);
				p *= fsign(p[j]);

				double w = - density * Conan::dot(x, q) * Conan::dot(x, p) * Conan::dot(x, n) / 6.;
				blub(X, w);
				blub(X - dx[j], -w);
				blub(X - dx[i], -w);
				blub(X - dx[i] - dx[j], w);

				p = Conan::cross(q, Point(dx[j]));
				p *= (Conan::dot(p, n) > 0 ? 1 : -1);
				w = - density * Conan::dot(x, q) * Conan::dot(x, p) * x[j] / 6.;
				blub(X, w);
				blub(X - dx[j], -w);
				blub(X - dx[i], -w);
				blub(X - dx[i] - dx[j], w);

				w = - density * x[k] * fsign(n[k]) * x[i] * x[j] / 3.;
				blub(X, w);
				blub(X - dx[j], -w);
				blub(X - dx[i], -w);
				blub(X - dx[i] - dx[j], w);
			});
		});
	}
}

