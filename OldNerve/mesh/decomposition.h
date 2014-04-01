#pragma once

#include <algorithm>
#include <vector>
#include "../conan/mvector.h"
#include "decom_data.h"

using Conan::mVector;

template <int R>
class Decomposition;

template <int R>
class Cell_iterator
{
	Decomposition<R> const &data;
	unsigned n;

	public:
		Cell_iterator(Decomposition<R> const &data_, unsigned n_):
			data(data_), n(n_) {}

		bool operator==(Cell_iterator const &other) const { return n == other.n; }
		bool operator!=(Cell_iterator const &other) const { return n != other.n; }
		Cell_iterator &operator++() { ++n; return *this; }

		bool contains(mVector<double, R> const &p) const { return data.cell_contains(n, p); }
		double unit_volume() const { return data.cell_unit_volume(n); }
		double volume() const { return data.cell_volume(n); }
};

template <>
class Decomposition<2>
{
	enum { R = 2 };
	static Decom_data_2 const dd_odd;
	static Decom_data_2 const dd_even;
	Decom_data_2 const &D;

	mVector<double, R> const *d;
	double d_cell_volume[2];

	public:
		typedef unsigned char byte;
		typedef std::vector<byte> Vec;

		double calc_volume(Vec const &v) const;
		bool cell_contains(byte c, mVector<double, 2> const &p) const;
		Vec const &cell_face(byte n, byte i) const { return D.E[D.C[n].edges[i]].vertices; }
		double cell_unit_volume(byte n) const { return D.C[n].unit_volume; }
		double cell_volume(byte n) const { return d_cell_volume[n]; }

		Cell_iterator<2> cells_begin() const { return Cell_iterator<2>(*this, 0); }
		Cell_iterator<2> cells_end() const { return Cell_iterator<2>(*this, 2); }

		Decomposition(int parity, mVector<double, 2> const *d_):
			D(parity == -1 ? dd_odd : dd_even)
		{
			d = d_;
			d_cell_volume[0] = fabs(calc_volume(D.C[0].vertices));
			d_cell_volume[1] = fabs(calc_volume(D.C[1].vertices));
		}

		template <typename Func>
		void for_each_intersection(Func f) const;

		template <typename Func>
		void loop_vefas_for_vertex(byte v_, Func f) const;
};

template <>
class Decomposition<3>
{
	enum { R = 3 };
	static Decom_data_3 const dd_odd; 
	static Decom_data_3 const dd_even; 
	Decom_data_3 const &D;

	mVector<double, 3> const *d;
	double d_cell_volume[5];

	mVector<double, 3> d_normals[16];

	public:
		typedef unsigned char byte;
		typedef std::vector<byte> Vec;
		typedef mVector<double, 3> Point;

		double calc_volume(Vec const &v) const;
		bool cell_contains(byte c, mVector<double, 3> const &p) const;
		double cell_unit_volume(byte n) const { return D.C[n].unit_volume; }
		double cell_volume(byte n) const { return d_cell_volume[n]; }
		Vec const &cell_face(byte n, byte i) const { return D.F[D.C[n].faces[i]].vertices; }

		mVector<double, 3> face_normal(Vec const &v) const;

		Cell_iterator<3> cells_begin() const { return Cell_iterator<3>(*this, 0); }
		Cell_iterator<3> cells_end() const { return Cell_iterator<3>(*this, 5); }

		Decomposition(int parity, mVector<double, 3> const *d_):
			D(parity == -1 ? dd_odd : dd_even)
		{
			d = d_;
			for (byte i = 0; i < 5; ++i)
			{
				d_cell_volume[i] = fabs(calc_volume(D.C[i].vertices));
			}

			for (byte i = 0; i < D.F.size(); ++i)
			{
				d_normals[i] = face_normal(D.F[i].vertices);
			}
		}

		template <typename Func>
		void for_each_intersection(Func f) const;

		template <typename Func>
		void loop_vefas_for_vertex(byte v_, Func f) const;
};

#include "decom_int.h"

