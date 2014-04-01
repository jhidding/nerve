#pragma once
#include <iostream>
#include <functional>
#include <algorithm>

#include "conan/boxconfig.h"
#include "conan/getopt.h"
#include "conan/cube.h"
#include "kdtree.h"
#include "floodfill.h"
#include "mesh/decomposition.h"
#include "mesh/remesh.h"
#include "mesh/so_cell.h"
#include "repeat.h"

template <int R>
double calculate_density(Conan::mVector<double, R> *f, cVector const &x, unsigned N, double L);

template <>
inline double calculate_density<2>(Conan::mVector<double, 2> *f, 
		cVector const &x, unsigned N, double L)
{
	Conan::mVector<double, 2> C[4];
	for (cVector dx(2, 1, 0); dx < (1 << 2); ++dx)
		C[dx] = Conan::mVector<double, 2>(
				modN(f[x + dx], x, x+dx, N, L));

	return 0.5 * (
		(C[1][0] - C[0][0]) * (C[2][1] - C[0][1]) - 
		(C[1][1] - C[0][1]) * (C[2][0] - C[0][0]) +
		(C[2][0] - C[3][0]) * (C[1][1] - C[3][1]) - 
		(C[2][1] - C[3][1]) * (C[1][0] - C[3][0]));
}

template <>
inline double calculate_density<3>(Conan::mVector<double, 3> *f,
		cVector const &x, unsigned N, double L)
{
	Conan::mVector<double, 3> C[8];
	for (cVector dx(3, 1, 0); dx < (1 << 3); ++dx)
		C[dx] = Conan::mVector<double, 3>(
				modN(f[x + dx], x, x+dx, N, L));

	int idx_parity = 1;
	for (unsigned k = 0; k < 3; ++k)
		idx_parity *= (x[k] % 1 == 0 ? 1 : -1);
	Decomposition<3> deco(idx_parity, C);
	double a = 0.0; for (unsigned i = 0; i < 5; ++i) a += deco.cell_volume(i);
	return a;
}

inline int random_colour(size_t seed)
{
	if (seed < 0) return -1;
	srand(seed);
	return rand() % 1024;
}

std::function<double ()> random_normal(double m, double s)
{
	std::shared_ptr<std::mt19937> engine(new std::mt19937(time(0)));
	std::shared_ptr<std::normal_distribution<>> dist(new std::normal_distribution<>(m, s));

	return [engine, dist] () { return (*dist)(*engine); };
}

template <int R>
std::function<mVector<double, R> ()> random_point(double m, double s)
{
	std::function<double ()> f = random_normal(m, s);
	return [f] () {
		mVector<double, R> x;
		std::generate(x.begin(), x.end(), f);
		return x;
	};
}

template <int R>
std::function<Standalone_cell<R> ()> random_tet()
{
	std::function<mVector<double, R> ()> f = random_point<R>(128, 64);

	return [f] () {
		std::vector<mVector<double, R>> v(R + 1);
		std::generate(v.begin(), v.end(), f);
		std::for_each(v.begin(), v.end(), 
			[&] (mVector<double, R> const &a)
		{
			std::cerr << a << " : ";
		});
		std::cerr << "\n";
		return Standalone_cell<R>(v);
	};
}

template <unsigned R>
void from_txt(std::istream &in, std::ostream &out)
{
	std::cerr << "Reading input data ... ";

	size_t N;
	in >> N;

	std::vector<Standalone_cell<R>> data;
	for (size_t q = 0; q < N; ++q)
	{
		std::vector<mVector<double, R>> triangle;

		for (unsigned i = 0; i < (R+1); ++i)
		{
			mVector<double, R> p;
			for (unsigned k = 0; k < R; ++k)
				in >> p[k];
			triangle.push_back(p / BoxConfig::scale());
		}

		data.push_back(Standalone_cell<R>(triangle));
	}

	std::cerr << "[done], read " << data.size() << " triangles.\n";

	std::cerr << "Plotting them on a mesh ... ";
	auto rho_m = Remesh<R>::from_cells(data.begin(), data.end());
	std::cerr << "[done]\n";

	size_t idx = 0;
	for (unsigned i = 0; i < BoxConfig::N(); ++i)
	{
		for (unsigned j = 0; j < BoxConfig::N(); ++j)
		{
			out << (*rho_m.first)[idx]  + (*rho_m.second)[idx] << " ";
			++idx;
		}
		out << "\n";
	}

	out << "\n\n";
	idx = 0;
	for (unsigned i = 0; i < BoxConfig::N(); ++i)
	{
		for (unsigned j = 0; j < BoxConfig::N(); ++j)
		{
			out << (*rho_m.first)[idx] << " ";
			++idx;
		}
		out << "\n";
	}
}

template <int R>
void nerve(std::istream &in, std::ostream &out)
{
//	uncomment as soon as needed
//	Conan::GetOpt &options = Conan::GetOpt::instance();

//	unsigned N = Conan::BoxConfig::N();
//	double L = Conan::BoxConfig::L();
//	size_t size = Conan::BoxConfig::size();
//	unsigned bits = Conan::BoxConfig::bits();

	// read particle positions
	std::cerr << "Reading input data ... ";
	Conan::Cube<Conan::mVector<double, R>> X;
	X.convert_from_file(in);
	std::cerr << "[done]\n";
/*
	// permutation cube
	Conan::Cube<size_t> P;
	size_t i = 0;
	std::generate(P.begin(), P.end(), [&i] () { return i++; });

	// build KdTree
	std::cerr << "Building kdTree ... ";
	kdTree::Tree<size_t, R> T(P.begin(), P.end(),
		[&X] (size_t i, int d) { return X[i][d]; });
	std::cerr << "[done]\n";

	// calculate Lagrangian density
	std::cerr << "Calculating density ... ";
	Conan::Cube<double> rho_q;
	for (cVector x(R, bits, 0); x < size; ++x)
	{
		rho_q[x] = calculate_density<R>(X.data(), x, N, L);
	}
	std::cerr << "[done]\n";

	// flood-fill algorithm to find tessellation on Lagrangian space
	std::cerr << "Finding 6d phase-space tessellation ... ";
	Conan::Cube<size_t> tessellation;
	std::fill(tessellation.begin(), tessellation.end(), 0);
	Conan::FloodFill<R> FF;

	i = 1;
	for (cVector x(R, bits, 0); x < size; ++x)
	{
		if (FF.unary_fill([&] (cVector const &y)
			{
				if (tessellation[y] != 0)
					return false;

				bool s = ((rho_q[x] * rho_q[y]) >= 0);
				if (s)
				{
					tessellation[y] = i;
					return true;
				}
				else
				{
					return false;
				}
			}, x))
		{
			++i;
		}
	}
	std::cerr << i << " tiles! ... [done]\n";
*/
	// eulerian density
	std::cerr << "Remeshing to calculate eulerian density ... ";
	auto rho_x = Remesh<R>::from_displacement(X.data());
	std::cerr << "[done]\n";
	rho_x->to_file(out);

	//Repeat<Standalone_cell<R>> tets(random_tet<R>(), 10);
	//auto rho_m = Remesh<R>::from_cells(tets.begin(), tets.end());
	//rho_m.first->to_file(out);
	//rho_m.second->to_file(out);

	//Conan::Cube<double> sum("sum");
	//for (size_t x = 0; x < BoxConfig::size(); ++x)
	//	sum[x] = (*rho_m.first)[x] + (*rho_m.second)[x];
	//sum.to_file(out);
	
	/*Conan::Cube<unsigned> stream_count;
	std::fill(stream_count.begin(), stream_count.end(), 0);
	*/

	/*
	for (cVector x(R, bits, 0); x < size; ++x)
	{
		std::cout << x[0] << " " << x[1] << " " << rho_q[x] << " "
			<< (*rho_x)[x] << " " 
			<< (*rho_m)[x].first << " "
			<< (*rho_m)[x].second << " "
			<< random_colour(tessellation[x]) << std::endl;
	}*/
}

