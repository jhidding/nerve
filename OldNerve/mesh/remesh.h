#pragma once

#include <functional>
#include <numeric>
#include <algorithm>
#include <utility>
#include "../progress.h"
#include "../conan/mvector.h"
#include "../conan/boxconfig.h"
#include "../conan/ndspace.h"
#include "../cvector/cvector.h"
#include "decomposition.h"

using Conan::mVector;
using Conan::BoxConfig;

template <typename T, int rank>
mVector<T, rank> modN(
	mVector<T, rank> const &v, cVector x1, cVector x2, 
	unsigned N, float L)
{
	mVector<T, rank> u = v;
	for (unsigned k = 0; k < rank; ++k)
	{
		if (int(x2[k]) - int(x1[k]) < 0)
		{
			u[k] += L;
		}
	}
	return u;
}

template <int R>
class Remesh_traits;

template <>
class Remesh_traits<2>
{
	public:
		enum { R = 2, tttp = 4, n_cells = 2 };

};

template <>
class Remesh_traits<3>
{
	public:
		enum { R = 3, tttp = 8, n_cells = 5 };
};

typedef std::shared_ptr<Conan::Cube<double>> dPtr;
typedef std::pair<dPtr, dPtr> dPair;

template <int R>
class Remesh: public Remesh_traits<R>
{
	public:
		typedef std::vector<mVector<double, R>> Vec;
		typedef Remesh_traits<R> Rt;

		static mVector<double, R> box_min(Vec const &d)
		{
			mVector<double, R> m = d.front();
			for (unsigned k = 0; k < R; ++k)
			{
				for (unsigned dx = 0; dx < d.size(); ++dx)
				{
					m[k] = std::min(m[k], d[dx][k]);
				}
			}
			return m;
		}

		static mVector<double, R> box_max(Vec const &d)
		{
			mVector<double, R> m = d.front();
			for (unsigned k = 0; k < R; ++k)
			{
				for (unsigned dx = 0; dx < d.size(); ++dx)
				{
					m[k] = std::max(m[k], d[dx][k]);
				}
			}
			return m;
		}

		template <typename T>
		static std::shared_ptr<Conan::Cube<double>> from_displacement(mVector<T, R> *src);

		template <typename Itr>
		static dPair from_cells(Itr begin, Itr end);
};

template <int R>
template <typename Itr>
dPair Remesh<R>::from_cells(Itr begin, Itr end)
{
	double Q = 1.0;
	dPair result(dPtr(new Conan::Cube<double>("grid")), dPtr(new Conan::Cube<double>("edges")));
	std::fill(result.first->begin(), result.first->end(), 0.0);
	std::fill(result.second->begin(), result.second->end(), 0.0);
//	std::shared_ptr<Conan::Cube<double>> result(new Conan::Cube<double>);
//	std::fill(result->begin(), result->end(), 0);

	// unit vectors at grid point
	mVector<int, R> I[Rt::tttp];
	for (cVector dx(R, 1, 0); dx < Rt::tttp; ++dx)
	{
		for (unsigned k = 0; k < R; ++k)
			I[dx][k] = (dx[k] == 0 ? -1 : 1);
	}

	size_t N = end - begin;
	Itr i = begin;
	#pragma omp parallel
	{
	#pragma omp for nowait
	for (size_t idx = 0; idx < N; ++idx)
	{
/*		std::for_each(i->points().begin(), i->points().end(),
			[&] (mVector<double, R> const &a)
		{
			std::cerr << a << " | ";
		}); std::cerr << "\n";

		std::cerr << "#";*/
		mVector<double, R> 	b_min = box_min(i->points()),		// calculate bounding box of mesh cell
					b_max = box_max(i->points());

		mVector<int, R> b_orig(b_min.floor());
		mVector<int, R> b_size = mVector<int, R>(b_max.ceil()) - b_orig;

//		std::cerr << " " << b_size << " ";
		Conan::NdSpace<R> box(b_size);
		for (auto j = box.begin(); j != box.end(); ++j)
		{
//			if ((*j)[0] == 0 and (*j)[1] == 0) std::cerr << ".";
			mVector<int, R> Y(*j + b_orig);
			if (i->contains(Y))
			{
				cVector y(R, BoxConfig::bits(), BoxConfig::M(Y));
				for (cVector dy(R, 1, 0); dy < Rt::tttp; ++dy)
				{
					#pragma omp critical
					{
					(*(result.first))[y - dy] += Q * Conan::trace(I[dy], Y);
					}
				}
			}
		}

		i->for_each_vefa([&] (mVector<int, R> const &p, double a)
		{
			#pragma omp critical 
			{
			(*(result.second))[BoxConfig::M(p)] += Q * a;
			}
		});

		++i;
	}
	}

	std::cerr << "\n";
	return result;
}

template <int R>
template <typename T>
std::shared_ptr<Conan::Cube<double>> Remesh<R>::from_displacement(mVector<T, R> *src)
{
	std::shared_ptr<Conan::Cube<double>> result(new Conan::Cube<double>);
	std::fill(result->begin(), result->end(), 0.0);

	// unit vectors at grid point
	mVector<int, R> I[Rt::tttp];
	for (cVector dx(R, 1, 0); dx < Rt::tttp; ++dx)
	{
		for (unsigned k = 0; k < R; ++k)
			I[dx][k] = (dx[k] == 0 ? 1 : -1) * (R % 2 == 0 ? 1 : -1);
	}

	// for all elements in the displacement field
	ProgressBar pb(BoxConfig::size(), "remeshing displacement field");
	//for (cVector x(R, BoxConfig::bits(), 0); x < BoxConfig::size(); ++x)
	#pragma omp parallel
	{
	#pragma omp for nowait
	for (size_t xxx = 0; xxx < BoxConfig::size(); ++xxx)
	{
		cVector x(R, BoxConfig::bits(), xxx);
		Vec 	d(Rt::tttp); 			// copy current working mesh cell;
		for (cVector dx(R, 1, 0); dx < Rt::tttp; ++dx)
			d[dx] = modN(src[x + dx], x, x+dx, BoxConfig::N(), BoxConfig::L()) / BoxConfig::scale();

		mVector<double, R> 	b_min = box_min(d),		// calculate bounding box of mesh cell
					b_max = box_max(d);

		mVector<int, R> b_orig(b_min.floor());
		mVector<int, R> b_size = mVector<int, R>(b_max.ceil()) - b_orig;

		// select cell decomposition, calculate all normals needed for further processing
		int idx_parity = 0;
		for (unsigned k = 0; k < R; ++k)
			idx_parity += x[k];
		idx_parity = (idx_parity % 2 ? 1 : -1);

		Decomposition<R> decom(idx_parity, &(d.front()));

		// the constant added depends on which part of a cubicle the grid point
		// ends up in, and how many of them. we need to check grid points against all
		// tetrahedra in the source cell. -> cost: box_size * 5
		Conan::NdSpace<R> box(b_size);
		for (auto i = box.begin(); i != box.end(); ++i)
		{
			mVector<int, R> Y(*i + b_orig);
			for (auto t = decom.cells_begin(); t != decom.cells_end(); ++t)
			{
				double density = t.unit_volume() / t.volume();
				if (t.contains(Y))
				{
					cVector y(R, BoxConfig::bits(), BoxConfig::M(Y));
					for (cVector dy(R, 1, 0); dy < Rt::tttp; ++dy)
					{
						#pragma omp critical
						{ (*result)[y - dy] += density * Conan::trace(I[dy], Y); }
					}
				}
			}
		}

		/*// loop vertices on source mesh -> cost: 120
		for (unsigned i = 0; i < Rt::tttp; ++i)
		{
			cVector y(R, BoxConfig::bits(), BoxConfig::M(d[i]));
		
		//	// calculates P.T * P.N * P.B voor each vertex-edge-face adjacency
			decom.loop_vefas_for_vertex(i, [&] (double v)
			{
				(*result)[y] += v;
			});
		}*/

		// wisdom 1: at each crossing an edge exits a cell and enters a cell
		// wisdom 2: each edge 'knows' which faces to involve
		// wisdom 3: vector T is the intersection with the face and the crossing plane,
		// 	for each combination of edge-face, this is known, by setting one
		// 	component of a precomputed vector to zero.
		// wisdom 4: the vector N into the cell, along the face is know because of (3)
		// wisdom 5: vector B is normal to the face and known
		// loop over edges on source mesh, trace crossings with integer points

		// cost: 
	
		// wisdom 1: each edge on the grid separates four cells
		// wisdom 2: six edges are involed, four of which are intersections of the face
		// 	with the grid.
		// wisdom 3: 	
		// loop over faces on source mesh, project along z-axis find xy integer points in triangle, repeat three times?
		//
		
		decom.for_each_intersection(
			[&] (mVector<int, R> const &p, double v)
		{
			#pragma omp critical
			{ (*result)[BoxConfig::M(p)] += v; }
		});

		pb.tic();
	}
	}
	pb.finish();

	return result;
}


