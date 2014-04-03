#pragma once
#include "mvector.hh"
#include "mdrange.hh"
#include "fourier.hh"

namespace System
{
	template <unsigned R>
	using GridSpace = System::Map<System::MdRange<R>, System::mVector<double, R>>;

	template <unsigned R>
	GridSpace<R> grid_space(unsigned N, double L)
	{
		typename GridSpace<R>::arg_type shape(N);
		System::MdRange<R> X(shape);

		return GridSpace<R>(X, [N, L] (typename GridSpace<R>::arg_type const &x)
		{ 
			mVector<double, R> y;
			for (unsigned i = 0; i < R; ++i)
				y[i] = x[i] * (L/N);
			return y; 
		});
	}

	template <unsigned R>
	class Box
	{
		bool			mPeriodic;
		unsigned 		mN;
		double 			mL, mRes, mRes2;
		std::vector<unsigned>	mShape;
		std::vector<size_t>	mStride;
		size_t			mSize;

		public:
			Fourier::KSpace<R> const K;
			GridSpace<R>	   const G;
			System::MdRange<R> const I;
			std::vector<mVector<int, R>> block;
			std::vector<mVector<int, R>> dx;

			Box(unsigned N_, double L_, bool p_ = true):
				mPeriodic(p_),
				mN(N_), mL(L_), mRes(L_/N_), mRes2(mRes*mRes),
				mShape(R, N_), mStride(R),
				mSize(ipow(N_, R)),
				K(Fourier::kspace<R>(N_, L_)),
				G(grid_space<R>(N_, L_)), I(mVector<int, R>(N_)),
				block(1 << R), dx(R, mVector<int, R>(0))
			{
				mStride[0] = 1;
				for (unsigned i = 1; i < R; ++i)
				{
					mStride[i] = mStride[i-1] * N_;
				}
				
				for (unsigned i = 0; i < R; ++i)
				{
					std::cerr << mStride[i] << " ";
					dx[i][i] = 1;
				}
				std::cerr << std::endl;

				MdRange<R> b(mVector<int, R>(2));
				for (unsigned i = 0; i < (1 << R); ++i)
				{
					block[i] = b[i];
				}
			}

			bool periodic() const { return mPeriodic; }
			double L() const { return mL; }
			unsigned N() const { return mN; }
			double scale() const { return mRes; }
			double scale2() const { return mRes2; }
			std::vector<unsigned> const &shape() const { return mShape; }
			size_t size() const { return mSize; }

			size_t idx(mVector<int, R> const &x) const
			{
				size_t a = 0;
				for (unsigned i = 0; i < R; ++i)
					a += modulus(x[i], mN) * mStride[i];
				return a;
			}
	};
}

