#include "cvector/cvector.h"
#include "mvector.h"

inline double fraction(double a)
{
	return a - floor(a);
}

inline int integer(double a)
{
	return static_cast<int>(floor(a));
}

template <typename T, int rank>
class Interpol
{
	T		*f;
	unsigned	bits;
	Metric		M;

	public:
		Interpol() {}
		Interpol(T *f_, bits_):
			f(f_), bits(bits_), M(rank, bits)
		{}

		T operator()(mVector<double, rank> const &x) const
		{
			mVector<double, rank> A[2];
			std::transform(x.begin(), x.end(), A[1].begin(), fraction);
			std::fill(A[0].begin(), A[0].end(), 1.0);
			A[0] -= A[1];

			mVector<int, rank> start;
			std::transform(x.begin(), x.end(), start.begin(), integer);
			cVector S(rank, bits, M(start));

			T v = 0;
			for (cVector i(rank, 1, 0); i < (1 << rank); ++i)
			{
				float z = 1;
				for (unsigned k = 0; k < rank; ++k) z *= A[i[k]][k];
				v += f[S + i] * z;
			}

			return v;
		}
};

