#pragma once

#include "../base/system.hh"

namespace Misc
{
	using System::Array;
	using System::mVector;
	using System::round_up;
	using System::round_down;

	template <unsigned R>
	using dVector = mVector<double, R>;

	template <unsigned R>
	using iVector = mVector<int, R>;

	template <unsigned R>
	using BoxPtr = ptr<Box<R>>;
}

