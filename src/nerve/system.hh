#pragma once

#include "../base/system.hh"

namespace Nerve
{
	using System::mVector;

	template <unsigned R>
	using dVector = mVector<double, R>;

	template <unsigned R>
	using iVector = mVector<int, R>;

	template <unsigned R>
	using BoxPtr = ptr<Box<R>>;
}

