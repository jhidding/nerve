#pragma once
#include "../base/system.hh"
#include "../base/box.hh"

#include "../misc/algorithms.hh"
#include "../misc/linalg.hh"

namespace Nerve
{
	using System::Array;
	using System::mVector;

	template <unsigned R>
	void from_txt(ptr<Box<R>> box, std::istream &fi, std::ostream &fo)
	{
		std::cerr << "Reading triangles from file.\n";

		Array<Cell<R>> data;
	}
}

