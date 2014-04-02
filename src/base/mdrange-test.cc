#ifdef UNITTEST
#include "unittest.hh"
#include "mdrange.hh"
#include "box.hh"
#include <iostream>

using namespace System;

Test::Unit mdrange_test("0013",
	"MdRange gives a multi-dimensional range of vectors, being "
	"indices into an array.",
	[] ()
{
	Box<3> b(4, 1.0);
	for (auto v : b.G)
	{
		std::cout << v << std::endl;
	}

	for (unsigned i = 0; i < 16; ++i)
	{
		std::cout << b.G[i] << std::endl;
	}
	return true;
});

Test::Unit box_test("0014",
	"test Box functionality",
	[] ()
{
	Box<2> b(8, 1.0);
	for (unsigned i = 6; i < 13; ++i)
	{
		auto x = b.I[i];
		std::cout << x << std::endl;
		for (auto &dx : b.block)
		{
			std::cout << "    " << x + dx << " -> " << b.idx(x + dx) << std::endl;
		}
	}
	return true;
});

#endif

