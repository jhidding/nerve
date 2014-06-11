#pragma once

#include "system.hh"

namespace Nerve
{
	template <R>
	void from_grid(BoxPtr<R> box, std::ifstream &fi, std::ofstream &fo)
	{
		auto X = System::load_from_file<dVector<R>>(fi, "pos");
	}
}

