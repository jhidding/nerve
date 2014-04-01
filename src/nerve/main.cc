#include "../base/system.hh"
#include "../base/format.hh"
#include "../base/box.hh"

#include "from_txt.hh"

#include <fstream>

using namespace System;
using Misc::format;

void cmd_nerve(int argc, char **argv)
{
	std::ostringstream ss;
	ss << time(NULL);
	std::string timed_seed = ss.str();

	Argv C = read_arguments(argc, argv,
		Option({0, "h", "help", "false",
			"print help on the use of this program."}),

		Option({Option::VALUED | Option::CHECK, "i", "id", date_string(),
			"identifier for filenames."}),

		Option({Option::VALUED | Option::CHECK, "d", "dim", "2",
			"number of dimensions in which to run the code. "
			"Usually this should by 2 or 3."}),

		Option({Option::VALUED | Option::CHECK, "N", "box", "256",
			"size of the mass-box."}),

		Option({Option::VALUED | Option::CHECK, "L", "size", "100",
			"physical size of the simulation in units of Mpc/h."}),

		Option({Option::VALUED | Option::CHECK, "", "txt", "",
			"takes a list of triangles and plots those on the specified box. "
			"Parameter gives name of file to be read. Output to <id>.nerve.conan. "}));

	if (C.get<bool>("help"))
	{
		std::cout << "Cosmic workset Conan, by Johan Hidding.\n\n";
		C.print(std::cout);
		exit(0);
	}

	if (C["txt"] != "")
	{
		std::ifstream fi(C["txt"]);
		std::ofstream fo(format(C["id"], ".nerve.conan"));

		switch (dim)
		{
			case 2: Nerve::from_txt<2>(std::make_shared<Box<2>(N, L), fi, fo);
			case 3: Nerve::from_txt<3>(std::make_shared<Box<3>(N, L), fi, fo);
		}
	}
}

Global<Command> _NERVE("nerve", cmd_nerve);

