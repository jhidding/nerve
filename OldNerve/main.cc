#include <iostream>
#include <cstdlib>

#include "conan/boxconfig.h"
#include "conan/getopt.h"
#include "conan/history.h"
#include "nerve.h"

using namespace Conan;
using namespace std;

void usage()
{
	cout << "usage: nerve -input <filename> -output <filename>\n";
}

int main(int argc, char **argv)
{
	GetOpt &options = GetOpt::initialize(argc, argv);

	if (options.has("-help"))
	{
		usage();
		exit(0);
	}

	if (not (options.has("-input") and options.has("-output")))
	{
		usage();
		exit(1);
	}

	string fn_input = options.get("-input");
	string fn_output = options.get("-output");

	if (options.has("-txt"))
	{	
		ifstream fi;
		fi.open(fn_input.c_str(), ios::in | ios::binary);

		unsigned 	rank = 2;
		unsigned	N;
		float		L;

		fi >> N;
		fi >> L;

		BoxConfig::initialize(rank, N, L);

		ofstream fo;
		fo.open(fn_output.c_str(), ios::out | ios::binary);
		
		from_txt<2>(fi, fo);

		fi.close();
		fo.close();
	}
	else
	{
		ifstream fi;
		fi.open(fn_input.c_str(), ios::in | ios::binary);

		Header 		H(fi);
		History 	I(fi);

		unsigned	rank 	= from_string<unsigned>(H["rank"]);
		unsigned	N	= from_string<unsigned>(H["N"]);
		float		L 	= from_string<float>(H["L"]);

		BoxConfig::initialize(rank, N, L);

		ofstream fo;
		fo.open(fn_output.c_str(), ios::out | ios::binary);

		H.to_file(fo);
		I.to_file(fo);

		try {
			switch (rank)
			{
				//case 1:	nerve<1>(fi, fo); break;
				case 2:	nerve<2>(fi, fo); break;
				case 3:	nerve<3>(fi, fo); break;
			}
		}
		catch (char const *msg)
		{
			cerr << msg << endl;
			exit(1);
		}

		fi.close();
		fo.close();
	}

	return 0;
}

