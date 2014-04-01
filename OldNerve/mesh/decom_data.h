#pragma once
#include <vector>
#include <map>
#include "misc.h"

struct Decom_data_2
{
	typedef unsigned char byte;
	typedef std::vector<byte> Vec;
	typedef std::map<byte, byte> Map;

	struct Cell_structure
	{
		double unit_volume;
		byte vertices[3];
	};

	class Cell
	{
		public:
		double unit_volume;
		Vec vertices;
		Vec edges;
	};

	class Edge: public Map
	{
		public:
		Vec vertices;
	};

	class Vertex: public Map
	{
		public:
		Vec cells;
	};

	std::vector<Cell> C;
	std::vector<Edge> E;
	std::vector<Vertex> V;

	Decom_data_2(Cell_structure data[2]);
};

struct Decom_data_3
{
	typedef unsigned char byte;
	typedef std::vector<byte> Vec;
	typedef std::map<byte, byte> Map;

	struct Cell_structure
	{
		double unit_volume;
		byte vertices[4];
	};

	class Cell
	{
		public:
		double unit_volume;
		Vec vertices;
		Vec edges;
		Vec faces;
	};

	class Face: public Map
	{
		public:
		Vec vertices;
		Vec edges;
	};

	class Edge: public Map
	{
		public:
		Vec vertices;
		Vec cells;
	};

	class Vertex: public Map
	{
		public:
		Vec faces;
		Vec cells;
	};

	std::vector<Cell> C;
	std::vector<Face> F;
	std::vector<Edge> E;
	std::vector<Vertex> V;

	Decom_data_3(Cell_structure data[5]);
};

