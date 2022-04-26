#pragma once
#include "pre.h"
#include "util.h"

namespace coudert {

using namespace std;
using namespace boost;

struct VertexInformation {
	VertexInformation() : id(0), color(0) {}
	typedef int NodeID;
	VertexInformation(NodeID _id) : id(_id), color(0) {}
	NodeID id;
	unsigned color;	// unused
	double centrality;
	bool clique;
	unsigned rank;	// rank in coloring order
};

typedef adjacency_list<vecS, vecS, undirectedS, VertexInformation> Graph;
typedef graph_traits<Graph>::vertex_descriptor Node;
typedef unsigned Color;

// read a graph in DIMACS standard format
// based on "Clique and coloring problems graph format", revision of May 8, 1993
inline Graph readDIMACS(ifstream& f) {
	char type;
	string s;
	Graph G;

	while (!f.eof()) {
		f >> type;
		if (f.eof()) break;
		switch (type) {
			case 'c':
				getline(f, s);
				// comment: skip
				break;
			case 'p':
				// allocate graph
				unsigned n, e;
				f >> s >> n >> e;
				for (unsigned i = 0; i < n; i++)
					add_vertex(VertexInformation(i), G);
				break;
			case 'n':
				// node description, ignore for now
				getline(f, s);
				break;
			case 'e':
				// edge description
				unsigned src, dst;
				f >> src >> dst;
				if (!is_adjacent(G, src - 1, dst - 1)) add_edge(src - 1, dst - 1, G);
				break;
			default:
				throw("File format error");
		}
	}
	return G;
}

class vertex_writer {
	const Graph& G;

 public:
	enum GraphType { none, queen };
	GraphType t;

	vertex_writer(const Graph& _G, GraphType _t = none) : G(_G), t(_t) {}
	template <typename VertexOrEdge>
	void operator()(ostream& out, const VertexOrEdge& v) const {
		out << "[";
		if (G[v].clique)
			out << "shape=box";
		else
			out << "shape=circle";
		out << ",label=\"" << G[v].rank
				<< "\",fontcolor=\"/white\",fillcolor=" << G[v].color;
		if (t == queen) {
			const unsigned scale = 60;
			unsigned n = sqrt(num_vertices(G));
			out << ",pos=\"" << scale * (G[v].id % n) << "," << scale * (G[v].id / n)
					<< "\"";
		}
		out << "]";
	}
};

struct graph_writer {
	const Graph& G;

 public:
	graph_writer(const Graph& _G) : G(_G) {}
	void operator()(std::ostream& out) const {
		out << "graph [outputorder=edgesfirst]" << std::endl;
		out << "node  [colorscheme=spectral11,style=filled]" << std::endl;
	}
};

vector<VertexInformation::NodeID> maxClique(const Graph& G);
pair<unsigned, unsigned> seqColor(Graph& G,
																	const vector<VertexInformation::NodeID>& C,
																	unsigned m, timer tm);

}	// namespace coudert
