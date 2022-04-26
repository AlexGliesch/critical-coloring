/**
 * \file gc.cpp
 *   \author Marcus Ritt <ritt@informatik.uni-tuebingen.de>
 *   \version $Id: gc.cpp 8854 2018-05-05 21:06:12Z ritt $
 *   \date Time-stamp: <2018-05-05 18:05:28 ritt>
 *
 * Coloring graphs after Coudert (1997). The coloring is a simple backtracking
 * using the DSATUR heuristics (select the vertex with most colors in its
 * neighborhood, break ties by largest number of uncolored neighbors).
 *
 * The main observation of the article is: most graphs are 1-perfect, i.e. the
 * size of maximum clique equals the chromatic number. This makes finding are
 * large maximum clique important, and he proposes an improved algorithm for
 * that.
 *
 * Modified by Alex Gliesch in May 2018.
 */
//#include "gc.H"
#include "coudert.h"

namespace coudert {

// forward declarations
vector<VertexInformation::NodeID> maxCliqueRec(
		const Graph& G, vector<VertexInformation::NodeID>& C,
		vector<VertexInformation::NodeID>& B, unsigned ub, unsigned level);

// remove the first vertex which is non-adjacent to v in G
// complexity O(V)
bool remove_first_non_adjacent(Graph& G, Node& v) {
	assert(v < num_vertices(G));

	graph_traits<Graph>::vertex_iterator vp, ve;
	for (tie(vp, ve) = vertices(G); vp != ve; vp++)
		if (*vp != v && !is_adjacent(G, v, *vp)) {
			clear_vertex(*vp, G);
			remove_vertex(*vp, G);
			return true;
		}
	return false;
}

// remove all vertices of G, except the neighborhood N(v) of node v
//
// a horrible operation; there must be a better way
// our problem is that every removal invalidates all node and edge descriptors
// so after each removal we have to find our "hook" again; this depends
// on the implementation with an adjacency_list and vecS,vecS
// the complexity is also a horrible O(V^2)
void keep_neighborhood(Graph& G, Node v) {
	VertexInformation::NodeID id = G[v].id;
	while (remove_first_non_adjacent(G, v)) {
		// now v is invalid: find it again. the dirty part: we know it's v or v-1.
		// even more dirty: G[v].id == id, even when v is invalid!
		if (G[v].id != id || !(v < num_vertices(G))) v--;
		assert(G[v].id == id && v < num_vertices(G));
	}
	clear_vertex(v, G);
	remove_vertex(v, G);
}

// determine the node of maximum degree of graph G
Node maxDegree(const Graph& G) {
	assert(num_vertices(G) > 0);
	graph_traits<Graph>::vertex_iterator vp, ve;
	Node r = 0;
	unsigned max_degree = 0;
	double max_central = 0;

	for (tie(vp, ve) = vertices(G); vp != ve; vp++) {
		if (degree(*vp, G) > max_degree ||
				(degree(*vp, G) == max_degree && G[*vp].centrality > max_central)) {
			max_degree = degree(*vp, G);
			max_central = G[*vp].centrality;
			r = *vp;
		}
	}
	return r;
}

// find the number of used colors in the neighborhood of v
pair<unsigned, unsigned> saturation(const Graph& G, const vector<Color>& color,
																		Node v) {
	set<Color> c;
	unsigned uncolored = 0;
	graph_traits<Graph>::adjacency_iterator up, ue;
	for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
		if (color[*up] != 0)
			c.insert(color[*up]);
		else
			uncolored++;
	return make_pair(c.size(), uncolored);
}

// determine the uncolored node of maximum saturation number of graph G
Node maxSatur(const Graph& G, vector<Color>& color) {
	assert(num_vertices(G) > 0);
	assert(num_vertices(G) <= color.size());

	graph_traits<Graph>::vertex_iterator vp, ve;
	Node r = 0;
	int max_satur = -1;
	unsigned max_uncolored = 0;
	double max_central = 0;

	for (tie(vp, ve) = vertices(G); vp != ve; vp++) {
		if (color[*vp] == 0) {
			int satur;
			unsigned uncolored;
			tie(satur, uncolored) = saturation(G, color, *vp);
			if (satur > max_satur ||
					(satur == max_satur && uncolored > max_uncolored) ||
					(satur == max_satur && uncolored == max_uncolored &&
					 G[*vp].centrality > max_central)) {
				max_satur = satur;
				max_uncolored = uncolored;
				max_central = G[*vp].centrality;
				r = *vp;
			}
		}
	}
	return r;
}

// sequential coloring of graph G without backtracking
// pre-condition: color[i]==0 in color array
// post-condition color[i]!=0 (all vertices colored)
// returns: number of colors used
unsigned color_H3(const Graph& G, vector<Color>& color) {
	assert(color.size() >= num_vertices(G));
	assert(all_of(color.begin(), color.begin() + num_vertices(G),
								[](unsigned c) { return c == 0; }));

	unsigned colored = 0;
	unsigned k = 0;

	// color all vertices sequentially without backtracking
	while (colored < num_vertices(G)) {
		// Node v = find(color.begin(), color.end(), 0u)-color.begin();
		// pick uncolored vertex using DSATUR heuristic, i.e. with maximal
		// saturation
		Node v = maxSatur(G, color);

		// try all colors
		for (unsigned c = 1; c <= k + 1; c++) {
			// check if c has no conflicts
			graph_traits<Graph>::adjacency_iterator up, ue;
			for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
				if (color[*up] == c) break;

			if (up != ue) continue;

			color[v] = c;
			colored++;
			k = max(c, k);
			break;
		}
	}

	// sanity checks and post-conditions
	assert(/*0 <= k && */ k <= num_vertices(G));
	assert(all_of(color.begin(), color.begin() + num_vertices(G),
								[](unsigned c) { return c != 0; }));
	return k;
}

unsigned color_H3(const Graph& G) {
	vector<Color> color(num_vertices(G));
	return color_H3(G, color);
}

// reduction rule B: a vertex whose degree+current clique size is smaller than
// best can't contribute
bool rule_B(Graph& G, int limit) {
	graph_traits<Graph>::vertex_iterator vp, ve;
	for (tie(vp, ve) = vertices(G); vp != ve; vp++)
		if (int(degree(*vp, G)) < limit) {
			clear_vertex(*vp, G);
			remove_vertex(*vp, G);
			return true;
		}
	return false;
}

// reduction rule C: a vertex with high degree must be part of the clique
bool rule_C(Graph& G, vector<VertexInformation::NodeID>& C) {
	graph_traits<Graph>::vertex_iterator vp, ve;
	for (tie(vp, ve) = vertices(G); vp != ve; vp++)
		if (degree(*vp, G) >= num_vertices(G) - 2) {
			C.push_back(G[*vp].id);
			keep_neighborhood(G, *vp);
			return true;
		}
	return false;
}

// reduction rule Q (Theorem 1)
bool rule_Q(Graph& G, const vector<Color>& color, const unsigned k,
						const int limit) {
	graph_traits<Graph>::vertex_iterator vp, ve;
	for (tie(vp, ve) = vertices(G); vp != ve; vp++)
		if (int(k) - int(saturation(G, color, *vp).first) > limit) {
			clear_vertex(*vp, G);
			remove_vertex(*vp, G);
			return true;
		}
	return false;
}

unsigned maxClique_backtracks;

// determine the maximum clique of graph G
vector<VertexInformation::NodeID> maxClique(const Graph& G) {
	TIME_BLOCK("coudert::maxClique");
	maxClique_backtracks = 0;
	vector<VertexInformation::NodeID> C, B;
	return maxCliqueRec(G, C, B, numeric_limits<unsigned>::max(), 0);
}

// recursive helper function to find maximum clique of graph G
// G: current rest graph
// C: current clique (vector of node IDs)
// B: best clique found so far (vector of node IDs)
// ub: current upper bound
vector<VertexInformation::NodeID> maxCliqueRec(
		const Graph& G, vector<VertexInformation::NodeID>& C_,
		vector<VertexInformation::NodeID>& B, unsigned ub, unsigned level) {
	Graph G_(G);
	vector<VertexInformation::NodeID> C(C_);

	unsigned k;
	vector<Color> color(num_vertices(G));	// for the q-rule
	do {
		// (1) all vertices done? C is our current best
		if (num_vertices(G_) == 0) return C;

		// (2) criterion (a)
		if (C.size() + num_vertices(G_) <= B.size()) return B;

		// (3) color G somehow to get an upper bound for the rest of the clique
		fill(color.begin(), color.end(), 0);
		k = color_H3(G_, color);

		// (4) update upper bound: addition to clique isn't larger than # colors k
		ub = min(ub, unsigned(C.size() + k));

		// cut: upper bound lower than current best? stop here
		if (ub <= B.size()) return B;

		// (5) apply reduction rules (b), (c), and q-pruning
	} while (rule_B(G_, int(B.size()) - int(C.size())) || rule_C(G_, C) ||
					 rule_Q(G_, color, k, int(C.size()) - int(B.size()) + int(k)));

	// Note: rule (d) costs too much; rule (e) has not been implemented

	unsigned Bold = B.size();

	// select vertex of maximum degree
	Node v = maxDegree(G_);

	// force v into clique
	{
		C.push_back(G_[v].id);
		Graph G1(G_);

		assert(num_vertices(G1) == num_vertices(G_));
		assert(G_[v].id == G1[v].id);	// ensure properties have been copied
		// let G1 be the graph induced by N(v); TBD: this is quick and _very_ dirty
		keep_neighborhood(G1, v);

		assert(num_vertices(G1) ==
					 degree(v, G_));	// ensure the resulting graph
														// equals the number of neighbors
		B = maxCliqueRec(G1, C, B, ub, level + 1);
		C.pop_back();
	}

	if (ub == B.size()) return B;

	// exclude v
	Graph G0(G_);
	clear_vertex(v, G0);
	remove_vertex(v, G0);
	B = maxCliqueRec(G0, C, B, ub, level + 1);
	// count # backtracks, definition: last child fails, i.e. no child contributed
	if (B.size() == Bold) maxClique_backtracks++;
	return B;
}

unsigned seqColorRec(Graph& G, unsigned k, unsigned B, const unsigned lb,
										 vector<Color>& color, unsigned colors, unsigned& back,
										 unsigned m, timer tm);

pair<unsigned, unsigned> seqColor(Graph& G,
																	const vector<VertexInformation::NodeID>& C,
																	unsigned m, timer tm) {
	TIME_BLOCK("coudert::seqColor");
	unsigned chi = 0;
	unsigned back = 0;

	vector<Color> color(num_vertices(G));

	// color the clique C
	for (unsigned i = 0; i < C.size(); i++)
		color[C[i]] = i + 1;

	// and go
	chi = seqColorRec(G, C.size(), num_vertices(G) + 1, C.size(), color, C.size(),
										back, m, tm);
	return make_pair(chi, back);
}

// partially colored graph G
// using k colors until now
// best solution until now uses B colors
// lb is a lower bound for the number of colors
// color is the color array, colors the number of colors used
// m: the algorithm stops at m colors
unsigned seqColorRec(Graph& G, unsigned k, unsigned B, const unsigned lb,
										 vector<Color>& color, unsigned colored, unsigned& back,
										 unsigned m, timer tm) {
	assert(num_vertices(G) == color.size());
	if (tm.timed_out()) throw logic_error("timed out");
	if (B <= m) return B;

	// test if G is entirely colored: then this is a solution with k colors
	if (num_vertices(G) == colored) {
		if (k < B) {
			// store coloring
			graph_traits<Graph>::vertex_iterator vp, ve;
			for (tie(vp, ve) = vertices(G); vp != ve; vp++)
				G[*vp].color = color[*vp];
		}
		return k;
	}

	// hit local lower bound: we are done, too.
	if (k >= B) return B;

	// pick any uncolored vertex; here: pick the first
	// TBD: Compare with DSATUR heuristic (see maxSatur)
	// unsigned v = find(color.begin(), color.end(), 0u)-color.begin();
	Node v = maxSatur(G, color);
	G[v].rank = colored + 1;
	// cout << "Chosen " << v << " had sat." << saturation(G,color,v).first << "
	// and " << saturation(G,color,v).second << " uncolored neighbors." << endl;

	// try all colors
	for (Color c = 1; c <= min(k + 1, B - 1); c++) {
		// check if c has no conflicts
		graph_traits<Graph>::adjacency_iterator up, ue;
		for (tie(up, ue) = adjacent_vertices(v, G); up != ue; up++)
			if (color[*up] == c) break;

		if (up != ue) continue;

		color[v] = c;
		// recursive call: B always decreases
		// if there is no coloring with fewer than best, the neighborhood check
		// fails always unsigned B0 = B;
		B = seqColorRec(G, max(c, k), B, lb, color, colored + 1, back, m, tm);
		// if (B<B0)
		//   cout << "L" << colored << " " << B << " " << B0 << " " << lb << endl;
		color[v] = 0;

		// hit some lower bound, or timed out: we are done
		if (tm.timed_out()) throw logic_error("timed out");
		if (lb == B || k >= B || B <= m) return B;
	}
	back++;
	return B;
}

// usage: gc <filename> <onlyclique>?
// int main(int argc, char *argv[]) {
//   assert(1<argc && argc<4);
//
//   std::chrono::system_clock::time_point start =
//   std::chrono::system_clock::now();
//
//   // (1) input the graph
//   ifstream f(argv[1]);
//   Graph G = readDIMACS(f);
//
//   // (2) compute centrality
//   brandes_betweenness_centrality(G,get(&VertexInformation::centrality,G));
//
//   // (3) find the maximum clique
//   vector<VertexInformation::NodeID> C = maxClique(G);
//   for (auto c : C)
//     G[c].clique=true;
//
//   // report: name, |V|, |E|, |C|
//   cout << filesystem::path(argv[1]).filename().string() << " " <<
//   num_vertices(G) << " " << num_edges(G) << " " << C.size()  << " " << flush;
//   if (argc==3) {
//     // report time, backtracks
//     cout <<
//     std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start).count()
//     << " " << maxClique_backtracks << " "; cout << "[";
//     copy(C.begin(),C.end(),ostream_iterator<VertexInformation::NodeID>(cout,"
//     ")); cout << "]" << endl;
//     // (3.1) write it, for control
//     ofstream fG(filesystem::path(argv[1]).stem().string()+".dot");
//     write_graphviz(fG, G, vertex_writer(G,vertex_writer::queen));
//     fG.close();
//     exit(0);
//   }
//
//   // (4) color the graph sequentially
//   pair<unsigned,unsigned> r = seqColor(G,C);
//   unsigned chi = r.first;
//   unsigned back= r.second;
//
//   // (5) report: colors, backtracks, time
//   cout << chi << " " << back << " " <<
//   std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now()-start).count()
//   << endl; ofstream fG(filesystem::path(argv[1]).stem().string()+".dot");
//   write_graphviz(fG, G, vertex_writer(G,vertex_writer::queen),
//   default_writer(), graph_writer(G)); fG.close();
// }

}	// namespace coudert
