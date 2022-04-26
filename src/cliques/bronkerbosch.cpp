#include "bronkerbosch.h"
#include "main.h"
#include "subgraph.h"

namespace {	// Local variables
vi res;
bool sure;
}	// namespace

// #TODO can probably pre-allocate R,P and X since they'll have at most n
// elements at all times. Especially R seems easy to deal with.
unsigned long long calls_to_bk_impl = 0;
void bk_impl(vi R, vi P, vi X, timer t, int max_k) {
	++calls_to_bk_impl;
	if (R.size() + P.size() <= res.size()) return;	// Simple lower bound
	if (P.empty() and X.empty()) {
		if (R.size() > res.size()) {
			// pr("max_clique: {}\n", R.size());
			res = move(R);
		}
		return;
	}
	if (t.timed_out()) {
		sure = false;
		return;
	}
	if (P.empty()) return;
	if ((int)res.size() >= max_k) return;
	assert(P.size() > 0);
	int u = P[0];	// Pivot. Better choice?

	for (int i = P.size() - 1; i >= 0; --i) {
		int v = P[i];
		if (not AM[u][v]) {
			vi Ra = R, Pa, Xa;
			Ra.push_back(v);
			for (int i : P)
				if (AM[i][v]) Pa.push_back(i);
			for (int i : X)
				if (AM[i][v]) Xa.push_back(i);
			bk_impl(move(Ra), move(Pa), move(Xa), t, max_k);
			X.push_back(v);
			// swap(P[i], P.back());
			// P.pop_back();
			P.erase(P.begin() + i);
		}
	}
}

vi get_degeneracy_ordering(subgraph g) {
	vi ord;
	while ((int)ord.size() != g.n) {
		int pos = -1;
		for (int i = 0; i < (int)g.ss.size(); ++i)
			if (pos == -1 or g.deg[g.ss[i]] < g.deg[g.ss[pos]]) pos = i;
		int v = g.ss[pos];
		ord.push_back(v);
		swap(g.ss[pos], g.ss.back());
		g.ss.pop_back();
		g.deg[v] = nli::infinity();
		for (int i : g.ss)
			g.deg[i] -= AM[i][v];
	}
	return move(ord);
}

pair<vi, bool> bron_kerbosch(const subgraph& h, timer t, int max_k) {
	TIME_BLOCK("bron_kerbosch");
	vi P = get_degeneracy_ordering(h), X;
	reverse(begin(P), end(P));
	sure = true;
	res.clear();

	while (P.size() and (int) res.size() < max_k) {
		if (t.timed_out()) {
			sure = false;
			break;
		}
		int i = P.back();
		// pr("ord {}, {} left\n", P[i], P.size());
		vi R, Pa, Xa;
		R.push_back(i);
		for (int j : P)
			if (AM[i][j]) Pa.push_back(j);
		for (int j : X)
			if (AM[i][j]) Xa.push_back(j);
		bk_impl(R, Pa, Xa, t, max_k);
		P.pop_back();
		X.push_back(i);
	}
	// pr("calls_to_bk_impl: {}\n", calls_to_bk_impl);
	return mp(move(res), sure);
}
