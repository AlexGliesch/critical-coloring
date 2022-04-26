/******************************************************************************/
//  This code implements the Hybrid Evolutionary Algorithm of Galinier and Hao.
//  The local search routines are based on the tabu search algorithm written
//  by Ivo Bloechliger http://rose.epfl.ch/~bloechli/coloring/
//  The remaining code was written by R. Lewis www.rhydLewis.eu
//
//	See: Lewis, R. (2015) A Guide to Graph Colouring: Algorithms and
// Applications. Berlin, Springer.
//       ISBN: 978-3-319-25728-0. http://www.springer.com/us/book/9783319257280
//
//	for further details
/******************************************************************************/#pragma once
#include <vector>
#include "../util.h"
#include "graph.h"

namespace hybridea {
using namespace std;
bool solIsOptimal(const vector<int>& sol, Graph& g, int k);
void makeAdjList(int** neighbors, Graph& g);
void freeAdjList(int**, int);
void replace(vector<vector<int>>& population, vector<int>& parents,
						 vector<int>& osp, vector<int>& popCosts, Graph& g, int oCost);

// #TODO remove "res" later, if we're not using it
int hea(Graph& g, timer tm, unsigned long long maxChecks = 100000000,
				int targetCols = 2, int randomSeed = 1, int popSize = 10,
				int maxIterations = 16, int verbose = 0, int constructiveAlg = 1,
				int xOverType = 1, bool measuringDiversity = false,
				std::vector<int>* res = nullptr);

}	 // namespace hybridea
