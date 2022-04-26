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
/******************************************************************************/#ifndef MAKESOLUTION_H
#define MAKESOLUTION_H

#include "graph.h"
#include <iostream>
#include <limits.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace hybridea {
using namespace std;

int generateInitialK(Graph& g, int alg, vector<int>& bestColouring);
void makeInitSolution(Graph& g, vector<int>& sol, int k, int verbose);

void prettyPrintSolution(vector<vector<int>>& candSol);
void checkSolution(vector<vector<int>>& candSol, Graph& g, int verbose);

} // namespace hybridea
#endif
