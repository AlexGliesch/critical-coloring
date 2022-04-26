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
/******************************************************************************/#ifndef TABU_INCLUDED
#define TABU_INCLUDED
#pragma GCC system_header

#include "graph.h"
#include <vector>

namespace hybridea {
using namespace std;

int tabu(Graph& g, vector<int>& c, int k, int maxIterations, int verbose,
         int** neighbors);
} // namespace hybridea
#endif
