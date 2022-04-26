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
/******************************************************************************/#ifndef MANIPULATEARRAYS_INCLUDED
#define MANIPULATEARRAYS_INCLUDED

#include "graph.h"
#include <iostream>
#include <vector>

namespace hybridea {
using namespace std;

void initializeArrays(int**& nodesByColor, int**& conflicts, int**& tabuStatus,
                      int*& nbcPosition, Graph& g, vector<int>& c, int k);

void freeArrays(int**& nodesByColor, int**& conflicts, int**& tabuStatus,
                int*& nbcPosition, int k, int n);

void moveNodeToColorForTabu(int bestNode, int bestColor, Graph& g,
                            vector<int>& c, int** nodesByColor, int** conflicts,
                            int* nbcPosition, int** neighbors,
                            int* nodesInConflict, int* confPosition,
                            int** tabuStatus, long totalIterations,
                            int tabuTenure);
} // namespace hybridea
#endif
