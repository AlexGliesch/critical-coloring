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
/******************************************************************************/#ifndef PETURB_H
#define PETURB_H

#include "graph.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <vector>

namespace hybridea {
using namespace std;
void doRandomPeturbation(vector<int>& osp, int k, Graph& g);
} // namespace hybridea
#endif
