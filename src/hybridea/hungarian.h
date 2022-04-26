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
/******************************************************************************/#ifndef HUNGARIAN_H
#define HUNGARIAN_H
#include <iostream>
#include <limits.h>
#include <stdlib.h>
#include <string>
#include <vector>

namespace hybridea {
using namespace std;

struct cell {
  int weight;
  bool lined, visible;
  // initialisation
  cell() {
    lined = false;
    visible = true;
    weight = -1;
  }
};

void doHungarian(int n, vector<vector<int>>& matrix, vector<int>& matching);
} // namespace hybridea
#endif // HUNGARIAN_H
