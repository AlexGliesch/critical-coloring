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
/******************************************************************************/#include "graph.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace hybridea {
using namespace std;

Graph::Graph() {
  matrix = NULL;
  n = matrixSize = 0;
  nbEdges = 0;
}

Graph::Graph(int m) {
  matrix = NULL;
  resize(m);
}

int* Graph::operator[](int index) {
  if (index < 0 || index >= this->n) {
    cout << "First node index out of range: " << index << "\n";
    exit(EXIT_FAILURE);
    // matrix[-1] = 0; // Make it crash.
  }
  return matrix + this->n * index;
}

void Graph::resize(int m) {
  if (matrix != NULL) {
    delete[] matrix;
  }
  if (m > 0) {
    n = m;
    matrixSize = m;
    nbEdges = 0;
    matrix = new int[m * m];
    for (int i = 0; i < m * m; i++) {
      matrix[i] = 0;
    }
  }
}

Graph::~Graph() { resize(0); }

} // namespace hybridea
