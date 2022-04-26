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
/******************************************************************************/#ifndef GraphIncluded
#define GraphIncluded

namespace hybridea {
class Graph {

public:
  Graph();
  Graph(int n);
  ~Graph();

  void resize(int n);

  int* matrix;
  int matrixSize;
  int n;       // number of nodes
  int nbEdges; // number of edges

  int* operator[](int index);
};

} // namespace hybridea
#endif
