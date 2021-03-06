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
/******************************************************************************/#include "inputgraph.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace hybridea {
using namespace std;

void inputDimacsGraph(Graph& g, char* file) {

  char c;
  char str[400];
  ifstream IN(file, ios::in);
  int line = 0;
  g.nbEdges = 0;
  int edges = -1;
  int blem = 1;
  int multiple = 0;
  while (!IN.eof()) {
    line++;
    IN.get(c);
    if (IN.eof()) break;
    switch (c) {
    case 'p':
      IN.get(c);
      IN.getline(str, 39, ' ');
      if (strcmp(str, "edge") && strcmp(str, "edges")) {
        cerr << "File " << file << " line " << line << ":\n";
        cerr << "Error reading 'p' line: no 'edge' keyword found.\n";
        cerr << "'" << str << "' found instead\n";
        exit(-1);
      }
      IN >> g.n;
      IN >> edges;
      blem = 0;
      g.resize(g.n);
      break;
    case 'n':
      if (blem) {
        cerr << "File " << file << " line " << line << ":\n";
        cerr << "Found 'n' line before a 'p' line.\n";
        exit(-1);
      }
      int node;
      IN >> node;
      if (node < 1 || node > g.n) {
        cerr << "File " << file << " line " << line << ":\n";
        cerr << "Node number " << node << " is out of range!\n";
        exit(-1);
      }
      node--;
      cout << "Tags (n Lines) not implemented in g object\n";
      break;
    case 'e':
      int node1, node2;
      IN >> node1 >> node2;
      if (node1 < 1 || node1 > g.n || node2 < 1 || node2 > g.n) {
        cerr << "File " << file << " line " << line << ":\n";
        cerr << "Node " << node1 << " or " << node2 << " is out of range!\n";
        exit(-1);
      }
      node1--;
      node2--;
      // Undirected graph
      // cout << "Read edge from " << node1 << " to " << node2 << endl;
      if (g[node1][node2] == 0) {
        g.nbEdges++;
      } else {
        multiple++;
        if (multiple < 5) {
          cerr << "Warning: in graph file " << file << " at line " << line
               << ": edge is defined more than once.\n";
          if (multiple == 4) {
            cerr << "  No more multiple edge warnings will be issued\n";
          }
        }
      }
      g[node1][node2] = 1;
      g[node2][node1] = 1;
      break;
    case 'd':
    case 'v':
    case 'x':
      cerr << "File " << file << " line " << line << ":\n";
      cerr << "'" << c << "' lines are not implemented yet...\n";
      IN.getline(str, 399, '\n');
      break;
    case 'c':
      IN.putback('c');
      IN.get(str, 399, '\n');
      break;
    default:
      cerr << "File " << file << " line " << line << ":\n";
      cerr << "'" << c << "' is an unknown line code\n";
      exit(-1);
    }
    IN.get(); // Kill the newline;
  }
  IN.close();
  if (multiple) {
    cerr << multiple << " multiple edges encountered\n";
  }
}

} // namespace hybridea
