/* This C++ source code of a DSatur-based backtracking algorithm for graph
  colouring has been adapted from the C code of Joseph Culberson and Dennis Papp
  which is available from Joseph Culberson's Coloring Page
  http://web.cs.ualberta.ca/~joe/Coloring/lclindex.html (Copyright (c) 1997
  Joseph Culberson. All rights reserved.)

  Points to note:
  * The program accepts as input a DIMACS formatted graph colouring test file
  * Other optional command line settings can be viewed by running the programm
  from the command line with no arguments
  * Some of the parameters in this implementation have been fixed. See the
  publication:

  See:
  Lewis, R. (2015) A Guide to Graph Colouring: Algorithms and Applications.
  Berlin, Springer. ISBN: 978-3-319-25728-0.
  http://www.springer.com/us/book/978331925728
  for further details
*/
#pragma once
#include "../util.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-result"
#pragma GCC diagnostic ignored "-Wreturn-type"
#include "bktdsat.h"
#include "colorrtns.h"
#include "graph.h"
#include "mysys.h"
#include <vector>

namespace btdsatur {
/* If res != nullptr, populates res with the n colors */
int colorsearch(int targetnumcolors, timer t, std::vector<int>* res = nullptr);
extern unsigned long long numConfChecks;
extern unsigned long long maxChecks;
extern int verbose;
} // namespace btdsatur
#pragma GCC diagnostic pop
