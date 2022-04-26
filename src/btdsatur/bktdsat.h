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
#ifndef BKTDSATDEF
#define BKTDSATDEF

#include "../util.h"
#include "colorrtns.h"

namespace btdsatur {
#define ENDLIST 65535

extern void bktdsat(popmembertype* m, int branch, colortype targetclr, int min,
                    int max, timer t);

/* use abacktrack version of dsatur, with
limited branching and forbidden branch ranges */

/* interrupt  control */
extern void cpulimbk();
extern void intfcbk();

} // namespace btdsatur
#endif
