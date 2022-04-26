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
#include "bkdmain.h"
#include <fstream>
#include <limits.h>
#include <map>
#include <unordered_map>

namespace btdsatur {

extern colortype bestcolor;

// This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft
// Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

// Global Variables
unsigned long long numConfChecks;
unsigned long long maxChecks;
int verbose;
////ofstream timeStream, checksStream;
clock_t startTime;

int colorsearch(int targetnumcolors, timer t, std::vector<int>* res) {
  numConfChecks = 0;
  partitionflag = 0;
  cheatflag = 0;
  popmembertype m;
  int i, branch, min, max;
  colortype targetcolor;
  //   char info[256];

  for (i = 0; i < order; i++)
    m.vc[i].vertex = i;
  computedeg();
  qsort((char*)m.vc, (int)order, sizeof(struct vrtxandclr), (compfunc)decdeg);

  // These paramters are fixed in this version (in Culberson's original C code
  // they can be altered)
  targetcolor = targetnumcolors;

  // ADDITIONAL PARAMETERS
  // Entering a branching factor of 0 causes the algorithm to behave like
  // DSATUR, essentially performing a sequential search. For larger values, when
  // the program backtracks to some point and takes an alternate branch, the
  // branching factor is reduced by one for the entire subtree. If the branch
  // factor at the root of a subtree is 0 then no further branching is allowed
  // in the subtree. If the branching factor is set very large, and no backtrack
  // region is excluded (see below) then this algorithm guarantees an optimal
  // coloring for all graphs (given excess time)
  branch = INT_MAX;

  // For some graphs improvements might only occur when branching is permitted
  // at certain leaves of the search tree. Entering a pair of numbers such as
  // min=30, max=280 on a graph of 300 vertices means that no branching will
  // occur at depths in that range. Here, min=max=1 meaning branching is
  // permitted at all levels
  min = 1;
  max = 1;

  // Now start the timer,
  startTime = clock();

  // And go to the backtracking algorithm itself
  bktdsat(&m, branch, targetcolor, min, max, t);
  //   getcolorinfo(&m);
  //   printinfo(&m);

  // Print the final lines to the output files for consistency
  //// timeStream<<"1\tX\n";
  //// checksStream<<"1\tX\n";

  // verify things and end
  //// verifycolor(&m);
  //// fileres(name, &m,info);
  //   cout << bestcolor << endl;
  if (res != nullptr) {
    res->resize(order);
    std::unordered_map<int,int> color_map;
    for (int i=0; i<order; ++i) {
      (*res)[m.vc[i].vertex] = m.vc[i].color;
      if (color_map.count(m.vc[i].color)==0) 
        color_map[m.vc[i].color]=color_map.size();      
    }
    for (int& i : *res)
      i = color_map[i];
  }
  return bestcolor;
}

int mainDSatur(int argc, char* argv[]) {

  if (argc <= 1) {
    cout << "Backtracking DSatur Algorithm for Graph Colouring\n\n"
         << "USAGE:\n"
         << "<InputFile>     (Required. File must be in DIMACS format)\n"
         << "-s <int>        (Stopping criteria expressed as number of "
            "constraint checks. Can be anything up to 9x10^18. DEFAULT = "
            "100,000,000.)\n"
         << "-r <int>        (Random seed. DEFAULT = 1)\n"
         << "-T <int>        (Target number of colours. Algorithm halts if "
            "this is reached. DEFAULT = 1.)\n"
         << "-v              (Verbosity. If present, output is sent to screen. "
            "If -v is repeated, more output is given.)\n"
         << "****\n";
    exit(1);
  }

  // The following variables are used for keeping track of the number of
  // constraint checks
  numConfChecks = 0;

  // Read in the graph and runtime parameters and/or set default values
  char filename[200];
  int seed = 1, targetNumCols = 2, i;
  verbose = 0;
  maxChecks = 100000000;
  for (i = 1; i < argc; i++) {
    if (strcmp("-T", argv[i]) == 0) {
      targetNumCols = atoi(argv[++i]);
    } else if (strcmp("-r", argv[i]) == 0) {
      seed = atoi(argv[++i]);
    } else if (strcmp("-s", argv[i]) == 0) {
      maxChecks = strtoull(argv[++i], NULL, 10);
    } else if (strcmp("-v", argv[i]) == 0) {
      verbose++;
    } else {
      strcpy(filename, argv[i]);
      // 			cout<<"Backtracking DSatur Algorithm using <"<<filename<<">\n\n";
    }
  }

  srand(seed);
  //   getgraph(filename);

  // Open the output streams
  //// timeStream.open("teffort.txt"); checksStream.open("ceffort.txt");
  ////	if (timeStream.fail() || checksStream.fail()){cout << "ERROR OPENING
  /// output FILE";exit(1);}

  // Produce some output
  if (verbose >= 1) cout << " COLS     CPU-TIME(ms)\tCHECKS" << endl;

  // This is the main algorithm bit
  startTime = clock();

  colorsearch(/*filename, */ targetNumCols, timer());

  // checksStream.close();
  // timeStream.close();

  return (0);
}
} // namespace btdsatur
