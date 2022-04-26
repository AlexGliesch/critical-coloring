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
/******************************************************************************/
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#include "main.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include "diversity.h"
#include "graph.h"
#include "inputgraph.h"
#include "makesolution.h"
#include "tabu.h"
#include "xover.h"

// This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft
// Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

namespace hybridea {

using namespace std;

// int ASS = INT_MIN;
// bool solIsOptimal(vector<int>& sol, Graph& g, int k);
// void makeAdjList(int** neighbors, Graph& g);
// void replace(vector<vector<int>>& population, vector<int>& parents,
//              vector<int>& osp, vector<int>& popCosts, Graph& g, int oCost);

unsigned long long numConfChecks;

//*********************************************************************
inline bool solIsOptimal(const vector<int>& sol, Graph& g, int k) {
  int i, j;
  for (i = 0; i < g.n; i++) {
    if (sol[i] < 1 || sol[i] > k) {
      cout << "Error: Colour of node " << i << " out of bounds " << endl;
      exit(1);
    }
  }

  for (i = 0; i < (g.n) - 1; i++) {
    for (j = i + 1; j < g.n; j++) {
      if (sol[i] == sol[j] && g[i][j]) return (false);
    }
  }
  // set<int> cols;
  // for (int i:sol) cols.insert(i);
  // cout <<"ok, solIsOptimal. colors: " << cols.size() << endl;
  // If we are here then we have established a solution with k or less colours
  return (true);
}

int hea(Graph& g, timer tm, unsigned long long maxChecks, int targetCols,
        int randomSeed, int popSize, int maxIterations, int verbose,
        int constructiveAlg, int xOverType, bool measuringDiversity,
        std::vector<int>* res) {
  int i, k;
  bool solFound = false, doKempeMutation = false;
  vector<int> parents;

  // This variable keeps count of the number of times information about the
  // instance is looked up
  numConfChecks = 0;

  // Set the number of parents in each crossover and decide if the Kempe
  // mutation is going to be used
  if (xOverType == 3)
    parents.resize(4);
  else
    parents.resize(2);
  if (xOverType == 2) doKempeMutation = true;

  // set tabucol limit
  maxIterations = maxIterations * g.n;
  if (targetCols < 2 || targetCols > g.n) targetCols = 2;

  // Now set up some output files
  // 	ofstream timeStream, confStream;
  // 	timeStream.open("teffort.txt"); confStream.open("ceffort.txt");
  // 	if (timeStream.fail() || confStream.fail()){cout << "ERROR OPENING output
  // FILE";exit(1);}

  // Do a check to see if we have the empty graph. If so, end immediately.
  if (g.nbEdges <= 0) {
    // 		confStream<<"1\t0\n0\tX\t0\n";
    // 		timeStream<<"1\t0\n0\tX\t0\n";
    if (verbose >= 1)
      cout << "Graph has no edges. Optimal solution is obviously using one "
              "colour. Exiting."
           << endl;
    // 		confStream.close();
    // 		timeStream.close();
    exit(1);
  }

  // Make the adjacency list structure
  static int** neighbors = nullptr;
  if (neighbors == nullptr) {
    neighbors = new int*[g.matrixSize];
    for (int i = 0; i < g.matrixSize; ++i)
      neighbors[i] = new int[g.matrixSize + 1];
  }
  makeAdjList(neighbors, g);

  /*int** neighbors = new int*[g.n];
  makeAdjList(neighbors, g);*/

  // Produce some output
  if (verbose >= 1) cout << " COLS     CPU-TIME\tCHECKS" << endl;

  // Seed and start timer
  clock_t clockStart = clock();
  srand(randomSeed);
  numConfChecks = 0;

  // Data structures used for population and offspring
  vector<vector<int>> population(popSize, vector<int>(g.n));
  vector<int> popCosts(popSize);
  vector<int> osp(g.n), bestColouring(g.n);

  auto getsol = [&](const vector<int>& c) {
    if (res != nullptr) {
      *res = c;
      std::unordered_map<int, int> color_map;
      for (int i = 0; i < g.n; ++i)
        if (color_map.count(c[i]) == 0) color_map[c[i]] = color_map.size();
      for (int& i : *res)
        i = color_map[i];
    }
  };

  // Generate the initial value for k using greedy or dsatur algorithm
  k = generateInitialK(g, constructiveAlg, bestColouring);
  getsol(bestColouring);

  //..and write the results to the output file
  int duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
  if (verbose >= 1)
    cout << setw(5) << k << setw(11) << duration << "ms\t" << numConfChecks
         << " (via constructive)" << endl;
  // 	confStream<<k<<"\t"<<numConfChecks<<"\n";
  // 	timeStream<<k<<"\t"<<duration<<"\n";
  if (k <= targetCols) {
    if (verbose >= 1)
      cout << "\nSolution with  <=" << targetCols
           << " colours has been found. Ending..." << endl;
    // 		confStream<<"1\t"<<"X"<<"\n";
    // 		timeStream<<"1\t"<<"X"<<"\n";
  }

  // MAIN ALGORITHM
  int bestK = k;
  k--;
  while (numConfChecks < maxChecks && k + 1 > targetCols) {
    solFound = false;

    // First build the population
    for (i = 0; i < popSize; i++) {
      // Build a solution using modified DSatur algorithm
      makeInitSolution(g, population[i], k, verbose);
      // Check to see whether this solution is alrerady optimal or if the cutoff
      // point has been reached. If so, we end
      if (solIsOptimal(population[i], g, k)) {
        solFound = true;
        for (int j = 0; j < g.n; j++)
          osp[j] = population[i][j];
        // getsol(osp);
        break;
      }
      if (numConfChecks >= maxChecks || tm.timed_out()) {
        for (int j = 0; j < g.n; j++)
          osp[j] = population[i][j];
        // getsol(osp);
        break;
      }
      // Improve each solution via tabu search and record their costs
      popCosts[i] = tabu(g, population[i], k, maxIterations, 0, neighbors);
      // Check to see whether this solution is now optimal or if the cuttoff
      // point is reached. If so, we end
      if (verbose >= 2)
        cout << "          -> Individual " << setw(4) << i
             << " constructed. Cost = " << popCosts[i] << endl;
      if (popCosts[i] == 0) {
        solFound = true;
        for (int j = 0; j < g.n; j++)
          osp[j] = population[i][j];
        // getsol(osp);
        break;
      }
      if (numConfChecks >= maxChecks || tm.timed_out()) {
        for (int j = 0; j < g.n; j++)
          osp[j] = population[i][j];
        // getsol(osp);
        break;
      }
    }

    // Now evolve the population
    int rIts = 0, oCost = 1, best = INT_MAX;
    while (numConfChecks < maxChecks && !tm.timed_out() && !solFound) {
      // Choose parents and perform crossover to produce a new offspring
      doCrossover(xOverType, osp, parents, g, k, population);

      // Improve the offspring via tabu search and record its cost
      oCost = tabu(g, osp, k, maxIterations, 0, neighbors);

      // Write osp over weaker parent and update popCosts
      replace(population, parents, osp, popCosts, g, oCost);

      if (verbose >= 2) {
        cout << "          -> Offspring " << setw(5) << rIts
             << " constructed. Cost = " << oCost;
        if (measuringDiversity)
          cout << "\tDiversity = " << measureDiversity(population, k);
        cout << endl;
      }

      rIts++;

      if (oCost < best) best = oCost;
      if (oCost == 0) solFound = true;
    }

    // Algorithm has finished at this k
    duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
    if (solFound) {
      bestK = min(bestK, k);
      if (verbose >= 1)
        cout << setw(5) << k << setw(11) << duration << "ms\t" << numConfChecks
             << endl;
      // 			confStream<<k<<"\t"<<numConfChecks<<"\n";
      // 			timeStream<<k<<"\t"<<duration<<"\n";
      // Copy the current solution as the best solution
      for (int i = 0; i < g.n; i++)
        bestColouring[i] = osp[i] - 1;
      getsol(osp);

      if (k <= targetCols) {
        if (verbose >= 1)
          cout << "\nSolution with  <=" << targetCols
               << " colours has been found. Ending..." << endl;
        break;
      }
    } else {
      if (verbose >= 1)
        cout << "\nRun limit exceeded. No solution using " << k
             << " colours was achieved (Checks = " << numConfChecks << ", "
             << duration << "ms)" << endl;
      // 			confStream<<k<<"\tX\t"<<numConfChecks<<"\n";
      // 			timeStream<<k<<"\tX\t"<<duration<<"\n";
    }

    k--;
  }
  // freeAdjList(neighbors, g.n);
  // 	ofstream solStrm;
  // 	solStrm.open("solution.txt");
  // 	solStrm<<g.n<<"\n";
  // 	for(int i=0;i<g.n;i++)solStrm<<bestColouring[i]<<"\n";
  // 	solStrm.close();
  //   cout << bestK << endl;
  return bestK;
}

void usage() {
  cout
      << "Hybrid EA for Graph Colouring\n\n"
      << "USAGE:\n"
      << "<InputFile>     (Required. File must be in DIMACS format)\n"
      << "-t <int>        (Time limit, in seconds. DEFAULT = 3600 seconds.)\n"
      << "-s <int>        (Stopping criteria expressed as number of constraint "
         "checks. Can be anything up to 9x10^18. DEFAULT = 100,000,000.)\n"
      << "-I <int>        (Number of iterations of tabu search per cycle. This "
         "figure is multiplied by the graph size |V|. DEFAULT = 16)\n"
      << "-r <int>        (Random seed. DEFAULT = 1)\n"
      << "-T <int>        (Target number of colours. Algorithm halts if this "
         "is reached. DEFAULT = 1.)\n"
      << "-v              (Verbosity. If present, output is sent to screen. If "
         "-v is repeated, more output is given.)\n"
      << "-p <int>        (Population Size. Should be 2 or more. DEFAULT = "
         "10)\n"
      << "-a <int>        (Choice of construction algorithm to determine "
         "initial value for k. DSsatur = 1, Greedy = 2. DEFAULT = 1.)\n"
      << "-x <int>        (Crossover Operator. 1 = GPX (2 parents)\n"
      << "                                     2 = GPX (2 parents + Kempe "
         "chain mutation)\n"
      << "                                     3 = MPX (4 parent crossover "
         "with q=2)\n"
      << "                                     4 = GGA (2 parents)\n"
      << "                                     5 = nPoint (2 parents)\n"
      << "                 DEFAULT = 1)\n"
      << "-d              (If present population diversity is measured after "
         "each crossover)"
      << "****\n";
  exit(1);
}

int hea_main(int argc, char** argv) {
  if (argc <= 1) {
    usage();
  }

  int timeLimitSeconds = 3600;

  Graph g;
  int popSize = 10, maxIterations = 16, verbose = 0, randomSeed = 1,
      constructiveAlg = 1, targetCols = 2, xOverType = 1;
  bool measuringDiversity = false;
  unsigned long long maxChecks = 100000000;

  for (int i = 1; i < argc; i++) {
    if (strcmp("-p", argv[i]) == 0) {
      popSize = atoi(argv[++i]);
      if (popSize < 2) {
        cout << "Error: PopSize should be >= 2\n";
        exit(1);
      }
    } else if (strcmp("-a", argv[i]) == 0) {
      constructiveAlg = atoi(argv[++i]);
    } else if (strcmp("-I", argv[i]) == 0) {
      maxIterations = atoi(argv[++i]);
    } else if (strcmp("-r", argv[i]) == 0) {
      randomSeed = atoi(argv[++i]);
    } else if (strcmp("-T", argv[i]) == 0) {
      targetCols = atoi(argv[++i]);
    } else if (strcmp("-v", argv[i]) == 0) {
      verbose++;
    } else if (strcmp("-s", argv[i]) == 0) {
      maxChecks = strtoull(argv[++i], NULL, 10);
    } else if (strcmp("-x", argv[i]) == 0) {
      xOverType = atoi(argv[++i]);
    } else if (strcmp("-d", argv[i]) == 0) {
      measuringDiversity = true;
    } else if (strcmp("-t", argv[i]) == 0) {
      timeLimitSeconds = atoi(argv[++i]);
    } else {
      // 			cout<<"Hybrid Evolutionary Algorithm using <"<<argv[i]<<">\n\n";
      inputDimacsGraph(g, argv[i]);
    }
  }

  hea(g, timer(timeLimitSeconds), maxChecks, targetCols, randomSeed, popSize,
      maxIterations, verbose, constructiveAlg, xOverType, measuringDiversity);

  return (0);
}

//*********************************************************************
void makeAdjList(int** neighbors, Graph& g) {
  // Makes the adjacency list corresponding to G
  // for (int i = 0; i < g.n; i++) {
  //  neighbors[i] = new int[g.n + 1];
  // neighbors[i][0] = 0;
  //}
  for (int i = 0; i < g.n; i++) {
    neighbors[i][0] = 0;
    for (int j = 0; j < g.n; j++) {
      if (g[i][j] && i != j) {
        neighbors[i][++neighbors[i][0]] = j;
      }
    }
  }
}

void freeAdjList(int** neighbors, int n) {
  for (int i = 0; i < n; ++i)
    delete[] neighbors[i];
  delete[] neighbors;
}

//*********************************************************************
void replace(vector<vector<int>>& population, vector<int>& parents,
             vector<int>& osp, vector<int>& popCosts, Graph& g, int oCost) {
  // Go through the parents and ID the worst one
  int toDie, i, max = INT_MIN;
  for (i = 0; i < parents.size(); i++) {
    if (popCosts[parents[i]] > max) {
      max = popCosts[parents[i]];
      toDie = parents[i];
    }
  }
  // Copy osp over the parent selected toDie
  population[toDie] = osp;
  popCosts[toDie] = oCost;
}
} // namespace hybridea
#pragma GCC diagnostic pop
