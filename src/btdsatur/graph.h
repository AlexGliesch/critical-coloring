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
#ifndef GRAPHDEFS
#define GRAPHDEFS

namespace btdsatur {
/*
!!!WARNING!!!!
MAXVERTEX must be divisible by 8
*/
#define MAXVERTEX 8000

#define SHIFT 3
#define MASK 7
#define ROWSIZE ((MAXVERTEX >> SHIFT) + 1)
#define BTDSATUR_GRAPHSIZE (MAXVERTEX * ROWSIZE)

/* Definitions useful for checking and setting edges */

/*                     ***NOTE***
set and clear are asymmetric - use setedge(i,j) setedge(j,i) etc.
*/
#define setedge(i, j) graph[((i)*ROWSIZE) + ((j) >> SHIFT)] |= (1 << ((j)&MASK))
#define clearedge(i, j)                                                        \
  graph[((i)*ROWSIZE) + ((j) >> SHIFT)] &= ~(1 << ((j)&MASK))

#define edge(i, j) (graph[((i)*ROWSIZE) + ((j) >> SHIFT)] & (1 << ((j)&MASK)))

/* for loops involving potential neighbors */
#define initnbr(x, i) (x) = graph + ((i)*ROWSIZE)
#define isnbr(x, i) (((x)[(i) >> SHIFT]) & (1 << ((i)&MASK)))

typedef int vertextype;
typedef unsigned char adjacencytype;

extern adjacencytype graph[BTDSATUR_GRAPHSIZE];
extern vertextype order;

/* CHEAT INFORMATION FROM GRAPH */
extern int partset[MAXVERTEX];
extern int partitionflag;
extern int partitionnumber;

extern int cheatflag;

extern void printgraph();
extern void getgraph(char a[]);

} // namespace btdsatur
#endif
