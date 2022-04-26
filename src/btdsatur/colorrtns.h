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
#ifndef COLORRTNDEF
#define COLORRTNDEF

#include "graph.h"

    /* COLOR STORAGE STRUCTURES */
    namespace btdsatur {
typedef unsigned short colortype;

struct clrinfo {
  colortype numcolors; /* number of colors used */
  int total;           /* sum of colors (weight)  used */
};
typedef struct clrinfo clrinfotype;

struct vrtxandclr {
  vertextype vertex;
  colortype color;
};
typedef struct vrtxandclr vrtxandclrtype;

struct popmember {
  clrinfotype clrdata;
  vrtxandclrtype vc[MAXVERTEX];
};
typedef struct popmember popmembertype;

/* COLOR MANIPULATION ROUTINES */
extern void printinfo(popmembertype* member);
/* print number of colors and color sum */

extern void printcoloring(popmembertype* member);
/* print the colors of vertices in order */

extern void printpurity(popmembertype* member);
/* the purity computation for closeness to hidden color */

extern void getcolorinfo(popmembertype* member);
/* computes the maximum and sum of colors in member */

extern void permute(popmembertype* member, vertextype first, vertextype last);
/* randomly permute the vertex ordering of member in the
range first to last -1 */

extern void trivial_color(popmembertype* m);
/* apply the colors 1 to order to the vertices of m,
whatever order they are in */

extern void verifycolor(popmembertype* m);
/* verify the coloring of m as to correctness and print
and appropriate message */

/* print a permutation member - for debugging mostly*/
extern void printperm(popmembertype* m);

extern void getacoloring(popmembertype* m, char* name, int* which);
/* open the indicated file and obtain the coloring asked for by the user */

/* PROPERTIES OF GRAPH */

extern int degseq[];

extern int decdeg(vrtxandclrtype* a, vrtxandclrtype* b);
/* comparison routine for decreasing sort by degree */

extern int incdeg(vrtxandclrtype* a, vrtxandclrtype* b);
/* comparison routine for increasing sort by degree */

extern void computedeg();
/* Compute degree sequence.  */

/* FINAL OUTPUT TO RESULTS FILE */
extern void fileres(/*char* name, */popmembertype* m/*, char* prgm*/);
/* name will be appended with .res, coloring data will be appended to
file name.res, prgm is name of program  */

extern void about(char* pgrmname);
} // namespace btdsatur
#endif
