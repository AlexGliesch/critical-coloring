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
#ifndef MYSYS
#define MYSYS 1
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

    namespace btdsatur {
using namespace std;

/* time information */
extern long seconds, microsecs;
extern struct rusage tmp;

extern int getrusage(int who, struct rusage* rusage);

#include <signal.h>
extern int setrlimit(), getrlimit();

/* useful to help link from output files to res files */
// extern int getpid();

// #ifndef linux
// extern void ungetc();
// 
// extern int printf();
// extern int fprintf();
// extern int scanf();
// extern int fscanf();
// 
// extern int fgetc();
// 
// extern int fread();
// extern void fflush();
// extern void fclose();
// extern void rewind();
// #endif

/* the following required to make qsort calling args
silent */
typedef int (*compfunc)(const void*, const void*);
} // namespace btdsatur
#endif
