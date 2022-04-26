/*
* A new heuristic for finding verifiable k-vertex-critical subgraphs
* 
* Copyright (c) 2022 Alex Gliesch, Marcus Ritt
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPY lRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
#include "subgraph.h"
#include "cons.h"
#include "ls.h"
subgraph& subgraph::update_all() {
  n = ss.size();
  deg.assign(::n, 0);
  for (int i = 0; i < ::n; ++i)
    for (int j = 0; j < n; ++j)
      if (inrange(ss[j], 0, ::n - 1)) deg[i] += AM[i][ss[j]];
  m = 0;
  for (int v : ss)
    if (inrange(v, 0, ::n - 1)) m += deg[v];
  m /= 2;
  return *this;
}
void subgraph::swap(int i, int v) {
  int old = ss[i];
  ss[i] = v;
  assert(inrange(old, 0, ::n - 1) and inrange(v, 0, ::n - 1));
  for (int j = 0; j < ::n; ++j)
    deg[j] = deg[j] - AM[old][j] + AM[v][j];
  m = 0;
  for (int v : ss)
    if (inrange(v, 0, ::n - 1)) {
      m += deg[v];
    }
  m /= 2;
}
void subgraph::set_ss(vi ss2) {
  ::swap(ss, ss2);
  update_all();
}
