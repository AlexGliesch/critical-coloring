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
#include "cons.h"
#include "cliques/mntshao.h"
#include "subgraph.h"
void add_partial(subgraph& h, int cur_sz, int target_sz, double alpha) {
  if (cur_sz >= target_sz) return;
  for (int i : h.ss)
    if (i != -1) h.deg[i] = -1;
  int j = cur_sz;
  while (j < target_sz) {
    auto bst = max_element(begin(h.deg), end(h.deg)) - begin(h.deg);
    reservoir_sampling rs;
    int c = bst;
    for (int i = 0; i < n; ++i)
      if (h.deg[i] != -1 and
          ff(double(h.deg[bst] - h.deg[i]) / double(h.deg[bst])) < ff(alpha) and
          rs.consider())
        c = i;
    h.ss[j++] = c;
    for (int i = 0; i < n; ++i)
      if (h.deg[i] != -1) {
        if (i == c)
          h.deg[i] = -1;
        else
          h.deg[i] += AM[i][c];
      }
  }
  h.update_all();
}
subgraph cons_add(int sz, int u, double alpha) {
  TIME_BLOCK("cons_add");
  subgraph h(sz);
  h.ss[0] = u;
  h.update_all();
  add_partial(h, 1, sz, alpha);
  return h;
}
subgraph cons_drop(int sz, double alpha) {
  TIME_BLOCK("cons_drop");
  subgraph h(ind_n);
  h.update_all();
  while ((int)h.ss.size() > sz) {
    int bst = -1;
    for (int i : h.ss)
      if (bst == -1 or h.deg[i] < h.deg[bst]) bst = i;
    reservoir_sampling rs;
    int c = bst;
    for (int i : h.ss)
      if (h.deg[i] != nli::infinity() and
          ff(double(h.deg[i] - h.deg[bst]) / h.deg[i]) < ff(alpha) and
          rs.consider())
        c = i;
    int pos = find(begin(h.ss), end(h.ss), c) - begin(h.ss);
    swap(h.ss[pos], h.ss.back());
    h.ss.pop_back();
    h.deg[c] = nli::infinity();
    for (int i : h.ss)
      if (h.deg[i] != nli::infinity() and i != c) h.deg[i] -= AM[i][c];
  }
  h.update_all();
  return h;
}
