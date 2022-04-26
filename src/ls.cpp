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
#include "ls.h"
#include "subgraph.h"
void ls_impl(subgraph& h, int tenure, int maxnonimpr, int pmin, int pmax,
             int pstep) {
  static vi Vmh;
  Vmh.clear();
  sort(begin(ind_n), end(ind_n));
  sort(begin(h.ss), end(h.ss));
  assert(is_sorted(begin(h.ss), end(h.ss)));
  set_difference(begin(ind_n), end(ind_n), begin(h.ss), end(h.ss),
                 back_inserter(Vmh));
  assert((int)Vmh.size() == n - h.n);
  tabu_list tabu(n, tenure);
  int nonimpr = 0;
  int pcur = pmin;
  bool improved = true;
  const bool skip_brooks_pruning = irace_its;
  while (improved) {
    static int i = -1, j = -1;
    if (not inrange(i, 0, (int)Vmh.size() - 1))
      i = rand_int(0, (int)Vmh.size() - 1);
    if (not inrange(j, 0, h.n - 1)) j = rand_int(0, h.n - 1);
    int bst_m = -1, bst_i = -1, bst_j = -1;
    reservoir_sampling rs;
    improved = false;
    for (uint ict = 0; ict < Vmh.size(); ++ict, i = (i + 1) % Vmh.size()) {
      if (improved) break;
      int v = Vmh[i];
      if (tabu.is_tabu(v)) continue;
      if (not skip_brooks_pruning and h.deg[v] <= k - 2) continue;
      for (int jct = 0; jct < h.n; ++jct, j = (j + 1) % h.n) {
        if (improved) break;
        int u = h.ss[j];
        if (not skip_brooks_pruning and h.deg[v] - AM[v][u] <= k - 2) continue;
        if (tabu.is_tabu(u)) continue;
        int m_try = h.m_cost_swap(j, v);
        if (bst_m == -1 or bst_m < m_try) {
          bst_m = m_try, bst_i = i, bst_j = j;
          rs.reset();
        } else if (bst_m == m_try and rs.consider()) {
          bst_m = m_try, bst_i = i, bst_j = j;
        }
        improved = bst_m > h.m;
      }
    }
    auto do_move = [&](int i, int j) {
      assert(not linear_in(h.ss, Vmh[i]));
      int old = h.ss[j];
      tabu.add(old);
      tabu.add(Vmh[i]);
      tabu.advance_iter();
      h.swap(j, Vmh[i]);
      Vmh[i] = old;
      ++stats::tot_ls_moves;
    };
    if (bst_m == -1) break;
    if (improved) {
      assert(not linear_in(h.ss, Vmh[bst_i]));
      assert(inrange(bst_j, 0, (int)h.ss.size() - 1));
      do_move(bst_i, bst_j);
      nonimpr = 0;
    } else {
      ++nonimpr;
      if (nonimpr < maxnonimpr) {
        do_move(bst_i, bst_j);
      } else {
        pcur = pcur + pstep;
        if (pcur > pmax) break;
        for (int p = 0; p < pcur; ++p)
          do_move(rand_int(0, (int)Vmh.size() - 1),
                  rand_int(0, (int)h.ss.size() - 1));
      }
    }
  }
}
void ls(subgraph& h) {
  TIME_BLOCK("ls");
  ls_impl(h,
          n * tenure_mult,
          max_nonimpr,
          pmin,
          pmax_mult * n,
          pstep
  );
}
