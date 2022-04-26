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
#pragma once       
#include "stats.h"
#include "util.h"
inline string input_filename;
inline string output_filename;
inline double time_limit_secs;
inline int k;
inline size_t random_seed;
inline int verb;
inline double mu;
inline int R;
inline double xi;
inline int max_iter;
inline double heu_secs;
inline double exact_secs;
inline string cons_alg;
inline double tenure_mult;
inline int max_nonimpr;
inline int pmin;
inline double pmax_mult;
inline int pstep;
inline double clique_alg_time_1st;
inline double cons_alpha;
inline double confirm_crit_timelimit;
inline bool do_confirm_criticality;
inline bool do_force_confirm;
inline bool print_summary_each_iter;
inline bool no_heuristic_coloring;
inline bool irace_its;
inline bool irace;
inline bool just_exact_coloring;
inline bool just_heuristic_coloring;
inline bool normal_run;
inline int imax;
inline bool no_exact_coloring;
inline bool no_fst_phase;
inline bool no_postproc;
inline bool no_dense_search;
inline bool just_pproc;
inline string instance_name;
inline vi ind_n;
inline vi ind_deg;
inline vi ind_deg_cum;
inline int n_ori, m_ori;
inline vvi AM_ori;
inline vvi AL_ori;
inline int n, m;
inline vvi AM;
inline vvi AL;
inline vi vmap;
inline vi best_fin;
inline bool best_fin_chroma_k = true;
inline bool best_fin_crit = false;
inline vi best_gen;
inline bool best_gen_sure = true;
inline bool did_postproc = false;
inline timer global_timer;
inline int global_iter_last_improve = 0;
inline void update_global_best(const vi& ss, bool chroma_k, bool crit) {
  if (((int)best_fin.size() == n and
       inrange((int)ss.size(), 1, (int)best_fin.size() - 1)) or
      mt(not chroma_k, not crit, ss.size()) <
          mt(not best_fin_chroma_k, not best_fin_crit, best_fin.size())) {
    stats::ttb = global_timer.elapsed_secs();
    stats::itb = stats::global_iter;
    stats::max_iter_diff =
        max(stats::max_iter_diff, stats::global_iter - global_iter_last_improve);
    global_iter_last_improve = stats::global_iter;
    best_fin = ss, best_fin_chroma_k = chroma_k, best_fin_crit = crit;
    if (verb >= 1)
      pr("Updated global best ({}), chroma_k: {}, crit: {}\n", ss.size(),
         best_fin_chroma_k, best_fin_crit);
    sort(begin(best_fin), end(best_fin));
  }
  if (print_summary_each_iter) {
    stats::time = global_timer.elapsed_secs();
    stats::print_summary();
  }
}
inline void sort_by_degree(vi& ss) {
  vi deg(n, 0);
  shuffle(begin(ind_n), end(ind_n), rng);
  for (uint i = 0; i < ss.size(); ++i) {
    for (uint j = 0; j < ss.size(); ++j)
      deg[ss[i]] += AM[ss[i]][ss[j]];
  }
  sort(begin(ss), end(ss),
       [&](int i, int j) { return tie(deg[i], ind_n[i]) > tie(deg[j], ind_n[j]); });
}
inline int count_edges(const vi& ss) {
  int m = 0;
  for (uint i = 0; i < ss.size(); ++i)
    for (uint j = i + 1; j < ss.size(); ++j)
      m += AM[ss[i]][ss[j]];
  return m;
}
