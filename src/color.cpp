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
#include "color.h"
#include "btdsatur/bkdmain.h"
#include "hybridea/main.h"
#include "main.h"
void setup_btdsatur(const vi& ss) {
  srand(rand_int(nli::min(), nli::max()));
  btdsatur::order = ss.size();
  btdsatur::maxChecks = 100000000000000LL;
  btdsatur::verbose = 0;
  memset(btdsatur::graph, 0, BTDSATUR_GRAPHSIZE);
  for (int i = 0; i < (int)ss.size(); ++i)
    for (int j = i + 1; j < (int)ss.size(); ++j)
      if (AM[ss[i]][ss[j]]) {
        btdsatur::setedge(i, j);
        btdsatur::setedge(j, i);
      }
}
hybridea::Graph& setup_hea(const vi& ss) {
  static hybridea::Graph g;
  if (g.n == 0) {
    g.resize(n);
  }
  g.n = ss.size();
  g.nbEdges = 0;
  fill(g.matrix, g.matrix + g.n * g.n, 0);
  for (int i = 0; i < (int)ss.size(); ++i)
    for (int j = i + 1; j < (int)ss.size(); ++j) {
      g[i][j] = g[j][i] = AM[ss[i]][ss[j]];
      ++g.nbEdges;
    }
  return g;
}
bool is_k_colorable_exact(int k, const vi& ss, timer t, vi* res) {
  setup_btdsatur(ss);
  if (t.secs_left() <= 0) return true;
  double tm = t.elapsed_secs();
  bool suc;
  try {
    if (res) res->clear();
    int colors = btdsatur::colorsearch(k, t, res);
    if (verb >= 3) pr("k: {}, |ss|: {}, colors: {}\n", k, ss.size(), colors);
    suc = colors <= k;
    stats::suc_ext_cals += suc;
    stats::unsuc_ext_cals += !suc;
  } catch (timeout_exception& e) {
    suc = true;
  }
  stats::color_time += t.elapsed_secs() - tm;
  return suc;
}
bool is_k_colorable_heuristic(int k, const vi& ss, timer t, vi* res) {
  auto& g = setup_hea(ss);
  if (t.secs_left() <= 0) return true;
  double tm = t.elapsed_secs();
  bool suc;
  try {
    int rnd_seed = rand_int(nli::min(), nli::max());
    unsigned long long num_checks = 100000000000000LL;
    if (res) res->clear();
    int colors = hybridea::hea(g, t, num_checks, k, rnd_seed, 10, 16, 0 ,
                               1, 1, false, res);
    if (verb >= 3) pr("k: {}, |ss|: {}, colors: {}\n", k, ss.size(), colors);
    suc = colors <= k;
    stats::suc_heu_cals += suc;
    stats::unsuc_heu_cals += !suc;
  } catch (timeout_exception& e) {
    suc = true;
  }
  stats::color_time += t.elapsed_secs() - tm;
  return suc;
}
bb check_colorability(int k, const vi& ss, timer t, vi* res) {
  TIME_BLOCK("check_colorability");
  ++stats::cals_to_coloring;
  static int min_size_exact_times_out = nli::max();
  if (not no_exact_coloring and min_size_exact_times_out > (int)ss.size()) {
    timer exact_timer(exact_secs, t);
    if (verb >= 3) pr("Running exact algorithm on size {}\n", ss.size());
    bool exactly_colorable = is_k_colorable_exact(k, ss, exact_timer, res);
    if (exact_timer.timed_out()) {
      min_size_exact_times_out = ss.size();
      if (verb >= 2)
        pr("Exact timed out, min_size_exact_times_out: {}\n",
           min_size_exact_times_out);
    } else {
      return mp(exactly_colorable, true);
    }
  }
  const int num_heu_reruns = 0;
  if (not no_heuristic_coloring)
    for (int j = 0; j < 1 + num_heu_reruns and not t.timed_out(); ++j)
      if (is_k_colorable_heuristic(k, ss, timer(heu_secs, t), res)) {
        if (j > 0) ++stats::heu_mistk;
        if (verb >= 3)
          pr("Actually, I was mistaken: it is indeed {}-colorable.\n", k);
        return mp(true, true);
      }
  return mp(t.timed_out(), false);
}
bb is_k_vcs(int k, const vi& ss, timer t, vi* res) {
  auto r = check_colorability(k - 1, ss, t, res);
  if (no_heuristic_coloring and r.second == false) {
    r.first = false;
  } else {
    r.first = not r.first;
  }
  return r;
}
int color_exactly(const vi& ss, timer t, int lb) {
  setup_btdsatur(ss);
  if (t.secs_left() <= 0) return true;
  try {
    int colors = btdsatur::colorsearch(lb, t, nullptr);
    return colors;
  } catch (timeout_exception& e) {
    return -1;
  }
}
int color_heuristically(const vi& ss, timer t, int lb) {
  auto& g = setup_hea(ss);
  try {
    int rnd_seed = rand_int(nli::min(), nli::max());
    unsigned long long num_checks = 100000000000000LL;
    int colors = hybridea::hea(g, t, num_checks, lb, rnd_seed, 10, 16, 0 ,
                               1, 1, false, nullptr);
    return colors;
  } catch (timeout_exception& e) {
    return -1;
  }
}
