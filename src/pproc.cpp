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
#include "pproc.h"
#include "color.h"
void mark_trivially_critical(const vi& ss, const vi& color, vb& crit, int k) {
  vi in_ss(n, -1);
  vi color_seen(k, -1);
  vi num_col_forced(k, 0);
  vi col_forced(k, -1);
  for (uint i = 0; i < ss.size(); ++i)
    in_ss[ss[i]] = i;
  for (int i = 0; i < (int)ss.size(); ++i) {
    ++num_col_forced[color[i]];
    for (int j : AL[ss[i]])
      if (in_ss[j] >= 0) color_seen[color[in_ss[j]]] = i;
    bool can_change = false;
    for (int j = 0; j < k; ++j)
      if (j != color[i] and color_seen[j] != i) {
        can_change = true;
        --num_col_forced[color[i]];
        break;
      }
    if (not can_change) col_forced[color[i]] = i;
  }
  int num_crits_found = 0;
  for (int i = 0; i < k; ++i)
    if (num_col_forced[i] == 1) {
      assert(inrange(col_forced[i], 0, (int)ss.size() - 1));
      assert(inrange(ss[col_forced[i]], 0, n - 1));
      if (not crit[ss[col_forced[i]]]) ++num_crits_found;
      crit[ss[col_forced[i]]] = true;
    }
  if (num_crits_found > 0) {
    stats::trivial_crits_found += num_crits_found;
    if (verb >= 3)
      pr("Marked {} new vertices as critical. Total: {}\n", num_crits_found,
         accumulate(begin(crit), end(crit), 0));
  }
}
int choose_v_sun(const vi& ss, const vb& crit) {
  int b = -1, b_score = 0;
  reservoir_sampling rs;
  for (int i = 0; i < (int)ss.size(); ++i) {
    assert(inrange(ss[i], 0, n - 1));
    if (not crit[ss[i]]) {
      int i_score = 0;
      for (int j : ss) {
        assert(inrange(j, 0, n - 1));
        if (AM[ss[i]][j]) i_score += (crit[j] ? m : 1);
      }
      if (b == -1 or i_score < b_score or (i_score == b_score and rs.consider())) {
        if (i_score != b_score) rs.reset();
        b = i, i_score = b_score;
      }
    }
  }
  return b;
}
int choose_lexi(const vi& ss, const vb& crit) {
  for (int i = 0; i < (int)ss.size(); ++i)
    if (not crit[ss[i]]) return i;
  return -1;
}
pair<bool, bool> reduce_subset_one_by_one(int k, vi& ss, bool chroma_k_bef, timer t) {
  vi color;
  vb surely_crit(n, false);
  auto color_ok = [&]() {
    if (color.size() != ss.size()) return false;
    if ((int)surely_crit.size() != n) return false;
    for (uint i = 0; i < color.size(); ++i)
      if (not inrange(color[i], 0, k - 1)) return false;
    return true;
  };
  check_colorability(k, ss, t, &color);
  if (color_ok()) mark_trivially_critical(ss, color, surely_crit, k);
  if (verb >= 2)
    pr("Trying to reduce subset of size {}, chroma_k_bef: {}...\n", ss.size(),
       chroma_k_bef);
  if ((int)ss.size() == k) return mp(true, true);
  bool chroma_k = chroma_k_bef;
  bool crit = true;
  while (not t.timed_out()) {
    int i = choose_v_sun(ss, surely_crit);
    if (i == -1) break;
    int v = ss[i];
    assert(not surely_crit[v]);
    swap(ss[i], ss.back());
    ss.pop_back();
    auto [colorable, sure] = check_colorability(k - 1, ss, t, &color);
    if (t.timed_out()) break;
    crit = crit and sure;
    if (colorable) {
      bool ok = color_ok();
      if (ok) mark_trivially_critical(ss, color, surely_crit, k - 1);
      surely_crit[v] = true;
      ss.push_back(v);
      swap(ss[i], ss.back());
      ++i;
      if (not ok) {
        ++stats::color_not_ok;
        continue;
      }
    } else {
      chroma_k = sure;
      if (verb >= 1)
        pr("Removed {} (chroma_k: {}), reduced ss size {}->{}\n", v, sure, ss.size() + 1,
           ss.size());
      update_global_best(ss, chroma_k, false);
    }
  }
  return mp(chroma_k, crit);
}
tuple<vi, bool, bool> pproc(const vi& in, bool chroma_gek, int k, timer t) {
  if (no_postproc) {
    return tuple(in, chroma_gek, false);
  }
  timer ptm;
  vi best = in, ss = in;
  auto [chroma_gek2, surely_crit] = reduce_subset_one_by_one(k, ss, chroma_gek, t);
  if (mp(not chroma_gek, ss.size()) < mp(not chroma_gek, best.size())) {
    best = move(ss);
    chroma_gek = chroma_gek2;
  }
  update_global_best(best, chroma_gek, surely_crit);
  stats::pp_reduced += int(in.size() - best.size());
  stats::pp_time += ptm.elapsed_secs();
  return mt(best, chroma_gek, surely_crit);
}
