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
#include "main.h"
#include "cliques/mcqd.h"
#include "cliques/mntshao.h"
#include "color.h"
#include "cons.h"
#include "hybridea/main.h"
#include "ls.h"
#include "pproc.h"
#include "readall.h"
#include "sm.h"
#include "stats.h"
#include "subgraph.h"
namespace {
bool exited = false;
vi cur_ss_gen;
vi random_bfs(int sz) {
  int start = rand_int(0, n - 1);
  vector<int> q;
  vb visited(n, false);
  q.push_back(start);
  visited[start] = true;
  vi ans;
  while (q.size()) {
    int ind = rand_int(0, (int)q.size() - 1);
    int v = q[ind];
    q[ind] = q.back();
    q.pop_back();
    ans.push_back(v);
    if ((int)ans.size() >= sz) return ans;
    for (int u : AL[v])
      if (not visited[u]) {
        visited[u] = true;
        q.push_back(u);
      }
  }
  return random_sample(sz,ind_n);
}
}
pair<vi, bool> find_k_vcs_fixed_size(int sz, int k, timer t) {
  if (sz == n) return mp(ind_n, true);
  ++stats::num_size_att;
  subgraph h;
  for (int iter = 1; iter <= R; ++iter) {
    if (no_dense_search) {
      vi ss = random_bfs(sz);
      auto [is_vcs, sure] = is_k_vcs(k, ss, t);
      if (is_vcs) {
        ++stats::suc_size_att;
        return mp(ss, sure);
      }
      continue;
    }
    int r = rand_int(0, 2 * m - 1);
    auto it = lower_bound(begin(ind_deg_cum), end(ind_deg_cum), r);
    int i = ind_deg[it - begin(ind_deg_cum)];
    if (t.timed_out()) break;
    ++stats::num_gen_subsets;
    if (cons_alg == "adddrop") {
      int lo = (k + 2) + (n - k + 2) * 0.25, hi = (k + 2) + (n - k + 2) * 0.75;
      if (sz < lo)
        h = cons_add(sz, i, cons_alpha);
      else if (sz > hi)
        h = cons_drop(sz, cons_alpha);
      else if (rand_int(lo, hi) < sz)
        h = cons_drop(sz, cons_alpha);
      else
        h = cons_add(sz, i, cons_alpha);
    } else if (cons_alg == "drop") {
      h = cons_drop(sz, cons_alpha);
    } else if (cons_alg == "add") {
      h = cons_add(sz, i, cons_alpha);
    }
    stats::tot_cons_edges += h.m;
    ls(h);
    ++stats::num_ls;
    sort(begin(h.ss), end(h.ss));
    auto [is_vcs, sure] = is_k_vcs(k, h.ss, t);
    if (is_vcs) {
      ++stats::suc_size_att;
      return mp(move(h.ss), sure);
    }
  }
  return mp(vi(), true);
}
pair<vi, bool> find_k_vcs(int k, timer t) {
  assert((int)ind_n.size() == n);
  for (int sz = k + 2, last = k; sz <= n and not t.timed_out();
       sz = min(n, int(sz * mu))) {
    if (verb >= 2) pr("Trying size = {}...\n", sz);
    vi ss;
    bool sure_chroma_k;
    tie(ss, sure_chroma_k) = find_k_vcs_fixed_size(sz, k, t);
    if (ss.size()) {
      if (verb >= 2) {
        pr("Found k-vcs of size = {} ({}sure_chroma_k)!\n", ss.size(),
           sure_chroma_k ? "" : "not ");
      }
      update_global_best(ss, sure_chroma_k, false);
      int lo = last, hi = sz - 1;
      while (sz != k and lo <= hi and not t.timed_out()) {
        int mid = (lo + hi) / 2;
        if (verb >= 2) pr("Trying size = {} (bs)...\n", mid);
        vi ss2;
        bool sure2;
        tie(ss2, sure2) = find_k_vcs_fixed_size(mid, k, t);
        if (ss2.size())
          hi = mid - 1;
        else
          lo = mid + 1;
        if (ss2.size() and ss2.size() < ss.size()) {
          swap(ss, ss2);
          sure_chroma_k = sure2;
          update_global_best(ss, sure_chroma_k, false);
          if (verb >= 2)
            pr("Reduced size to |ss| = {} ({}sure_chroma_k)!\n", ss.size(),
               sure_chroma_k ? "" : "not ");
        }
      }
      return mp(ss, sure_chroma_k);
    }
    last = sz;
  }
  return mp(ind_n, true);
}
pair<vi, bool> our_algorithm(int k, timer t) {
  unordered_map<int, pair<vi, bool>> iter_cache;
  unordered_map<vi, pair<vi, bool>> gen_cache;
  const bool use_gen_cache = false;
  using stats::global_iter;
  for (global_iter = 1;
       global_iter <= max_iter and not global_timer.timed_out() and not stats::infeas;
       ++global_iter) {
    vi ss;
    bool chroma_k = true;
    if (global_iter == 1) {
      if (no_fst_phase) {
        ss = best_gen = best_fin = ind_n;
        chroma_k = true;
        stats::gen_best_1st_phase = best_gen.size();
        stats::fin_best_1st_phase = best_fin.size();
      } else {
        TIME_BLOCK("iter = 1");
        tie(ss, chroma_k) = find_k_vcs(k, t);
        stats::gen_best_1st_phase = best_gen.size();
        stats::fin_best_1st_phase = best_fin.size();
        if ((int)ss.size() == n) break;
      }
    } else {
      if (global_iter - global_iter_last_improve > imax) break;
      TIME_BLOCK("iter >= 1");
      const int ub = min(n, (int)ceil(xi * (double)best_gen.size()));
      if (inrange(ub, k + 2, n)) {
        for (int sz = ub, i = 1; sz >= k + 2; --sz, ++i) {
          if (verb >= 2) pr("Trying size {} (iteration {}.{})...\n", sz, global_iter, i);
          auto it = iter_cache.find(sz);
          if (it != iter_cache.end()) {
            assert((int)it->second.first.size() == sz);
            tie(ss, chroma_k) = it->second;
            if (verb >= 2) pr("Already cached! ({}chroma_k)!\n", chroma_k ? "" : "not ");
          } else {
            vi ss2;
            bool sure2;
            tie(ss2, sure2) = find_k_vcs_fixed_size(sz, k, t);
            if (ss2.empty()) break;
            iter_cache[sz] = mp(ss2, sure2);
            swap(ss2, ss), swap(sure2, chroma_k);
            if (verb >= 2 and ss.size())
              pr("Found k-vcs of size = {} ({}chroma_k)!\n", ss.size(),
                 chroma_k ? "" : "not ");
          }
        }
      } else
        tie(ss, chroma_k) = find_k_vcs(k, t);
    }
    if (ss.empty()) {
      if (verb >= 1)
        pr(">> Global iteration #{}: gen.: empty (b: {}), fin.: empty (b: {}), "
           "not sure\n",
           global_iter, best_gen.size(), best_fin.size());
      continue;
    }
    iter_cache.erase((int)ss.size());
    ++stats::num_gen_subsets;
    stats::tot_size_gen += ss.size();
    sort(begin(ss), end(ss));
    cur_ss_gen = ss;
    bool crit = false;
    if (best_gen.empty() or (int) best_gen.size() == n or
        pair(not chroma_k, ss.size()) < pair(not best_gen_sure, best_gen.size()))
      best_gen = ss, best_gen_sure = chroma_k;
    update_global_best(ss, chroma_k, false);
    if (use_gen_cache) {
      auto it = gen_cache.find(ss);
      if (it != gen_cache.end()) {
        if (verb >= 2) pr("Gen.~duplicate, skipping...\n");
        stats::tot_size_fin += it->second.first.size();
        tie(ss, chroma_k) = it->second;
        goto print_iter;
      }
    }
    if (not ss.empty()) {
      did_postproc = true;
      auto [ss_n, sure_n, crit_n] =
          pproc(ss, chroma_k and ((int)ss.size() != n), k, t);
      if (use_gen_cache) gen_cache[ss] = mp(ss_n, sure_n);
      if (mt(not sure_n, not crit_n, ss_n.size()) <
          mt(not chroma_k, not crit, ss.size())) {
        ss = move(ss_n), chroma_k = sure_n, crit = crit_n;
      }
      ++stats::num_pp;
    }
    stats::tot_size_fin += ss.size();
    update_global_best(ss, chroma_k, crit);
  print_iter:
    if (verb >= 1)
      pr(">> Global iteration #{}: gen.: {} (b: {}), fin.: {} (b: {}), "
         "chroma_k: {}, crit: {}\n",
         global_iter, cur_ss_gen.size(), best_gen.size(), ss.size(), best_fin.size(),
         chroma_k, crit);
    if ((int)best_fin.size() == k + 2) break;
  }
  return mp(best_fin, best_fin_chroma_k);
}
void check_kvcs_init() {
  auto is = is_k_colorable_heuristic(k - 1, ind_n, timer(15.0, global_timer));
  if (is or global_timer.timed_out()) {
    stats::infeas = true;
    if (global_timer.timed_out()) {
      best_fin_chroma_k = false;
      if (verb >= 1)
        pr("\nTimed out before testing whether input graph is k-VCS. Stop.\n", k);
    } else {
      best_fin_chroma_k = true;
      if (verb >= 1) pr("\nInput graph does not contain a k-VCS. Stop.\n", k);
    }
    exit(EXIT_SUCCESS);
  }
}
void check_criticality() {
  if ((int)best_fin.size() == n) best_fin_chroma_k = true;
  if ((int)best_fin.size() == k) best_fin_chroma_k = best_fin_crit = true;
  if (best_fin_chroma_k) stats::confirmed_chroma_early = stats::confirmed_chroma = true;
  stats::confirmed_crit_early = best_fin_crit;
  if (stats::confirmed_chroma and best_fin_crit) stats::confirmed_crit = true;
  if (do_force_confirm || (do_confirm_criticality and (int) best_fin.size() >= k)) {
    if (verb >= 1) {
      pr("\n>> Best solution has size {}, will now use {}s to confirm if "
         "it's "
         "really critical...\n",
         best_fin.size(), confirm_crit_timelimit);
    }
    timer confirm_chroma_timer(confirm_crit_timelimit);
    if (do_force_confirm || not best_fin_chroma_k) {
      stats::confirmed_chroma =
          not is_k_colorable_exact(k - 1, best_fin, confirm_chroma_timer);
    }
    stats::confirm_chroma_time = confirm_chroma_timer.elapsed_secs();
    pr("(confirmed chroma: {} in {}s)\n", stats::confirmed_chroma,
       stats::confirm_chroma_time);
    stats::confirmed_crit = best_fin_crit;
    timer confirm_crit_timer(confirm_crit_timelimit);
    if (do_force_confirm || (stats::confirmed_chroma and not best_fin_crit)) {
      stats::confirmed_crit = true;
      for (int i : best_fin) {
        if (confirm_crit_timer.timed_out()) {
          stats::confirmed_crit = false;
          break;
        }
        vi v = best_fin;
        remove(begin(v), end(v), i);
        v.pop_back();
        const double tl_each = max(0.5, double(best_fin.size()) / confirm_crit_timelimit);
        timer tmr(tl_each);
        bool is_kc = is_k_colorable_exact(k - 1, v, tmr);
        stats::confirmed_crit = stats::confirmed_crit and is_kc and not tmr.timed_out();
        if (not stats::confirmed_crit) break;
      }
    }
    stats::confirm_crit_time = confirm_crit_timer.elapsed_secs();
    if (verb >= 1) {
      pr(">> Solution {}critical, and \n", stats::confirmed_crit ? "is " : "may not be ");
      pr("{} chromatic number >= k.\n", stats::confirmed_chroma ? "has" : "may not have");
      pr(">> Time used: {} + {} seconds.\n", stats::confirm_chroma_time,
         stats::confirm_crit_time);
    }
  }
}
void output_to_file() {
  if (not output_filename.empty() and not best_fin.empty()) {
    string s;
    for (uint i = 0; i < best_fin.size(); ++i) {
      if (i != 0) s += " ";
      if (inrange(best_fin[i], 0, (int)vmap.size() - 1)) {
        s += to_string(vmap[best_fin[i]]);
      }
    }
    ofstream of(output_filename);
    if (of.good()) {
      of << s;
      of.close();
    }
  }
}
void do_irace_its() {
  exited = true;
  int res = 0;
  const double targets[] = {0.1, 0.25, 0.5};
  for (double target : targets) {
    k = ceil(target * n);
    int best_m = 0;
    timer t(time_limit_secs);
    for (int iter = 1; not t.timed_out(); ++iter) {
      subgraph h;
      int r = rand_int(0, 2 * m - 1);
      auto it = lower_bound(begin(ind_deg_cum), end(ind_deg_cum), r);
      int i = ind_deg[it - begin(ind_deg_cum)];
      if (t.timed_out()) break;
      int lo = n * 0.25, hi = n * 0.75;
      if (k < lo)
        h = cons_add(k, i, cons_alpha);
      else if (k > hi)
        h = cons_drop(k, cons_alpha);
      else if (rand_int(lo, hi) < k)
        h = cons_drop(k, cons_alpha);
      else
        h = cons_add(k, i, cons_alpha);
      ls(h);
      best_m = max(best_m, h.m);
    }
    if (verb >= 1) print("Target: {}, value: {}\n", k, -best_m);
    res -= best_m;
  }
  print("{}", res);
}
void do_just_exact_coloring() {
  int lb = max(k, (int)size(hao_mnts_max_clique(ind_n, k, timer(clique_alg_time_1st))));
  int colors = color_exactly(ind_n, global_timer, lb);
  pr("{} {} {}\n", instance_name, global_timer.elapsed_secs(), colors);
}
void do_just_heuristic_coloring() {
  int lb = max(k, (int)size(hao_mnts_max_clique(ind_n, k, timer(clique_alg_time_1st))));
  int colors = color_heuristically(ind_n, global_timer, lb);
  pr("{} {} {}\n", instance_name, global_timer.elapsed_secs(), colors);
}
void exit_fun(int sig) {
  if (not exited) {
    exited = true;
    if (last_exit_code == EXIT_SUCCESS) {
      stats::time = global_timer.elapsed_secs();
      if (sig != SIGINT)
        check_criticality();
      output_to_file();
      stats::print_stats();
    }
  }
  exit(last_exit_code);
}
void exit_fun_2() {
  if (not exited) exit_fun(0);
}
vi mcqd_max_clique(const vi& ind, timer t) {
  int n = ind.size();
  bool** conn = new bool*[n];
  for (int i = 0; i < n; ++i)
    conn[i] = new bool[n];
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      conn[i][j] = AM[ind[i]][ind[j]];
  Maxclique m(conn, n, t);
  int* qmax;
  int qsize;
  m.mcq(qmax, qsize);
  vi clique(qsize);
  for (int i = 0; i < qsize; ++i)
    clique[i] = qmax[i];
  delete[] qmax;
  for (int i = 0; i < n; ++i)
    delete[] conn[i];
  delete[] conn;
  return clique;
}
int main(int argc, char** argv) {
  cmd_line(argc, argv);
  if (normal_run) {
    signal(SIGINT, exit_fun);
    atexit(exit_fun_2);
    for (int i = 0; i < argc; ++i)
      pr("{} ", argv[i]);
    pr("\n");
  }
  bool do_preprocess = normal_run;
  read_dimacs(input_filename, do_preprocess);
  if (irace_its) {
    do_irace_its();
    return EXIT_SUCCESS;
  }
  if (just_exact_coloring) {
    do_just_exact_coloring();
    return EXIT_SUCCESS;
  }
  if (just_heuristic_coloring) {
    do_just_heuristic_coloring();
    return EXIT_SUCCESS;
  }
  best_fin = best_gen = ind_n;
  best_fin_chroma_k = true;
  global_timer.reset(time_limit_secs);
  if (verb >= 1) pr("\n");
  vi clique;
  try {
    clique = hao_mnts_max_clique(ind_n, k, timer(clique_alg_time_1st));
  } catch (std::exception& e) {
    if (verb >= 1) pr("Exception on max-clique: {}\n", e.what());
    clique.assign(1, 0);
  };
  if (verb >= 1) pr("Maximum clique of size {} found.\n", size(clique));
  stats::cliq_1st_size = clique.size();
  if (stats::cliq_1st_size >= k) {
    while ((int)clique.size() > k)
      clique.pop_back();
    tie(best_fin, best_fin_chroma_k) = mp(clique, true);
  } else {
    if (just_pproc) {
      shuffle(begin(ind_n), end(ind_n), rng);
      pproc(ind_n, true, k, global_timer);
    } else {
      our_algorithm(k, global_timer);
    }
  }
}
