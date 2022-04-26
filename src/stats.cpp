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
#include "stats.h"
#include "cliques/mntshao.h"
#include "main.h"
#include "sm.h"
#include "subgraph.h"
namespace stats {
subgraph best_sg;
void print_summary() {
  best_sg = subgraph(best_fin);
  pr("summary_line ");
  pr("instance={} ", instance_name);
  pr("k={} ", k);
  pr("n_red={} ", n);
  pr("m_red={} ", m);
  pr("n={} ", n_ori);
  pr("m={} ", m_ori);
  pr("heu_timelimit={} ", heu_secs);
  pr("exact_timelimit={} ", exact_secs);
  pr("do_confirm_crit={} ", (int)do_confirm_criticality);
  pr("confirm_timelimit={} ", confirm_crit_timelimit);
  pr("R={} ", R);
  pr("reltg_size={} ", xi);
  pr("do_heu={} ", int(not no_heuristic_coloring));
  pr("ss={} ", best_fin.size());
  pr("ss_gen_1st={} ", gen_best_1st_phase);
  pr("ss_fin_1st={} ", fin_best_1st_phase);
  pr("chroma_k_early={} ", (int)confirmed_chroma_early);
  pr("chroma_k={} ", (int)confirmed_chroma);
  pr("crit_early={} ", (int)confirmed_crit_early);
  pr("crit={} ", (int)confirmed_crit);
  pr("edges={} ", best_sg.m);
  pr("time={} ", time);
  pr("time_pp={} ", pp_time);
  pr("time_confirm_chroma={} ", confirm_chroma_time);
  pr("time_confirm_crit={} ", confirm_crit_time);
  pr("time_color={} ", color_time - confirm_crit_time - confirm_chroma_time);
  pr("iter={} ", global_iter);
  pr("ttb={} ", ttb);
  pr("itb={} ", itb);
  pr("max_iter_diff={} ", max_iter_diff);
  pr("ss_fin_avg={} ", divOrNA(tot_size_fin, num_gen_subsets));
  pr("ss_gen_avg={} ", divOrNA(tot_size_gen, num_gen_subsets));
  pr("ss_gen_best={} ", best_gen.size() ? to_string(best_gen.size()) : "NA");
  pr("num_size_att={} ", num_size_att);
  pr("heu_mistakes={} ", heu_mistk);
  pr("num_pp={} ", num_pp);
  pr("avg_pp_reduced={} ", divOrNA(pp_reduced, num_pp));
  pr("avg_pp_skipped={} ", divOrNA(trivial_crits_found, num_pp));
  pr("calls_to_coloring={} ", cals_to_coloring);
  pr("num_gen_subsets={} ", num_gen_subsets);
  pr("infeas={} ", (int)infeas);
  pr("clique_start={} ", cliq_1st_size);
  pr("color_not_ok={} ", color_not_ok);
  pr("seed={} ", random_seed);
  pr("\n");
}
void print_stats() {
  if (not do_print) return;
  best_sg = subgraph(best_fin);
  if (verb >= 1) {
    pr("\n");
    pr("-- Time: {}\n", global_timer.elapsed_secs());
    pr("-- Time post-processing: {}\n", pp_time);
    pr("-- Iter: {}\n", global_iter);
    pr("-- Generated: {}, avg.: {}\n", best_gen.size(),
       divOrNA(tot_size_gen, num_gen_subsets));
    pr("-- Subset size: {}, avg.: {}\n", best_fin.size(),
       divOrNA(tot_size_fin, num_gen_subsets));
    if (best_sg.ss.size() and best_sg.n == (int)best_sg.ss.size()) {
      pr("-- Edges: {}/{} (density: {})\n", best_sg.m, (best_sg.n * best_sg.n - 1) / 2,
         best_sg.m / double((best_sg.n * best_sg.n - 1) / 2.0));
      pr("-- Chroma k: {}{}\n", confirmed_chroma ? "yes" : "not sure",
         confirmed_chroma_early ? " (early)" : "");
      pr("-- Critical: {}{}\n", confirmed_crit ? "yes" : "not sure",
         confirmed_crit_early ? " (early)" : "");
      pr("-- t.t.b.: {}\n", ttb);
      pr("-- i.t.b.: {}\n", itb);
      pr("-- Calls to p.p. skipped: {}/{}\n", global_iter - num_pp, global_iter);
      pr("-- P.p. reduced (avg.): {}\n", divOrNA(pp_reduced, num_pp));
      pr("-- Total subsets generated: {}\n", n_subsets);
      pr("-- Calls to coloring: {}\n", cals_to_coloring);
      pr("-- Skips cons.+LS cache: {}%\n",
         100.0 * skip_cons_ls_cache / double(n_subsets));
      pr("-- Skips trivially critical (p.p.): {}\n", trivial_crits_found);
      pr("-- Successful size attempts: {}\n", suc_size_att);
    }
    pr("-- Clique (1st): {}\n", cliq_1st_size);
    pr("-- Infeas: {}\n", infeas);
    pr("-- Seed: {}\n", random_seed);
    pr("\n");
  }
#if 0
 if (verb >= 1) {
  print_timed_blocks();
  pr("\n");
 }
#endif
  if (irace) {
    do_print_summary = false;
    pr("{}", (1 - int(confirmed_crit)) * n + best_fin.size());
  }
  if (do_print_summary) print_summary();
}
}
