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
#include "util.h"
namespace stats {
inline double pp_time;
inline int num_pp = 0;
inline int heu_mistk = 0;
inline int pp_reduced = 0;
inline int suc_ext_cals = 0;
inline int unsuc_ext_cals = 0;
inline int suc_heu_cals = 0;
inline int unsuc_heu_cals = 0;
inline int cals_to_coloring = 0;
inline int n_subsets =
    0;
inline int skip_cons_ls_cache =
    0;
inline int num_ls = 0;
inline int tot_ls_moves = 0;
inline int64_t tot_cons_edges =
    0;
inline int num_size_att = 0;
inline int suc_size_att =
    0;
inline int cliq_1st_size = 0;
inline int gen_rep = 0;
inline double tot_size_gen =
    0;
inline int num_gen_subsets = 0;
inline double tot_size_fin = 0;
inline int trivial_crits_found =
    0;
inline int global_iter = 0;
inline double ttb;
inline int itb;
inline bool confirmed_chroma_early =
    false;
inline bool confirmed_crit_early = false;
inline bool confirmed_chroma =
    false;
inline bool confirmed_crit =
    false;
inline double confirm_crit_time = 0.0;
inline double confirm_chroma_time = 0.0;
inline double time = 0.0;
inline double color_time = 0.0;
inline bool infeas = false;
inline int max_iter_diff = 0;
inline int gen_best_1st_phase = 0;
inline int fin_best_1st_phase = 0;
inline int color_not_ok = 0;
void print_stats();
void print_summary();
inline bool do_print = true;
inline bool do_print_summary = true;
}
