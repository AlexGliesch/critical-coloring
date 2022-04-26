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
#include "../../util.h"

inline vvi AM, AL;
inline int n, m, k;

// Parameters
inline bool do_exact_coloring = false; // Whether to attempt an exact method
inline double heu_secs;   // Time limit (seconds) of the heuristic method
inline double exact_secs; // Time limit (seconds) of the exact method
inline int num_heu_reruns =
    0; // Number of times to call the heuristic method each time
inline double confirm_crit_timelimit;
inline bool
    do_confirm_criticality; // If enabled, and if criticality was not proven
                            // during the execution of the algorithm, at the
                            // end an exact coloring algorithm (BTDSatur) is
                            // executed on the best subgraph obtained with a
                            // time limit of "confirmcrittime" seconds, where
                            // confirmcrittime is a parameter.
inline double time_limit_secs;
inline string input_filename;
inline string
    output_filename; // If set, will output the best solution found to given
                     // file as a list of 0-based vertex indices.
inline size_t random_seed;
inline int R;               // Maximum number of perturbations
inline double gama, lambda; // For perturbation, see readme.txt
inline bool multistart;
inline vi global_best; // Best known solution
inline int global_best_m = -1;
inline bool global_best_sure = false;
inline timer global_timer;

// Statistics
namespace stats {
inline int num_r = 0;
inline double color_time = 0.0;
inline int64_t tot_pert = 0;  // Total perturbation strength
inline int num_pert = 0;      // Number of perturbations
inline int64_t num_flips = 0; // Number of flips between subsets
inline int64_t num_bt = 0;    // Number of backtracks
inline double ttb = 0.0;      // Time to best
inline int itb = 0;           // Iterations to best
inline int heu_mistk = 0; // Number of mistakes in the heuristic; only detected
                          // if we are doing multiple reruns
inline int suc_ext_cals = 0;   // Number of calls to the exact coloring method,
                               // where it successfully found a k-coloring
inline int unsuc_ext_cals = 0; // Number of calls to the exact coloring method,
                               // where it failed to find a k-coloring
inline int suc_heu_cals = 0;   // See above, for the heuristic coloring method
inline int unsuc_heu_cals = 0; // See above, for the heuristic coloring method
inline int cals_to_coloring = 0; // Total number of calls to coloring
inline bool confirmed_chroma_early =
    false; // Set to true if we were able to confirm chroma without
           // running an exact algorithm post-hoc, meaning we solved colorings
           // exactly during the execution of the algorithm.
inline bool confirmed_chroma =
    false; // This is set to true if the generated solution's chroma is
           // confirmed to be equal to k, by an exact algorithm.
inline bool confirmed_crit =
    false; // Set to true if the generated solution is confirmedly critical, by
           // attempting to remove each of its vertices and then running an
           // exact algorithm.
inline double confirm_time = 0.0; // Time used to confirm criticality
inline double time = 0.0;   // Running time of the algorithm, excluding the time
                            // needed to confirm criticality (if enabled)
inline bool infeas = false; // #TODO document
inline int repl_count = 0;  //#TODO document
} // namespace stats

inline bool do_print_summary = true;

// Dummies, not used currently
inline int verb = 0;

// Whether we're printing intermediate results at R=20
inline bool print_at_r20 = false;
