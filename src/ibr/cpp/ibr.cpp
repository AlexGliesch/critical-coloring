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
#include "ibr.h"
#include "color.h"

enum vt { // Vertex classification of Sun et al.
  A = 1,  // critical
  B = 2,  // unknown
  C = 4   // not critical
};

using vvt = vector<vt>;

// Returns whether vertex sets 'which' are a k-VCS or not
bb is_k_vcs(const vvt& v, int which, const timer& t) {
  // TIME_BLOCK("is_k_vcs");
  vi g;
  for (int i = 0; i < n; ++i)
    if (v[i] & which) g.push_back(i);
  if ((int)g.size() < k) return mp(false, true); // surely not critical
  return is_k_vcs(k, g, t);
}

// Perturbation, Section 3.5
// Using Glover et al.'s scheme, see readme.txt
void perturbation(vvt& v, vi& flip_freq, const vi& elite_freq, int r) {
  // TIME_BLOCK("perturbation");
  int max_freq = *max_element(begin(flip_freq), end(flip_freq));
  vector<pair<int, double>> score_index;
  for (int i = 0; i < n; ++i)
    if (v[i] == A)
      score_index.push_back(
          mp((elite_freq[i] * (r - elite_freq[i])) / double(r * r) +
                 ((1 - flip_freq[i]) / double(max_freq)),
             i));
  int Asize = score_index.size();
  sort(begin(score_index), end(score_index), greater<pair<int, double>>());

  int gam = max(1, int(Asize * gama));
  double den = 0.0;
  for (int j = 0; j < gam; ++j)
    den += pow(1 + j, -lambda);
  [[maybe_unused]] int num_pert = 0;

  for (int j = 0; j < gam; ++j) {
    double Pj = pow(1 + j, -lambda) / den;
    if (rand_double(0.0, 1.0) <= Pj) {
      ++num_pert;
      int i = score_index[j].second;
      assert(inrange(i, 0, n - 1) and v[i] == A);
      v[i] = B;
      ++flip_freq[i];
    }
  }
  ++stats::num_pert;
  stats::tot_pert += num_pert;
  if (verb >= 1) pr("Perturbed {} vertices\n", num_pert);
}

// Update sets v, see Section 3.6
void update_set(vvt& v, vi& flip_freq, timer t) {
  // TIME_BLOCK("update_set");
  if (not is_k_vcs(v, A, t).first) {
    if (verb >= 1) pr("update_set: move all C to B\n");
    for (int i = 0; i < n; ++i)
      if (v[i] == C) { // move all from C to B
        v[i] = B;
        ++flip_freq[i];
      }
  } else {
    // if (verb >= 1) pr("update_set: move all B to C; A to B\n");
    for (int i = 0; i < n; ++i)
      if (v[i] == B) { // move all from B to C
        v[i] = C;
        ++flip_freq[i];
      } else if (v[i] == A) { // move all from A to B
        v[i] = B;
        ++flip_freq[i];
      }
  }
}

// Choose vertex with smallest neighborhood weight W; see Section 3.4
int choose_remove(const vvt& v) {
  // TIME_BLOCK("choose_remove");
  int r = -1, sr = -1;
  for (int i = 0; i < n; ++i)
    if (v[i] == B) {
      int si = 0;
      for (int j : AL[i])
        si += (v[j] == A ? m : (v[j] == B ? 1 : 0));
      if (r == -1 or sr > si) r = i, sr = si;
    }
  return r;
}

// Update best known solution so far
void update_global_best(const vvt& v, int which, int r, bool sure) {
  // TIME_BLOCK("update_global_best");
  int S = 0;
  for (int i = 0; i < n; ++i)
    S += bool(v[i] & which);
  if (global_best.empty() or
      mp((int)global_best.size(), not global_best_sure) > mp(S, not sure)) {
    if (verb >= 1)
      pr("New global best: {}->{}, sure: {}\n", global_best.size(), S, sure);
    stats::ttb = global_timer.elapsed_secs();
    stats::itb = r;
    global_best.clear();
    for (int i = 0; i < n; ++i)
      if (v[i] & which) global_best.push_back(i);
    global_best_m = 0;
    for (uint i = 0; i < global_best.size(); ++i)
      for (uint j = i + 1; j < global_best.size(); ++j)
        global_best_m += AM[global_best[i]][global_best[j]];
    global_best_sure = sure;
  }
}

// Main iterated backtracking algorithm
void ibr(timer t) {
  vvt v(n, B);         // v: which group each vertex belongs to
  vi flip_freq(n, 0);  // flip frequency, see Section 3.5
  vi elite_freq(n, 0); // elite frequency, see Section 3.5
  int r = 0;           // r is the number of performed perturbations
  while (r < R and not t.timed_out()) {
    bool isvcs, sure;
    ++stats::num_r;
    stats::num_flips = accumulate(begin(flip_freq), end(flip_freq), 0);
    if (verb >= 1) pr("Iteration #{}...\n", r + 1);
    stack<int> moved_to_C; // Used for backtracking
    for (int i = 0; i < n; ++i)
      if (v[i] == C) moved_to_C.push(i);
    while (true) {
      tie(isvcs, sure) = is_k_vcs(v, A, t);
      if (isvcs or t.timed_out()) {
        break;
      }
      int i = choose_remove(v); // Vertex removal, Section 3.2
      if (i == -1) break;
      assert(v[i] == B);
      v[i] = C;
      ++flip_freq[i];
      moved_to_C.push(i);
      // if (verb >= 1) pr("Moving i={} from B to C\n", i);
      tie(isvcs, sure) = is_k_vcs(v, A | B, t);
      if (t.timed_out()) {
        break;
      }
      if (not isvcs) {
        // if (verb >= 1) pr("AuB not VCS anymore! Moving i={} to A\n", i);
        v[i] = A;
        ++flip_freq[i];
        moved_to_C.pop();
        // Backtracking, Section 3.4: move back in the order they were put in
        while (true) {
          tie(isvcs, sure) = is_k_vcs(v, A | B, t);
          if (isvcs or t.timed_out() or moved_to_C.empty()) {
            break;
          }
          assert(moved_to_C.size() > 0);
          int l = moved_to_C.top();
          moved_to_C.pop();
          assert(v[l] == C);
          v[l] = B;
          ++flip_freq[l];
          ++stats::num_bt;
          // if (verb >= 1)
          // pr("AuB still not VCS! Backtracking l={} from C to B\n", l);
        }
      }
      update_global_best(v, A | B, r, sure);
      [[maybe_unused]] int A_size = count(begin(v), end(v), A),
                           C_size = moved_to_C.size(),
                           B_size = n - A_size - C_size;
      if (verb >= 1) pr("A: {}, B: {}, C: {}\n", A_size, B_size, C_size);
    }
    ++r;
    if (not t.timed_out()) {
      // if (verb >= 1) pr("Found k-VCS of size {}!\n", count(begin(v), end(v),
      // A));
      update_global_best(v, A | B, r, sure);
      if (r < R) {
        // assert(is_k_vcs(v, A, timer()));
        for (int i = 0; i < n; ++i)
          elite_freq[i] += (v[i] == A);
        perturbation(v, flip_freq, elite_freq, r);
        update_set(v, flip_freq, t);
      }
    }
  }
}
