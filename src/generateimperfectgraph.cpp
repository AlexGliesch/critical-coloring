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
#include "btdsatur/bkdmain.h"
#include "cliques/mntshao.h"
#include "hybridea/main.h"
#include "util.h"
int n, k, seed, cl;
double d;
vvi AM;
double clique_tl = 5;
double exa_tl = 5;
double heu_tl = 5;
int cols_heu() {
  static hybridea::Graph g;
  if (g.n == 0) {
    g.resize(n);
  }
  g.n = n;
  g.nbEdges = 0;
  fill(g.matrix, g.matrix + g.n * g.n, 0);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j) {
      g[i][j] = g[j][i] = AM[i][j];
      ++g.nbEdges;
    }
  int rnd_seed = rand_int(nli::min(), nli::max());
  int colors = hybridea::hea(g, timer(heu_tl), 100000000000000LL, cl, rnd_seed,
                             10, 16, 0 , 1, 1, false, nullptr);
  return colors;
}
pair<bool, int> cols_exa() {
  srand(rand_int(nli::min(), nli::max()));
  btdsatur::order = n;
  btdsatur::maxChecks = 100000000000000LL;
  btdsatur::verbose = 0;
  memset(btdsatur::graph, 0, BTDSATUR_GRAPHSIZE);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      if (AM[i][j]) {
        btdsatur::setedge(i, j);
        btdsatur::setedge(j, i);
      }
  int colors = 0;
  bool timeout = false;
  try {
    colors = btdsatur::colorsearch(1, timer(exa_tl), nullptr);
  } catch (timeout_exception& e) {
    timeout = true;
  }
  return make_pair(not timeout, colors);
}
int clique() {
  static vi indn;
  if (indn.empty()) {
    indn.resize(n);
    iota(begin(indn), end(indn), 0);
  }
  return hao_mnts_max_clique(indn, n, timer(clique_tl)).size();
}
int main(int argc, char** argv) {
  if (argc != 5) {
    pr("Usage: ./generateimperfectgraph num_vertices density seed out_file. "
       "Incorrect "
       "parameter types will not be checked against.\n");
    exit(EXIT_SUCCESS);
  }
  n = atoi(argv[1]);
  d = atof(argv[2]);
  seed = atoi(argv[3]);
  string out_file = argv[4];
  rng.seed(seed);
  AM.assign(n, vi(n, 0));
  timer t;
  int iter = 0;
  int num_edges = int((n * (n - 1) * d) / 2.0);
  while (t.elapsed_secs() < 3600) {
    ++iter;
    vii edges;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j)
        edges.emplace_back(i, j);
    set<int> chosen;
    while ((int)chosen.size() != num_edges)
      chosen.insert(rand_int(0, (int)edges.size() - 1));
    for (int i : chosen) {
      AM[edges[i].first][edges[i].second] =
          AM[edges[i].second][edges[i].first] = 1;
    }
    cl = clique();
    auto [exa_ok, exa_cols] = cols_exa();
    if (exa_ok) {
      k = exa_cols;
    } else {
      k = cols_heu();
    }
    if (k != cl) {
      assert(k > cl);
      ofstream of(out_file);
      of << format("c Random graph, n {}, d {}, s {}, clique {}, k {}\n", n, d,
                   seed, cl, k);
      of << format("p edge {} {}\n", n, num_edges);
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
          if (AM[i][j]) of << format("e {} {}\n", i + 1, j + 1);
      pr("{} {}\n", out_file, k);
      of.close();
      return 0;
    }
  }
  pr("Timed out, could not generate a graph with parameters n {}, d {}, s "
     "{}.\n",
     n, d, seed);
}
