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
#include "readall.h"
#include "main.h"
#include "stats.h"
#include <boost/any.hpp>
void cmd_line(int argc, char** argv) {
  namespace po = boost::program_options;
  po::variables_map vm;
  po::options_description desc("Finding critical sub-graphs for graph coloring");
  desc.add_options()("help,h", "Show help menu.");
  desc.add_options()("in,i", po::value<string>(&input_filename)->required(),
                     "Input filename.");
  desc.add_options()("out,o", po::value<string>(&output_filename)->default_value(""),
                     "If set, will output the best subgraph found (only if the "
                     "chromatic number is confirmed) to the given filename, as "
                     "a list of 0-based vertex indices.");
  desc.add_options()("k,k", po::value<int>(&k)->required(), "Number of colors.");
  desc.add_options()("time,t", po::value<double>(&time_limit_secs)->default_value(600),
                     "Time limit, in seconds.");
  desc.add_options()("heutime", po::value<double>(&heu_secs)->default_value(0.5),
                     "Time limit of each HEA call, in seconds.");
  desc.add_options()("exacttime", po::value<double>(&exact_secs)->default_value(1.0),
                     "Time limit of each BTDSatur call, in seconds.");
  desc.add_options()("cliquetime",
                     po::value<double>(&clique_alg_time_1st)->default_value(2.0),
                     "Time limit of the initial MN/TS call, in seconds.");
  desc.add_options()("iter", po::value<int>(&max_iter)->default_value(nli::max()),
                     "Maximum number of iterations (=1: only do first phase).");
  desc.add_options()("mu", po::value<double>(&mu)->default_value(1.5),
                     "Subset size increase multiplicative.");
  desc.add_options()("R", po::value<int>(&R)->default_value(200),
                     "Number of subgraphs to generate at each size attempt.");
  desc.add_options()(
      "xi", po::value<double>(&xi)->default_value(1.08),
      "Multistart iterations other than the first one will use "
      "xi*sb as a target subset size, where \'sb\' is the size of "
      "the best generated subset. Use a value < 1.0 to disable. If disabled, "
      "the algorithm looks for a new best size at each iteration.");
  desc.add_options()("imax", po::value<int>(&imax)->default_value(nli::max()),
                     "Maximum non-improving iterations of the second phase.");
  desc.add_options()("consalg", po::value<string>(&cons_alg)->default_value("adddrop"),
                     "Constructive algorithm, in [add,drop,adddrop].");
  desc.add_options()("alpha", po::value<double>(&cons_alpha)->default_value(0.1),
                     "Alpha parameter of the constructive algorithm.");
  desc.add_options()("tenure", po::value<double>(&tenure_mult)->default_value(0.1),
                     "Tabu tenure.");
  desc.add_options()("imaxits", po::value<int>(&max_nonimpr)->default_value(10000),
                     "Maximum non-improving iterations of the iterated tabu "
                     "search heuristic.");
  desc.add_options()("pmin", po::value<int>(&pmin)->default_value(1),
                     "p-min (VNS part).");
  desc.add_options()("pmax", po::value<double>(&pmax_mult)->default_value(1.0),
                     "p-max (VNS part), relative to n.");
  desc.add_options()("pstep", po::value<int>(&pstep)->default_value(20),
                     "p-step (VNS part).");
  desc.add_options()("noheu",
                     po::bool_switch(&no_heuristic_coloring)->default_value(false),
                     "Do not run a heuristic algorithm for coloring.");
  desc.add_options()("noexact",
                     po::bool_switch(&no_exact_coloring)->default_value(false),
                     "Do not run an exact algorithm for coloring.");
  desc.add_options()("nopproc", po::bool_switch(&no_postproc)->default_value(false),
                     "Do not run post-processing on generated subgraphs.");
  desc.add_options()("nofstphase", po::bool_switch(&no_fst_phase)->default_value(false),
                     "Do not run the first phase; start second phase from full graph");
  desc.add_options()(
      "nodensesearch", po::bool_switch(&no_dense_search)->default_value(false),
      "Generate subgraphs using random samples instead of dense search.");
  desc.add_options()(
    "justpproc", po::bool_switch(&just_pproc)->default_value(false),
    "Only do one round of post-processing, and not the full algorithm.");
  desc.add_options()(
      "confirmcrit", po::bool_switch(&do_confirm_criticality)->default_value(false),
      "If enabled, and if criticality was not proven during the execution of "
      "the algorithm, at the end an exact coloring algorithm (BTDSatur) is "
      "executed on the best subgraph obtained with a time limit of "
      "\"confirmcrittime\" seconds, where confirmcrittime is a parameter.");
  desc.add_options()("confirmtime",
                     po::value<double>(&confirm_crit_timelimit)->default_value(60),
                     "Time limit on BTDSatur to confirm criticality, at the end.");
  desc.add_options()("forceconfirm",
                     po::bool_switch(&do_force_confirm)->default_value(false),
                     "Run the post-hoc check even if chromaticity/criticality were "
                     "confirmed by the base heuristic.");
  options_counter verbosec;
  desc.add_options()("verbose,v", po::value(&verbosec)->zero_tokens(),
                     "Verbosity. If present, output is sent to screen. If -v "
                     "is repeated, more output is given.");
  desc.add_options()("summaryeachiter",
                     po::bool_switch(&print_summary_each_iter)->default_value(false),
                     "Print a summary line after each iteration.");
  desc.add_options()("seed,s", po::value<size_t>(&random_seed)->default_value(0),
                     "Random seed. If 0, a random value will be used.");
  desc.add_options()(
      "iraceits", po::bool_switch(&irace_its)->default_value(false),
      "If set, will run the iterated tabu search algorithm for "
      "finding the subgraph of size U of maximum density. Assumes input graph "
      "is in a special format which specifies target density; see the "
      "read_dimacs function. Outputs a single number (the number of edges).");
  desc.add_options()(
      "irace", po::bool_switch(&irace)->default_value(false),
      "If set, will output a single number X=(1-c)n+s, where c\\in {0,1} is "
      "whether the result was confirmedly critical, s is the resulting "
      "subgraph size, and n is the number of vertices in the input graph."
      " Requires the input graph to be in a special format which specifies "
      "target k, see the read_dimacs function.");
  desc.add_options()(
      "justexactcolor", po::bool_switch(&just_exact_coloring)->default_value(false),
      "If set, will run an exact coloring algorithm with the given time limit, "
      "and stop. Outputs the chromatic number found, or -1, if the exact "
      "method timed out. Uses the given k as a starting lower bound.");
  desc.add_options()(
      "justheuristiccolor",
      po::bool_switch(&just_heuristic_coloring)->default_value(false),
      "If set, will run a heuristic coloring algorithm with the given time "
      "limit, and stop. Outputs the best upper bound found. Uses the given k "
      "as a lower bound.");
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      stats::do_print_summary = false;
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }
    if (cons_alg != "add" and cons_alg != "drop" and cons_alg != "adddrop") {
      throw po::validation_error(po::validation_error::invalid_option_value, "consalg",
                                 cons_alg);
    }
    if (random_seed == 0) random_seed = unique_random_seed();
    rng.seed(random_seed);
    normal_run =
        not(irace_its or irace or just_exact_coloring or just_heuristic_coloring);
    instance_name = boost::filesystem::path(input_filename).stem().filename().string();
    global_timer = timer(time_limit_secs);
    verb = max(0, min(3, verbosec.count));
    if (verb >= 1) {
      pr("--instance: {}\n", instance_name);
      pr("--k: {}\n", k);
      pr("--time limit: {}\n", time_limit_secs);
      pr("--random seed: {}\n", random_seed);
    }
    if (no_exact_coloring and no_heuristic_coloring) {
      pr("Error: options --noheu and --noexact cannot be enabled at the same time.");
      exit(EXIT_SUCCESS);
    }
  } catch (po::error& e) {
    print("Error parsing options: {}.\n", e.what());
    stats::do_print = false;
    exit(EXIT_SUCCESS);
  }
}
void read_dimacs(const string& filename, bool do_preprocess) {
  ifstream f(filename);
  if (not f) {
    print("Could not open input file.\n");
    exit(EXIT_FAILURE);
  }
  string buf, tmp;
  int num_edges = 0;
  bool invalid = false;
  while (not invalid and getline(f, buf)) {
    boost::trim_left(buf);
    if (buf.empty()) continue;
    if (buf[0] == 'c') {
      if ((irace_its and boost::starts_with(buf, "c k ")) or
          (irace and boost::starts_with(buf, "c Random graph,"))) {
        vs s;
        boost::split(s, buf, boost::is_any_of(" "));
        k = stoi(s.back());
      }
      continue;
    } else if (buf[0] == 'p') {
      stringstream ss(buf);
      ss >> tmp >> tmp >> n >> m;
      if (boost::to_lower_copy(tmp) != "edge" and boost::to_lower_copy(tmp) != "col") {
        invalid = true;
        break;
      }
      AL.resize(n);
      AM.assign(n, vi(n, 0));
    } else if (buf[0] == 'e') {
      stringstream ss(buf);
      int v1, v2;
      ss >> tmp >> v1 >> v2;
      assert(inrange(v1, 1, n) and inrange(v2, 1, n));
      --v1, --v2;
      if (AM[v1][v2] == 0) {
        assert(AM[v2][v1] == 0);
        AL[v1].push_back(v2);
        AL[v2].push_back(v1);
        AM[v1][v2] = AM[v2][v1] = 1;
        ++num_edges;
      }
    } else
      invalid = true;
  }
  if (invalid) {
    print("Error: instance {} does not seem to be in DIMACS col format.\n", filename);
    print("Line: {}\n", buf);
    exit(EXIT_FAILURE);
  }
  m = num_edges;
  if (verb >= 1) pr("--n: {}, m: {}\n", n, m);
  vmap.resize(n);
  iota(begin(vmap), end(vmap), 0);
  n_ori = n, m_ori = m, AM_ori = AM, AL_ori = AL;
  if (do_preprocess) preprocess();
  ind_n.resize(n);
  iota(begin(ind_n), end(ind_n), 0);
  ind_deg = ind_n;
  sort_by_degree(ind_deg);
  ind_deg_cum = vi(n, 0);
  ind_deg_cum[0] = AL[ind_deg[0]].size();
  for (int i = 1; i < n; ++i)
    ind_deg_cum[i] = ind_deg_cum[i - 1] + AL[ind_deg[i]].size();
}
void preprocess() {
  while (true) {
    auto n_bef = n, m_bef = m;
    auto AM_bef = move(AM);
    auto AL_bef = move(AL);
    auto vmap_bef = move(vmap);
    assert(vmap.empty());
    vi vs;
    for (int i = 0; i < n_bef; ++i)
      if ((int)AL_bef[i].size() >= k - 1) {
        vs.push_back(i);
        vmap.push_back(vmap_bef[i]);
      }
    n = vs.size();
    m = 0;
    AL.resize(n);
    AM.assign(n, vi(n, 0));
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) {
        if (AM_bef[vs[i]][vs[j]]) {
          AM[i][j] = AM[j][i] = 1;
          AL[i].push_back(j);
          AL[j].push_back(i);
          ++m;
        }
      }
    if (verb >= 1)
      pr("Preprocess, n,m before: {},{}; n,m after: {},{}\n", n_bef, m_bef, n, m);
    if (n_bef == n and m_bef == m) break;
  }
  if (n == 0) {
    stats::infeas = true;
    if (verb >= 1) {
      pr("Graph is not a k-VCS, by pre-processing.\n");
    }
    exit(EXIT_SUCCESS);
  }
}
