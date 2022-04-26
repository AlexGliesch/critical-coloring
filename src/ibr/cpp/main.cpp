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
#include "color.h"
#include "ibr.h"
#include "../../fmt/ranges.h"

// Local
namespace {
bool exited = false;
bool started = false;
} // namespace

void read_cmd_line(int argc, char** argv) {
  namespace po = boost::program_options;
  po::variables_map vm;
  po::options_description desc(
      "Reimplementation of Sun et. al's algorithm \"Iterated backtrack removal "
      "search for finding k-vertex-critical subgraphs\" (doi: "
      "10.1007/s10732-017-9358-5)");
  desc.add_options()("help,h", "Show help menu.");
  // Verbosity
  desc.add_options()("verbose,v", "Verbose.");
  // Input filename
  desc.add_options()("in,i", po::value<string>(&input_filename)->required(),
                     "Input filename.");
  // Output filename
  desc.add_options()("out,o",
                     po::value<string>(&output_filename)->default_value(""),
                     "If set, will output the best solution found to given "
                     "file as a list of 0-based vertex indices.");
  // k
  desc.add_options()("k,k", po::value<int>(&k)->required(),
                     "Number of colors.");
  // Time
  desc.add_options()("time,t",
                     po::value<double>(&time_limit_secs)->default_value(600),
                     "Time limit, in seconds.");
  // Heuristic time
  desc.add_options()("heutime",
                     po::value<double>(&heu_secs)->default_value(0.5),
                     "Time limit to each HEA call, in seconds.");
  // Exact time
  desc.add_options()("exacttime",
                     po::value<double>(&exact_secs)->default_value(1.0),
                     "Time limit to each BTDSatur call, in seconds.");
  // Don't run exact coloring
  desc.add_options()("exact",
                     po::bool_switch(&do_exact_coloring)->default_value(false),
                     "Do not run exact algorithm for coloring.");
  // Lambda
  desc.add_options()(
      "lambda", po::value<double>(&lambda)->default_value(1.2),
      "Lambda parameter; see Glover et al. (doi: 10.1007/s10288-009-0115-y).");
  // Lambda
  desc.add_options()(
      "gamma", po::value<double>(&gama)->default_value(0.25),
      "Gamma parameter; see Glover et al. (doi: 10.1007/s10288-009-0115-y).");
  // Number of iterations
  desc.add_options()("maxpert,R", po::value<int>(&R)->default_value(nli::max()),
                     "Maximum number of perturbations.");
  // Use a multistart strategy
  desc.add_options()(
      "multistart", po::bool_switch(&multistart)->default_value(false),
      "Restart algorithm after R perturbations, instead of terminating.");

  // Do confirm criticality?
  desc.add_options()(
      "confirmcrit",
      po::bool_switch(&do_confirm_criticality)->default_value(false),
      "If enabled, and if criticality was not proven during the execution of "
      "the algorithm, at the end an exact coloring algorithm (BTDSatur) is "
      "executed on the best subgraph obtained with a time limit of "
      "\"confirmcrittime\" seconds, where confirmcrittime is a parameter.");
  // Confirm time limit
  desc.add_options()(
      "confirmtime",
      po::value<double>(&confirm_crit_timelimit)->default_value(600),
      "Time limit on BTDSatur to confirm criticality, at the end.");
  // Seed
  desc.add_options()("seed,s",
                     po::value<size_t>(&random_seed)->default_value(0),
                     "Random seed. If 0, a random value will be used.");

  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
      cout << desc << endl;
      exit(EXIT_SUCCESS);
    } else {
      po::notify(vm);
    }
  } catch (po::error& e) {
    print("Error parsing command line: {}.\n{}\n", e.what(), desc);
    exit(EXIT_SUCCESS);
  }
  global_timer = timer(time_limit_secs);
  if (random_seed == 0) random_seed = unique_random_seed();
  rng.seed(random_seed);
  verb = vm.count("verbose");
}

void read_dimacs(const string& filename) {
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
    if (buf[0] == 'c')
      continue; // comment
    else if (buf[0] == 'p') {
      stringstream ss(buf);
      ss >> tmp >> tmp >> n >> m;
      if (boost::to_lower_copy(tmp) != "edge" and
          boost::to_lower_copy(tmp) != "col") {
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
    print("Error: instance does not seem to be in DIMACS col format.\n");
    print("Line: {}\n", buf);
    exit(EXIT_FAILURE);
  }
  m = num_edges;
}

void print_stats() {
  if (not started) return;
  // print_timed_blocks();

  string instance_name =
      boost::filesystem::path(input_filename).stem().filename().string();
  if (verb >= 1) {
    pr("\n");
    pr("-- Time: {}\n", global_timer.elapsed_secs());
    pr("-- Vertices: {}\n", global_best.size());
    pr("-- Edges: {}\n", global_best_m);
    pr("-- Witness: {}\n", global_best);
    pr("-- Coloring time: {}\n", stats::color_time);
    pr("-- Calls to color: {}\n", stats::cals_to_coloring);
    pr("-- Iter.: {}\n", stats::num_r);
    pr("-- t.t.b.: {}\n", stats::ttb);
    pr("-- i.t.b.: {}\n", stats::itb + 1);
    pr("-- Flips: {}\n", stats::num_flips);
    pr("-- Backtracks: {}\n", stats::num_bt);
    pr("-- Avg. pert.: {}\n", divOrNA(stats::tot_pert, stats::num_pert));
    pr("-- Infeas: {}\n", stats::infeas);
    pr("-- Seed: {}\n", random_seed);
    pr("\n");
  }

  // Summary
  if (do_print_summary) {
    pr("summary_line ");
    pr("instance={} ", instance_name);
    pr("k={} ", k); // colors?
    pr("n={} ", n);
    pr("m={} ", m);
    pr("heu_timelimit={} ", heu_secs);
    pr("exact_timelimit={} ", exact_secs);
    pr("do_confirm_crit={} ", (int)do_confirm_criticality);
    pr("confirm_timelimit={} ", confirm_crit_timelimit);
    pr("do_exact={} ", int(do_exact_coloring));
    pr("ss={} ", global_best.size());
    pr("chroma_k_early={} ", (int)stats::confirmed_chroma_early);
    pr("chroma_k={} ", (int)stats::confirmed_chroma);
    pr("crit={} ", (int)stats::confirmed_crit);
    pr("edges={} ", global_best_m);
    pr("num_perturb={} ", stats::num_r);
    pr("repl={} ", stats::repl_count);
    pr("time={} ", stats::time);
    pr("time_confirm={} ", stats::confirm_time);
    pr("time_color={} ", stats::color_time - stats::confirm_time);
    pr("calls_to_coloring={} ", stats::cals_to_coloring);
    // pr("exact_calls={} ", stats::suc_ext_cals + stats::unsuc_ext_cals);
    // pr("heu_calls={} ", stats::suc_heu_cals + stats::unsuc_heu_cals);
    pr("heu_mistakes={} ", stats::heu_mistk);
    pr("avg_perturb_str={} ", divOrNA(stats::tot_pert, stats::num_pert));
    pr("ttb={} ", stats::ttb);
    pr("itb={} ", stats::itb + 1);
    pr("num_flips={} ", stats::num_flips);
    pr("num_bt={} ", stats::num_bt);
    pr("infeas={} ", (int)stats::infeas);
    pr("seed={} ", random_seed);
    pr("\n");
  }
}

// Check if input graph is a k-VCS
void check_kvcs_init() {
  exact_secs *= 15, heu_secs *= 15; 
  vi ind_n(n);
  iota(begin(ind_n), end(ind_n), 0);

  auto is = is_k_vcs(k, ind_n, global_timer);
  exact_secs /= 15, heu_secs /= 15;
  if (not is.first or global_timer.timed_out()) {
    stats::infeas = true;
    if (global_timer.timed_out()) {
      global_best_sure = false;
      pr("\nTimed out before testing whether input graph is k-VCS. Stop.\n", k);
    } else {
      global_best_sure = is.second;
      pr("\nInput graph does not contain a k-VCS. Stop.\n", k);
    }
    exit(EXIT_SUCCESS);
  }
}

void check_criticality() {
  if (do_confirm_criticality and (int) global_best.size() >= k) {
    if (verb >= 1) {
      pr("\n>> Best solution has size {}, will now use {}s to confirm if "
         "it's "
         "really critical...\n",
         global_best.size(), confirm_crit_timelimit);
    }
    stats::confirmed_chroma =
        ((int)global_best.size() == n) // Assume k is already a lower bound
        or not is_k_colorable_exact(k - 1, global_best,
                                    timer(confirm_crit_timelimit));
    stats::confirmed_crit = false;
    if (stats::confirmed_chroma) {
      stats::confirmed_crit = true;
      for (int i : global_best) {
        vi v = global_best;
        remove(begin(v), end(v), i);
        v.pop_back();
        const double tl_each =
            max(0.5, double(global_best.size()) / confirm_crit_timelimit);
        bool is_kc = is_k_colorable_exact(k - 1, v, timer(tl_each));
        stats::confirmed_crit = stats::confirmed_crit and is_kc;
      }
    }
    stats::confirm_time = global_timer.elapsed_secs() - stats::time;
    if (verb >= 1) {
      pr(">> Solution {}critical, and \n",
         stats::confirmed_crit ? "is " : "may not be ");
      pr("{} chromatic number >= k.\n",
         stats::confirmed_chroma ? "has" : "may not have");
      pr(">> Time used: {}s.\n", stats::confirm_time);
    }
  }
}

void output_to_file() {
  if (not output_filename.empty() and not global_best.empty()) {
    ofstream of(output_filename);
    if (of.good()) of << format("{}", global_best);
  }
}

void exit_fun(int) {
  if (not exited) {
    exited = true;
    if (last_exit_code == EXIT_SUCCESS) {
      stats::time = global_timer.elapsed_secs();
      if (not do_exact_coloring) global_best_sure = false;
      check_criticality();
      output_to_file();
      print_stats();
    }
  }
  exit(last_exit_code);
}

void exit_fun_2() {
  if (not exited) exit_fun(0);
}

int main(int argc, char** argv) {
  signal(SIGINT, exit_fun);
  atexit(exit_fun_2);
  read_cmd_line(argc, argv);
  read_dimacs(input_filename);
  for (int i = 0; i < argc; ++i)
    pr("{} ", argv[i]);
  pr("\n");
  started = true;
  // For now, let's assume all input instances are already k-VCSs.
  //check_kvcs_init(); 
  while (not global_timer.timed_out()) {
    if (verb >= 1) pr("Repl. #{}\n", ++stats::repl_count);
    ibr(global_timer);
    if (not multistart) break;
  }
}
