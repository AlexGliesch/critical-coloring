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
#include "pre.h"
#include "random.h"
inline int last_exit_code = EXIT_SUCCESS;
#define exit(x) \
  (exit)(last_exit_code = (x))
using uint = unsigned int;
using ii = pair<int, int>;
using bb = pair<bool, bool>;
using vi = vector<int>;
using vb = vector<bool>;
using nli = numeric_limits<int>;
using nld = numeric_limits<double>;
using vii = vector<ii>;
using vvi = vector<vi>;
using vvb = vector<vb>;
using vs = vector<string>;
using vd = vector<double>;
#define pr fmt::print
#define mp make_pair
#define mt make_tuple
using fmt::format;
using fmt::print;
struct ff {
  static constexpr double EPS = 1e-7;
  static bool fEq(double a, double b) { return abs(a - b) < EPS; }
  ff() = default;
  ff(double v) : v(v) {}
  bool operator==(const ff& o) const { return fEq(v, o.v); }
  bool operator!=(const ff& o) const { return not fEq(v, o.v); }
  bool operator<=(const ff& o) const { return v < o.v or fEq(v, o.v); }
  bool operator>=(const ff& o) const { return v > o.v or fEq(v, o.v); }
  bool operator<(const ff& o) const { return v + EPS < o.v; }
  bool operator>(const ff& o) const { return o < *this; }
  operator double() const { return v; }
  static void fixZero(double& v) {
    if (ff(v) == ff(0.0)) v = 0.0;
  }
  double v;
};
inline auto ffuple(double f) { return tuple(ff(f)); }
template <typename... Ts> auto ffuple(double f, Ts... v) {
  static_assert((is_same<Ts, double>::value && ...));
  return tuple_cat(tuple(ff(f)), ffuple(v...));
}
template <typename ContainerType>
void insert_to_sorted_vector(ContainerType& v,
                             const typename ContainerType::value_type& t,
                             bool allowDuplicates) {
  auto ub = upper_bound(begin(v), end(v), t);
  if (allowDuplicates or ub == begin(v) or *(ub - 1) != t) v.insert(ub, t);
}
template <typename ContainerType>
void erase_from_sorted_vector(ContainerType& v,
                              const typename ContainerType::value_type& t) {
  auto it = lower_bound(begin(v), end(v), t);
  if (it != end(v) and *it == t) v.erase(it);
}
template <typename T> class IsStreamable {
  template <typename U> static auto t(const U* u) -> decltype(std::cout << *u);
  static auto t(...) -> std::false_type;
public:
  enum { value = !std::is_same<decltype(t((T*)nullptr)), std::false_type>::value };
};
template <typename T>
typename enable_if<!is_void<decltype(begin(T()))>::value && !IsStreamable<T>::value,
                   ostream&>::type
operator<<(ostream& o, const T& v) {
  for (auto it = begin(v); it != end(v); ++it) {
    if (next(it) != end(v))
      o << *it << ' ';
    else
      o << *it;
  }
  return o;
}
template <typename T1, typename T2>
ostream& operator<<(ostream& o, const std::pair<T1, T2>& p) {
  o << "(" << p.first << ", " << p.second << ")";
  return o;
}
template <typename ContainerType, typename ValueType>
size_t index(const ContainerType& c, const ValueType& v) {
  return find(begin(c), end(c), v) - begin(c);
}
template <typename ContainerType, typename ValueType>
bool linear_in(const ContainerType& c, const ValueType& v) {
  return index(c, v) != size(c);
}
class timer {
  using clock = chrono::steady_clock;
  using timepoint = clock::time_point;
  timepoint tpstart;
  double tmlim;
public:
  timer(double time_lim_secs = nld::max()) { reset(time_lim_secs); }
  timer(double time_lim_secs, const timer& parent) {
    reset(min(time_lim_secs, parent.secs_left()));
  }
  void reset(double time_lim_secs = nld::max()) {
    tmlim = time_lim_secs, tpstart = clock::now();
  }
  double elapsed_secs() const {
    return chrono::duration_cast<chrono::duration<double>>(clock::now() - tpstart)
        .count();
  }
  double secs_left() const { return tmlim - elapsed_secs(); }
  bool timed_out() const { return elapsed_secs() >= tmlim; }
};
struct timeout_exception {
  timer t;
  timeout_exception(timer t) : t(t) {}
  const char* what() const {
    return format("Timeout exception (timer {}/{})", t.elapsed_secs(),
                  t.elapsed_secs() + t.secs_left())
        .c_str();
  }
};
inline double interp(double val, double min1, double max1, double min2, double max2) {
  return min2 + ((max2 - min2) * (val - min1)) / double(max1 - min1);
}
template <typename V> inline bool inrange(V v, V lo, V hi) { return v >= lo and v <= hi; }
inline bool is_unique(auto v) {
  sort(begin(v), end(v));
  return adjacent_find(begin(v), end(v)) == end(v);
}
namespace std {
template <typename T1, typename T2> struct hash<pair<T1, T2>> {
  size_t operator()(const pair<T1, T2>& p) const {
    size_t seed(0);
    boost::hash_combine(seed, p.first);
    boost::hash_combine(seed, p.second);
    return seed;
  }
};
template <typename T> struct hash<std::vector<T>> {
  size_t operator()(const vector<T>& p) const {
    size_t seed(0);
    for (auto& t : p)
      boost::hash_combine(seed, t);
    return seed;
  }
};
}
#define USE_TIMED_BLOCKS 
#ifdef NDEBUG
#undef USE_TIMED_BLOCKS
#endif
#ifdef USE_TIMED_BLOCKS
inline unordered_map<string, double> timedBlocks;
#endif
struct TimedBlock {
#ifdef USE_TIMED_BLOCKS
  TimedBlock(const string& name) : name(name) {}
  ~TimedBlock() { timedBlocks[name] += tm.elapsed_secs(); }
  string name;
  timer tm;
#else
  TimedBlock(const string&) {}
#endif
};
#define COMBINE1(X,Y) X ##Y
#define COMBINE(X,Y) COMBINE1(X, Y)
#define TIME_BLOCK(s) TimedBlock COMBINE(tbDummy, __LINE__)(s)
inline void print_timed_blocks() {
#ifdef USE_TIMED_BLOCKS
  vector<pair<double, string>> tbS;
  for (const auto& p : timedBlocks)
    tbS.emplace_back(-p.second, p.first);
  sort(begin(tbS), end(tbS));
  if (!timedBlocks.empty()) {
    vs s;
    size_t maxSize = 0;
    for (const auto& p : tbS) {
      s.push_back(fmt::format("{}: {}", p.second, -p.first));
      maxSize = std::max(maxSize, s.back().size());
    }
    fmt::print("{}\n| Timed blocks:{}|\n", string(maxSize + 4, '='),
               string(maxSize - 12, ' '));
    for (const auto& i : s)
      fmt::print("| {}{} |\n", i, string(maxSize - i.size(), ' '));
    fmt::print("{}\n", string(maxSize + 4, '='));
  }
#endif
}
struct tabu_list {
  tabu_list(int size, int tenure) : tabu(size, -tenure), iter(1), ten(tenure), sz(size) {}
  bool is_tabu(int u) const { return not tabu.empty() and tabu[u] + ten > iter; }
  void add(int u) { tabu[u] = iter; }
  void remove(int u) { tabu[u] = -ten; }
  void advance_iter() { ++iter; }
  int since(int u) const { return is_tabu(u) ? tabu[u] : nli::max(); }
  void reset() { tabu.assign(sz, -ten), iter = 1; }
  vi tabu;
  int iter, ten, sz;
};
inline string valOrNA(bool yes, double val) { return yes ? format("{}", val) : "NA"; }
inline string divOrNA(double num, double den) { return valOrNA(den, num / den); }
struct options_counter {
  int count = 0;
};
inline void validate(boost::any& v, std::vector<std::string> const&, options_counter*,
                     long) {
  if (v.empty())
    v = options_counter{1};
  else
    ++boost::any_cast<options_counter&>(v).count;
}
