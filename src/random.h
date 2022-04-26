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
#include <unistd.h>
#include <ctime>
#include <random>
using namespace std;
inline mt19937 rng;
inline int rand_int(int from, int to) {
 static uniform_int_distribution<int> d;
 return d(rng, decltype(d)::param_type{from, to});
}
inline double rand_double(double from, double to) {
 static std::uniform_real_distribution<double> d;
 return d(rng, decltype(d)::param_type{from, to});
}
inline bool rand_bool() { return (bool)rand_int(0, 1); }
inline size_t unique_random_seed() {
 size_t a = (size_t)clock(), b = (size_t)time(nullptr), c = (size_t)getpid();
 a = (a - b - c) ^ (c >> 13);
 b = (b - c - a) ^ (a << 8);
 c = (c - a - b) ^ (b >> 13);
 a = (a - b - c) ^ (c >> 12);
 b = (b - c - a) ^ (a << 16);
 c = (c - a - b) ^ (b >> 5);
 a = (a - b - c) ^ (c >> 3);
 b = (b - c - a) ^ (a << 10);
 c = (c - a - b) ^ (b >> 15);
 return c;
}
template <typename T>
inline void random_choice(const vector<T>& v, vector<T>& result, int k) {
 int n = v.size();
 result.resize(k);
 for (int i = 0; i < k; ++i)
  result[i] = v[i];
 for (int i = k + 1; i < n; ++i) {
  int j = rand_int(0, i);
  if (j < k) result[j] = v[i];
 }
}
struct reservoir_sampling {
 bool consider() { return (1.0 / ++num) >= rand_double(0.0, 1.0); }
 void reset(double num = 0.0) { this->num = num; }
 double num = 0.0;
};
template<typename T>
std::vector<T> random_sample(uint k, const std::vector<T> v) {
  std::vector<int> indices(k);
  iota(begin(indices), end(indices), 0);
  shuffle(begin(indices), end(indices), rng);
  std::vector<T> res;
  for (int i : indices) res.push_back(v[i]);
  return res;
}
