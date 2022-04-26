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
#if 0
#include "sm.h"
#include "main.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
vd get_eigen(vd& data, int sz) {
  assert((int)data.size() == sz * sz);
  vd eigen;
  gsl_matrix_view m = gsl_matrix_view_array(data.data(), sz, sz);
  gsl_vector_complex* eval = gsl_vector_complex_alloc(sz);
  gsl_matrix_complex* evec = gsl_matrix_complex_alloc(sz, sz);
  gsl_eigen_nonsymmv_workspace* w = gsl_eigen_nonsymmv_alloc(sz);
  gsl_eigen_nonsymmv(&m.matrix, eval, evec, w);
  gsl_eigen_nonsymmv_free(w);
  gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
  for (int i = 0; i < 4; ++i) {
    gsl_complex eval_i = gsl_vector_complex_get(eval, i);
    if (ff(GSL_IMAG(eval_i)) == ff(0.0)) eigen.push_back(GSL_REAL(eval_i));
  }
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  return move(eigen);
}
double energy(const vi& ss) {
  TIME_BLOCK("energy");
  int sz = ss.size();
  vd adj_mat(sz * sz);
  for (int i = 0; i < sz; ++i)
    for (int j = 0; j < sz; ++j)
      adj_mat[i * sz + j] = AM[ss[i]][ss[j]];
  vd eig = get_eigen(adj_mat, sz);
  double en = accumulate(begin(eig), end(eig), 0.0) / double(eig.size());
  return en;
}
double algebraic_connectivity(const vi& ss) {
  TIME_BLOCK("algebraic_connectivity");
  int sz = ss.size();
  vd lap_mat(sz * sz);
  for (int i = 0; i < sz; ++i)
    for (int j = 0; j < sz; ++j)
      lap_mat[i * sz + j] = (i == j) ? size(AL[ss[i]]) : -AM[ss[i]][ss[j]];
  vd eig = get_eigen(lap_mat, sz);
  sort(begin(eig), end(eig));
  if (eig.size() <= 1)
    return nli::min();
  else
    return eig[1];
}
#endif
