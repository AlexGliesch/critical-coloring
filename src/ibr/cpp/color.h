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

// Whether the subgraph induced by ss \in V is k-colorable. If t times out,
// the function returns true.
bool is_k_colorable_exact(int k, const vi& ss, timer t, vi* res = nullptr);

// Whether the subgraph induced by ss \in V is heuristically k-colorable. If t
// times out, returns true.
bool is_k_colorable_heuristic(int k, const vi& ss, timer t, vi* res = nullptr);

// Whether the subgraph induced by ss \in V is k-colorable. Return value:
// .first: whether it is k-colorable or not, and .second: whether we are sure of
// this result. If a k-coloring is found, we are always sure of it. Otherwise,
// if a k-coloring cannot be found, it could be because the exact method failed
// (here, we are sure of the result), or multiple calls to the heuristic method
// with a time limit failed (in this case, we cannot be sure one does not
// exist). If the method times out, it a pair (true, false).
bb check_colorability(int k, const vi& ss, timer t, vi* res = nullptr);

// Simple wrapper around check_colorability, to be more readable; a subset is a
// k-vcs if it is not k-1 colorable.
bb is_k_vcs(int k, const vi& ss, timer t, vi* res = nullptr);