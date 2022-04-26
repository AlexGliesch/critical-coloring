/*  Multi-neighborhood tabu search for the maximum weight clique problem
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * This program demonstrates the use of the heuristic algorithm for solving
 * the problem of the maximum weight clique problem. For more information about
 * this algorithm, please visit the website
 * http://www.info.univ-angers.fr/pub/hao/mnts_clique.html or email to:
 * wu@info-univ.angers.fr.
 */
#pragma once
#include <vector>

struct timer;

// Returns the maximum clique of size <= k in the vertex set s, wrt. global AM.
// If omega(s) >= k the algorithm returns a clique of size k or greater. Timer t
// specifies a time limit.
std::vector<int> hao_mnts_max_clique(const std::vector<int>& s, int k,
                                     const timer& t);
