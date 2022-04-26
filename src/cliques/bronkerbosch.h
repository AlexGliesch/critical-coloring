#pragma once
#include "subgraph.h"
#include "util.h"

// Returns .first: the maximal clique found by Bron-Kerbosch under the time
// limit t, and .second: whether timer t timed out (meaning .first is optimal).
// Also stops as soon as a clique size >= max_k is found.
pair<vi, bool> bron_kerbosch(const subgraph& h, timer t,
														 int max_k = nli::infinity());
