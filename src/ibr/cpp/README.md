C++ reimplementation of Sun et. al's algorithm "Iterated backtrack removal search for finding k-vertex-critical subgraphs" (doi: 10.1007/s10732-017-9358-5).

Inside directory "ibr/jl" you can find a partial reimplementation in Julia.

Some notes:

* The perturbation step (Section 3.5 in the paper) is unclear in the original paper, deterministic (when it probably wasn't meant to be), and with ill-defined parameter values, namely eta and p. For this reason, in this implementation I chose to use the similar approach of Glover et al. (doi: 10.1007/s10288-009-0115-y) which selects a number of items (I used \lambda=|A|/2) probabilistically to be moved from A to B. See Glover's paper for details.

* We use HEA* and BTDSatur to solve colorings.
