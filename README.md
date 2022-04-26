# Supplementary material for *A new heuristic for finding verifiable k-vertex-critical subgraphs*

This repository holds the source code, detailed result tables and instances for the following paper:

> **A new heuristic for finding verifiable k-vertex-critical subgraphs**<br>
> Alex Gliesch, Marcus Ritt<br>
> Journal of Heuristics 28, pp. 61–91 <br>
> https://doi.org/10.1007/s10732-021-09487-9

**Bibtex**

```bibtex
@article{GlieschRitt/2022,
  title     = {A new heuristic for finding verifiable k-vertex-critical subgraphs},
  author    = {Gliesch, Alex and Ritt, Marcus},
  journal   = {Journal of Heuristics},
  volume    = {28},
  pages     = {61--91},
  year      = {2022},
  publisher = {Springer},
  doi       = {10.1007/s10732-021-09487-9}
}
```

Do use the reference above if you use this material in your research.

## Running the code 

1. Unpack the instances in `instances.tar.gz`.
2. Compile the code under `src` using `make`. Make sure you have a C++17 compatible compiler.
3. Run using `./critcol -i {instance} -k {numColors} -t {timeLimit} -o {outFile}`. Add `-v` for increased verbosity. For more options, see `--help`. 

Our reimplementation of Sun et al. (2017)'s IBR algorithm is in `src/ibr`, see README there.

## Instance generator for imperfect graphs (with high probability) 

1. Compile the code under `src` using `make`. 
2. Run `./generateimperfectgraph {numVertices} {density} {seed} {outFile}`.
3. The first line of the output DIMACS file is a commentary containing the (heuristic) chromatic and clique numbers for the generated graph.