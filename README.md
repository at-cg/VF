VF
========================================================================

This is an open-source prototype software implementation of our algorithms for variant filtering in variation graphs. Its purpose is to reduce graph size while preserving subsequent mapping accuracy. The graph size reduction naturally benefits read-to-graph mapping algorithms. The intuition of the proposed graph reduction framework is following. Consider a complete variation graph constructed from a set of given haplotypes. Any substring of a haplotype has a corresponding path in the complete variation graph. Not including some variants will introduce errors in the corresponding paths. If the number of such errors is matched with the error tolerance built into sequence-to-graph mapping algorithms, the same identical paths can still be found. Mathematically, this problem is phrased in terms of minimizing variation graph size subject to preserving paths of length α with at most δ differences. See our [preprint](https://www.doi.org) for more details. 

## Dependencies
- A C++ compiler with c++11 support, e.g., GNU g++ (version 5+)
- [vcftools](https://vcftools.github.io/)
- [Gurobi](https://www.gurobi.com)
- [clipp](https://github.com/muellan/clipp)

## Installation
The above dependencies can be handled by running script `dependencies.sh`.
```sh
git clone https://github.com/NedaTavakoli/VF.git
cd VF
./dependencies.sh
make
```

After a successful compilation, expect three executables `greedy_snp`, `lp_snp`, `greedy_sv_indels` in a directory named `build`. 

## Usage
The three executables implement different algorithms to achieve variant graph size reduction, but they all have a similar interface.
```
SYNOPSIS
        greedy_snp       -a <alpha> -d <delta> -vcf <file> -chr <id>
        lp_snp           -a <alpha> -d <delta> -vcf <file> -chr <id>
        greedy_sv_indels -a <alpha> -d <delta> -vcf <file> -chr <id>


OPTIONS
        <alpha>     path length in variation graph (e.g., 500)
        <delta>     differences allowed (e.g., 10)
        <file>      uncompressed vcf file (something.vcf)
        <id>        chromosome id (e.g., 1 or chr1), make it consistent with vcf file
```

A few [example runs](examples) are made available for user's reference. In practice, α should be a function of read lengths whereas δ is determined based on sequencing errors and error-tolerance of read-to-graph mapping algorithms. NOTE: At runtime, `lp_snp` executable might complain if you don't have a valid Gurobi license file. It is straight-forward and free to get one for academic use [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic).

## Benchmark

We evaluated the magnitude of graph reduction achieved in human chromosome variation graphs using VF (v1.0) with multiple α and δ parameter values corresponding to short and long-read resequencing characteristics. When our algorithm is run with parameter settings amenable to long-read mapping (α=10 kbp, δ=1000), 99.99% SNPs and 69% indel structural variants could be safely excluded from human chromosome 1 variation graph.

## Publications

- **Chirag Jain, Neda Tavakoli and Srinivas Aluru**. "[A Theoretical Framework for Variant Selection in Graph Representation of Genomes](https://www.doi.org)". *bioRxiv 2020*.
