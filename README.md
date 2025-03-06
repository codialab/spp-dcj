[![Anaconda-Server Badge](https://anaconda.org/bioconda/spp-dcj/badges/version.svg)](https://anaconda.org/bioconda/spp-dcj) [![Anaconda-Server Badge](https://anaconda.org/bioconda/spp-dcj/badges/platforms.svg)](https://anaconda.org/bioconda/spp-dcj) [![Anaconda-Server Badge](https://anaconda.org/bioconda/spp-dcj/badges/license.svg)](https://anaconda.org/bioconda/spp-dcj)

# SPP-DCJ

## Introduction
The Small Parsimony Problem (SPP) aims at finding the gene orders at internal nodes of a given phylogenetic tree such that the overall genome rearrangement distance along the tree branches is minimized. This problem is intractable in most genome rearrangement models, especially when gene duplication and loss are considered.
`SPP_DCJ` is _Integer Linear Program_-based algorithm to solve the SPP for natural genomes, _i.e._, genomes that contain conserved, unique, and duplicated markers. The evolutionary model that we consider is the _DCJ-indel model_ that includes the Double-Cut and Join rearrangement operation and the insertion and deletion of genome segments. 

`SPP_DCJ` is an extension of [DING](https://gitlab.ub.uni-bielefeld.de/gi/ding).

## Installation

### From bioconda channel

Make sure you have [conda](https://conda.io)/[mamba](https://anaconda.org/conda-forge/mamba) installed!

```
mamba install -c conda-forge -c bioconda spp-dcj
```

### From repository

```
# clone repository
git clone https://github.com/codialab/spp-dcj.git 

# install with pip
pip install spp-dcj/
```

## Dependencies
- [Gurobi](https://www.gurobi.com/products/gurobi-optimizer/) 
- [Python 3](https://www.python.org/downloads/)
- Python 3 libraries:
    - Matplotlib
    - NetworkX

## How to run 

The following steps show howto run `SPP_DCJ` with Gurobi.

0. `SPP_DCJ` requires (i) a given phylogeny and (ii) a table with candidate adjacencies for all genomes corresponding to nodes of the given phylogeny. The candidate adjacency table has the following columns:

    ```#Species     Gene_1  Ext_1   Species Gene_2  Ext_2   Weight```

    The phylogeny must be given as table format:

    ```#Child	Parent```

1. Generate ILP (`-a` is a parameter of the objective function that is set here to 1):
    
    ```spp-dcj ilp -a 1 -m example/idmap.txt example/tree.txt example/adjacencies.txt > example/example.ilp```

2. Compute heuristic solution for warm starting the solver (optional)

    ```spp-dcj heuristic -a 1 example/tree.txt example/adjacencies.txt example/idmap.txt > example/example_init.sol```

3. Run ILP (you may omit the `InputFile` argument in case you don't want to warm-start the solver)
    
    ```gurobi_cl InputFile=example/example_init.sol ResultFile=example/example.sol example/example.ilp```

4. Extract adjacencies from solution

    ```spp-dcj sol2adj example/example.sol example/idmap.txt > example/resolved_adjacencies.txt```

5. Visualize candidate and predicted adjacencies

    ```spp-dcj draw -i example/resolved_adjacencies.txt example/adjacencies.txt > example/adjacencies.pdf```

## Example
The code above is included in the directory in `example`.


