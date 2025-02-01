#!/usr/bin/env bash
spp-dcj ilp -a 1 -m idmap.txt tree.txt adjacencies.txt > example.ilp
spp-dcj heuristic -a 1 tree.txt adjacencies.txt idmap.txt > example_init.sol
gurobi_cl InputFile=example_init.sol ResultFile=example.sol example.ilp 
spp-dcj sol2adj example.sol idmap.txt > resolved_adjacencies.txt
