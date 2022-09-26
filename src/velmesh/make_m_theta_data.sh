#!/bin/sh

# velmesh
parallel ~/velmesh -m {1} -t {2} {3} -S    -L {4} :::: deltaMs :::: thetas ::: -N -M :::: list_of_lists
parallel ~/velmesh -m {1} -t {2} {3} -S -c -L {4} :::: deltaMs :::: thetas ::: -N -M :::: list_of_lists

# coord
#parallel ~/exact_sim.exe -m {1} -t {2} -V {3} -N -S :::: deltaMs :::: thetas :::: velocity.list
