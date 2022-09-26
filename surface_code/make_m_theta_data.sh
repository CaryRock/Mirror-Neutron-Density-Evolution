#!/bin/sh

parallel ~/exact_sim.exe -m {1} -t {2} -N -S :::: deltaMs :::: thetas

#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-8 0.5E-8 9.5E-8) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-7 0.5E-7 9.5E-7) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-6 0.5E-6 9.5E-6) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-5 0.5E-5 9.5E-5) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-4 0.5E-4 9.5E-4) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-3 0.5E-3 9.5E-3) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-2 0.5E-2 9.5E-2) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E-1 0.5E-1 9.5E-1) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E00 0.5E00 9.5E00) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E01 0.5E01 9.5E01) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E02 0.5E02 9.5E02) ::: $(seq 0.0 0.001 0.78)
#parallel ~/exact_sim.exe -m {1} -t {2} ::: $(seq 1.0E03 0.5E03 9.5E03) ::: $(seq 0.0 0.001 0.78)
