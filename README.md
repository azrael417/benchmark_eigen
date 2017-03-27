# benchmark_eigen
This is an Eigensolver benchmark for DFT calculations.

Use the Makefile to compile (requires intel compiler) and run it with:

eigen.x nprow npcol N mode,

where nprow are the number of processes per row, npcol the number of processes per column, N the rank of the NxN matrix being diagonalized and mode is the function called for diagonalization. Currently, only PDSYEVD is supported.
