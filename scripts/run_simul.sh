#!/bin/bash -l
#$ -cwd
#$ -N simul
#$ -o logs/simul.log
#$ -e logs/simul.err

module -f unload compilers mpi gcc-libs
module load r/recommended

echo $n_batches

R --no-save < /home/ucakble/Projects/multibergm-methods/scripts/run_simul.R > logs/simul_${n_nets}_${n_nodes}.log
