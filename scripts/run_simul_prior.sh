#!/bin/bash -l
#$ -cwd
#$ -N simul
#$ -o logs/simul.log
#$ -e logs/simul.err
#$ -v "decay=${DECAY}"

module -f unload compilers mpi gcc-libs
module load r/recommended

echo $n_batches

R --no-save < /home/ucakble/Projects/multibergm-methods/scripts/run_simul_prior.R > logs/simul_prior${n_nets}_${cov_scale}.log
