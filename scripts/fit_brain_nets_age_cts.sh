#!/bin/bash -l
#$ -cwd
#$ -N brain
#$ -o logs/brain.log
#$ -e logs/brain.err

module -f unload compilers mpi gcc-libs
module load r/recommended

echo $n_batches

R --no-save < /home/ucakble/Projects/multibergm-methods/scripts/fit_brain_nets_age_cts.R > logs/brain_cts_${n_nets}.log
