#!/bin/bash -l
#$ -cwd
#$ -N simul
#$ -o logs/simul_2grp.log
#$ -e logs/simul_2grp.err

module -f unload compilers mpi gcc-libs
module load r/recommended

echo $n_batches

R --no-save < /home/ucakble/Projects/multibergm-methods/scripts/run_simul_two_group.R > logs/simul_2grp_${n_nets}.log
