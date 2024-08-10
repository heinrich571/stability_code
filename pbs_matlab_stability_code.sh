#!/bin/sh
#PBS -N stagnation_point_flow_stability_tests_batch_run
#PBS -q karp_q
#PBS -M dan.heinrich@campus.technion.ac.il
#PBS -mbea
#PBS -l select=1:ncpus=160:mpiprocs=160
#PBS -l place=scatter

# working directory:
cd $PBS_O_WORKDIR

/usr/local/bin/matlab -nodesktop  < main.m  > output
 
