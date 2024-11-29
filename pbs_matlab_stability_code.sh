#!/bin/sh
#PBS -N stagnation_point_flow_stability_tests_batch_run
#PBS -q karp_q
#PBS -M dan.heinrich@campus.technion.ac.il
#PBS -mbea
#PBS -l select=1:ncpus=160:mpiprocs=80
#PBS -l place=scatter

# working directory:
cd $PBS_O_WORKDIR

echo running on nodes
cat $PBS_NODEFILE | uniq

matlab -nodesktop  < batch_run.m  > output

echo done
 
