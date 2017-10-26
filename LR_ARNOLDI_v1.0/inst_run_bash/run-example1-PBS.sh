#!/bin/bash
#PBS -l nodes=4:ppn=8:bwgrid
#PBS -l walltime=23:59:00
cd $PBS_O_WORKDIR
source ~/.bashrc
echo '======================= All ===================================='
qstat -f $PBS_JOBID
echo '======================= Hosts ===================================='
qstat -n $PBS_JOBID
echo '======================= Run ===================================='
pwd
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=6
#export PSC_OMP_AFFINITY=TRUE
echo "MPIEXEC IS FROM"
which mpiexec
time mpiexec  -n 4 -npernode 1 -x KMP_STACKSIZE=16m  -x OMP_NUM_THREADS=6 -x MKL_NUM_THREADS=1  ./boson_MCTDHB_intel  > out.out 2>err.out
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo '======================= Done  ===================================='
