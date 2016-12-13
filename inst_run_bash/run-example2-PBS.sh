#!/bin/sh
#PBS -N MCTDHB-2D-N4
#PBS -l nodes=24:ppn=6,walltime=16:00:00
echo Working directory is $PBS_O_WORKDIR >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Running on host `hostname` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Time is `date` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Directory is `pwd` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME 
echo This jobs runs on the following processors: >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME 
NODES=`cat $PBS_NODEFILE`
echo $NODES >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
cd $PBS_O_WORKDIR
mpirun -np 24 -npernode 1 -x OMP_NUM_THREADS=8 -x KMP_AFFINITY=verbose,compact,1,0 -x KMP_STACKSIZE=16777216 boson_MCTDHB_intel 

#cp input.l1800 input.in
#mpirun -np 12 -npernode 1 -x OMP_NUM_THREADS=6 -x KMP_AFFINITY=verbose,compact -x KMP_STACKSIZE=16777216 boson_MCTDHB_intel
#mkdir l1800
#cp *bin l1800
