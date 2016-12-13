#!/bin/sh
#PBS -N MCTDHB-1D-N2-M8-PROP
#PBS -V 
#PBS -l nodes=1:ppn=32
#PBS -l walltime=24:00:00
echo Working directory is $PBS_O_WORKDIR >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Running on host `hostname` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Time is `date` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
echo Directory is `pwd` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME 
echo This jobs runs on the following processors: >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME 
NODES=`cat $PBS_NODEFILE`
echo $NODES >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME
cd $PBS_O_WORKDIR
echo Directory is `pwd` >> $PBS_O_WORKDIR/environment_$PBS_JOBID_$PBS_JOBNAME 
export OMP_NUM_THREADS=8
aprun -n 4 -d 8 $PBS_O_WORKDIR/boson_MCTDHB_intel >> OUTPUT_$PBS_JOBNAME

