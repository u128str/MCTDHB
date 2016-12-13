#!/bin/bash
#$ -N MCTDHB
#$ -l mem=5000,disk=5000
#$ -cwd
#$ -o job_name.out
#$ -e job_name.err
#$ -pe shmem 16
./boson_MCTDHB_intel
