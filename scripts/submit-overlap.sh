#! /bin/bash
#PBS -q regular
#PBS -l mppwidth=128
#PBS -l walltime=0:40:00
#PBS -A m1628
#PBS -M tua01780@temple.edu
#PBS -N PBE-OH-overlap
#PBS -m eab



cd $PBS_O_WORKDIR'/config292440'
../overlap-create.sh 3 64 65 94 186 187 197
cd $PBS_O_WORKDIR'/config292600'
../overlap-create.sh 1 12 65 67 77 190 191

