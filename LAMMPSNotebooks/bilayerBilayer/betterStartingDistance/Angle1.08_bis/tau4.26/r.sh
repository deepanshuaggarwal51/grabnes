#!/bin/sh
##PBS -l nodes=2:mem128:ppn=20
##PBS -l nodes=2:mem512:ppn=20
#PBS -l nodes=1:ppn=20
#PBS -q mem128
#PBS -N ChernFloating00.0002

##ulimit -s unlimited
#NPROCS=`wc -l < $PBS_NODEFILE`

hostname
date

#cd $PBS_O_WORKDIR

date

export OMP_NUM_THREADS=1

mpirun /home/nleconte/newCodes/lammps/build/lmp -in lammps.in

