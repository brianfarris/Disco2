#!/bin/sh
#DISCO_hot_sub.sh
#Torque script to submit MPI program to run DISCO with Cooling parfile

#Torque directives
#PBS -N DISCO_hot_sub.sh
#PBS -W group_list=hpcastro
#PBS -l nodes=2:ppn=16,walltime=00:04:00:00,mem=2000mb
#PBS -M djd2134@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o /hpc/astrostats/astro/users/djd2134/Howdy/Disco2/hpcout
#PBS -e /hpc/astrostats/astro/users/djd2134/Howdy/Disco2/hpcerr

mpirun -n 32 bin/disco parfiles/Cooling_Tests.par > DISCO_hot_out
