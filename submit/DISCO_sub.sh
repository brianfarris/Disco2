#!/bin/sh
#DISCO_sub.sh
#Torque script to submit MPI program to run DISCO with milos parfile

#Torque directives
#PBS -N DISCO_sub.sh
#PBS -W group_list=hpcastro
#PBS -l nodes=4:ppn=8,walltime=00:12:00:00,mem=4000mb
#PBS -M djd2134@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o /hpc/astro/users/djd2134/codes/FLASH32/TRQout
#PBS -e /hpc/astro/users/djd2134/codes/FLASH32/TRQout/err

mpirun -n 32 bin/disco parfiles/milos_macfadyen.par > DISCO_out
