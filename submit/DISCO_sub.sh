#!/bin/sh
#DISCO_sub.sh
#Torque script to submit MPI program to run DISCO with milos parfile

#Torque directives
#PBS -N DISCO_sub.sh
#PBS -W group_list=yetiastro
#PBS -l nodes=4:ppn=8,walltime=00:06:00:00,mem=32000mb
#PBS -M djd2134@columbia.edu
#PBS -m abe
#PBS -V

#set output and error directories
#PBS -o /vega/astro/users/djd2134/Disco2/hpcout
#PBS -e /vega/astro/users/djd2134/Disco2/hpcerr

mpirun -n 32 bin/disco parfiles/milos_macfadyen.par > DISCO_out
