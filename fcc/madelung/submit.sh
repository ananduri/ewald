#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 03:59:00

cd /home/ananduri/ewald/fcc/madelung

prod/a $1 $2 $3 $4 
