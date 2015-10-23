#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 03:59:00

cd /home/ananduri/ewald/fcc

prod/a $1 1 10 6 $2 $3
