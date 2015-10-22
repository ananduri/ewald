#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 72:59:00

cd /home/ananduri/ewald/fcc

prod/a $1 5 4 18
