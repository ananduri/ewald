#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 00:59:00

cd /home/ananduri/ewald/fcc

prod/a 2 $1 2 40
