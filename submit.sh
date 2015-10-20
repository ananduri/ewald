#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -n 1
#SBATCH -t 08:59:00

cd /home/ananduri/ewald

prod/a $1 3 3 10
