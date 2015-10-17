#!/bin/bash
#SBATCH -J ewald
#SBATCH -p New
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 23:59:00

cd /home/ananduri/ewald

#for i in {3..6}
#do
#	./a $i 3 3 10
#done

./a 8 3 3 10
