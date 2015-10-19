#!/bin/bash
#SBATCH -J ewald
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 23:59:00

cd /home/ananduri/ewald

for i in {1..101}
do
	./a 1 $i 50 50 >> out.txt
done
