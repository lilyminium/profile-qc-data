#!/bin/bash
#SBATCH -J benchmark
#SBATCH -p standard
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mem=64GB
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out

date
hostname

source ~/.bashrc
conda activate yammbs

# Copy the force field files to the current directory


#python benchmark-forcefield.py         \
#	-i "../offqcdata/data/yammbs-unsafe/optimizations.sqlite" \
#	-np 32

python benchmark-forcefield.py         \
        -i "../offqcdata/data/yammbs-unsafe/torsiondrives.sqlite" \
        -np 32
date
date
