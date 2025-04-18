#!/bin/bash
#SBATCH -J convert
#SBATCH -p standard
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --mem=32GB
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out

date
hostname

source ~/.bashrc
conda activate yammbs

# Copy the force field files to the current directory
#python convert-to-yammbs.py

# python convert-optimization-to-yammbs-unsafe.py


date


python convert-torsiondrive-to-yammbs-unsafe.py


date
