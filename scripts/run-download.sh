#!/bin/bash
#SBATCH -J download
#SBATCH -p standard
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
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
python download-qcsubmit.py


date
