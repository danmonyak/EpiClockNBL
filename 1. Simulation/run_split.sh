#!/bin/bash
#SBATCH --job-name=run_split
#SBATCH -o output_files/run_split-%j.out
#SBATCH -e error_files/run_split-%j.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64G

python -u run_split.py $1
