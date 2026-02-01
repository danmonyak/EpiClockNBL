#!/bin/bash
#SBATCH --job-name=run_base_ensemble
#SBATCH -o run_base_ensemble.out
#SBATCH -e run_base_ensemble.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64G

python -u run_base_ensemble.py
