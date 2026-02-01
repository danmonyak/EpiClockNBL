#!/bin/bash
#SBATCH --job-name=move_and_compress
#SBATCH -o move_and_compress.out
#SBATCH -e move_and_compress.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=64G

rm -r send/*
cp -r 90_sites_NB_split_* send/
rm -r send/90_sites_NB_split_base/splits
zip -r send.zip send
