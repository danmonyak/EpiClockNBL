#!/bin/bash
#SBATCH --job-name=job_scheduler_run_split
#SBATCH -o job_scheduler_run_split.out
#SBATCH -e job_scheduler_run_split.err
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G

n_jobs=50

while read LINE; do
	while [ $(squ | wc -l) -gt $n_jobs ]; do
		sleep 1
	done
    sbatch run_split.sh $LINE
done < 90_sites_NB_split_base/split_jobs.txt
