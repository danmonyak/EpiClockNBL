"""
run_ensemble.py
=======
Author - Daniel Monyak
12-27-25
=======

Source code for running 
    
"""

import numpy as np
import sys
import os
from time import time
from math import floor
import simulation as sim
from EpiClockInvasiveBRCA.src.simulation_util import *

###########################################################################
################################ Arguments ################################
###########################################################################
split_limit = int(1e7)
n_split = 200
# split_limit = int(1e5)
# n_split = 10

sigma = 0.05842
digits = 1
delta_t_raw = 1/(2-sigma)
delta_t = floor(delta_t_raw * 10**digits) / 10**digits
total_time_years_param = 1
###########################################################################
###########################################################################


prog_params_list = [
    {'output_dir':'90_sites_NB_split_base', 'n_CpGs_each':30,
               # 'flip_rate':(0.0042 / .17),
               'flip_rate':0.0075,
               'death_rate':0.11158, 'nyears':1, 'seed':0},
    {'output_dir':'90_sites_NB_split_base', 'n_CpGs_each':30,
               # 'flip_rate':(0.0042 / .17),
               'prolif_rate':delta_t,
               'flip_rate':1e-3,
               'death_rate':(1-sigma) * delta_t,
               'nyears':total_time_years_param / delta_t,
               'delta_t':delta_t,
               'seed':0}
]

prog_params = prog_params_list[1]

print(f'Running simulation with the following parameters: {prog_params}')

# Constant parameters
FLIP_RATE = 0.002
PROLIF_RATE = 0.17

if 'flip_rate' not in prog_params:
    prog_params['flip_rate'] = FLIP_RATE
if 'prolif_rate' not in prog_params:
    prog_params['prolif_rate'] = PROLIF_RATE

os.makedirs(prog_params['output_dir'], exist_ok=True)

split_jobs_filepath = os.path.join(prog_params['output_dir'], 'split_jobs.txt')
if os.path.exists(split_jobs_filepath):
    os.remove(split_jobs_filepath)


init_params = {
    'flip_rate':prog_params['flip_rate'], # flip rate per cell division per allele
    'prolif_rate':prog_params['prolif_rate'], # cell divisions per day
    'death_rate':prog_params['death_rate'], # cell deaths per day,
    'init_site_state_counts':[prog_params['n_CpGs_each'], prog_params['n_CpGs_each'], 0, prog_params['n_CpGs_each']],
    }
gen = np.random.default_rng(prog_params['seed'])
ensmbl = sim.Ensemble(init_params, gen,
                      split=True, split_limit=split_limit, n_split=n_split,
                      save_info = {
                          'saveAsDirectory_parent_outdir':os.path.join(prog_params['output_dir'], 'splits'),
                          'split_jobs_filepath':split_jobs_filepath
                      }
                     )

total_days = int(prog_params['nyears'] * 365)

beta_list = []   # hold beta values over time
n_cells_list = [] # hold # cells over time

beta_values_outfilepath = os.path.join(prog_params['output_dir'], 'beta_values.txt')
n_cells_outfilepath = os.path.join(prog_params['output_dir'], 'n_cells.txt')

if os.path.exists(beta_values_outfilepath):
    os.remove(beta_values_outfilepath)
if os.path.exists(n_cells_outfilepath):
    os.remove(n_cells_outfilepath)

total_before = time()
i = k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
while (not ensmbl.atCapacity()) and (i < total_days+1):
    if k == 1e9:     # don't loop forever
        sys.exit()

    n_cells = ensmbl.getNumCells()
    print(f'Day {i} -------- {n_cells} cells -------- {(time() - total_before) / 60:.1f} minutes')
    
    # Print progress
    if (i > 0) and (i % 50 == 0):
        print('Flushing results to files...')

        writeBetaValues(beta_list, beta_values_outfilepath)
        writeNcells(n_cells_list, n_cells_outfilepath)

    # Attach day to data
    beta_values_day = np.concatenate([[i], ensmbl.getBetaValues()])
    beta_list.append(beta_values_day)

    # Attach day to data
    n_cells_day = [i, ensmbl.getNumCells()]
    n_cells_list.append(n_cells_day)

    ##################################
    ############ Pass day ############
    ##################################
    response = ensmbl.passDay()
    
    if response and (response['result'] == 'success'): # passDay was successful
        i += 1
    elif response and (response['result'] == 'split'):
        ens_splits = response['data']
        print(f'Split into {ensmbl.n_split} ensembles.')
        print('#'*50)
        # for i, ens in enumerate(ens_splits):
        #     print(f'Writing ensemble {i}...')
        #     split_name = f'split_{i}'
        #     ens.saveAsDirectory(parent_outdir=os.path.join(prog_params['output_dir'], 'splits'),
        #                         outdir=split_name)
        #     writeLine(split_jobs_filepath, split_name)
        print('All ensembles written')
        break
    # had to restart (e.g. all cells died)
    else:
        i = 0

        resetObjsDeleteFiles(
            ensmbl,
            [beta_list, n_cells_list],
            [beta_values_outfilepath, n_cells_outfilepath]
        )

    k += 1


after_total = time()
print(f'Total time: {after_total - total_before:.1f} seconds')

writeBetaValues(beta_list, beta_values_outfilepath)
writeNcells(n_cells_list, n_cells_outfilepath)

