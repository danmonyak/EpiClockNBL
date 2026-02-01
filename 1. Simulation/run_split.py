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
import simulation as sim
from math import floor
from EpiClockInvasiveBRCA.src.simulation_util import *

###########################################################################
################################ Arguments ################################
###########################################################################
if len(sys.argv) == 1:
    sys.exit('Enter index of split...')
split_i = sys.argv[1]

split_name = f'{split_i}'

base_output_dir = '90_sites_NB_split_base'
target_cell_count_limit = int(55e5)
# target_nyears = 1.025

#########
sigma = 0.05842
digits = 1
delta_t_raw = 1/(2-sigma)
delta_t = floor(delta_t_raw * 10**digits) / 10**digits
#########
target_nyears = 1.025 / delta_t
###########################################################################
###########################################################################


prog_params = {
    'output_dir':os.path.join('90_sites_NB_split_splitOutputs', split_name),
    'nyears':target_nyears
}

print(f'Running simulation with the following parameters: {prog_params}')


os.makedirs(prog_params['output_dir'], exist_ok=True)

ensmbl = sim.loadFromDirectory(indir=os.path.join(base_output_dir, 'splits', split_name))
ensmbl.split = False

total_days = int(prog_params['nyears'] * 365)

beta_list = []   # hold beta values over time
n_cells_list = [] # hold # cells over time

beta_values_outfilepath = os.path.join(prog_params['output_dir'], 'beta_values.txt')
n_cells_outfilepath = os.path.join(prog_params['output_dir'], 'n_cells.txt')

if os.path.exists(beta_values_outfilepath):
    os.remove(beta_values_outfilepath)
if os.path.exists(n_cells_outfilepath):
    os.remove(n_cells_outfilepath)

split_jobs_filepath = os.path.join(prog_params['output_dir'], 'split_jobs.txt')
if os.path.exists(split_jobs_filepath):
    os.remove(split_jobs_filepath)

total_before = time()
i = ensmbl.day
k = 0    # i is day, k is iteration (i \neq j iff the simulation restarts)
while (not ensmbl.atCapacity()):
    if k == 1e9:     # don't loop forever
        sys.exit()

    n_cells = ensmbl.getNumCells()
    print(f'Day {i} -------- {n_cells} cells -------- {(time() - total_before) / 60:.1f} minutes')

    if (n_cells >= target_cell_count_limit) and (i >= (total_days + 1)):
        print(f'Surpassed {target_cell_count_limit:,} cells and surpassed {total_days} days.')
        print('Ending simulation.')
        break
    
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
    elif response and (response['result'] != 'success'):
        raise Exception('Tried to split...')
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

